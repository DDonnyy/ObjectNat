from typing import Literal

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from osmnx import graph_to_gdfs
from pyproj import CRS
from scipy.spatial import KDTree
from shapely import Point
from shapely.ops import unary_union

from objectnat import config

logger = config.logger


def get_accessibility_isochrones(
    points: gpd.GeoDataFrame,
    weight_value: float,
    weight_type: Literal["time_min", "length_meter"],
    graph_nx: nx.Graph,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Calculate accessibility isochrones from a gpd.GeoDataFrame based on the provided city graph.

    Isochrones represent areas that can be reached from a given point within a specific time or distance,
    using a graph that contains road and transport network data.

    Parameters
    ----------
    points : gpd.GeoDataFrame
        A GeoDataFrame containing the geometry from which accessibility isochrones should be calculated.
        The CRS of this GeoDataFrame must match the CRS of the provided graph.
    weight_value : float
        The maximum distance or time threshold for calculating isochrones.
    weight_type : Literal["time_min", "length_meter"]
        The type of weight to use for distance calculations. Either time in minutes ("time_min") or distance
        in meters ("length_meter").
    graph_nx : nx.Graph
        A NetworkX graph representing the city network.
        The graph must contain the appropriate CRS and, for time-based isochrones, a speed attribute.

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        A tuple containing:
        - isochrones : GeoDataFrame with the calculated isochrone geometries.
        - public transport stops (if applicable) : GeoDataFrame with public transport stops within the isochrone, or None if not applicable.
        - public transport routes (if applicable) : GeoDataFrame with public transport routes within the isochrone, or None if not applicable.

    Examples
    --------
    >>> from iduedu import get_intermodal_graph
    >>> graph = get_intermodal_graph(polygon=my_territory_polygon)
    >>> points = gpd.GeoDataFrame(geometry=[Point(30.33, 59.95)], crs=4326).to_crs(graph.graph['crs'])
    >>> isochrones, pt_stops, pt_routes = get_accessibility_isochrones(points, weight_value=15, weight_type="time_min", graph_nx=my_graph)

    """

    assert points.crs == CRS.from_epsg(
        graph_nx.graph["crs"]
    ), f'CRS mismatch , points.crs = {points.crs.to_epsg()}, graph["crs"] = {graph_nx.graph["crs"]}'

    nodes_with_data = list(graph_nx.nodes(data=True))
    logger.info("Calculating isochrones distances...")
    coordinates = np.array([(data["x"], data["y"]) for node, data in nodes_with_data])
    tree = KDTree(coordinates)

    target_coord = [(p.x, p.y) for p in points.representative_point()]
    distances, indices = tree.query(target_coord)

    nearest_nodes = [nodes_with_data[idx][0] for idx in indices]
    del nodes_with_data
    dist_nearest = pd.DataFrame(data=distances, index=nearest_nodes, columns=["dist"])
    speed = 0
    if graph_nx.graph["type"] in ["walk", "intermodal"] and weight_type == "time_min":
        try:
            speed = graph_nx.graph["walk_speed"]
        except KeyError:
            logger.warning("There is no walk_speed in graph, set to the default speed - 83.33 m/min")
            speed = 83.33
        dist_nearest = dist_nearest / speed
    elif weight_type == "time_min":
        speed = 20 * 1000 / 60
        dist_nearest = dist_nearest / speed

    if (dist_nearest > weight_value).all().all():
        raise RuntimeError(
            "The point(s) lie further from the graph than weight_value, it's impossible to "
            "construct isochrones. Check the coordinates of the point(s)/their projection"
        )

    data = []
    for source in nearest_nodes:
        dist, path = nx.single_source_dijkstra(graph_nx, source, weight=weight_type, cutoff=weight_value)
        for target_node, way in path.items():
            source = way[0]
            distance = dist.get(target_node, np.nan)
            data.append((source, target_node, distance))
    del dist
    dist_matrix = pd.DataFrame(data, columns=["source", "destination", "distance"])
    del data
    dist_matrix = dist_matrix.pivot_table(index="source", columns="destination", values="distance", sort=False)

    dist_matrix = dist_matrix.add(dist_nearest.dist, axis=0)
    dist_matrix = dist_matrix.mask(dist_matrix >= weight_value, np.nan)
    dist_matrix.dropna(how="all", inplace=True)

    results = []
    logger.info("Building isochrones geometry...")
    for _, row in dist_matrix.iterrows():
        geometry = []
        for node_to, value in row.items():
            if not pd.isna(value):
                node = graph_nx.nodes[node_to]
                point = Point(node["x"], node["y"])
                geometry.append(
                    point.buffer(round((weight_value - value) * speed * 0.8, 2))
                    if weight_type == "time_min"
                    else point.buffer(round((weight_value - value) * 0.8, 2))
                )
        geometry = unary_union(geometry)
        results.append(geometry)

    isochrones = gpd.GeoDataFrame(data=points, geometry=results, crs=graph_nx.graph["crs"])
    isochrones["weight_type"] = weight_type
    isochrones["weight_value"] = weight_value

    isochrones_subgraph = graph_nx.subgraph(dist_matrix.columns)
    nodes = pd.DataFrame.from_dict(dict(isochrones_subgraph.nodes(data=True)), orient="index")
    if "desc" in nodes.columns and "stop" in nodes["desc"].unique():
        pt_nodes = nodes[nodes["desc"] == "stop"]
        nodes, edges = graph_to_gdfs(isochrones_subgraph.subgraph(pt_nodes.index))
        nodes.reset_index(drop=True, inplace=True)
        nodes = nodes[["desc", "route", "geometry"]]
        edges.reset_index(drop=True, inplace=True)
        edges = edges[["type", "route", "geometry"]]
        return isochrones, nodes, edges
    return isochrones, None, None
