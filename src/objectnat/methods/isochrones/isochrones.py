from typing import Literal

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from pyproj.exceptions import CRSError
from shapely import Point
from shapely.ops import unary_union

from objectnat import config
from objectnat.methods.utils.graph_utils import graph_to_gdf, remove_weakly_connected_nodes, get_closest_nodes_from_gdf

logger = config.logger


def get_accessibility_isochrones(
    points: gpd.GeoDataFrame,
    weight_value: float,
    weight_type: Literal["time_min", "length_meter"],
    nx_graph: nx.Graph,
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
    nx_graph : nx.Graph
        A NetworkX graph representing the city network.
        The graph must contain the appropriate CRS and, for time-based isochrones, a speed attribute.

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        A tuple containing:
        - isochrones :
            GeoDataFrame with the calculated isochrone geometries.
        - public transport stops (if applicable) :
            GeoDataFrame with public transport stops within the isochrone, or None if not applicable.
        - public transport routes (if applicable) :
            GeoDataFrame with public transport routes within the isochrone, or None if not applicable.

    Examples
    --------
    >>> from iduedu import get_intermodal_graph
    >>> graph = get_intermodal_graph(polygon=my_territory_polygon)
    >>> points = gpd.GeoDataFrame(geometry=[Point(30.33, 59.95)], crs=4326).to_crs(graph.graph['crs'])
    >>> isochrones, stops, routes = get_accessibility_isochrones(points,15,weight_type="time_min", nx_graph=graph)

    """
    if weight_value <= 0:
        raise ValueError("Weight value must be greater than 0")
    if weight_type not in ["time_min", "length_meter"]:
        raise UserWarning("Weight type should be either 'time_min' or 'length_meter'")

    try:
        local_crs = nx_graph.graph["crs"]
    except KeyError as exc:
        raise ValueError("Graph does not have crs attribute") from exc
    try:
        points = points.to_crs(local_crs)
    except CRSError as e:
        raise CRSError(f"Graph crs ({local_crs}) has invalid format.") from e

    nx_graph = remove_weakly_connected_nodes(nx_graph)
    distances, nearest_nodes = get_closest_nodes_from_gdf(points, nx_graph)
    points["nearest_node"] = nearest_nodes

    dist_nearest = pd.DataFrame(data=distances, index=nearest_nodes, columns=["dist"]).drop_duplicates()
    speed = 0
    if nx_graph.graph["type"] in ["walk", "intermodal"] and weight_type == "time_min":
        try:
            speed = nx_graph.graph["walk_speed"]
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

    data = {}
    for source in nearest_nodes:
        dist = nx.single_source_dijkstra_path_length(nx_graph, source, weight=weight_type, cutoff=weight_value)
        data.update({source: dist})
    dist_matrix = pd.DataFrame.from_dict(data, orient="index")
    dist_matrix = dist_matrix.add(dist_nearest.dist, axis=0)
    dist_matrix = dist_matrix.mask(dist_matrix > weight_value, np.nan)
    dist_matrix.dropna(how="all", inplace=True)
    dist_matrix.dropna(how="all", axis=1, inplace=True)

    logger.info("Building isochrones geometry...")
    subgraph = nx_graph.subgraph(dist_matrix.columns.to_list())
    nodes = pd.DataFrame.from_dict(dict(subgraph.nodes(data=True)), orient="index")
    nodes["geometry"] = nodes.apply(lambda x: Point(x["x"], x["y"]), axis=1)
    nodes = gpd.GeoDataFrame(nodes, geometry="geometry", crs=local_crs)

    results = []
    for source in dist_matrix.index:
        buffers = (weight_value - dist_matrix.loc[source]) * 0.8
        if weight_type == "time_min":
            buffers = buffers * speed
        buffers = nodes.merge(buffers, left_index=True, right_index=True)
        buffers.geometry = buffers.geometry.buffer(buffers[source],resolution=8)
        results.append(buffers.union_all())

    isochrones = gpd.GeoDataFrame(geometry=results, index=dist_matrix.index, crs=local_crs)
    isochrones = (
        points.drop(columns="geometry")
        .merge(isochrones, left_on="nearest_node", right_index=True, how="left")
        .drop(columns="nearest_node")
    )
    isochrones = gpd.GeoDataFrame(isochrones, geometry="geometry", crs=local_crs)
    isochrones["weight_type"] = weight_type
    isochrones["weight_value"] = weight_value

    if "desc" in nodes.columns and "stop" in nodes["desc"].unique():
        pt_nodes = nodes[nodes["desc"] == "stop"]
        edges = graph_to_gdf(subgraph.subgraph(pt_nodes.index),nodes=False)
        pt_nodes.reset_index(drop=True, inplace=True)
        pt_nodes = pt_nodes[["desc", "route", "geometry"]]
        edges.reset_index(drop=True, inplace=True)
        edges = edges[["type", "route", "geometry"]]
        return isochrones, pt_nodes, edges
    return isochrones, None, None


def get_accessibility_isochrones_by_roads(
    points: gpd.GeoDataFrame,
    weight_value: float,
    weight_type: Literal["time_min", "length_meter"],
    nx_graph: nx.Graph,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    if weight_value <= 0:
        raise ValueError("Weight value must be greater than 0")
    if weight_type not in ["time_min", "length_meter"]:
        raise UserWarning("Weight type should be either 'time_min' or 'length_meter'")

    try:
        local_crs = nx_graph.graph["crs"]
    except KeyError as exc:
        raise ValueError("Graph does not have crs attribute") from exc
    try:
        points = points.to_crs(local_crs)
    except CRSError as e:
        raise CRSError(f"Graph crs ({local_crs}) has invalid format.") from e

    nx_graph = remove_weakly_connected_nodes(nx_graph)
    distances, nearest_nodes = get_closest_nodes_from_gdf(points, nx_graph)
    points["nearest_node"] = nearest_nodes

    dist_nearest = pd.DataFrame(data=distances, index=nearest_nodes, columns=["dist"]).drop_duplicates()
    speed = 0
    if nx_graph.graph["type"] in ["walk", "intermodal"] and weight_type == "time_min":
        try:
            speed = nx_graph.graph["walk_speed"]
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

    data = {}
    for source in nearest_nodes:
        dist = nx.single_source_dijkstra_path_length(nx_graph, source, weight=weight_type, cutoff=weight_value)
        data.update({source: dist})
    dist_matrix = pd.DataFrame.from_dict(data, orient="index")
    dist_matrix = dist_matrix.add(dist_nearest.dist, axis=0)
    dist_matrix = dist_matrix.mask(dist_matrix >= weight_value, np.nan)
    dist_matrix.dropna(how="all", inplace=True)
    dist_matrix.dropna(how="all", axis=1, inplace=True)

    return isochrones, None, None
