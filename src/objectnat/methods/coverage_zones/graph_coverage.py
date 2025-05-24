from typing import Literal

import geopandas as gpd
import networkx as nx
import pandas as pd
from pyproj.exceptions import CRSError
from shapely import Point, concave_hull

from objectnat.methods.utils.graph_utils import get_closest_nodes_from_gdf, reverse_graph


def get_graph_coverage(
    gdf_to: gpd.GeoDataFrame,
    nx_graph: nx.Graph,
    weight_type: Literal["time_min", "length_meter"],
    weight_value_cutoff: float = None,
    zone: gpd.GeoDataFrame = None,
):
    """
    Calculate coverage zones from source points through a graph network using Dijkstra's algorithm
    and Voronoi diagrams.

    The function works by:
    1. Finding nearest graph nodes for each input point
    2. Calculating all reachable nodes within cutoff distance using Dijkstra
    3. Creating Voronoi polygons around graph nodes
    4. Combining reachability information with Voronoi cells
    5. Clipping results to specified zone boundary

    Parameters
    ----------
    gdf_to : gpd.GeoDataFrame
        Source points to which coverage is calculated.
    nx_graph : nx.Graph
        NetworkX graph representing the transportation network.
    weight_type : Literal["time_min", "length_meter"]
        Edge attribute to use as weight for path calculations.
    weight_value_cutoff : float, optional
        Maximum weight value for path calculations (e.g., max travel time/distance).
    zone : gpd.GeoDataFrame, optional
        Boundary polygon to clip the resulting coverage zones. If None, concave hull of reachable nodes will be used.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with coverage zones polygons, each associated with its source point, returns in the same CRS as
        original gdf_from.

    Notes
    -----
    - The graph must have a valid CRS attribute in its graph properties
    - MultiGraph/MultiDiGraph inputs will be converted to simple Graph/DiGraph

    Examples
    --------
    >>> from iduedu import get_intermodal_graph # pip install iduedu to get OSM city network graph
    >>> points = gpd.read_file('points.geojson')
    >>> graph = get_intermodal_graph(osm_id=1114252)
    >>> coverage = get_graph_coverage(points, graph, "time_min", 15)
    """
    original_crs = gdf_to.crs
    try:
        local_crs = nx_graph.graph["crs"]
    except KeyError as exc:
        raise ValueError("Graph does not have crs attribute") from exc

    try:
        points = gdf_to.copy()
        points.to_crs(local_crs, inplace=True)
    except CRSError as e:
        raise CRSError(f"Graph crs ({local_crs}) has invalid format.") from e

    nx_graph, reversed_graph = reverse_graph(nx_graph, weight_type)

    points.geometry = points.representative_point()

    _, nearest_nodes = get_closest_nodes_from_gdf(points, nx_graph)

    points["nearest_node"] = nearest_nodes

    nearest_paths = nx.multi_source_dijkstra_path(
        reversed_graph, nearest_nodes, weight=weight_type, cutoff=weight_value_cutoff
    )
    reachable_nodes = list(nearest_paths.keys())
    graph_points = pd.DataFrame(
        data=[{"node": node, "geometry": Point(data["x"], data["y"])} for node, data in nx_graph.nodes(data=True)]
    ).set_index("node")
    nearest_nodes = pd.DataFrame(
        data=[path[0] for path in nearest_paths.values()], index=reachable_nodes, columns=["node_to"]
    )
    graph_nodes_gdf = gpd.GeoDataFrame(
        graph_points.merge(nearest_nodes, left_index=True, right_index=True, how="left"),
        geometry="geometry",
        crs=local_crs,
    )
    graph_nodes_gdf["node_to"] = graph_nodes_gdf["node_to"].fillna("non_reachable")
    voronois = gpd.GeoDataFrame(geometry=graph_nodes_gdf.voronoi_polygons(), crs=local_crs)
    graph_nodes_gdf = graph_nodes_gdf[graph_nodes_gdf["node_to"] != "non_reachable"]
    zone_coverages = voronois.sjoin(graph_nodes_gdf).dissolve(by="node_to").reset_index().drop(columns=["node"])
    zone_coverages = zone_coverages.merge(
        points.drop(columns="geometry"), left_on="node_to", right_on="nearest_node", how="inner"
    ).reset_index(drop=True)
    zone_coverages.drop(columns=["node_to", "nearest_node"], inplace=True)
    if zone is None:
        zone = concave_hull(graph_nodes_gdf[~graph_nodes_gdf["node_to"].isna()].union_all(), ratio=0.5)
    else:
        zone = zone.to_crs(local_crs)
    return zone_coverages.clip(zone).to_crs(original_crs)
