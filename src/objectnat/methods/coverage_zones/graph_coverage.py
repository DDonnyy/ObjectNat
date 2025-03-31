from typing import Literal

import geopandas as gpd
import networkx as nx
import pandas as pd
from pyproj.exceptions import CRSError
from shapely import Point, concave_hull

from objectnat.methods.utils.graph_utils import get_closest_nodes_from_gdf


def get_graph_coverage(
    gdf_from: gpd.GeoDataFrame,
    nx_graph: nx.Graph,
    weight_type: Literal["time_min", "length_meter"],
    weight_value_cutoff: float = None,
    zone: gpd.GeoDataFrame = None,
):
    try:
        local_crs = nx_graph.graph["crs"]
    except KeyError as exc:
        raise ValueError("Graph does not have crs attribute") from exc

    try:
        points = gdf_from.copy()
        points.to_crs(local_crs, inplace=True)
    except CRSError as e:
        raise CRSError(f"Graph crs ({local_crs}) has invalid format.") from e

    if type(nx_graph) is nx.MultiDiGraph:
        nx_graph = nx.DiGraph(nx_graph)
    if type(nx_graph) is nx.MultiGraph:
        nx_graph = nx.Graph(nx_graph)

    sparse_matrix = nx.to_scipy_sparse_array(nx_graph, weight=weight_type)
    transposed_matrix = sparse_matrix.transpose()
    reversed_graph = nx.from_scipy_sparse_array(
        transposed_matrix, edge_attribute=weight_type, create_using=type(nx_graph)
    )

    points.geometry = points.representative_point()

    distances, nearest_service_nodes = get_closest_nodes_from_gdf(points, nx_graph)

    nearest_paths = nx.multi_source_dijkstra_path(
        reversed_graph, nearest_service_nodes, weight=weight_type, cutoff=weight_value_cutoff
    )
    graph_points = gpd.GeoDataFrame(
        geometry=[
            Point(data["x"], data["y"])
            for node, data in list(nx_graph.nodes(data=True))
            if node in nearest_service_nodes
        ],
        crs=local_crs,
    )
    nearest_service_nodes = pd.DataFrame(
        data=[path[0] for path in nearest_paths.values()], index=nearest_paths.keys(), columns=["node_from"]
    )
    graph_nodes_gdf = gpd.GeoDataFrame(
        nearest_service_nodes.merge(graph_points, left_index=True, right_index=True), geometry="geometry", crs=local_crs
    )
    voronois = gpd.GeoDataFrame(geometry=graph_nodes_gdf.voronoi_polygons(), crs=local_crs)
    merged = voronois.sjoin(graph_nodes_gdf).dissolve(by="node_from").reset_index()
    if zone is None:
        zone = concave_hull(graph_nodes_gdf.union_all(), ratio=0.5)
    else:
        zone = zone.to_crs(local_crs)
    return merged.clip(zone)
