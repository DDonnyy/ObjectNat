from typing import Literal

import networkx as nx
import geopandas as gpd
import numpy as np
from scipy.spatial import KDTree
import pandas as pd
from shapely import Point, concave_hull


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

    if type(nx_graph) is nx.MultiDiGraph:
        nx_graph = nx.DiGraph(nx_graph)
    if type(nx_graph) is nx.MultiGraph:
        nx_graph = nx.Graph(nx_graph)
    sparse_matrix = nx.to_scipy_sparse_array(nx_graph, weight="time_min")
    transposed_matrix = sparse_matrix.transpose()
    reversed_graph = nx.from_scipy_sparse_array(
        transposed_matrix, edge_attribute=weight_type, create_using=type(nx_graph)
    )
    points = gdf_from.copy()
    points.to_crs(local_crs, inplace=True)
    points.geometry = points.representative_point()
    nodes_with_data = list(nx_graph.nodes(data=True))
    coordinates = np.array([(data["x"], data["y"]) for node, data in nodes_with_data])
    tree = KDTree(coordinates)
    target_coord = [(p.x, p.y) for p in points.representative_point()]
    distances, indices = tree.query(target_coord)
    nearest_paths = nx.multi_source_dijkstra_path(
        reversed_graph, indices.tolist(), weight=weight_type, cutoff=weight_value_cutoff
    )

    graph_points = gpd.GeoDataFrame(geometry=[Point(x, y) for x, y in coordinates], crs=local_crs)
    nearest_service_node = pd.DataFrame(
        data=[path[0] for path in nearest_paths.values()], index=nearest_paths.keys(), columns=["node_from"]
    )
    graph_nodes_gdf = gpd.GeoDataFrame(
        nearest_service_node.merge(graph_points, left_index=True, right_index=True), geometry="geometry", crs=local_crs
    )
    voronois = gpd.GeoDataFrame(geometry=graph_nodes_gdf.voronoi_polygons(), crs=local_crs)
    merged = voronois.sjoin(graph_nodes_gdf).dissolve(by="node_from").reset_index()
    if zone is None:
        zone = concave_hull(graph_nodes_gdf.union_all(), ratio=0.5)
    else:
        zone = zone.to_crs(local_crs)
    return merged.clip(zone)