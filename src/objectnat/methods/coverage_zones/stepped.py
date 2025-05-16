from typing import Literal

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from pyproj.exceptions import CRSError
from shapely import Point, concave_hull

from objectnat.methods.utils.graph_utils import get_closest_nodes_from_gdf, remove_weakly_connected_nodes


def get_stepped_graph_coverage(
    gdf_to: gpd.GeoDataFrame,
    nx_graph: nx.Graph,
    weight_type: Literal["time_min", "length_meter"],
    step_type: Literal["voronoi", "separate"],
    weight_value_cutoff: float = None,
    zone: gpd.GeoDataFrame = None,
    step: float = None,

):

    if step is None:
        if weight_type == "length_meter":
            step = 100
        else:
            step = 1

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

    if nx_graph.is_multigraph():
        if nx_graph.is_directed():
            nx_graph = nx.DiGraph(nx_graph)
        else:
            nx_graph = nx.Graph(nx_graph)

    nx_graph = remove_weakly_connected_nodes(nx_graph)

    mapping = {old_label: new_label for new_label, old_label in enumerate(nx_graph.nodes())}
    nx_graph = nx.relabel_nodes(nx_graph, mapping)

    sparse_matrix = nx.to_scipy_sparse_array(nx_graph, weight=weight_type)
    transposed_matrix = sparse_matrix.transpose()
    reversed_graph = nx.from_scipy_sparse_array(
        transposed_matrix, edge_attribute=weight_type, create_using=type(nx_graph)
    )

    points.geometry = points.representative_point()

    distances, nearest_nodes = get_closest_nodes_from_gdf(points, nx_graph)

    points["nearest_node"] = nearest_nodes
    points["distance"] = distances

    dist = nx.multi_source_dijkstra_path_length(
        reversed_graph, nearest_nodes, weight=weight_type, cutoff=weight_value_cutoff
    )

    graph_points = pd.DataFrame(
        data=[{"node": node, "geometry": Point(data["x"], data["y"])} for node, data in nx_graph.nodes(data=True)]
    )

    nearest_nodes = pd.DataFrame.from_dict(dist, orient="index", columns=["dist"]).reset_index()

    graph_nodes_gdf = gpd.GeoDataFrame(
        graph_points.merge(nearest_nodes, left_on="node", right_on="index", how="left").reset_index(drop=True),
        geometry="geometry",
        crs=local_crs,
    )
    graph_nodes_gdf["dist"] = np.minimum(np.ceil(graph_nodes_gdf["dist"] / step) * step, weight_value_cutoff)
    voronois = gpd.GeoDataFrame(geometry=graph_nodes_gdf.voronoi_polygons(), crs=local_crs)

    zone_coverages = voronois.sjoin(graph_nodes_gdf).dissolve(by="dist", as_index=False,dropna=False)

    if zone is None:
        zone = concave_hull(graph_nodes_gdf[~graph_nodes_gdf["node_to"].isna()].union_all(), ratio=0.5)
    else:
        zone = zone.to_crs(local_crs)
    return zone_coverages.clip(zone).to_crs(original_crs)
