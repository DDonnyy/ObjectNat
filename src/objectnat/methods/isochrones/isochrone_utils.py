from typing import Literal

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from pyproj.exceptions import CRSError

from objectnat import config
from objectnat.methods.utils.graph_utils import get_closest_nodes_from_gdf, remove_weakly_connected_nodes

logger = config.logger


def _validate_inputs(
    points: gpd.GeoDataFrame, weight_value: float, weight_type: Literal["time_min", "length_meter"], nx_graph: nx.Graph
) -> tuple[str, str]:
    """Validate common inputs for accessibility functions."""
    if weight_value <= 0:
        raise ValueError("Weight value must be greater than 0")
    if weight_type not in ["time_min", "length_meter"]:
        raise UserWarning("Weight type should be either 'time_min' or 'length_meter'")

    try:
        local_crs = nx_graph.graph["crs"]
    except KeyError as exc:
        raise ValueError("Graph does not have crs attribute") from exc
    try:
        graph_type = nx_graph.graph["type"]
    except KeyError as exc:
        raise ValueError("Graph does not have type attribute") from exc

    try:
        points.to_crs(local_crs, inplace=True)
    except CRSError as e:
        raise CRSError(f"Graph crs ({local_crs}) has invalid format.") from e

    return local_crs, graph_type


def _prepare_graph_and_nodes(
    points: gpd.GeoDataFrame, nx_graph: nx.Graph, graph_type: str, weight_type: str, weight_value: float
) -> tuple[nx.Graph, gpd.GeoDataFrame, pd.DataFrame, float]:
    """Prepare graph and calculate nearest nodes with distances."""
    nx_graph = remove_weakly_connected_nodes(nx_graph)
    distances, nearest_nodes = get_closest_nodes_from_gdf(points, nx_graph)
    points["nearest_node"] = nearest_nodes

    dist_nearest = pd.DataFrame(data=distances, index=nearest_nodes, columns=["dist"]).drop_duplicates()

    # Calculate speed adjustment if needed
    speed = 0
    if graph_type in ["walk", "intermodal"] and weight_type == "time_min":
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

    return nx_graph, points, dist_nearest, speed


def _process_pt_data(
    nodes: gpd.GeoDataFrame, edges: gpd.GeoDataFrame, graph_type: str
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame] | tuple[None, None]:
    """Process public transport data if available."""
    if "desc" in nodes.columns and "stop" in nodes["desc"].unique():
        pt_nodes = nodes[nodes["desc"] == "stop"]
        if graph_type == "intermodal":
            edges = edges[~edges["type"].isin(["walk", "boarding"])]
        pt_nodes = pt_nodes[["desc", "route", "geometry"]]
        edges = edges[["type", "route", "geometry"]]
        return pt_nodes, edges
    return None, None


def _calculate_distance_matrix(
    nx_graph: nx.Graph,
    nearest_nodes: np.ndarray,
    weight_type: str,
    weight_value: float,
    dist_nearest: pd.DataFrame,
) -> tuple[pd.DataFrame, nx.Graph]:
    """Calculate distance matrix from nearest nodes."""

    data = {}
    for source in nearest_nodes:
        dist = nx.single_source_dijkstra_path_length(nx_graph, source, weight=weight_type, cutoff=weight_value)
        data.update({source: dist})

    dist_matrix = pd.DataFrame.from_dict(data, orient="index")
    dist_matrix = dist_matrix.add(dist_nearest.dist, axis=0)
    dist_matrix = dist_matrix.mask(dist_matrix > weight_value, np.nan)
    dist_matrix.dropna(how="all", inplace=True)
    dist_matrix.dropna(how="all", axis=1, inplace=True)

    subgraph = nx_graph.subgraph(dist_matrix.columns.to_list())

    return dist_matrix, subgraph


def _create_isochrones_gdf(
    points: gpd.GeoDataFrame,
    results: list,
    dist_matrix: pd.DataFrame,
    local_crs: str,
    weight_type: str,
    weight_value: float,
) -> gpd.GeoDataFrame:
    """Create final isochrones GeoDataFrame."""
    isochrones = gpd.GeoDataFrame(geometry=results, index=dist_matrix.index, crs=local_crs)
    isochrones = (
        points.drop(columns="geometry")
        .merge(isochrones, left_on="nearest_node", right_index=True, how="left")
        .drop(columns="nearest_node")
    )
    isochrones = gpd.GeoDataFrame(isochrones, geometry="geometry", crs=local_crs)
    isochrones["weight_type"] = weight_type
    isochrones["weight_value"] = weight_value
    return isochrones
