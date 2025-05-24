from typing import Literal

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from pyproj.exceptions import CRSError
from shapely import Point, concave_hull

from objectnat.methods.isochrones.isochrone_utils import create_separated_dist_polygons
from objectnat.methods.utils.graph_utils import get_closest_nodes_from_gdf, reverse_graph


def get_stepped_graph_coverage(
    gdf_to: gpd.GeoDataFrame,
    nx_graph: nx.Graph,
    weight_type: Literal["time_min", "length_meter"],
    step_type: Literal["voronoi", "separate"],
    weight_value_cutoff: float = None,
    zone: gpd.GeoDataFrame = None,
    step: float = None,
):
    """
    Calculate stepped coverage zones from source points through a graph network using Dijkstra's algorithm
    and Voronoi-based or buffer-based isochrone steps.

    This function combines graph-based accessibility with stepped isochrone logic. It:
    1. Finds nearest graph nodes for each input point
    2. Computes reachability for increasing weights (e.g. time or distance) in defined steps
    3. Generates Voronoi-based or separate buffer zones around network nodes
    4. Aggregates zones into stepped coverage layers
    5. Optionally clips results to a boundary zone

    Parameters
    ----------
    gdf_to : gpd.GeoDataFrame
        Source points from which stepped coverage is calculated.
    nx_graph : nx.Graph
        NetworkX graph representing the transportation network.
    weight_type : Literal["time_min", "length_meter"]
        Type of edge weight to use for path calculation:
        - "time_min": Edge travel time in minutes
        - "length_meter": Edge length in meters
    step_type : Literal["voronoi", "separate"]
        Method for generating stepped zones:
        - "voronoi": Stepped zones based on Voronoi polygons around graph nodes
        - "separate": Independent buffer zones per step
    weight_value_cutoff : float, optional
        Maximum weight value (e.g., max travel time or distance) to limit the coverage extent.
    zone : gpd.GeoDataFrame, optional
        Optional boundary polygon to clip resulting stepped zones. If None, concave hull of reachable area is used.
    step : float, optional
        Step interval for coverage zone construction. Defaults to:
        - 100 meters for distance-based weight
        - 1 minute for time-based weight

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with polygons representing stepped coverage zones for each input point, annotated by step range.

    Notes
    -----
    - Input graph must have a valid CRS defined.
    - MultiGraph or MultiDiGraph inputs will be simplified.
    - Designed for accessibility and spatial equity analyses over multimodal networks.

    Examples
    --------
    >>> from iduedu import get_intermodal_graph
    >>> points = gpd.read_file('destinations.geojson')
    >>> graph = get_intermodal_graph(osm_id=1114252)
    >>> stepped_coverage = get_stepped_graph_coverage(
    ...     points, graph, "time_min", step_type="voronoi", weight_value_cutoff=30, step=5
    ... )
    >>> # Using buffer-style zones
    >>> stepped_separate = get_stepped_graph_coverage(
    ...     points, graph, "length_meter", step_type="separate", weight_value_cutoff=1000, step=200
    ... )
    """
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

    nx_graph, reversed_graph = reverse_graph(nx_graph, weight_type)

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
    graph_nodes_gdf.drop(columns=["index", "node"], inplace=True)
    if weight_value_cutoff is None:
        weight_value_cutoff = max(nearest_nodes["dist"])
    if step_type == "voronoi":
        graph_nodes_gdf["dist"] = np.minimum(np.ceil(graph_nodes_gdf["dist"] / step) * step, weight_value_cutoff)
        voronois = gpd.GeoDataFrame(geometry=graph_nodes_gdf.voronoi_polygons(), crs=local_crs)
        zone_coverages = voronois.sjoin(graph_nodes_gdf).dissolve(by="dist", as_index=False, dropna=False)
        zone_coverages = zone_coverages[["dist", "geometry"]].explode(ignore_index=True)
        if zone is None:
            zone = concave_hull(graph_nodes_gdf[~graph_nodes_gdf["node_to"].isna()].union_all(), ratio=0.5)
        else:
            zone = zone.to_crs(local_crs)
        zone_coverages = zone_coverages.clip(zone).to_crs(original_crs)
    else:  # step_type == 'separate':
        speed = 83.33  # TODO HARDCODED WALK SPEED
        weight_value = weight_value_cutoff
        zone_coverages = create_separated_dist_polygons(graph_nodes_gdf, weight_value, weight_type, step, speed)
        if zone is not None:
            zone = zone.to_crs(local_crs)
            zone_coverages = zone_coverages.clip(zone).to_crs(original_crs)
    return zone_coverages
