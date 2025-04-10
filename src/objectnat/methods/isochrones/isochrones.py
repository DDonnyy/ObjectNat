from typing import Literal

import geopandas as gpd
import networkx as nx
import numpy as np
from shapely.ops import polygonize

from objectnat import config
from objectnat.methods.isochrones.isochrone_utils import (
    _calculate_distance_matrix,
    _create_isochrones_gdf,
    _prepare_graph_and_nodes,
    _process_pt_data,
    _validate_inputs,
)
from objectnat.methods.utils.geom_utils import polygons_to_multilinestring, remove_inner_geom
from objectnat.methods.utils.graph_utils import graph_to_gdf

logger = config.logger


def get_accessibility_isochrone_stepped(
    isochrone_type: Literal["radius", "ways", "separate"],
    point: gpd.GeoDataFrame,
    weight_value: float,
    weight_type: Literal["time_min", "length_meter"],
    nx_graph: nx.Graph,
    step: float = None,
    **kwargs,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Calculate stepped accessibility isochrones for a single point with specified intervals.

    Parameters
    ----------
    isochrone_type : Literal["radius", "ways", "separate"]
        Visualization method for stepped isochrones:
        - "radius": Voronoi-based in circular buffers
        - "ways": Voronoi-based in road network polygons
        - "separate": Circular buffers for each step
    point : gpd.GeoDataFrame
        Single source point for isochrone calculation (uses first geometry if multiple provided).
    weight_value : float
        Maximum travel time (minutes) or distance (meters) threshold.
    weight_type : Literal["time_min", "length_meter"]
        Type of weight calculation:
        - "time_min": Time-based in minutes
        - "length_meter": Distance-based in meters
    nx_graph : nx.Graph
        NetworkX graph representing the transportation network.
    step : float, optional
        Interval between isochrone steps. Defaults to:
        - 100 meters for distance-based
        - 1 minute for time-based
    **kwargs
        Additional buffer parameters:
        - buffer_factor: Size multiplier for buffers (default: 0.7)
        - road_buffer_size: Buffer size for road edges in meters (default: 5)

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        Tuple containing:
        - stepped_isochrones: GeoDataFrame with stepped polygons and distance/time attributes
        - pt_stops: Public transport stops within isochrones (if available)
        - pt_routes: Public transport routes within isochrones (if available)

    Examples
    --------
    >>> from iduedu import get_intermodal_graph # pip install iduedu to get OSM city network graph
    >>> graph = get_intermodal_graph(polygon=my_territory_polygon)
    >>> point = gpd.GeoDataFrame(geometry=[Point(30.33, 59.95)], crs=4326)
    >>> # Stepped radius isochrones with 5-minute intervals
    >>> radius_stepped, stops, _ = get_accessibility_isochrone_stepped(
    ...     "radius", point, 30, "time_min", graph, step=5
    ... )
    >>> # Stepped road isochrones with 200m intervals
    >>> ways_stepped, _, routes = get_accessibility_isochrone_stepped(
    ...     "ways", point, 1000, "length_meter", graph, step=200
    ... )
    >>> # Voronoi-based stepped isochrones
    >>> separate_stepped, stops, _ = get_accessibility_isochrone_stepped(
    ...     "separate", point, 15, "time_min", graph
    ... )
    """
    buffer_params = {
        "buffer_factor": 0.7,
        "road_buffer_size": 5,
    }

    buffer_params.update(kwargs)
    original_crs = point.crs
    point = point.copy()
    if len(point) > 1:
        logger.warning(
            f"This method processes only single point. The GeoDataFrame contains {len(point)} points - "
            "only the first geometry will be used for isochrone calculation. "
        )
        point = point.iloc[[0]]

    local_crs, graph_type = _validate_inputs(point, weight_value, weight_type, nx_graph)

    if step is None:
        if weight_type == "length_meter":
            step = 100
        else:
            step = 1
    nx_graph, points, dist_nearest, speed = _prepare_graph_and_nodes(
        point, nx_graph, graph_type, weight_type, weight_value
    )

    dist_matrix, subgraph = _calculate_distance_matrix(
        nx_graph, points["nearest_node"].values, weight_type, weight_value, dist_nearest
    )

    logger.info("Building isochrones geometry...")
    nodes, edges = graph_to_gdf(subgraph)
    nodes.loc[dist_matrix.columns, "dist"] = dist_matrix.iloc[0]
    steps = np.arange(0, weight_value + step, step)
    if steps[-1] > weight_value:
        steps[-1] = weight_value  # Ensure last step doesn't exceed weight_value

    if isochrone_type == "separate":
        for i in range(len(steps) - 1):
            min_dist = steps[i]
            max_dist = steps[i + 1]
            nodes_in_step = nodes["dist"].between(min_dist, max_dist, inclusive="left")
            nodes_in_step = nodes_in_step[nodes_in_step].index
            if not nodes_in_step.empty:
                buffer_size = (max_dist - nodes.loc[nodes_in_step, "dist"]) * 0.7
                if weight_type == "time_min":
                    buffer_size = buffer_size * speed
                nodes.loc[nodes_in_step, "buffer_size"] = buffer_size
        nodes.geometry = nodes.geometry.buffer(nodes["buffer_size"])
        nodes["dist"] = np.round(nodes["dist"], 0)
        nodes = nodes.dissolve(by="dist", as_index=False)
        polygons = gpd.GeoDataFrame(
            geometry=list(polygonize(nodes.geometry.apply(polygons_to_multilinestring).union_all())),
            crs=local_crs,
        )
        polygons_points = polygons.copy()
        polygons_points.geometry = polygons.representative_point()

        stepped_iso = polygons_points.sjoin(nodes, predicate="within").reset_index()
        stepped_iso = stepped_iso.groupby("index").agg({"dist": "mean"})
        stepped_iso["geometry"] = polygons
        stepped_iso = gpd.GeoDataFrame(stepped_iso, geometry="geometry", crs=local_crs).reset_index(drop=True)
    else:
        if isochrone_type == "radius":
            isochrone_geoms = _build_radius_isochrones(
                dist_matrix, weight_value, weight_type, speed, nodes, buffer_params["buffer_factor"]
            )
        else:  # isochrone_type == 'ways':
            if graph_type in ["intermodal", "walk"]:
                isochrone_edges = edges[edges["type"] == "walk"]
            else:
                isochrone_edges = edges.copy()
            all_isochrones_edges = isochrone_edges.buffer(buffer_params["road_buffer_size"], resolution=1).union_all()
            all_isochrones_edges = gpd.GeoDataFrame(geometry=[all_isochrones_edges], crs=local_crs)
            isochrone_geoms = _build_ways_isochrones(
                dist_matrix=dist_matrix,
                weight_value=weight_value,
                weight_type=weight_type,
                speed=speed,
                nodes=nodes,
                all_isochrones_edges=all_isochrones_edges,
                buffer_factor=buffer_params["buffer_factor"],
            )
        nodes = nodes.clip(isochrone_geoms[0], keep_geom_type=True)
        nodes["dist"] = np.minimum(np.ceil(nodes["dist"] / step) * step, weight_value)
        voronois = gpd.GeoDataFrame(geometry=nodes.voronoi_polygons(), crs=local_crs)
        stepped_iso = (
            voronois.sjoin(nodes[["dist", "geometry"]]).dissolve(by="dist", as_index=False).drop(columns="index_right")
        )
        stepped_iso = stepped_iso.clip(isochrone_geoms[0], keep_geom_type=True)

    pt_nodes, pt_edges = _process_pt_data(nodes, edges, graph_type)
    if pt_nodes is not None:
        pt_nodes.to_crs(original_crs, inplace=True)
    if pt_edges is not None:
        pt_edges.to_crs(original_crs, inplace=True)
    return stepped_iso.to_crs(original_crs), pt_nodes, pt_edges


def get_accessibility_isochrones(
    isochrone_type: Literal["radius", "ways"],
    points: gpd.GeoDataFrame,
    weight_value: float,
    weight_type: Literal["time_min", "length_meter"],
    nx_graph: nx.Graph,
    **kwargs,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Calculate accessibility isochrones from input points based on the provided city graph.

    Supports two types of isochrones:
    - 'radius': Circular buffer-based isochrones
    - 'ways': Road network-based isochrones

    Parameters
    ----------
    isochrone_type : Literal["radius", "ways"]
        Type of isochrone to calculate:
        - "radius": Creates circular buffers around reachable nodes
        - "ways": Creates polygons based on reachable road network
    points : gpd.GeoDataFrame
        GeoDataFrame containing source points for isochrone calculation.
    weight_value : float
        Maximum travel time (minutes) or distance (meters) threshold.
    weight_type : Literal["time_min", "length_meter"]
        Type of weight calculation:
        - "time_min": Time-based accessibility in minutes
        - "length_meter": Distance-based accessibility in meters
    nx_graph : nx.Graph
        NetworkX graph representing the transportation network.
        Must contain CRS and speed attributes for time calculations.
    **kwargs
        Additional buffer parameters:
        - buffer_factor: Size multiplier for buffers (default: 0.7)
        - road_buffer_size: Buffer size for road edges in meters (default: 5)

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        Tuple containing:
        - isochrones: GeoDataFrame with calculated isochrone polygons
        - pt_stops: Public transport stops within isochrones (if available)
        - pt_routes: Public transport routes within isochrones (if available)

    Examples
    --------
    >>> from iduedu import get_intermodal_graph # pip install iduedu to get OSM city network graph
    >>> graph = get_intermodal_graph(polygon=my_territory_polygon)
    >>> points = gpd.GeoDataFrame(geometry=[Point(30.33, 59.95)], crs=4326)
    >>> # Radius isochrones
    >>> radius_iso, stops, routes = get_accessibility_isochrones(
    ...     "radius", points, 15, "time_min", graph, buffer_factor=0.8
    ... )
    >>> # Road network isochrones
    >>> ways_iso, stops, routes = get_accessibility_isochrones(
    ...     "ways", points, 1000, "length_meter", graph, road_buffer_size=7
    ... )
    """

    buffer_params = {
        "buffer_factor": 0.7,
        "road_buffer_size": 5,
    }
    original_crs = points.crs
    buffer_params.update(kwargs)

    points = points.copy()
    local_crs, graph_type = _validate_inputs(points, weight_value, weight_type, nx_graph)

    nx_graph, points, dist_nearest, speed = _prepare_graph_and_nodes(
        points, nx_graph, graph_type, weight_type, weight_value
    )

    weight_cutoff = (
        weight_value + (100 if weight_type == "length_meter" else 1) if isochrone_type == "ways" else weight_value
    )

    dist_matrix, subgraph = _calculate_distance_matrix(
        nx_graph, points["nearest_node"].values, weight_type, weight_cutoff, dist_nearest
    )

    logger.info("Building isochrones geometry...")
    nodes, edges = graph_to_gdf(subgraph)
    if isochrone_type == "radius":
        isochrone_geoms = _build_radius_isochrones(
            dist_matrix, weight_value, weight_type, speed, nodes, buffer_params["buffer_factor"]
        )
    else:  # isochrone_type == 'ways':
        if graph_type in ["intermodal", "walk"]:
            isochrone_edges = edges[edges["type"] == "walk"]
        else:
            isochrone_edges = edges.copy()
        all_isochrones_edges = isochrone_edges.buffer(buffer_params["road_buffer_size"], resolution=1).union_all()
        all_isochrones_edges = gpd.GeoDataFrame(geometry=[all_isochrones_edges], crs=local_crs)
        isochrone_geoms = _build_ways_isochrones(
            dist_matrix=dist_matrix,
            weight_value=weight_value,
            weight_type=weight_type,
            speed=speed,
            nodes=nodes,
            all_isochrones_edges=all_isochrones_edges,
            buffer_factor=buffer_params["buffer_factor"],
        )
    isochrones = _create_isochrones_gdf(points, isochrone_geoms, dist_matrix, local_crs, weight_type, weight_value)
    pt_nodes, pt_edges = _process_pt_data(nodes, edges, graph_type)
    if pt_nodes is not None:
        pt_nodes.to_crs(original_crs, inplace=True)
    if pt_edges is not None:
        pt_edges.to_crs(original_crs, inplace=True)
    return isochrones.to_crs(original_crs), pt_nodes, pt_edges


def _build_radius_isochrones(dist_matrix, weight_value, weight_type, speed, nodes, buffer_factor):
    results = []
    for source in dist_matrix.index:
        buffers = (weight_value - dist_matrix.loc[source]) * buffer_factor
        if weight_type == "time_min":
            buffers = buffers * speed
        buffers = nodes.merge(buffers, left_index=True, right_index=True)
        buffers.geometry = buffers.geometry.buffer(buffers[source], resolution=8)
        results.append(buffers.union_all())
    return results


def _build_ways_isochrones(dist_matrix, weight_value, weight_type, speed, nodes, all_isochrones_edges, buffer_factor):
    results = []
    for source in dist_matrix.index:
        reachable_nodes = dist_matrix.loc[source]
        reachable_nodes = reachable_nodes[reachable_nodes <= weight_value]
        reachable_nodes = (weight_value - reachable_nodes) * buffer_factor
        if weight_type == "time_min":
            reachable_nodes = reachable_nodes * speed
        reachable_nodes = nodes.merge(reachable_nodes, left_index=True, right_index=True)
        clip_zone = reachable_nodes.buffer(reachable_nodes[source], resolution=4).union_all()

        isochrone_edges = all_isochrones_edges.clip(clip_zone, keep_geom_type=True).explode(ignore_index=True)
        geom_to_keep = isochrone_edges.sjoin(reachable_nodes, how="inner").index.unique()
        isochrone = remove_inner_geom(isochrone_edges.loc[geom_to_keep].union_all())
        results.append(isochrone)
    return results
