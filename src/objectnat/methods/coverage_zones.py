from typing import Literal

import geopandas as gpd
import networkx as nx
import pandas as pd
from dongraphio import GraphType

from .isochrones import get_accessibility_isochrones


def get_radius_zone_coverage(services: gpd.GeoDataFrame, radius: int) -> gpd.GeoDataFrame:
    """
    Create a buffer zone with a defined radius around each service location.

    Parameters
    ----------
    services : gpd.GeoDataFrame
        GeoDataFrame containing the service locations.
    radius : int
        The radius for the buffer in meters.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with the buffer zones around each service location.

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point

    >>> # Create a sample GeoDataFrame for services
    >>> services = gpd.read_file('services.geojson')

    >>> # Define the radius
    >>> radius = 50

    >>> # Get radius zone coverage
    >>> radius_zones = get_radius_zone_coverage(services, radius)
    >>> print(radius_zones)
    """
    services["geometry"] = services["geometry"].buffer(radius)
    return services


def get_isochrone_zone_coverage(
    services: gpd.GeoDataFrame,
    weight_type: Literal["time_min", "length_meter"],
    weight_value: int,
    city_graph: nx.Graph,
    graph_type: list[GraphType],
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Create isochrones for each service location based on travel time/distance.

    Parameters
    ----------
    services : gpd.GeoDataFrame
        GeoDataFrame containing the service locations.
    weight_type : str
        Type of weight used for calculating isochrones, either "time_min" or "length_meter".
    weight_value : int
        The value of the weight, representing time in minutes or distance in meters.
    city_graph : nx.Graph
        The graph representing the city's transportation network.
    graph_type : list[GraphType]
        List of graph types to be used for isochrone calculations.

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        The calculated isochrone zones, optionally including routes and stops.

    Examples
    --------
    >>> import networkx as nx
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point
    >>> from dongraphio import GraphType

    >>> # Create a sample city graph or download it from osm with get_intermodal_graph_from_osm()
    >>> city_graph = nx.MultiDiGraph()

    >>> # Create a sample GeoDataFrame for services
    >>> services = gpd.read_file('services.geojson')

    >>> # Define parameters
    >>> weight_type = "time_min"
    >>> weight_value = 10
    >>> graph_type = [GraphType.PUBLIC_TRANSPORT, GraphType.WALK]

    >>> # Get isochrone zone coverage
    >>> isochrone_zones = get_isochrone_zone_coverage(services, weight_type, weight_value, city_graph, graph_type)
    >>> isochrone_zones[0] # represent isochrones geodataframe
    """
    assert services.crs == city_graph.graph["crs"], "CRS not match"
    points = services.geometry.representative_point()
    iso, routes, stops = get_accessibility_isochrones(
        points, graph_type, weight_value, weight_type, city_graph, points.crs.to_epsg()
    )
    services_ = services.copy()
    iso = gpd.GeoDataFrame(
        pd.concat([iso.drop(columns=["point", "point_number"]), services_.drop(columns=["geometry"])], axis=1)
    )
    return iso, routes, stops
