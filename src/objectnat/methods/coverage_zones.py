from typing import Literal

import geopandas as gpd
import networkx as nx

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
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Create isochrones for each service location based on travel time/distance.

    Parameters
    ----------
    services : gpd.GeoDataFrame
        Layer containing the service locations.
    weight_type : str
        Type of weight used for calculating isochrones, either "time_min" or "length_meter".
    weight_value : int
        The value of the weight, representing time in minutes or distance in meters.
    city_graph : nx.Graph
        The graph representing the city's transportation network.

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        The calculated isochrone zones, optionally including routes and stops.

    Examples
    --------
    >>> import networkx as nx
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point
     >>> from iduedu import get_intermodal_graph

    >>> # Create a sample city graph with get_intermodal_graph()
    >>> graph = get_intermodal_graph(polygon=my_territory_polygon)

    >>> # Create a sample GeoDataFrame for services
    >>> services = gpd.read_file('services.geojson')

    >>> # Define parameters
    >>> weight_type = "time_min"
    >>> weight_value = 10

    >>> # Get isochrone zone coverage
    >>> isochrones, pt_stops, pt_routes = get_isochrone_zone_coverage(services, weight_type, weight_value, city_graph)
    """
    iso, routes, stops = get_accessibility_isochrones(services, weight_value, weight_type, city_graph)
    return iso, routes, stops
