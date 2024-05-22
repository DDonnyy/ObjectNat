import geopandas as gpd
from dongraphio import GraphType


def get_radius_zone_coverage(services: gpd.GeoDataFrame, radius: int) -> gpd.GeoDataFrame:
    """
    Creates a buffer with a defined radius for each service.

    Parameters
    ---------
    service_type: str
        The type of services to run the method on.
    radius: int, optional
        The radius for the buffer.

    Example
    -------
    Get coverage zones for schools with radius of 50 meters.
        CityMetricsMethods.Coverage_Zones(city_model).get_radius_zone(service_type='schools', radius=50)
    """

    services["geometry"] = services["geometry"].buffer(radius)

    return services


def get_isochrone_zone_coverage(services: gpd.GeoDataFrame, weight_type: str, weight_value: int, graph_type: list[GraphType]):
    """
    Creates an isochrone with defined way of transportation and time to travel (in minutes) for each service.
    The method calls Accessibility_Isochrones_v2.get_isochrone.

    Parameters
    ---------
    service_type: str
        The type of services to run the method on.
    travel_type: str
        From Accessibility_Isochrones_v2.
        One of the given ways of transportation: "public_transport", "walk" or "drive".
    weight_value: int
        From Accessibility_Isochrones_v2.
        Minutes to travel.

    Returns
    -------
    FeatureCollection

    Example
    -------
    Get coverage zones for dentistries using 10 mins pedestrian-ways isochrone.
        CityMetricsMethods.Coverage_Zones(city_model)._get_isochrone_zone(
            service_type='dentists', travel_type='walk', weight_value = 10)
    """
    points = services.geometry.representative_point()
    isochrone_res =

    return isochrone_res
