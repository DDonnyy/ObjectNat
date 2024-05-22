import geopandas as gpd
import networkx as nx
from dongraphio import DonGraphio, GraphType
from shapely import Point


def get_accessibility_isochrones(
    graph_type: list[GraphType],
    points: Point | gpd.GeoSeries,
    weight_value: int,
    weight_type: str,
    city_graph: nx.MultiDiGraph,
    city_crs: int,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Calculate accessibility isochrones based on provided city graph for selected graph_type enum

    Args:
        graph_type (List[GraphType]): List of graph types to calculate isochrones for.
        weight_value (int): Weight value for the accessibility calculations.
        weight_type (str): The type of the weight, could be only "time_min" or "length_meter".
        city_graph (nx.MultiDiGraph): The graph representing the city.
        city_crs (int): The CRS (Coordinate Reference System) for the city.

    Returns:
        Tuple[gpd.GeoDataFrame, Optional[gpd.GeoDataFrame], Optional[gpd.GeoDataFrame]]:
        The accessibility isochrones, optional routes and stops.

    """
    dongraphio = DonGraphio(city_crs)
    dongraphio.set_graph(city_graph)
    return dongraphio.get_accessibility_isochrones(graph_type, points, weight_value, weight_type)
