import geopandas as gpd
import networkx as nx
from dongraphio import DonGraphio, GraphType
from shapely import Point


def get_accessibility_isochrones(
    points: Point | gpd.GeoSeries,
    graph_type: list[GraphType],
    weight_value: int,
    weight_type: str,
    city_graph: nx.MultiDiGraph,
    city_crs: int,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Calculate accessibility isochrones based on the provided city graph for the selected graph_type from point(s).

    Parameters
    ----------
    points : Point or gpd.GeoSeries
        Points from which the isochrones will be calculated.
    graph_type : list[GraphType]
        List of graph types to calculate isochrones for.
    weight_value : int
        Weight value for the accessibility calculations.
    weight_type : str
        The type of the weight, could be either "time_min" or "length_meter".
    city_graph : nx.MultiDiGraph
        The graph representing the city.
    city_crs : int
        The CRS (Coordinate Reference System) for the city.

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        A tuple containing the accessibility isochrones, and optionally the routes and stops.

    Examples
    --------
    >>> import networkx as nx
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point
    >>> from dongraphio import GraphType

    >>> # Create a sample city graph or download it from osm with get_intermodal_graph_from_osm()
    >>> city_graph = nx.MultiDiGraph()

    >>> # Define parameters
    >>> graph_type = [GraphType.PUBLIC_TRANSPORT, GraphType.WALK]
    >>> points = gpd.GeoSeries([Point(0, 0)])
    >>> weight_value = 15
    >>> weight_type = "time_min"
    >>> city_crs = 4326 # Should be the same with CRS of the city graph

    >>> # Calculate isochrones
    >>> isochrones, routes, stops = get_accessibility_isochrones(
    ...     graph_type, points, weight_value, weight_type, city_graph, city_crs
    ... )

    >>> print(isochrones)
    >>> print(routes)
    >>> print(stops)
    """
    dongraphio = DonGraphio(city_crs)
    dongraphio.set_graph(city_graph)
    return dongraphio.get_accessibility_isochrones(graph_type, points, weight_value, weight_type)
