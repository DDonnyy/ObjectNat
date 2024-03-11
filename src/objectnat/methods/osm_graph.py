import networkx as nx
from dongraphio import DonGraphio


def get_intermodal_graph_from_osm(
    city_osm_id: int, keep_city_boundary: bool = True, city_crs: int = 3857, dongraphio: DonGraphio = None
) -> nx.MultiDiGraph:
    """
    Generate an intermodal graph from OpenStreetMap data for the specified city.

    Args:
        city_osm_id (int): The OpenStreetMap ID of the city.
        keep_city_boundary (bool, optional): Flag to indicate whether to keep the city boundary. Defaults to True.
        city_crs (int, optional): The Coordinate Reference System (CRS) for the city. Defaults to 3857.
        dongraphio (DonGraphio, optional): An instance of DonGraphio for handling the graph. Defaults to None.

    Returns:
        nx.MultiDiGraph: The intermodal graph generated from OpenStreetMap data.
    """
    if dongraphio:
        return dongraphio.get_intermodal_graph_from_osm(city_osm_id, keep_city_boundary)
    dongraphio = DonGraphio(city_crs)
    return dongraphio.get_intermodal_graph_from_osm(city_osm_id, keep_city_boundary)
