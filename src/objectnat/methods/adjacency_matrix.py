import geopandas as gpd
import networkx as nx
import pandas as pd
from dongraphio import DonGraphio
from loguru import logger
from pydantic import ValidationError


def get_adjacency_matrix(
    buildings_from: gpd.GeoDataFrame,
    services_to: gpd.GeoDataFrame,
    weight: str,
    city_crs: int | None = None,
    nx_graph: nx.MultiDiGraph | None = None,
    dongraphio: DonGraphio | None = None,
    graph_type=None,
) -> pd.DataFrame:
    """
    Get the adjacency matrix for the specified city graph, buildings, and services.

    Args:
        nx_graph (nx.Graph): The networkx graph.
        buildings_from (gpd.GeoDataFrame): GeoDataFrame representing buildings to build matrix from.
        services_to (gpd.GeoDataFrame): GeoDataFrame representing services to build matrix to.
        weight (str): The weight attribute, could be only "time_min" or "length_meter".

    Returns:
        pd.DataFrame: The adjacency matrix.
    """
    try:
        if dongraphio:
            return dongraphio.get_adjacency_matrix(buildings_from, services_to, weight, graph_type=graph_type)

        dongraphio = DonGraphio(city_crs)
        dongraphio.set_graph(nx_graph)
        return dongraphio.get_adjacency_matrix(buildings_from, services_to, weight, graph_type=graph_type)
    except ValidationError as e:
        logger.error("Function get_adjacency_matrix() missing 'weight' argument")
        raise e
