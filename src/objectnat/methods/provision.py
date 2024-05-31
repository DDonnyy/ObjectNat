import geopandas as gpd
import networkx as nx
import pandas as pd
from dongraphio import DonGraphio
from loguru import logger
from provisio import get_service_provision
from provisio.provisio import CityProvision
from provisio.provisio_exceptions import CapacityKeyError, DemandKeyError

from .demands import get_demands


class NoWeightAdjacencyException(RuntimeError):
    pass


class NoOsmIdException(RuntimeError):
    pass


def get_provision(
    buildings: gpd.GeoDataFrame,
    services: gpd.GeoDataFrame,
    threshold: int,
    adjacency_matrix: pd.DataFrame | None = None,
    calculation_type: str = "gravity",
    demand_normative: float | None = None,
    population: int | None = None,
    city_graph: nx.MultiDiGraph | None = None,
    city_crs: int | None = 3857,
    city_osm_id: int | None = None,
    weight_adjacency_matrix: str | None = None,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Calculate the provision based on the specified buildings demands and services capacity.

    Args:
        buildings (gpd.GeoDataFrame): GeoDataFrame representing buildings.
            - Must contain either "living_area" or "storeys_count" columns if buildings are not evacuated;
            - Must contain "population" column if buildings are evacuated;
            - Must contain "demand" column if demands are calculated.

        services (gpd.GeoDataFrame): GeoDataFrame representing services.
            - Must contain "capacity" attribute indicating the capacity of each service.

        threshold (int): The threshold value for calculating provision.

        adjacency_matrix (pd.DataFrame): The adjacency matrix for the specified graph, buildings, and services.
            If None, a city_graph or city_osm_id is needed to calculate the matrix.

        calculation_type (str): The type of calculation to use for provision, can only be "linear" or "gravity".

        demand_normative (float): The normative value for calculating demands.
            Required if demands are not present in buildings.

        population (int): Total population of provided buildings data for resettlement if buildings are not evacuated.

        city_graph (nx.MultiDiGraph): The graph representing the city.
            Needed if adjacency_matrix is not provided.

        city_crs (int): The CRS (Coordinate Reference System) for the city. Default is 3857.

        city_osm_id (int): The OpenStreetMap ID of the city.
            Needed if adjacency_matrix and city_graph are not provided.

        weight_adjacency_matrix (str): The weight attribute for adjacency matrix calculation.
            Can only be "time_min" or "length_meter".

    Returns:
        Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
        The calculated provision for buildings, services and links.
    """
    buildings.to_crs(city_crs, inplace=True)
    services.to_crs(city_crs, inplace=True)

    calculate_matrix = False
    calculate_graph = False
    calculate_demands = False
    dngraph = DonGraphio(city_crs=city_crs)

    if adjacency_matrix is None:
        logger.warning("The adjacency matrix is not provided.")
        if not weight_adjacency_matrix:
            raise NoWeightAdjacencyException("No weight type ('time_min' or 'length_meter') was provided.")
        if not city_graph:
            logger.warning("The graph is not provided, attempting to load data from OSM.")
            if not city_osm_id:
                raise NoOsmIdException("osm_id of the city is not provided, unable to retrieve data.")
            calculate_graph = True
        calculate_matrix = True
        adjacency_matrix = pd.DataFrame(data=0, index=[], columns=buildings.index.astype(int).tolist())

    try:
        cp = CityProvision(
            services=services,
            demanded_buildings=buildings,
            adjacency_matrix=adjacency_matrix,
            threshold=threshold,
            calculation_type=calculation_type,
        )
        buildings = cp.demanded_buildings
        services = cp.services
    except DemandKeyError as demand_error:
        logger.warning("The 'demand' column is missing in the provided building data, attempting to calculate values.")
        if demand_normative is None:
            raise ValueError(
                "Unable to calculate demand and provision accordingly, 'demand_normative' value is not specified."
            ) from demand_error
        calculate_demands = True
    except CapacityKeyError as capacity_error:
        raise ValueError(
            "The 'capacity; column is missing in provided services data, unable to calculate provision."
        ) from capacity_error

    if calculate_demands:

        buildings = get_demands(
            buildings_with_people=buildings, normative=demand_normative, population=population, do_restoration=True
        )

    if calculate_graph:
        city_graph = dngraph.get_intermodal_graph_from_osm(city_osm_id)

    dngraph.set_graph(city_graph)

    if calculate_matrix:
        adjacency_matrix = dngraph.get_adjacency_matrix(buildings, services, weight_adjacency_matrix)

    return get_service_provision(
        services=services,
        adjacency_matrix=adjacency_matrix,
        demanded_buildings=buildings,
        threshold=threshold,
        calculation_type=calculation_type,
    )
