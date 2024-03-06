from typing import Optional, Tuple

import geopandas as gpd
import networkx as nx
import pandas as pd
import population_restorator.balancer.houses as b_build
import population_restorator.balancer.territories as b_terr
from dongraphio import DonGraphio, GraphType
from loguru import logger
from provisio import demands_from_buildings_by_normative, get_service_provision
from provisio.provisio import CityProvision
from provisio.provisio_exceptions import DemandKeyError, CapacityKeyError

from pydantic import ValidationError


def get_balanced_buildings(
    living_buildings: gpd.GeoDataFrame,
    population: int,
    territories: gpd.GeoDataFrame | None = None,
) -> gpd.GeoDataFrame:
    """
    Balance city population into living buildings according to their spot area * `storeys_count`.

    Args:
        living_buildings (gpd.GeoDataFrame): (multi)polygons of buildings where people live. Must contain either
            "storeys_count" or "living_area" attributes.

        population (int): Total city population

        territories (gpd.GeoDataFrame, optional): If more detailed population is known, it can be set via this parameter
            in "population" field. In case territories are given, only buildings inside them will be populated.

    Returns:
        gpd.GeoDataFrame: The balanced living buildings.

    Raises:
        ValueError: If neither city nor population values are set, or if living buildings are not set,
            or if living buildings are missing the `storeys_count` attribute, or if buildings contain no valid value.
    """

    if population is None:
        raise ValueError("Neither city nor population values are set")
    if living_buildings.shape[0] == 0:
        raise ValueError("Living buildings are not set")
    if "living_area" not in living_buildings.columns:
        if "storeys_count" not in living_buildings.columns:
            raise ValueError("Living buildings are missing `storeys_count` attribute")
    logger.debug(f"Evacuating {population} residents into the provided building")
    living_buildings = living_buildings.copy()
    if "living_area" not in living_buildings.columns:
        living_buildings["living_area"] = living_buildings.area * living_buildings["storeys_count"]
    living_buildings = living_buildings[living_buildings["living_area"].notna()]
    if living_buildings.shape[0] == 0:
        raise ValueError("Buildings contain no valid value")
    if territories is None:
        city_territory = b_terr.Territory("city", population, None, living_buildings)
    else:
        inner_territories = [
            b_terr.Territory(
                f"terr_{i}",
                terr.get("population"),
                None,
                living_buildings[terr.contains(living_buildings.centroid)].copy(),
            )
            for i, (_, terr) in enumerate(territories.iterrows())
        ]
        city_territory = b_terr.Territory("city", population, inner_territories, None)

    b_terr.balance_territories(city_territory)
    b_build.balance_houses(city_territory)

    houses = city_territory.get_all_houses()
    return houses


def get_demands(
    buildings_with_people: gpd.GeoDataFrame,
    normative: float,
    do_restoration: bool = False,
    population: int | None = None,
) -> gpd.GeoDataFrame:
    """
    Calculate demands according to the normative for buildings with people.

    Args:
        buildings_with_people (gpd.GeoDataFrame): GeoDataFrame representing buildings with people.

        normative (float): The normative value for calculating demands.

        do_restoration (bool, optional): Flag to indicate whether to attempt restoration if the 'population'
            attribute is missing. Defaults to False.

        population (int, optional): Total population of provided buildings data for restoration. Defaults to None.

    Returns:
        gpd.GeoDataFrame: The buildings with calculated demands.
    """
    try:
        logger.debug(f"Calculating demands according to the normative {normative} ...")
        buildings = demands_from_buildings_by_normative(buildings_with_people, normative)
        buildings["demand"] = buildings["demand"].round().astype(int)
    except KeyError as e:
        if do_restoration:
            logger.warning(
                "Buildings are not evacuated, the 'population' attribute is missing, attempting to evacuate..."
            )
            buildings = get_balanced_buildings(buildings_with_people, population)
            buildings = get_demands(buildings, normative)
        else:
            raise e
    return buildings


def get_accessibility_isochrones(
    graph_type: list[GraphType],
    x_from: float,
    y_from: float,
    weight_value: int,
    weight_type: str,
    city_graph: nx.MultiDiGraph,
    city_crs: int,
) -> Tuple[gpd.GeoDataFrame, Optional[gpd.GeoDataFrame], Optional[gpd.GeoDataFrame]]:
    """
    Calculate accessibility isochrones based on provided city graph for selected graph_type enum

    Args:
        graph_type (List[GraphType]): List of graph types to calculate isochrones for.
        x_from (float): X coordinate of the starting point in the corresponding coordinate system.
        y_from (float): Y coordinate of the starting point in the corresponding coordinate system.
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
    return dongraphio.get_accessibility_isochrones(graph_type, x_from, y_from, weight_value, weight_type)


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


def get_adjacency_matrix(
    buildings_from: gpd.GeoDataFrame,
    services_to: gpd.GeoDataFrame,
    weight: str,
    city_crs: int | None = None,
    nx_graph: nx.MultiDiGraph | None = None,
    dongraphio: DonGraphio | None = None,
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
            return dongraphio.get_adjacency_matrix(buildings_from, services_to, weight)

        dongraphio = DonGraphio(city_crs)
        dongraphio.set_graph(nx_graph)
        return dongraphio.get_adjacency_matrix(buildings_from, services_to, weight)
    except ValidationError as e:
        logger.error("Function get_adjacency_matrix() missing 'weight' argument")
        raise e


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
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
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

    if not adjacency_matrix:
        logger.warning("The adjacency matrix is not provided.")
        if not weight_adjacency_matrix:
            raise RuntimeError("No weight type ('time_min' or 'length_meter') was provided.")
        if not city_graph:
            logger.warning("The graph is not provided, attempting to load data from OSM.")
            if not city_osm_id:
                raise RuntimeError("osm_id of the city is not provided, unable to retrieve data.")
            calculate_graph = True
        calculate_matrix = True
        adjacency_matrix = pd.DataFrame(data=0, index=[], columns=buildings.index.astype(int).tolist())

    try:
        CityProvision(
            services=services,
            demanded_buildings=buildings,
            adjacency_matrix=adjacency_matrix,
            threshold=threshold,
            calculation_type=calculation_type,
        )
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

    return get_service_provision(services, adjacency_matrix, buildings, threshold, calculation_type)
