import geopandas as gpd
from loguru import logger
from provisio import demands_from_buildings_by_normative

from .balanced_buildings import get_balanced_buildings


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
