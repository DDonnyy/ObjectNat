import geopandas as gpd
import population_restorator.balancer.houses as b_build
import population_restorator.balancer.territories as b_terr
from loguru import logger


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
        ValueError: If population is missing, or if living buildings are not set,
            or if living buildings are missing the `storeys_count` attribute, or if buildings contain no valid value.
    """

    if population is None:
        raise ValueError("Population value is missing")
    if living_buildings.shape[0] == 0:
        raise ValueError("Living buildings are not set")
    if "living_area" not in living_buildings.columns and "storeys_count" not in living_buildings.columns:
        raise ValueError("Living buildings are missing `storeys_count` attribute")
    logger.debug(f"Evacuating {population} residents into the provided building")
    indexes = living_buildings.index
    living_buildings = living_buildings.copy()
    living_buildings.reset_index(drop=True, inplace=True)
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
    houses.set_index(indexes, drop=True, inplace=True)
    return houses
