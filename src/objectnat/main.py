from dongraphio import *
from provisio import demands_from_buildings_by_normative, get_service_provision
import geopandas as gpd

import population_restorator.balancer.territories as b_terr
import population_restorator.balancer.houses as b_build


def balance_population(
    living_buildings: gpd.GeoDataFrame,
    population: int | None = None,
    territories: gpd.GeoDataFrame | None = None,
) -> gpd.GeoDataFrame:
    """Balance city population into living buildings according to their spot area * `storeys_count`.

    :param living_buildings:gpd.GeoDataFrame:  (multi)polygons of buildings where people live.
    Must contain either "storeys_count" or "living_area" attributes.

    :param  population:int: total city population.

    :param territories:gpd.GeoDataFrame: if more detailed population is known, it can be set via this parameter
    in "population" field. In case territories are given, only buildings inside them will be populated.

    :returns gpd.GeoDataFrame:
    """
    if population is None:
        raise ValueError("Neither city nor population values are set")
    if living_buildings.shape[0] == 0:
        raise ValueError("Living buildings are not set")
    if "storeys_count" not in living_buildings.columns:
        raise ValueError("Living buildings are missing `storeys_count` attribute")
    living_buildings = living_buildings.copy()
    living_buildings["living_area"] = living_buildings.area * living_buildings["storeys_count"]
    living_buildings = living_buildings[living_buildings["living_area"].notna()]
    if living_buildings.shape[0] == 0:
        raise ValueError("No buildings contain valid ")
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

demands_from_buildings_by_normative()