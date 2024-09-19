from typing import Tuple

import geopandas as gpd
import numpy as np
import pandas as pd

from .city_provision import CityProvision


def get_service_provision(
    buildings: gpd.GeoDataFrame,
    adjacency_matrix: pd.DataFrame,
    services: gpd.GeoDataFrame,
    threshold: int,
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Calculate load from buildings with demands on the given services using the distances matrix between them.

    Args:
        services (gpd.GeoDataFrame): GeoDataFrame of services
        adjacency_matrix (pd.DataFrame): DataFrame representing the adjacency matrix
        buildings (gpd.GeoDataFrame): GeoDataFrame of demanded buildings
        threshold (int): Threshold value
    Returns:
        Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]: Tuple of GeoDataFrames representing provision
        buildings, provision services, and provision links
    """
    provision_buildings, provision_services, provision_links = CityProvision(
        services=services,
        demanded_buildings=buildings,
        adjacency_matrix=adjacency_matrix,
        threshold=threshold,
    ).get_provisions()
    return provision_buildings, provision_services, provision_links


def clip_provision(
    buildings: gpd.GeoDataFrame, services: gpd.GeoDataFrame, links: gpd.GeoDataFrame, selection_zone: gpd.GeoDataFrame
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:

    assert (
        selection_zone.crs == buildings.crs == services.crs == links.crs
    ), f"CRS mismatch: buildings_crs:{buildings.crs}, links_crs:{links.crs} , services_crs:{services.crs}, selection_zone_crs:{selection_zone.crs}"

    s = buildings.intersects(selection_zone.unary_union)
    buildings = buildings.loc[s[s].index]
    links = links[links["building_index"].isin(buildings.index.tolist())]
    services_to_keep = set(links["service_index"].tolist())
    services.drop(list(set(services.index.tolist()) - services_to_keep), inplace=True)
    return buildings, services, links


def recalculate_links(
    buildings: gpd.GeoDataFrame, services: gpd.GeoDataFrame, links: gpd.GeoDataFrame, new_max_dist: float
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    buildings = buildings.copy()
    services = services.copy()
    links = links.copy()

    max_dist = links["distance"].max()
    assert new_max_dist <= max_dist, "New distance exceeds max links distance"

    links_to_recalculate = links[links["distance"] > new_max_dist]
    links_to_keep = links[links["distance"] <= new_max_dist]

    free_demand = links_to_recalculate.groupby("building_index").agg({"demand": list, "distance": list})
    free_demand["distance"] = free_demand.apply(
        lambda x: sum((x1 * x2) for x1, x2 in zip(x.demand, x.distance)), axis=1
    )
    free_demand["demand"] = free_demand["demand"].apply(sum)
    free_demand = free_demand.reindex(buildings.index, fill_value=0)
    new_sum_time = (buildings["supplyed_demands_within"] + buildings["supplyed_demands_without"]) * buildings[
        "avg_dist"
    ] - free_demand["distance"]

    buildings["demand_left"] = buildings["demand_left"] + free_demand["demand"]
    buildings["supplyed_demands_without"] = buildings["supplyed_demands_without"] - free_demand["demand"]
    buildings["avg_dist"] = new_sum_time / (
        buildings["supplyed_demands_without"] + buildings["supplyed_demands_within"]
    )
    buildings["avg_dist"] = buildings.apply(
        lambda x: np.nan if (x["demand"] == x["demand_left"]) else round(x["avg_dist"], 2), axis=1
    )

    free_capacity = links_to_recalculate.groupby("service_index").agg({"demand": "sum"})
    free_capacity = free_capacity.reindex(services.index, fill_value=0)
    services["capacity_left"] = services["capacity_left"] + free_capacity["demand"]
    services["carried_capacity_without"] = services["carried_capacity_without"] - free_capacity["demand"]
    services["service_load"] = services["service_load"] - free_capacity["demand"]

    return buildings, services, links_to_keep
