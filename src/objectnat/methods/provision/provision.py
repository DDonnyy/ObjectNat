from typing import Literal, Tuple

import geopandas as gpd
import pandas as pd

from .city_provision import CityProvision


def get_service_provision(
    buildings: gpd.GeoDataFrame,
    adjacency_matrix: pd.DataFrame,
    services: gpd.GeoDataFrame,
    threshold: int,
    calculation_type: Literal["gravity", "linear"] = "gravity",
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Calculate load from buildings with demands on the given services using the distances matrix between them.

    Args:
        services (gpd.GeoDataFrame): GeoDataFrame of services
        adjacency_matrix (pd.DataFrame): DataFrame representing the adjacency matrix
        buildings (gpd.GeoDataFrame): GeoDataFrame of demanded buildings
        threshold (int): Threshold value
        calculation_type (str): Calculation type for provision, might be "gravity" or "linear"
    Returns:
        Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]: Tuple of GeoDataFrames representing provision
        buildings, provision services, and provision links
    """
    provision_buildings, provision_services, provision_links = CityProvision(
        services=services,
        demanded_buildings=buildings,
        adjacency_matrix=adjacency_matrix,
        threshold=threshold,
        calculation_type=calculation_type,
    ).get_provisions()
    return provision_buildings, provision_services, provision_links


def is_shown(
    buildings: gpd.GeoDataFrame, services: gpd.GeoDataFrame, links: gpd.GeoDataFrame, selection_zone: gpd.GeoDataFrame
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    buildings.reset_index(inplace=True)
    buildings = gpd.overlay(buildings, selection_zone, how="intersection")
    buildings.set_index("index", inplace=True)
    links = links[links["building_index"].isin(buildings.index.tolist())]
    services_to_keep = set(links["service_index"].tolist())
    services.drop(list(set(services.index.tolist()) - services_to_keep), inplace=True)
    return buildings, services, links
