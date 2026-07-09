from dataclasses import dataclass

import geopandas as gpd
import pandas as pd
from shapely import LineString

from objectnat import config

logger = config.logger


@dataclass(slots=True)
class ProvisionResult:
    """Result of service provision calculation."""

    flow: pd.DataFrame
    demand_rows: pd.DataFrame
    capacity_rows: pd.DataFrame
    distance_matrix: pd.DataFrame
    threshold: float


def _validate_provision_result(provision_result: ProvisionResult) -> None:
    if not isinstance(provision_result, ProvisionResult):
        raise TypeError(f"provision_result must be ProvisionResult, got {type(provision_result).__name__}")
    if provision_result.flow.index.has_duplicates:
        raise ValueError("provision_result.flow.index must be unique")
    if provision_result.flow.columns.has_duplicates:
        raise ValueError("provision_result.flow.columns must be unique")
    if provision_result.demand_rows.index.has_duplicates:
        raise ValueError("provision_result.demand_rows.index must be unique")
    if provision_result.capacity_rows.index.has_duplicates:
        raise ValueError("provision_result.capacity_rows.index must be unique")
    if not provision_result.flow.index.equals(provision_result.distance_matrix.index):
        raise ValueError("provision_result.flow.index must match provision_result.distance_matrix.index")
    if not provision_result.flow.columns.equals(provision_result.distance_matrix.columns):
        raise ValueError("provision_result.flow.columns must match provision_result.distance_matrix.columns")


def _drop_overwritten_columns(df: gpd.GeoDataFrame | pd.DataFrame, columns_to_join: pd.Index, df_name: str):
    overwritten_columns = df.columns.intersection(columns_to_join)
    if overwritten_columns.empty:
        return df

    logger.warning(
        f"{df_name} already contains provision columns that will be overwritten: {overwritten_columns.tolist()}"
    )
    return df.drop(columns=overwritten_columns)


def get_provision_buildings(
    buildings_gdf: gpd.GeoDataFrame | pd.DataFrame,
    provision_result: ProvisionResult,
) -> gpd.GeoDataFrame | pd.DataFrame:
    """Join provision demand metrics to buildings."""
    if not isinstance(buildings_gdf, (pd.DataFrame, gpd.GeoDataFrame)):
        raise TypeError(f"buildings_gdf must be DataFrame or GeoDataFrame, got {type(buildings_gdf).__name__}")
    _validate_provision_result(provision_result)
    if buildings_gdf.index.has_duplicates:
        raise ValueError("buildings_gdf.index must be unique")

    buildings_gdf = _drop_overwritten_columns(buildings_gdf, provision_result.demand_rows.columns, "buildings_gdf")
    return buildings_gdf.join(provision_result.demand_rows, how="left")


def get_provision_services(
    services_gdf: gpd.GeoDataFrame | pd.DataFrame,
    provision_result: ProvisionResult,
) -> gpd.GeoDataFrame | pd.DataFrame:
    """Join provision capacity metrics to services."""
    if not isinstance(services_gdf, (pd.DataFrame, gpd.GeoDataFrame)):
        raise TypeError(f"services_gdf must be DataFrame or GeoDataFrame, got {type(services_gdf).__name__}")
    _validate_provision_result(provision_result)
    if services_gdf.index.has_duplicates:
        raise ValueError("services_gdf.index must be unique")

    services_gdf = _drop_overwritten_columns(services_gdf, provision_result.capacity_rows.columns, "services_gdf")
    return services_gdf.join(provision_result.capacity_rows, how="left")


def get_provision_links(
    buildings_gdf: gpd.GeoDataFrame,
    services_gdf: gpd.GeoDataFrame,
    provision_result: ProvisionResult,
) -> gpd.GeoDataFrame:
    """Build building-service link geometries from non-zero provision flows."""
    if not isinstance(buildings_gdf, gpd.GeoDataFrame):
        raise TypeError(f"buildings_gdf must be GeoDataFrame, got {type(buildings_gdf).__name__}")
    if not isinstance(services_gdf, gpd.GeoDataFrame):
        raise TypeError(f"services_gdf must be GeoDataFrame, got {type(services_gdf).__name__}")
    _validate_provision_result(provision_result)

    if buildings_gdf.index.has_duplicates:
        raise ValueError("buildings_gdf.index must be unique")
    if services_gdf.index.has_duplicates:
        raise ValueError("services_gdf.index must be unique")

    flow_coo = provision_result.flow.sparse.to_coo()
    if flow_coo.nnz == 0:
        logger.warning(
            "Unable to create provision links - no demand could be matched with service locations. "
            "This is likely because either demand or capacity is zero, or no services are reachable."
        )
        return gpd.GeoDataFrame()

    services_for_links = services_gdf
    if buildings_gdf.crs is not None and services_gdf.crs is not None and buildings_gdf.crs != services_gdf.crs:
        services_for_links = services_gdf.to_crs(buildings_gdf.crs)

    building_points = buildings_gdf.representative_point()
    service_points = services_for_links.representative_point()
    building_labels = provision_result.flow.index.to_numpy()[flow_coo.row]
    service_labels = provision_result.flow.columns.to_numpy()[flow_coo.col]
    distance_values = provision_result.distance_matrix.to_numpy(dtype=float)[flow_coo.row, flow_coo.col]

    links = pd.DataFrame(
        {
            "building_index": building_labels,
            "service_index": service_labels,
            "demand": flow_coo.data.astype(int),
            "distance": distance_values.round(2),
        }
    )
    links = links[
        (links["building_index"].isin(buildings_gdf.index)) & (links["service_index"].isin(services_gdf.index))
    ]
    if links.empty:
        return gpd.GeoDataFrame()

    links["geometry"] = [
        LineString([building_points.loc[building_index], service_points.loc[service_index]])
        for building_index, service_index in zip(links["building_index"], links["service_index"])
    ]
    return gpd.GeoDataFrame(links, geometry="geometry", crs=buildings_gdf.crs)
