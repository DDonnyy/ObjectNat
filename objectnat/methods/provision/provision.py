import geopandas as gpd
import numpy as np
import pandas as pd

from objectnat import config

from .provision_model import Provision

logger = config.logger


def get_service_provision(
    buildings: gpd.GeoDataFrame,
    adjacency_matrix: pd.DataFrame,
    services: gpd.GeoDataFrame,
    threshold: int,
    buildings_demand_column: str = "demand",
    services_capacity_column: str = "capacity",
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Compute service provision between demand locations (buildings) and service facilities.

    The function implements a **gravity-based allocation model**: service capacity is
    distributed across nearby demand points with weights that **decay with the square
    of distance (or generalized cost)**. Closer buildings receive proportionally
    higher shares of the available capacity.

    Args:
        buildings (gpd.GeoDataFrame):
            GeoDataFrame of **demand locations** (e.g., residential buildings).
            Must include a numeric column with demand values
            (see ``buildings_demand_column``).
        adjacency_matrix (pd.DataFrame):
            A rectangular DataFrame representing **OD (origin–destination) costs**
            between ``buildings`` (rows) and ``services`` (columns).
            Units must match ``threshold`` (e.g., minutes or meters).
            Missing or infinite values (``NaN`` or ``inf``) are treated as **unreachable**.
            The row index must match ``buildings.index`` and column index must
            match ``services.index``.
        services (gpd.GeoDataFrame):
            GeoDataFrame of **service facilities** (e.g., schools, clinics).
            Must include a numeric column with service capacity
            (see ``services_capacity_column``).
        threshold (int):
            Maximum allowed cost value for assignment.
            Any OD entry **greater than this threshold** is considered unreachable.
            Units are the same as in ``adjacency_matrix``.
        buildings_demand_column (str):
            Column name of building demand values. Default is ``"demand"``.
        services_capacity_column (str):
            Column name of service capacity values. Default is ``"capacity"``.

    Returns:
        Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
            A tuple of three GeoDataFrames:

            - **buildings**: input buildings with updated provision metrics.
            - **services**: input services with updated load and capacity metrics.
            - **links**: building–service links within the threshold, containing
              allocated demand shares and distances/costs based on the gravity model.

    Notes:
        - The model is **gravity-based**, with cost weights decaying by the **square of distance**.
        - Unreachable OD pairs (``NaN`` or ``inf``) are ignored.
        - The function does not perform routing; it expects a precomputed OD matrix.
    """
    buildings = buildings.copy()
    services = services.copy()
    adjacency_matrix = adjacency_matrix.copy()
    buildings["demand"] = buildings[buildings_demand_column]
    services["capacity"] = services[services_capacity_column]

    provision_buildings, provision_services, provision_links = Provision(
        services=services,
        demanded_buildings=buildings,
        adjacency_matrix=adjacency_matrix,
        threshold=threshold,
    ).run()
    return provision_buildings, provision_services, provision_links


def clip_provision(
    buildings: gpd.GeoDataFrame, services: gpd.GeoDataFrame, links: gpd.GeoDataFrame, selection_zone: gpd.GeoDataFrame
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Clip service provision results to a specific geographic boundary.

    Keeps only:
      * buildings that intersect the ``selection_zone``;
      * links that connect to the kept buildings;
      * services referenced by those links.

    Args:
        buildings:
            GeoDataFrame of buildings **after** :func:`get_service_provision`.
        services:
            GeoDataFrame of services **after** :func:`get_service_provision`.
        links:
            GeoDataFrame of building–service links from
            :func:`get_service_provision`. Must include indices or columns
            to match buildings and services.
        selection_zone:
            GeoDataFrame (polygon or multipolygon) defining the clipping area.

    Returns:
        Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
            The filtered subsets of buildings, services, and links
            that fall inside the specified zone.

    Notes:
        - The function performs **spatial filtering only**.
          It does **not** recompute or redistribute demand/supply.
    """
    assert selection_zone.crs == buildings.crs == services.crs == links.crs, (
        f"CRS mismatch: buildings_crs:{buildings.crs}, "
        f"links_crs:{links.crs} , "
        f"services_crs:{services.crs}, "
        f"selection_zone_crs:{selection_zone.crs}"
    )
    buildings = buildings.copy()
    links = links.copy()
    services = services.copy()

    s = buildings.intersects(selection_zone.union_all())
    buildings = buildings.loc[s[s].index]
    links = links[links["building_index"].isin(buildings.index.tolist())]
    services_to_keep = set(links["service_index"].tolist())
    services.drop(list(set(services.index.tolist()) - services_to_keep), inplace=True)
    return buildings, services, links


def recalculate_links(
    buildings: gpd.GeoDataFrame, services: gpd.GeoDataFrame, links: gpd.GeoDataFrame, new_max_dist: float
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Recalculate provision aggregates after tightening the accessibility threshold.

    Removes all links whose cost (distance or time) exceeds ``new_max_dist``, then
    updates demand and capacity aggregates accordingly. This is done **without
    redistributing** removed demand to alternative services.

    Args:
        buildings:
            GeoDataFrame of buildings after :func:`get_service_provision`.
            Expected to include provision-related fields such as demand, demand_left,
            supplied demand, and average distance/cost.

        services:
            GeoDataFrame of services after :func:`get_service_provision`, with
            fields describing remaining capacity and service load.

        links:
            GeoDataFrame of building–service links containing at least:

            - ``building_index``
            - ``service_index``
            - ``distance`` (or time cost, in the same units as ``new_max_dist``)
            - ``demand`` (assigned portion)

        new_max_dist:
            New maximum allowed cost value (same units as OD/threshold).
            Links with cost **greater than** this value are removed.

    Returns:
        tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
            - **buildings**: updated aggregate demand metrics and recalculated
              average cost.
            - **services**: updated load and capacity fields after freeing excess capacity.
            - **links**: subset of links that remain within the new threshold.

    Notes:
        - If no links exceed ``new_max_dist``, the function logs a warning
          and returns the original inputs unchanged.
        - Average cost values are recomputed based on remaining links.
          If a building has no remaining assigned demand, ``avg_dist`` becomes ``NaN``.
        - Removed demand is **not reallocated** to other services.
    """
    buildings = buildings.copy()
    services = services.copy()
    links = links.copy()

    links_to_recalculate = links[links["distance"] > new_max_dist]
    if len(links_to_recalculate) == 0:
        logger.warning("To clip distance exceeds max links distance, returning full provision")
        return buildings, services, links

    links_to_keep = links[links["distance"] <= new_max_dist]
    free_demand = links_to_recalculate.groupby("building_index").agg({"demand": list, "distance": list})
    free_demand["distance"] = free_demand.apply(
        lambda x: sum((x1 * x2) for x1, x2 in zip(x.demand, x.distance)), axis=1
    )
    free_demand["demand"] = free_demand["demand"].apply(sum)
    free_demand = free_demand.reindex(buildings.index, fill_value=0)
    new_sum_time = (buildings["supplied_demands_within"] + buildings["supplied_demands_without"]) * buildings[
        "avg_dist"
    ] - free_demand["distance"]

    buildings["demand_left"] = buildings["demand_left"] + free_demand["demand"]
    buildings["supplied_demands_without"] = buildings["supplied_demands_without"] - free_demand["demand"]
    buildings["avg_dist"] = new_sum_time / (
        buildings["supplied_demands_without"] + buildings["supplied_demands_within"]
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
