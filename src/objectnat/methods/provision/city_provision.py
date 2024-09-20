# pylint: disable=singleton-comparison
from typing import Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import LineString

from objectnat import config

from .provision_exceptions import CapacityKeyError, DemandKeyError

logger = config.logger

from pandarallel import pandarallel


class CityProvision:
    """
    Represents the logic for city provision calculations using a gravity or linear model.

    Args:
        services (InstanceOf[gpd.GeoDataFrame]): GeoDataFrame representing the services available in the city.
        demanded_buildings (InstanceOf[gpd.GeoDataFrame]): GeoDataFrame representing the buildings with demands for services.
        adjacency_matrix (InstanceOf[pd.DataFrame]): DataFrame representing the adjacency matrix between buildings.
        threshold (int): Threshold value for the provision calculations.
        calculation_type (str, optional): Type of calculation ("gravity" or "linear"). Defaults to "gravity".

    Returns:
        CityProvision: The CityProvision object.

    Raises: KeyError: If the 'demand' column is missing in the provided 'demanded_buildings' GeoDataFrame,
    or if the 'capacity' column is missing in the provided 'services' GeoDataFrame. ValueError: If the 'capacity'
    column in 'services' or 'demand' column  'demanded_buildings' GeoDataFrame has no valid value.
    """

    _destination_matrix = None

    def __init__(
        self,
        services: gpd.GeoDataFrame,
        demanded_buildings: gpd.GeoDataFrame,
        adjacency_matrix: pd.DataFrame,
        threshold: int,
    ):
        self.services = self.ensure_services(services)
        self.demanded_buildings = self.ensure_buildings(demanded_buildings)
        self.adjacency_matrix = self.delete_useless_matrix_rows(adjacency_matrix.copy(), demanded_buildings, services)
        self.threshold = threshold
        self.check_crs(self.demanded_buildings, self.services)
        pandarallel.initialize(progress_bar=False, verbose=0)

    @staticmethod
    def ensure_buildings(v: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        if "demand" not in v.columns:
            raise DemandKeyError
        v = v.copy()
        v["demand_left"] = v["demand"]
        return v

    @staticmethod
    def ensure_services(v: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        if "capacity" not in v.columns:
            raise CapacityKeyError
        v = v.copy()
        v["capacity_left"] = v["capacity"]
        return v

    @staticmethod
    def check_crs(demanded_buildings, services):
        assert (
            demanded_buildings.crs == services.crs
        ), f"\nThe CRS in the provided geodataframes are different.\nBuildings CRS:{demanded_buildings.crs}\nServices CRS:{services.crs} \n"

    @staticmethod
    def delete_useless_matrix_rows(adjacency_matrix, demanded_buildings, services):
        adjacency_matrix.index = adjacency_matrix.index.astype(int)

        builds_indexes = set(demanded_buildings.index.astype(int).tolist())
        rows = set(adjacency_matrix.index.astype(int).tolist())
        dif = rows ^ builds_indexes
        adjacency_matrix.drop(index=(list(dif)), axis=0, inplace=True)

        service_indexes = set(services.index.astype(int).tolist())
        columns = set(adjacency_matrix.columns.astype(int).tolist())
        dif = columns ^ service_indexes
        adjacency_matrix.drop(columns=(list(dif)), axis=0, inplace=True)
        return adjacency_matrix

    def get_provisions(self) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
        self._destination_matrix = pd.DataFrame(
            0,
            index=self.adjacency_matrix.columns,
            columns=self.adjacency_matrix.index,
        )
        self.adjacency_matrix = self.adjacency_matrix.transpose()
        logger.debug(
            "Calculating provision from {} services to {} buildings.",
            len(self.services),
            len(self.demanded_buildings),
        )
        self.adjacency_matrix = self.adjacency_matrix.where(self.adjacency_matrix <= self.threshold * 3, np.inf)

        self._destination_matrix = self._provision_loop_gravity(
            self.demanded_buildings.copy(),
            self.services.copy(),
            self.adjacency_matrix.copy() + 1,
            self.threshold,
            self._destination_matrix.copy(),
        )
        _additional_options(
            self.demanded_buildings,
            self.services,
            self.adjacency_matrix,
            self._destination_matrix,
            self.threshold,
        )

        # self.demanded_buildings['provision_value'] = self.demanded_buildings['provision_value'].fillna(0)
        # self.services = self.services.fillna(0)

        return (
            self.demanded_buildings,
            self.services,
            _calc_links(
                self._destination_matrix,
                self.services,
                self.demanded_buildings,
                self.adjacency_matrix,
            ),
        )

    def _provision_loop_gravity(
        self,
        houses_table: gpd.GeoDataFrame,
        services_table: gpd.GeoDataFrame,
        distance_matrix: pd.DataFrame,
        selection_range,
        destination_matrix: pd.DataFrame,
        best_houses=0.9,
    ):
        def apply_function_based_on_size(df, func, axis, threshold=500):
            if len(df) > threshold:
                return df.parallel_apply(func, axis=axis)
            else:
                return df.apply(func, axis=axis)

        def _calculate_flows_y(loc):
            import numpy as np
            import pandas as pd

            c = services_table.loc[loc.name]["capacity_left"]
            p = 1 / loc / loc
            p = p / p.sum()
            threshold = p.quantile(best_houses)
            p = p[p >= threshold]
            p = p / p.sum()
            if p.sum() == 0:
                return loc
            rng = np.random.default_rng(seed=0)
            r = pd.Series(0, p.index)
            choice = np.unique(rng.choice(p.index, int(c), p=p.values), return_counts=True)
            choice = r.add(pd.Series(choice[1], choice[0]), fill_value=0)

            return choice

        def _balance_flows_to_demands(loc):
            import numpy as np
            import pandas as pd

            d = houses_table.loc[loc.name]["demand_left"]
            loc = loc[loc > 0]
            if loc.sum() > 0:
                p = loc / loc.sum()
                rng = np.random.default_rng(seed=0)
                r = pd.Series(0, p.index)
                choice = np.unique(rng.choice(p.index, int(d), p=p.values), return_counts=True)
                choice = r.add(pd.Series(choice[1], choice[0]), fill_value=0)
                choice = pd.Series(
                    data=np.minimum(loc.sort_index().values, choice.sort_index().values),
                    index=loc.sort_index().index,
                )
                return choice
            return loc

        temp_destination_matrix = apply_function_based_on_size(
            distance_matrix, lambda x: _calculate_flows_y(x[x <= selection_range]), 1
        )

        temp_destination_matrix = temp_destination_matrix.fillna(0)

        temp_destination_matrix = apply_function_based_on_size(temp_destination_matrix, _balance_flows_to_demands, 0)

        temp_destination_matrix = temp_destination_matrix.fillna(0)
        destination_matrix = destination_matrix.add(temp_destination_matrix, fill_value=0)

        axis_1 = destination_matrix.sum(axis=1)
        axis_0 = destination_matrix.sum(axis=0)

        services_table["capacity_left"] = services_table["capacity"].subtract(axis_1, fill_value=0)
        houses_table["demand_left"] = houses_table["demand"].subtract(axis_0, fill_value=0)

        distance_matrix = distance_matrix.drop(
            index=services_table[services_table["capacity_left"] == 0].index.values,
            columns=houses_table[houses_table["demand_left"] == 0].index.values,
            errors="ignore",
        )

        distance_matrix = distance_matrix.loc[~(distance_matrix == np.inf).all(axis=1)]
        distance_matrix = distance_matrix.loc[:, ~(distance_matrix == np.inf).all(axis=0)]

        selection_range += selection_range

        if best_houses > 0.1:
            best_houses -= 0.1
            if best_houses <= 0.1:
                best_houses = 0
        if len(distance_matrix.columns) > 0 and len(distance_matrix.index) > 0:
            return self._provision_loop_gravity(
                houses_table, services_table, distance_matrix, selection_range, destination_matrix, best_houses
            )
        return destination_matrix


def _calc_links(
    destination_matrix: pd.DataFrame,
    services: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame,
    distance_matrix: pd.DataFrame,
):
    def subfunc(loc):
        try:
            return [
                {
                    "building_index": int(k),
                    "demand": int(v),
                    "service_index": int(loc.name),
                }
                for k, v in loc.to_dict().items()
            ]
        except:
            return np.NaN

    def subfunc_geom(loc):
        return LineString(
            (
                buildings_["geometry"][loc["building_index"]],
                services_["geometry"][loc["service_index"]],
            )
        )

    buildings_ = buildings.copy()
    services_ = services.copy()
    buildings_.geometry = buildings_.representative_point()
    services_.geometry = services_.representative_point()
    flat_matrix = destination_matrix.transpose().apply(lambda x: subfunc(x[x > 0]), result_type="reduce")

    distribution_links = gpd.GeoDataFrame(data=[item for sublist in list(flat_matrix) for item in sublist])

    distribution_links["distance"] = distribution_links.apply(
        lambda x: distance_matrix.loc[x["service_index"]][x["building_index"]],
        axis=1,
        result_type="reduce",
    )

    sel = distribution_links["building_index"].isin(buildings_.index.values) & distribution_links["service_index"].isin(
        services_.index.values
    )
    sel = distribution_links.loc[sel[sel].index.values]
    distribution_links = distribution_links.set_geometry(sel.apply(subfunc_geom, axis=1)).set_crs(buildings_.crs)
    distribution_links["distance"] = distribution_links["distance"].astype(float).round(2)
    return distribution_links


def _additional_options(
    buildings,
    services,
    matrix,
    destination_matrix,
    normative_distance,
):
    buildings["avg_dist"] = 0
    buildings["supplyed_demands_within"] = 0
    buildings["supplyed_demands_without"] = 0
    services["carried_capacity_within"] = 0
    services["carried_capacity_without"] = 0
    for i in range(len(destination_matrix)):
        loc = destination_matrix.iloc[i]
        distances_all = matrix.loc[loc.name]
        distances = distances_all[distances_all <= normative_distance]
        s = matrix.loc[loc.name] <= normative_distance
        within = loc[s]
        without = loc[~s]
        within = within[within > 0]
        without = without[without > 0]
        buildings["avg_dist"] = (
            buildings["avg_dist"]
            .add(distances.multiply(within, fill_value=0), fill_value=0)
            .add(distances_all.multiply(without, fill_value=0), fill_value=0)
        )
        buildings["demand_left"] = buildings["demand_left"].sub(within.add(without, fill_value=0), fill_value=0)
        buildings["supplyed_demands_within"] = buildings["supplyed_demands_within"].add(within, fill_value=0)
        buildings["supplyed_demands_without"] = buildings["supplyed_demands_without"].add(without, fill_value=0)

        services.at[loc.name, "capacity_left"] = (
            services.at[loc.name, "capacity_left"] - within.add(without, fill_value=0).sum()
        )
        services.at[loc.name, "carried_capacity_within"] = (
            services.at[loc.name, "carried_capacity_within"] + within.sum()
        )
        services.at[loc.name, "carried_capacity_without"] = (
            services.at[loc.name, "carried_capacity_without"] + without.sum()
        )
    buildings["avg_dist"] = (buildings["avg_dist"] / (buildings["demand"] - buildings["demand_left"])).astype(
        np.float32
    )
    buildings["avg_dist"] = buildings.apply(
        lambda x: np.nan if (x["demand"] == x["demand_left"]) else round(x["avg_dist"], 2), axis=1
    )
    buildings["provison_value"] = (buildings["supplyed_demands_within"] / buildings["demand"]).astype(float).round(2)
    services["service_load"] = (services["capacity"] - services["capacity_left"]).astype(np.uint16)
    buildings["supplyed_demands_within"] = buildings["supplyed_demands_within"].astype(np.uint16)
    buildings["supplyed_demands_without"] = buildings["supplyed_demands_without"].astype(np.uint16)
    services["carried_capacity_within"] = services["carried_capacity_within"].astype(np.uint16)
    services["carried_capacity_without"] = services["carried_capacity_without"].astype(np.uint16)
