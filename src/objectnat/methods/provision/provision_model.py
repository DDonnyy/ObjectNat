# pylint: disable=singleton-comparison
from typing import Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from shapely import LineString

from objectnat import config

from .provision_exceptions import CapacityKeyError, DemandKeyError

logger = config.logger


class Provision:
    """
    Represents the logic for city provision calculations using a gravity or linear model.

    Args:
        services (InstanceOf[gpd.GeoDataFrame]): GeoDataFrame representing the services available in the city.
        demanded_buildings (InstanceOf[gpd.GeoDataFrame]): GeoDataFrame representing the buildings with demands for services.
        adjacency_matrix (InstanceOf[pd.DataFrame]): DataFrame representing the adjacency matrix between buildings.
        threshold (int): Threshold value for the provision calculations.
        calculation_type (str, optional): Type of calculation ("gravity" or "linear"). Defaults to "gravity".

    Returns:
        Provision: The CityProvision object.

    Raises: KeyError: If the 'demand' column is missing in the provided 'demanded_buildings' GeoDataFrame,
    or if the 'capacity' column is missing in the provided 'services' GeoDataFrame. ValueError: If the 'capacity'
    column in 'services' or 'demand' column  'demanded_buildings' GeoDataFrame has no valid value.
    """

    destination_matrix = None

    def __init__(
        self,
        services: gpd.GeoDataFrame,
        demanded_buildings: gpd.GeoDataFrame,
        adjacency_matrix: pd.DataFrame,
        threshold: int,
    ):
        self.services = self.ensure_services(services.copy())
        self.demanded_buildings = self.ensure_buildings(demanded_buildings.copy())
        self.adjacency_matrix = self.delete_useless_matrix_rows_columns(
            adjacency_matrix.copy(), demanded_buildings, services
        ).copy()
        self.threshold = threshold
        self.check_crs(self.demanded_buildings, self.services)
        pandarallel.initialize(progress_bar=False, verbose=0, use_memory_fs=config.pandarallel_use_file_system)

    @staticmethod
    def ensure_buildings(v: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        if "demand" not in v.columns:
            raise DemandKeyError
        v["demand_left"] = v["demand"]
        return v

    @staticmethod
    def ensure_services(v: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        if "capacity" not in v.columns:
            raise CapacityKeyError
        v["capacity_left"] = v["capacity"]
        return v

    @staticmethod
    def check_crs(demanded_buildings, services):
        assert (
            demanded_buildings.crs == services.crs
        ), f"\nThe CRS in the provided geodataframes are different.\nBuildings CRS:{demanded_buildings.crs}\nServices CRS:{services.crs} \n"

    @staticmethod
    def delete_useless_matrix_rows_columns(adjacency_matrix, demanded_buildings, services):
        adjacency_matrix.index = adjacency_matrix.index.astype(int)

        builds_indexes = set(demanded_buildings.index.astype(int).tolist())
        rows = set(adjacency_matrix.index.astype(int).tolist())
        dif = rows ^ builds_indexes
        adjacency_matrix.drop(index=(list(dif)), axis=0, inplace=True)

        service_indexes = set(services.index.astype(int).tolist())
        columns = set(adjacency_matrix.columns.astype(int).tolist())
        dif = columns ^ service_indexes
        adjacency_matrix.drop(columns=(list(dif)), axis=0, inplace=True)
        return adjacency_matrix.transpose()

    def run(self) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:

        def apply_function_based_on_size(df, func, axis, threshold=100):
            if len(df) > threshold:
                return df.parallel_apply(func, axis=axis)
            return df.apply(func, axis=axis)

        def calculate_flows_y(loc):
            import numpy as np  # pylint: disable=redefined-outer-name,reimported,import-outside-toplevel
            import pandas as pd  # pylint: disable=redefined-outer-name,reimported,import-outside-toplevel

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

        def balance_flows_to_demands(loc):
            import numpy as np  # pylint: disable=redefined-outer-name,reimported,import-outside-toplevel
            import pandas as pd  # pylint: disable=redefined-outer-name,reimported,import-outside-toplevel

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

        logger.debug(
            f"Calculating provision from {len(self.services)} services to {len(self.demanded_buildings)} buildings."
        )

        distance_matrix = self.adjacency_matrix
        destination_matrix = pd.DataFrame(
            0,
            index=distance_matrix.index,
            columns=distance_matrix.columns,
            dtype=int,
        )
        distance_matrix = distance_matrix.where(distance_matrix <= self.threshold * 3, np.inf)

        houses_table = self.demanded_buildings[["demand", "demand_left"]].copy()
        services_table = self.services[["capacity", "capacity_left"]].copy()
        distance_matrix = distance_matrix.drop(
            index=services_table[services_table["capacity_left"] == 0].index.values,
            columns=houses_table[houses_table["demand_left"] == 0].index.values,
            errors="ignore",
        )
        distance_matrix = distance_matrix.loc[~(distance_matrix == np.inf).all(axis=1)]
        distance_matrix = distance_matrix.loc[:, ~(distance_matrix == np.inf).all(axis=0)]

        distance_matrix = distance_matrix + 1
        selection_range = (self.threshold + 1) / 2
        best_houses = 0.9
        while len(distance_matrix.columns) > 0 and len(distance_matrix.index) > 0:
            objects_n = sum(distance_matrix.shape)
            logger.debug(
                f"Matrix shape: {distance_matrix.shape},"
                f" Total objects: {objects_n},"
                f" Selection range: {selection_range},"
                f" Best houses: {best_houses}"
            )

            temp_destination_matrix = apply_function_based_on_size(
                distance_matrix, lambda x: calculate_flows_y(x[x <= selection_range]), 1
            )

            temp_destination_matrix = temp_destination_matrix.fillna(0)
            temp_destination_matrix = apply_function_based_on_size(temp_destination_matrix, balance_flows_to_demands, 0)
            temp_destination_matrix = temp_destination_matrix.fillna(0)
            temp_destination_matrix_aligned = temp_destination_matrix.reindex(
                index=destination_matrix.index, columns=destination_matrix.columns, fill_value=0
            )
            del temp_destination_matrix
            destination_matrix_np = destination_matrix.to_numpy()
            temp_destination_matrix_np = temp_destination_matrix_aligned.to_numpy()
            del temp_destination_matrix_aligned
            destination_matrix = pd.DataFrame(
                destination_matrix_np + temp_destination_matrix_np,
                index=destination_matrix.index,
                columns=destination_matrix.columns,
            )
            del destination_matrix_np, temp_destination_matrix_np
            axis_1 = destination_matrix.sum(axis=1).astype(int)
            axis_0 = destination_matrix.sum(axis=0).astype(int)

            services_table["capacity_left"] = services_table["capacity"].subtract(axis_1, fill_value=0)
            houses_table["demand_left"] = houses_table["demand"].subtract(axis_0, fill_value=0)
            del axis_1, axis_0
            distance_matrix = distance_matrix.drop(
                index=services_table[services_table["capacity_left"] == 0].index.values,
                columns=houses_table[houses_table["demand_left"] == 0].index.values,
                errors="ignore",
            )
            distance_matrix = distance_matrix.loc[~(distance_matrix == np.inf).all(axis=1)]
            distance_matrix = distance_matrix.loc[:, ~(distance_matrix == np.inf).all(axis=0)]

            selection_range *= 1.5
            if best_houses <= 0.1:
                best_houses = 0
            else:
                objects_n_new = sum(distance_matrix.shape)
                best_houses = objects_n_new / (objects_n / best_houses)

        logger.debug("Done!")
        del distance_matrix, houses_table, services_table
        self.destination_matrix = destination_matrix

        _additional_options(
            self.demanded_buildings,
            self.services,
            self.adjacency_matrix,
            self.destination_matrix,
            self.threshold,
        )

        return (
            self.demanded_buildings,
            self.services,
            _calc_links(
                self.destination_matrix,
                self.services,
                self.demanded_buildings,
                self.adjacency_matrix,
            ),
        )


def _calc_links(
    destination_matrix: pd.DataFrame,
    services: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame,
    distance_matrix: pd.DataFrame,
):
    buildings_ = buildings.copy()
    services_ = services.copy()
    buildings_.geometry = buildings_.representative_point()
    services_.geometry = services_.representative_point()

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
                buildings_.geometry[loc["building_index"]],
                services_.geometry[loc["service_index"]],
            )
        )

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
    for i, loc in destination_matrix.iterrows():
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
    buildings["min_dist"] = matrix.min(axis=0).replace(np.inf, None)
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
    logger.debug("Done adding additional options")
