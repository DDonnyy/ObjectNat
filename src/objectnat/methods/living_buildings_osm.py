import geopandas as gpd
import numpy as np
import osmnx as ox
import pandas as pd
from iduedu import get_boundary
from shapely import MultiPolygon, Polygon

from objectnat import config

from ..utils import get_utm_crs_for_4326_gdf

logger = config.logger


def eval_is_living(row: gpd.GeoSeries):
    """
    Determine if a building is used for residential purposes based on its attributes.

    Parameters
    ----------
    row : gpd.GeoSeries
        A GeoSeries representing a row in a GeoDataFrame, containing building attributes.

    Returns
    -------
    bool
        A boolean indicating whether the building is used for residential purposes.

    Examples
    --------
    >>> buildings = download_buildings(osm_territory_id=421007)
    >>> buildings['is_living'] = buildings.apply(eval_is_living, axis=1)
    """
    if row["building"] in (
        "apartments",
        "house",
        "residential",
        "detached",
        "dormitory",
        "semidetached_house",
        "bungalow",
        "cabin",
        "farm",
    ):
        return True
    else:
        return False


def eval_population(source: gpd.GeoDataFrame, population_column: str, area_per_person: float = 33):
    """
    Estimate the population of buildings in a GeoDataFrame based on their attributes.

    Parameters
    ----------
    source : gpd.GeoDataFrame
        A GeoDataFrame containing building geometries and attributes.
    population_column : str
        The name of the column where the estimated population will be stored.
    area_per_person : float
        The standart living space per person im m², (default is 33)
    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with an added column for estimated population.

    Raises
    ------
    RuntimeError
        If the 'building:levels' column is not present in the provided GeoDataFrame.

    Examples
    --------
    >>> source = gpd.read_file('buildings.shp')
    >>> source['is_living'] = source.apply(eval_is_living, axis=1)
    >>> population_df = eval_population(source, 'approximate_pop')
    """
    if "building:levels" not in source.columns:
        raise RuntimeError("No 'building:levels' column in provided GeoDataFrame")
    df = source.copy()
    local_utm_crs = get_utm_crs_for_4326_gdf(source.to_crs(4326))
    df["area"] = df.to_crs(local_utm_crs.to_epsg()).geometry.area.astype(float)
    df["building:levels_is_real"] = df["building:levels"].apply(lambda x: False if pd.isna(x) else True)
    df["building:levels"] = df["building:levels"].fillna(1)
    df["building:levels"] = pd.to_numeric(df["building:levels"], errors="coerce")
    df = df.dropna(subset=["building:levels"])
    df["building:levels"] = df["building:levels"].astype(int)
    df[population_column] = np.nan
    df.loc[df["is_living"] == 1, population_column] = (
        df[df["is_living"] == 1]
        .apply(
            lambda row: (
                3
                if ((row["area"] <= 400) & (row["building:levels"] <= 2))
                else (row["building:levels"] * row["area"] * 0.8 / area_per_person)
            ),
            axis=1,
        )
        .round(0)
        .astype(int)
    )
    return df


def download_buildings(
    osm_territory_id: int | None = None,
    osm_territory_name: str | None = None,
    terr_polygon: Polygon | MultiPolygon | None = None,
    is_living_column: str = "is_living",
    population_column: str = "approximate_pop",
    area_per_person: float = 33,
) -> gpd.GeoDataFrame | None:
    """
    Download building geometries and evaluate 'is_living' and 'population' attributes for a specified territory from OpenStreetMap.

    Parameters
    ----------
    osm_territory_id : int, optional
        The OpenStreetMap ID of the territory to download buildings for.
    osm_territory_name : str, optional
        The name of the territory to download buildings for.
    terr_polygon : Polygon or MultiPolygon, optional
        A Polygon or MultiPolygon geometry defining the territory to download buildings for.
    is_living_column : str, optional
        The name of the column indicating whether a building is residential (default is "is_living").
    population_column : str, optional
        The name of the column for storing estimated population (default is "approximate_pop").
    area_per_person : float
        The standart living space per person im m², (default is 33)
    Returns
    -------
    gpd.GeoDataFrame or None
        A GeoDataFrame containing building geometries and attributes, or None if no buildings are found or an error occurs.

    Examples
    --------
    >>> buildings_df = download_buildings(osm_territory_name="Saint-Petersburg, Russia")
    >>> buildings_df.head()
    """
    polygon = get_boundary(osm_territory_id, osm_territory_name, terr_polygon)

    logger.debug("Downloading buildings from OpenStreetMap and counting population...")
    buildings = ox.features_from_polygon(polygon, tags={"building": True})
    if not buildings.empty:
        buildings = buildings.loc[
            (buildings["geometry"].geom_type == "Polygon") | (buildings["geometry"].geom_type == "MultiPolygon")
        ]
    if buildings.empty:
        logger.warning(f"There are no buildings in the specified territory. Output GeoDataFrame is empty.")
        return buildings
    else:
        buildings[is_living_column] = buildings.apply(eval_is_living, axis=1)
        buildings = eval_population(buildings, population_column, area_per_person)
        buildings.reset_index(drop=True, inplace=True)
        logger.debug("Done!")
        return buildings[
            [
                "building",
                "addr:street",
                "addr:housenumber",
                "amenity",
                "area",
                "name",
                "building:levels",
                "leisure",
                "design:year",
                is_living_column,
                "building:levels_is_real",
                population_column,
                "geometry",
            ]
        ]
