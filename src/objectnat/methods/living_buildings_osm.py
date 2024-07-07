import geopandas as gpd
import osm2geojson
import osmnx as ox
import pandas as pd
import requests
from loguru import logger
from shapely import MultiPolygon, Polygon

from ..utils import get_utm_crs_for_4326_gdf

DEFAULT_OVERPASS_URL = "http://overpass-api.de/api/interpreter"


def get_terr_polygon_osm_name(territory_name: str) -> Polygon | MultiPolygon:
    """
    Retrieve the polygon geometry of a specified territory using its name from OpenStreetMap.

    Parameters
    ----------
    territory_name : str
        The name of the territory to retrieve the polygon for.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the polygon geometry of the specified territory.

    Examples
    --------
    >>> territory_name = "Saint-Petersburg, Russia"
    >>> polygon = get_terr_polygon_osm_name(territory_name)
    """
    logger.info(f"Retrieving polygon geometry for '{territory_name}'")
    place = ox.geocode_to_gdf(territory_name)
    polygon = place.geometry.values[0]
    return polygon.unary_union


def get_terr_polygon_osm_id(osm_id: int) -> Polygon | MultiPolygon:
    """
    Retrieve the polygon geometry of a specified territory using its OSM ID from OpenStreetMap.

    Parameters
    ----------
    osm_id : int
        The OpenStreetMap ID of the territory to retrieve the polygon for.

    Returns
    -------
    Polygon | MultiPolygon
        A Polygon or MultiPolygon geometry of the specified territory.

    Examples
    --------
    >>> osm_id = 421007
    >>> polygon = get_terr_polygon_osm_id(osm_id)
    """
    overpass_query = f"""
                [out:json];
                (
                    relation({osm_id});
                );
                out geom;
                """
    logger.info(f"Retrieving polygon geometry for osm id '{osm_id}'")
    result = requests.get(DEFAULT_OVERPASS_URL, params={"data": overpass_query}, timeout=500)
    json_result = result.json()
    json_result = osm2geojson.json2geojson(json_result)
    json_result = gpd.GeoDataFrame.from_features(json_result["features"]).set_crs(4326)
    return json_result.geometry.unary_union


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
    df[population_column] = 0.0
    df.loc[df["is_living"] == 1, population_column] = df[df["is_living"] == 1].apply(
        lambda row: (
            3
            if ((row["area"] <= 400) & (row["building:levels"] <= 2))
            else (row["building:levels"] * row["area"] * 0.8 / area_per_person)
        ),
        axis=1,
    )
    df[population_column] = df[population_column].fillna(0).round(0).astype(int)
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
    if osm_territory_id is not None:
        polygon = get_terr_polygon_osm_id(osm_territory_id)
        return download_buildings(
            terr_polygon=polygon,
            area_per_person=area_per_person,
            is_living_column=is_living_column,
            population_column=population_column,
        )

    if osm_territory_name is not None:
        polygon = get_terr_polygon_osm_name(osm_territory_name)
        return download_buildings(
            terr_polygon=polygon,
            area_per_person=area_per_person,
            is_living_column=is_living_column,
            population_column=population_column,
        )

    logger.info("Downloading buildings from OpenStreetMap and counting population...")
    buildings = ox.features_from_polygon(terr_polygon, tags={"building": True})
    print(buildings.crs)
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
        logger.info("Done!")
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
                "is_living",
                "building:levels_is_real",
                "approximate_pop",
                "geometry",
            ]
        ]
