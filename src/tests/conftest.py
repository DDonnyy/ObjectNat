import pickle
from pathlib import Path

import pytest
import os
import geopandas as gpd
import pandas as pd
from iduedu import get_intermodal_graph, get_boundary,config
from shapely import Point

logger = config.logger

path_to_data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../examples/examples_data/")


@pytest.fixture(scope="session")
def buildings_data():
    data_path = os.path.join(path_to_data, "buildings.parquet")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Файл {data_path} не найден!")
    buildings_data = gpd.read_parquet(data_path)
    buildings_data.index = buildings_data.index.astype(int)
    return buildings_data


@pytest.fixture(scope="session")
def services_data():
    data_path = os.path.join(path_to_data, "services.parquet")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Файл {data_path} не найден!")
    services_data = gpd.read_parquet(data_path)
    services_data.index = services_data.index.astype(int)
    return services_data


@pytest.fixture(scope="session")
def matrix_time_data():
    data_path = os.path.join(path_to_data, "matrix_time.parquet")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Файл {data_path} не найден!")
    matrix_time_data = pd.read_parquet(data_path)
    matrix_time_data.index = matrix_time_data.index.astype(int)
    matrix_time_data.columns = matrix_time_data.columns.astype(int)
    return matrix_time_data


@pytest.fixture(scope="session")
def trees_data():
    data_path = os.path.join(path_to_data, "trees.parquet")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Файл {data_path} не найден!")
    return gpd.read_parquet(data_path)


@pytest.fixture(scope="session")
def boundary_osm_1114252():
    return get_boundary(osm_id=1114252)


@pytest.fixture(scope="session")
def intermodal_osm_1114252(boundary_osm_1114252):
    cache_dir = Path("cache")
    cache_dir.mkdir(exist_ok=True)
    cache_file = cache_dir / "graph.pickle"
    if cache_file.exists():
        logger.info(f'Found cached graph {cache_file}')
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    graph = get_intermodal_graph(polygon=boundary_osm_1114252, clip_by_bounds=True)
    with open(cache_file, "wb") as f:
        logger.info(f'Saving graph to {cache_file}')
        pickle.dump(graph, f)
    return graph


@pytest.fixture(scope="session")
def gdf_1point():
    return gpd.GeoDataFrame(geometry=[Point(30.27060176, 59.93546846)], crs=4326)


@pytest.fixture(scope="session")
def gdf_3points():
    points = [Point(30.27060176, 59.93546846), Point(30.295606, 59.9439234), Point(30.2312112, 59.9482336)]
    return gpd.GeoDataFrame(geometry=[points], crs=4326)
