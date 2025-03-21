import pytest
import os
import geopandas as gpd
import pandas as pd


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
