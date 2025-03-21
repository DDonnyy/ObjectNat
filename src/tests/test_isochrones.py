import time

import pytest
from iduedu import get_intermodal_graph
from shapely import Point
from shapely.ops import unary_union
import geopandas as gpd

from objectnat import get_accessibility_isochrones


@pytest.fixture(scope="module")
def intermodal_graph(buildings_data, services_data):
    time.sleep(0.5)
    polygon = unary_union(
        [buildings_data.to_crs(4326).geometry.to_list() + services_data.to_crs(4326).geometry.to_list()]
    ).convex_hull.buffer(0.001)
    print("\n Downloading intermodal graph for bounds \n")
    return get_intermodal_graph(polygon=polygon)


def test_1point_isochrone(intermodal_graph, buildings_data):
    build_point = buildings_data.sample(n=1)
    isochrones, stops, routes = get_accessibility_isochrones(
        points=build_point,
        weight_type="time_min",
        weight_value=15,
        graph_nx=intermodal_graph
    )

    assert isochrones is not None
    assert len(isochrones) == 1


def test_5points_isochrone(intermodal_graph, buildings_data):
    build_point = buildings_data.sample(n=5)
    isochrones, stops, routes = get_accessibility_isochrones(
        points=build_point,
        weight_type="time_min",
        weight_value=15,
        graph_nx=intermodal_graph
    )

    assert isochrones is not None
    assert len(isochrones) == 5
