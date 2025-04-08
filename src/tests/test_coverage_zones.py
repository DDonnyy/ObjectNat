import pytest
from pyproj.exceptions import CRSError
from objectnat import get_graph_coverage, get_radius_coverage
import geopandas as gpd


def test_graph_time_min(services_data, intermodal_osm_1114252, boundary_osm_1114252):
    zone = gpd.GeoDataFrame(geometry=[boundary_osm_1114252], crs=4326)
    result = get_graph_coverage(
        gdf_from=services_data,
        nx_graph=intermodal_osm_1114252,
        weight_type="time_min",
        weight_value_cutoff=10,
        zone=zone,
    )
    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) == len(services_data)


def test_graph_length_meter(services_data, intermodal_osm_1114252):
    result = get_graph_coverage(
        gdf_from=services_data, nx_graph=intermodal_osm_1114252, weight_type="length_meter", weight_value_cutoff=600
    )
    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) == len(services_data)


def test_graph_same_crs(services_data, intermodal_osm_1114252):
    services_data = services_data.to_crs(3857)
    result = get_graph_coverage(
        gdf_from=services_data, nx_graph=intermodal_osm_1114252, weight_type="length_meter", weight_value_cutoff=600
    )
    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) == len(services_data)
    assert result.crs == services_data.crs


def test_wrong_graph_crs(services_data, intermodal_osm_1114252):
    wrong_graph = intermodal_osm_1114252.copy()
    wrong_graph.graph["crs"] = "Wrong CRS"
    with pytest.raises(CRSError) as _:
        _ = get_graph_coverage(
            gdf_from=services_data, nx_graph=wrong_graph, weight_type="length_meter", weight_value_cutoff=600
        )
    wrong_graph.graph = {}
    with pytest.raises(ValueError) as _:
        _ = get_graph_coverage(
            gdf_from=services_data, nx_graph=wrong_graph, weight_type="length_meter", weight_value_cutoff=600
        )


def test_radius_coverage(services_data):
    result = get_radius_coverage(services_data, radius=1000)
    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) == len(services_data)
    assert result.crs == services_data.crs