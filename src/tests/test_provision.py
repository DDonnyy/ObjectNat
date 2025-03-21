import pytest

from objectnat import get_service_provision, recalculate_links, clip_provision
import numpy as np


@pytest.fixture(scope="module")
def basic_provision(buildings_data, services_data, matrix_time_data):
    build_prov, services_prov, links_prov = get_service_provision(
        buildings=buildings_data,
        services=services_data,
        adjacency_matrix=matrix_time_data,
        threshold=10,
    )
    return build_prov, services_prov, links_prov


def test_get_service_provision(basic_provision):
    build_prov, services_prov, links_prov = basic_provision

    assert build_prov is not None
    assert services_prov is not None
    assert links_prov is not None

    assert np.isin(["service_load", "capacity_left"], services_prov.columns).all()
    assert np.isin(["min_dist", "avg_dist", "provision_value"], build_prov.columns).all()
    assert np.isin(["distance", "demand"], links_prov.columns).all()

    assert not build_prov.empty
    assert not services_prov.empty
    assert not links_prov.empty


def test_recalculate_links(basic_provision):
    build_prov, services_prov, links_prov = basic_provision

    build_prov2, services_prov2, links_prov2 = recalculate_links(build_prov, services_prov, links_prov, 15)

    assert len(build_prov) == len(build_prov2)
    assert len(services_prov) == len(services_prov2)
    assert (links_prov2["distance"] <= 15).all()


def test_clip_links(basic_provision):
    build_prov, services_prov, links_prov = basic_provision

    to_clip_gdf = build_prov.iloc[:20].copy()
    to_clip_gdf['geometry'] = to_clip_gdf['geometry'].buffer(500)

    build_prov_clipped, services_prov_clipped, links_prov_clipped = clip_provision(
        build_prov, services_prov, links_prov, to_clip_gdf
    )

    assert build_prov_clipped is not None
    assert services_prov_clipped is not None
    assert links_prov_clipped is not None
