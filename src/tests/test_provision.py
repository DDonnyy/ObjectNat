from .fixtures import buildings_data,services_data,matrix_time_data
from objectnat import get_service_provision
import numpy as np

def test_get_service_provision(buildings_data, services_data, matrix_time_data):

    build_prov, services_prov, links_prov = get_service_provision(
        buildings=buildings_data,
        services=services_data,
        adjacency_matrix=matrix_time_data,
        threshold=10,
    )

    assert build_prov is not None
    assert services_prov is not None
    assert links_prov is not None

    assert np.isin(['service_load', 'capacity_left'], services_prov.columns).all()
    assert np.isin(['min_dist','avg_dist','provision_value'], build_prov.columns).all()
    assert np.isin(['distance','demand'], links_prov.columns).all()

    assert not build_prov.empty
    assert not services_prov.empty
    assert not links_prov.empty


from objectnat import recalculate_links

def test_recalculate_links(buildings_data, services_data, matrix_time_data):

    build_prov, services_prov, links_prov = get_service_provision(
        buildings=buildings_data,
        services=services_data,
        adjacency_matrix=matrix_time_data,
        threshold=10,
    )

    build_prov2, services_prov2, links_prov2 = recalculate_links(
        build_prov, services_prov, links_prov, 15
    )

