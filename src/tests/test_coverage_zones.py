from objectnat import get_isochrone_zone_coverage, get_radius_zone_coverage


def test_service_coverage_isochrones(intermodal_graph, services_data):
    build_point = buildings_data.sample(n=5)
    isochrones, stops, routes = get_accessibility_isochrones(
        points=build_point, weight_type="time_min", weight_value=15, graph_nx=intermodal_graph
    )

    assert isochrones is not None
    assert len(isochrones) == 5
