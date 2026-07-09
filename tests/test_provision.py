import os

import numpy as np
import pytest
from matplotlib import pyplot as plt

from objectnat import (
    ProvisionResult,
    clip_provision,
    get_provision_buildings,
    get_provision_links,
    get_provision_services,
    get_service_provision,
    recalculate_links,
)
from tests.conftest import output_dir


@pytest.fixture(scope="module")
def basic_provision_result(buildings_data, services_data, matrix_time_data):
    return get_service_provision(
        buildings=buildings_data, services=services_data, distance_matrix=matrix_time_data, threshold=10
    )


def materialize_provision(buildings_data, services_data, provision_result):
    return (
        get_provision_buildings(buildings_data, provision_result),
        get_provision_services(services_data, provision_result),
        get_provision_links(buildings_data, services_data, provision_result),
    )


def test_no_demand(buildings_data, services_data, matrix_time_data):
    buildings_data = buildings_data.copy()
    buildings_data["demand"] = 0
    provision_result = get_service_provision(
        buildings=buildings_data, services=services_data, distance_matrix=matrix_time_data, threshold=10
    )
    links_prov = get_provision_links(buildings_data, services_data, provision_result)
    assert links_prov.empty


def test_no_capacity(buildings_data, services_data, matrix_time_data):
    services_data = services_data.copy()
    services_data["capacity"] = 0
    provision_result = get_service_provision(
        buildings=buildings_data, services=services_data, distance_matrix=matrix_time_data, threshold=10
    )
    links_prov = get_provision_links(buildings_data, services_data, provision_result)
    assert links_prov.empty


def test_get_service_provision(basic_provision_result, buildings_data, services_data):
    assert isinstance(basic_provision_result, ProvisionResult)
    assert basic_provision_result.flow.sparse.to_coo().nnz > 0
    build_prov, services_prov, links_prov = materialize_provision(buildings_data, services_data, basic_provision_result)

    assert build_prov is not None
    assert services_prov is not None
    assert links_prov is not None

    assert np.isin(["service_load", "capacity_left"], services_prov.columns).all()
    assert np.isin(["min_dist", "avg_dist", "provision_value"], build_prov.columns).all()
    assert np.isin(["distance", "demand"], links_prov.columns).all()

    assert not build_prov.empty
    assert not services_prov.empty
    assert not links_prov.empty

    visualize_provision(buildings_data, build_prov, services_prov, links_prov, filename_suffix="initial")


def test_recalculate_links(basic_provision_result, buildings_data, services_data):
    threshold = 10
    provision_result2 = recalculate_links(basic_provision_result, threshold)
    build_prov = get_provision_buildings(buildings_data, basic_provision_result)
    services_prov = get_provision_services(services_data, basic_provision_result)
    build_prov2, services_prov2, links_prov2 = materialize_provision(buildings_data, services_data, provision_result2)

    assert len(build_prov) == len(build_prov2)
    assert len(services_prov) == len(services_prov2)
    assert (links_prov2["distance"] <= threshold).all()

    visualize_provision(
        buildings_data,
        build_prov2,
        services_prov2,
        links_prov2,
        title_suffix=f"(Recalculated with threshold={threshold})",
        filename_suffix="recalculated",
    )


def test_clip_links(basic_provision_result, buildings_data, services_data):
    build_prov, services_prov, links_prov = materialize_provision(buildings_data, services_data, basic_provision_result)

    to_clip_gdf = build_prov.iloc[:20].copy()
    to_clip_gdf["geometry"] = to_clip_gdf["geometry"].buffer(500)

    build_prov_clipped, services_prov_clipped, links_prov_clipped = clip_provision(
        build_prov, services_prov, links_prov, to_clip_gdf
    )

    assert build_prov_clipped is not None
    assert services_prov_clipped is not None
    assert links_prov_clipped is not None

    visualize_provision(
        buildings_data,
        build_prov_clipped,
        services_prov_clipped,
        links_prov_clipped,
        title_suffix="(Clipped by buildings)",
        filename_suffix="clipped",
    )


def visualize_provision(initial_buildings, build_prov, services_prov, links_prov, title_suffix="", filename_suffix=""):
    local_crs = initial_buildings.estimate_utm_crs()
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    initial_buildings = initial_buildings.to_crs(local_crs)
    build_prov = build_prov.to_crs(local_crs)
    build_prov.geometry = build_prov.geometry.buffer(10, resolution=4)
    services_prov = services_prov.to_crs(local_crs)
    links_prov = links_prov.to_crs(local_crs)

    minx, miny, maxx, maxy = initial_buildings.total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    initial_buildings.plot(ax=ax, edgecolor="gray", facecolor="none", linewidth=0.2)
    build_prov.plot(
        ax=ax,
        column="avg_dist",
        cmap="RdYlGn_r",
        alpha=0.8,
        label="Buildings",
        legend=True,
        legend_kwds={"label": "average distance in building to chosen services(meters)", "shrink": 0.5},
    )
    links_prov.plot(ax=ax, column="service_index", cmap="prism", linewidth=0.15, alpha=0.2, label="Links")
    services_prov.plot(ax=ax, color="red", markersize=10, label="Services")

    ax.set_title(f"Service provision {title_suffix}")
    ax.legend()
    ax.set_axis_off()

    output_path = os.path.join(output_dir, f"service_provision_{filename_suffix}.png")
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()

    return output_path
