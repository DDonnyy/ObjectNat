import os

import pytest
from matplotlib import pyplot as plt
from pyproj.exceptions import CRSError

from objectnat import get_accessibility_isochrone_stepped, get_accessibility_isochrones
from tests.conftest import output_dir


def test_1point_isochrone_radius(intermodal_osm_1114252, gdf_1point, buildings_data):
    weight_value = 15
    isochrones, stops, routes = get_accessibility_isochrones(
        isochrone_type="radius",
        points=gdf_1point,
        weight_value=weight_value,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
    )
    assert isochrones is not None
    assert len(isochrones) == 1
    visualize_isochrones(
        isochrones,
        gdf_1point,
        routes,
        buildings_data,
        title_suffix=f"(radius mode, {weight_value} minutes)",
        filename_suffix=f"radius_{weight_value}_min",
    )


def test_1point_isochrone_ways(intermodal_osm_1114252, gdf_1point, buildings_data):
    gdf_1point = gdf_1point.to_crs(4326)
    weight_value = 15
    isochrones, stops, routes = get_accessibility_isochrones(
        isochrone_type="ways",
        points=gdf_1point,
        weight_value=weight_value,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
    )
    assert isochrones is not None
    assert len(isochrones) == 1
    assert isochrones.crs == gdf_1point.crs
    assert stops.crs == gdf_1point.crs
    assert routes.crs == gdf_1point.crs
    visualize_isochrones(
        isochrones,
        gdf_1point,
        routes,
        buildings_data,
        title_suffix=f"(ways mode, {weight_value} minutes)",
        filename_suffix=f"ways_{weight_value}_min",
    )


def test_3point_isochrone_radius(intermodal_osm_1114252, gdf_3points):
    isochrones, stops, routes = get_accessibility_isochrones(
        isochrone_type="radius",
        points=gdf_3points,
        weight_value=15,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
    )
    assert isochrones is not None
    assert len(isochrones) == 3


def test_3point_isochrone_ways(intermodal_osm_1114252, gdf_3points):
    isochrones, stops, routes = get_accessibility_isochrones(
        isochrone_type="ways",
        points=gdf_3points,
        weight_value=5,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
    )
    assert isochrones is not None
    assert len(isochrones) == 3


def test_wrong_graph_crs(intermodal_osm_1114252, gdf_1point):
    wrong_graph = intermodal_osm_1114252.copy()
    wrong_graph.graph["crs"] = "Wrong CRS"
    with pytest.raises(CRSError) as _:
        _ = get_accessibility_isochrones(
            isochrone_type="ways",
            points=gdf_1point,
            weight_value=15,
            weight_type="time_min",
            nx_graph=wrong_graph,
        )
    wrong_graph.graph = {}
    with pytest.raises(ValueError) as _:
        _ = get_accessibility_isochrones(
            isochrone_type="ways",
            points=gdf_1point,
            weight_value=15,
            weight_type="time_min",
            nx_graph=wrong_graph,
        )


def test_isochrone_stepped_radius(intermodal_osm_1114252, gdf_1point, buildings_data):
    weight_value = 15
    stepped_iso, stops, routes = get_accessibility_isochrone_stepped(
        isochrone_type="radius",
        point=gdf_1point,
        weight_value=15,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
        step=3,
    )
    assert stepped_iso is not None
    assert len(stepped_iso) == 5

    visualize_stepped_isochrones(
        stepped_iso,
        gdf_1point,
        routes,
        buildings_data,
        title_suffix=f"(radius mode, {weight_value} minutes)",
        filename_suffix=f"stepped_radius_{weight_value}_min",
    )


def test_isochrone_stepped_ways(intermodal_osm_1114252, gdf_1point, buildings_data):
    weight_value = 15
    stepped_iso, stops, routes = get_accessibility_isochrone_stepped(
        isochrone_type="ways",
        point=gdf_1point,
        weight_value=15,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
        step=3,
    )
    assert stepped_iso is not None
    assert len(stepped_iso) == 5

    visualize_stepped_isochrones(
        stepped_iso,
        gdf_1point,
        routes,
        buildings_data,
        title_suffix=f"(ways mode, {weight_value} minutes)",
        filename_suffix=f"stepped_ways_{weight_value}_min",
    )


def test_isochrone_stepped_separate(intermodal_osm_1114252, gdf_1point, buildings_data):
    weight_value = 15
    stepped_iso, stops, routes = get_accessibility_isochrone_stepped(
        isochrone_type="separate",
        point=gdf_1point,
        weight_value=15,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
        step=3,
    )
    assert stepped_iso is not None
    visualize_stepped_isochrones(
        stepped_iso,
        gdf_1point,
        routes,
        buildings_data,
        title_suffix=f"(separate mode, {weight_value} minutes)",
        filename_suffix=f"stepped_separate_{weight_value}_min",
    )


def test_multipoint_in_stepped(intermodal_osm_1114252, gdf_3points):
    stepped_iso, stops, routes = get_accessibility_isochrone_stepped(
        isochrone_type="radius",
        point=gdf_3points,
        weight_value=15,
        weight_type="time_min",
        nx_graph=intermodal_osm_1114252,
        step=3,
    )
    assert stepped_iso is not None
    assert len(stepped_iso) == 5


def visualize_isochrones(isochrones, point_from, routes, buildings_data, title_suffix="", filename_suffix=""):
    local_crs = buildings_data.estimate_utm_crs()

    fig, ax = plt.subplots(figsize=(10, 10))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    isochrones = isochrones.to_crs(local_crs)
    buildings_data = buildings_data.to_crs(local_crs)
    routes = routes.to_crs(local_crs)
    point_from = point_from.to_crs(local_crs)

    minx, miny, maxx, maxy = buildings_data.total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    buildings_data.plot(ax=ax, edgecolor="gray", facecolor="none", linewidth=0.5)
    isochrones.plot(
        ax=ax,
        alpha=0.9,
        color="#d1e5f0",  # Светло-голубая заливка
        edgecolor="#2166ac",  # Темно-синяя граница
        linewidth=0.8,
        label="Isochrones",
    )
    routes.plot(ax=ax, column="type", linewidth=0.5, label="Public transport routes")
    point_from.plot(ax=ax, color="red", markersize=50, label="Start point")
    ax.set_title(f"Isochrone {title_suffix}")
    ax.legend()
    ax.set_axis_off()

    output_path = os.path.join(output_dir, f"isochrone_{filename_suffix}.png")
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()


def visualize_stepped_isochrones(
    stepped_isochrones, point_from, routes, buildings_data, title_suffix="", filename_suffix=""
):

    local_crs = buildings_data.estimate_utm_crs()

    stepped_isochrones = stepped_isochrones.to_crs(local_crs)
    buildings_data = buildings_data.to_crs(local_crs)
    routes = routes.to_crs(local_crs)
    point_from = point_from.to_crs(local_crs)

    fig, ax = plt.subplots(figsize=(12, 10))

    minx, miny, maxx, maxy = buildings_data.total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    buildings_data.plot(ax=ax, edgecolor="gray", facecolor="none", linewidth=0.5)

    stepped_isochrones.plot(
        ax=ax,
        column="dist",
        cmap="viridis",
        alpha=0.7,
        edgecolor="black",
        linewidth=0.2,
        legend=True,
        legend_kwds={"label": "Distance (meters)", "shrink": 0.5},
        label="Stepped isochrone",
    )
    routes.plot(ax=ax, column="type", linewidth=0.5, label="Public transport routes")
    point_from.plot(ax=ax, color="red", markersize=50, label="Start point")

    ax.set_title(f"Stepped isochrone {title_suffix}")
    ax.set_axis_off()

    output_path = os.path.join(output_dir, f"stepped_isochrone_{filename_suffix}.png")
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()
