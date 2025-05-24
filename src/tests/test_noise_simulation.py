import os

import geopandas as gpd
import pytest
from matplotlib import pyplot as plt
from shapely import Point

from objectnat import config, simulate_noise
from objectnat.methods.noise import InvalidStepError
from objectnat.methods.noise.noise_reduce import get_air_resist_ratio
from tests.conftest import output_dir

logger = config.logger


def test_basic_functionality(gdf_1point, buildings_data, trees_data):
    gdf_1point = gdf_1point.to_crs(4326)
    source_noise_db = 90
    target_noise_db = 40
    reflection_n = 3
    geometric_mean_freq_hz = 2000
    result = simulate_noise(
        source_points=gdf_1point,
        obstacles=buildings_data,
        source_noise_db=source_noise_db,
        geometric_mean_freq_hz=2000,
        standart_absorb_ratio=0.05,
        trees=trees_data,
        tree_resolution=4,
        air_temperature=20,
        target_noise_db=target_noise_db,
        db_sim_step=1,
        reflection_n=reflection_n,
        dead_area_r=5,
    )

    assert isinstance(result, gpd.GeoDataFrame)
    assert not result.empty
    assert "geometry" in result.columns
    assert "noise_level" in result.columns
    assert result["noise_level"].max() <= source_noise_db
    assert result["noise_level"].min() >= target_noise_db
    assert result.crs == gdf_1point.crs

    plot_simulation_result(
        result_gdf=result,
        source_points=gdf_1point,
        buildings=buildings_data,
        trees=trees_data,
        source_db_range=(source_noise_db, target_noise_db),
        reflection_n=reflection_n,
        frequency_desc=geometric_mean_freq_hz,
        output_filename="noise_simulation_1point",
    )


def test_multiple_sources(buildings_data, trees_data):
    p1 = Point(30.27060176, 59.93546846)
    p2 = Point(30.27303864, 59.9362777)
    p3 = Point(30.26804078, 59.93474246)
    gdf = gpd.GeoDataFrame(
        {
            "source_noise_db": [85, 90, 95],
            "geometric_mean_freq_hz": [500, 1000, 2000],
            "geometry": [p1, p2, p3],
        },
        crs=4326,
    )

    target_noise_db = 50
    reflection_n = 1
    result = simulate_noise(
        source_points=gdf,
        obstacles=buildings_data,
        standart_absorb_ratio=0.05,
        trees=trees_data,
        tree_resolution=1,
        air_temperature=20,
        target_noise_db=target_noise_db,
        db_sim_step=1,
        reflection_n=reflection_n,
        dead_area_r=5,
    )

    assert isinstance(result, gpd.GeoDataFrame)
    assert not result.empty
    assert "geometry" in result.columns
    assert "noise_level" in result.columns
    assert result["noise_level"].max() <= gdf["source_noise_db"].max()
    assert result["noise_level"].min() >= 40
    assert result.crs == gdf.crs

    plot_simulation_result(
        result_gdf=result,
        source_points=gdf,
        buildings=buildings_data,
        trees=trees_data,
        source_db_range=(gdf["source_noise_db"].max(), target_noise_db),
        reflection_n=reflection_n,
        frequency_desc="Mixed",
        output_filename="noise_simulation_3points",
    )


def test_wrong_step(gdf_1point, buildings_data):
    gdf_1point = gdf_1point.to_crs(4326)
    with pytest.raises(InvalidStepError) as _:
        _ = simulate_noise(
            source_points=gdf_1point,
            obstacles=buildings_data,
            source_noise_db=90,
            geometric_mean_freq_hz=2000,
            db_sim_step=4,
        )


def test_wrong_db_value(gdf_1point, buildings_data):
    gdf_1point = gdf_1point.to_crs(4326)
    with pytest.raises(ValueError) as _:
        _ = simulate_noise(
            source_points=gdf_1point,
            obstacles=buildings_data,
            source_noise_db=350,
            geometric_mean_freq_hz=2000,
            db_sim_step=4,
        )


def test_out_of_range_values(gdf_1point, buildings_data):
    out_of_range_hz = 10
    out_of_range_temperature = 40
    in_middle_hz = 1500
    in_middle_temperature = 15
    res = get_air_resist_ratio(out_of_range_temperature, out_of_range_hz, True)
    logger.info(f"Out of range result: {res}")
    res = get_air_resist_ratio(in_middle_temperature, in_middle_hz, True)
    logger.info(f"Between values result: {res}")
    res = get_air_resist_ratio(in_middle_temperature, 2000, True)
    logger.info(f"Between values result: {res}")
    res = get_air_resist_ratio(10, in_middle_hz, True)
    logger.info(f"Between values result: {res}")


def plot_simulation_result(
    result_gdf,
    source_points,
    buildings,
    trees,
    source_db_range,
    reflection_n,
    frequency_desc,
    output_filename,
):

    source_db, target_db = source_db_range
    local_crs = result_gdf.estimate_utm_crs()
    fig, ax = plt.subplots(figsize=(12, 10))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    minx, miny, maxx, maxy = result_gdf.to_crs(local_crs).total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    result_gdf.to_crs(local_crs).plot(
        ax=ax,
        column="noise_level",
        cmap="plasma",
        legend=True,
        alpha=0.8,
        edgecolor="white",
        linewidth=0.1,
        vmin=target_db,
        vmax=source_db,
        legend_kwds={"label": "Noise level, Decibels", "shrink": 0.5},
    )

    buildings.to_crs(local_crs).plot(ax=ax, facecolor="gray", edgecolor="black", linewidth=0.5, label="Buildings")
    trees.to_crs(local_crs).plot(ax=ax, edgecolor="green", facecolor="none", linewidth=1.5, label="Trees")
    source_points.to_crs(local_crs).plot(ax=ax, color="red", markersize=10, label="Noise sources")

    ax.set_title(
        f"Noise propagation {source_db}dB -> {target_db}dB\n"
        f"Frequency: {frequency_desc}, Reflection count: {reflection_n}, Temperature: 20°C"
    )
    ax.legend()
    ax.set_axis_off()

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()
