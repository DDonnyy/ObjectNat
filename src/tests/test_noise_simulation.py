import os

import geopandas as gpd
import pytest
from matplotlib import pyplot as plt

from objectnat import config, simulate_noise
from objectnat.methods.noise import InvalidStepError
from objectnat.methods.noise.noise_reduce import get_air_resist_ratio
from tests.conftest import output_dir

logger = config.logger


def test_basic_functionality(gdf_1point, buildings_data, trees_data):
    gdf_1point = gdf_1point.to_crs(4326)
    source_noise_db = 90
    target_noise_db = 40
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
        reflection_n=3,
        dead_area_r=5,
    )

    assert isinstance(result, gpd.GeoDataFrame)
    assert not result.empty
    assert "geometry" in result.columns
    assert "noise_level" in result.columns
    assert result["noise_level"].max() <= source_noise_db
    assert result["noise_level"].min() >= target_noise_db
    assert result.crs == gdf_1point.crs

    local_crs = result.estimate_utm_crs()
    fig, ax = plt.subplots(figsize=(12, 10))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    minx, miny, maxx, maxy = result.to_crs(local_crs).total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    result.to_crs(local_crs).plot(
        ax=ax,
        column="noise_level",
        cmap="plasma",
        legend=True,
        alpha=1,
        edgecolor="white",
        linewidth=0.1,
        vmin=target_noise_db,
        vmax=source_noise_db,
        label="Noise level (dB)",
    )

    buildings_data.to_crs(local_crs).plot(ax=ax, facecolor="gray", edgecolor="black", linewidth=0.5, label="Buildings")
    trees_data.to_crs(local_crs).plot(ax=ax, edgecolor="green", facecolor="none", linewidth=2, label="Trees")
    gdf_1point.to_crs(local_crs).plot(ax=ax, color="red", markersize=10, label="Noise source")

    ax.set_title(
        f"Noise propagation {source_noise_db}dB -> {target_noise_db}dB\n"
        f"Frequency: 2000Hz, Reflection count: 3, Temperature: 20C"
    )
    ax.legend()
    ax.set_axis_off()

    output_path = os.path.join(output_dir, "noise_simulation_test_result.png")
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()


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
