import os

from matplotlib import pyplot as plt

from objectnat import get_clusters_polygon
from tests.conftest import output_dir


def test_and_visualize_clusters(buildings_data):
    min_dist = 70
    min_points = 2
    buildings_data = buildings_data.to_crs(4326)
    clusters, buildings_clustered = get_clusters_polygon(buildings_data, min_dist=min_dist, min_point=min_points)

    assert clusters.crs == buildings_data.crs
    assert buildings_clustered.crs == buildings_data.crs

    clusters = clusters[~clusters["outlier"]]
    buildings_clustered = buildings_clustered[~buildings_clustered["outlier"]]

    fig, ax = plt.subplots(figsize=(12, 10))

    local_crs = buildings_data.estimate_utm_crs()
    buildings_data = buildings_data.to_crs(local_crs)
    clusters.to_crs(local_crs, inplace=True)
    buildings_clustered.to_crs(local_crs, inplace=True)

    minx, miny, maxx, maxy = buildings_data.total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    clusters.plot(ax=ax, column="cluster", cmap="prism", alpha=0.4, edgecolor="black", linewidth=1, categorical=True)

    buildings_clustered.plot(
        ax=ax,
        column="cluster",
        cmap="prism",
        categorical=True,
        markersize=20,
        alpha=0.8,
    )

    ax.set_title(f"Building clusterization\n" f"Min distance to cluster: {min_dist}м, method: HDBSCAN")
    ax.set_axis_off()
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "building_clusters.png")
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()
