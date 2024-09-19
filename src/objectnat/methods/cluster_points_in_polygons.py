from typing import Literal

import geopandas as gpd
import pandas as pd
from sklearn.cluster import DBSCAN, HDBSCAN

from objectnat import config

logger = config.logger


def _get_cluster(services_select, min_dist, min_point, method):
    services_coords = pd.DataFrame(
        {"x": services_select.geometry.representative_point().x, "y": services_select.geometry.representative_point().y}
    )
    if method == "DBSCAN":
        db = DBSCAN(eps=min_dist, min_samples=min_point).fit(services_coords.to_numpy())
    else:
        db = HDBSCAN(min_cluster_size=min_point, cluster_selection_epsilon=min_dist).fit(services_coords.to_numpy())
    services_select["cluster"] = db.labels_
    return services_select


def _get_service_ratio(loc):
    all_services = loc.shape[0]
    loc["service_code"] = loc["service_code"].astype(str)
    services_count = loc.groupby("service_code").size()
    return (services_count / all_services).round(2)


def get_clusters_polygon(
    points: gpd.GeoDataFrame,
    min_dist: float | int = 100,
    min_point: int = 5,
    method: Literal["DBSCAN", "HDBSCAN"] = "HDBSCAN",
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Generate cluster polygons for given points based on a specified minimum distance and minimum points per cluster.
    Optionally, calculate the relative ratio between types of services within the clusters.

    Parameters
    ----------
    points : gpd.GeoDataFrame
        GeoDataFrame containing the points to be clustered.
        Must include a 'service_code' column for service ratio calculations.
    min_dist : float | int, optional
        Minimum distance between points to be considered part of the same cluster. Defaults to 100.
    min_point : int, optional
        Minimum number of points required to form a cluster. Defaults to 5.
    method : Literal["DBSCAN", "HDBSCAN"], optional
        The clustering method to use. Must be either "DBSCAN" or "HDBSCAN". Defaults to "HDBSCAN".

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        A tuple containing the clustered polygons GeoDataFrame and the original points GeoDataFrame with cluster labels.

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point

    >>> points = gpd.GeoDataFrame({
    ...     'geometry': [Point(0, 0), Point(1, 1), Point(2, 2)],
    ...     'service_code': [1, 1, 2]
    ... }, crs=4326)

    >>> clusters, services = get_clusters_polygon(points, min_dist=50, min_point=2)
    """
    if method not in ["DBSCAN", "HDBSCAN"]:
        raise ValueError("Method must be either 'DBSCAN' or 'HDBSCAN'")

    services_select = _get_cluster(points, min_dist, min_point, method)

    if "service_code" not in points.columns:
        logger.warning(
            "No 'service_code' column in provided GeoDataFrame, cluster polygons will be without relative ratio."
        )
        points["service_code"] = 1

    services_normal = services_select[services_select["cluster"] != -1]
    services_outlier = services_select[services_select["cluster"] == -1]

    if len(services_normal) > 0:
        cluster_service = services_normal.groupby("cluster", group_keys=True).apply(_get_service_ratio)
        if isinstance(cluster_service, pd.Series):
            cluster_service = cluster_service.unstack(level=1, fill_value=0)

        polygons_normal = services_normal.dissolve("cluster").concave_hull(ratio=0.7, allow_holes=False)
        df_clusters_normal = pd.concat([cluster_service, polygons_normal.rename("geometry")], axis=1)
        cluster_normal = df_clusters_normal.index.max()
    else:
        df_clusters_normal = None
        cluster_normal = 0

    if len(services_outlier) > 0:
        clusters_outlier = cluster_normal + 1
        new_clusters = list(range(clusters_outlier, clusters_outlier + len(services_outlier)))
        services_outlier.loc[:, "cluster"] = new_clusters

        cluster_service = services_outlier.groupby("cluster", group_keys=True).apply(_get_service_ratio)
        if isinstance(cluster_service, pd.Series):
            cluster_service = cluster_service.unstack(level=1, fill_value=0)

        df_clusters_outlier = cluster_service.join(services_outlier.set_index("cluster")["geometry"])
    else:
        services_outlier = None
        df_clusters_outlier = None

    df_clusters = pd.concat([df_clusters_normal, df_clusters_outlier]).fillna(0).set_geometry("geometry")
    df_clusters["geometry"] = df_clusters["geometry"].buffer(min_dist / 2)
    df_clusters = df_clusters.rename(columns={"index": "cluster_id"})

    services = pd.concat([services_normal, services_outlier])

    return df_clusters, services
