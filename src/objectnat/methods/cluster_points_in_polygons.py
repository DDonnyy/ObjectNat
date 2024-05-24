import contextlib
import os

import geopandas as gpd
import pandas as pd
import pca
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import DBSCAN


def _get_cluster(services_select, min_dist, min_point):
    services_coords = pd.DataFrame({"x": services_select.geometry.representative_point().x,
                                    "y": services_select.geometry.representative_point().y})
    #clusterization = linkage(services_coords.to_numpy(), method="ward")
    #services_select["cluster"] = fcluster(clusterization, t=condition_value, criterion=condition).astype(int)

    eps = min_dist
    min_samples = min_point

    db = DBSCAN(eps=eps, min_samples=min_samples).fit(services_coords.to_numpy())

    #services_select_2 = services_select.copy()
    services_select["cluster"] = db.labels_

    return services_select


#def _find_dense_groups(geom, n_std):
#    if len(geom) > 1:
#        X = pd.DataFrame({"x": geom.representative_point().x, "y": geom.representative_point().y})
#        X = X.to_numpy()
#
#        # supress pca lib output
#        with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
#            outlier = pca.spe_dmodx(X, n_std=n_std)[0]["y_bool_spe"]
#        return pd.Series(data=outlier.values, index=geom.index)
#    else:
#        return pd.Series(data=True, index=geom.index)


def _get_service_ratio(loc):
    all_services = loc["id"].count()
    loc['service_code'] = loc['service_code'].astype(str)
    services_count = loc.groupby("service_code")["id"].count()
    return (services_count / all_services).round(2)


def get_clusters_polygon(points: gpd.GeoDataFrame, min_dist=100, min_point=5):
    services_select = _get_cluster(points, min_dist, min_point)

    # Find outliers of clusters and exclude it
    #outlier = services_select.groupby("cluster", group_keys=False)["geometry"].apply(
    #    lambda x: _find_dense_groups(x, n_std)
    #)
    #services_outlier = services_select.loc[outlier]

    services_normal = services_select[services_select['cluster'] != -1]
    services_outlier = services_select[services_select['cluster'] == -1]
    if len(services_normal) > 0:
        cluster_service = services_normal.groupby(["cluster"], group_keys=True).apply(
            lambda x: _get_service_ratio(x)
        )
        if isinstance(cluster_service, pd.Series):
            cluster_service = cluster_service.unstack(level=1, fill_value=0)

        # Get MultiPoint from cluster Points and make polygon
        polygons_normal = services_normal.dissolve("cluster").convex_hull
        df_clusters_normal = pd.concat([cluster_service, polygons_normal.rename("geometry")], axis=1)
        cluster_normal = df_clusters_normal.index.max()

    # Select outliers
    if any(services_outlier):
        # Reindex clusters
        clusters_outlier = cluster_normal + 1
        new_clusters = [c for c in range(clusters_outlier, clusters_outlier + len(services_outlier))]
        services_outlier.loc[:, "cluster"] = new_clusters
        cluster_service = services_outlier.groupby(["cluster"], group_keys=True).apply(
            lambda x: _get_service_ratio(x)
        )
        if isinstance(cluster_service, pd.Series):
            cluster_service = cluster_service.unstack(level=1, fill_value=0)
        df_clusters_outlier = cluster_service.join(services_outlier.set_index("cluster")["geometry"])
    else:
        services_outlier = None
        df_clusters_outlier = None

    df_clusters = pd.concat([df_clusters_normal, df_clusters_outlier]).fillna(0).set_geometry("geometry")
    df_clusters = df_clusters.fillna(0).set_geometry("geometry")
    df_clusters["geometry"] = df_clusters["geometry"].buffer(50, join_style=3)
    df_clusters = df_clusters.rename(columns={"index": "cluster_id"})

    services = pd.concat([services_normal, services_outlier])

    return df_clusters, services
