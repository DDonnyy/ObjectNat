# simplified version
import geopandas as gpd
import pandas as pd
from shapely.ops import polygonize, unary_union

from objectnat.methods.noise.noise_reduce import dist_to_target_db
from objectnat.methods.utils.geom_utils import (
    distribute_points_on_linestrings,
    distribute_points_on_polygons,
    polygons_to_multilinestring,
)
from objectnat.methods.visibility.visibility_analysis import get_visibility_accurate

MAX_DB_VALUE = 194


def calculate_simplified_noise_frame(
    noise_sources: gpd.GeoDataFrame, obstacles: gpd.GeoDataFrame, air_temperature, **kwargs
) -> gpd.GeoDataFrame:
    target_noise_db = kwargs.get("target_noise_db", 40)
    db_sim_step = kwargs.get("db_sim_step", 5)
    linestring_point_radius = kwargs.get("linestring_point_radius", 20)
    polygon_point_radius = kwargs.get("polygon_point_radius", 10)

    required_columns = ["source_noise_db", "geometric_mean_freq_hz"]
    for col in required_columns:
        if col not in noise_sources.columns:
            raise ValueError(f"'{col}' column is missing in provided GeoDataFrame")
        if noise_sources[col].isnull().any():
            raise ValueError(f"Column '{col}' contains missing (NaN) values")
    if (noise_sources["source_noise_db"] > MAX_DB_VALUE).any():
        raise ValueError(
            f"One or more values in 'source_noise_db' column exceed the physical limit of {MAX_DB_VALUE} dB."
        )
    original_crs = noise_sources.crs
    if len(obstacles) > 0:
        obstacles = obstacles.copy()
        obstacles.geometry = obstacles.geometry.simplify(tolerance=1)
        local_crs = obstacles.estimate_utm_crs()
        obstacles.to_crs(local_crs, inplace=True)
        noise_sources.to_crs(local_crs, inplace=True)
    else:
        local_crs = noise_sources.estimate_utm_crs()
        noise_sources.to_crs(local_crs, inplace=True)
        noise_sources.reset_index(drop=True)

    noise_sources = noise_sources.explode(ignore_index=True)
    noise_sources["geom_type"] = noise_sources.geom_type

    grouped_sources = noise_sources.groupby(by=["source_noise_db", "geometric_mean_freq_hz", "geom_type"])

    frame_result = []

    for (source_db, freq_hz, geom_type), group_gdf in grouped_sources:
        # calculating layer dist and db values
        dist_db = [(0, source_db)]
        cur_db = source_db - db_sim_step
        max_dist = 0
        while cur_db > target_noise_db - db_sim_step:
            if cur_db - db_sim_step < target_noise_db:
                cur_db = target_noise_db
            max_dist = dist_to_target_db(source_db, cur_db, freq_hz, air_temperature)
            dist_db.append((max_dist, cur_db))
            cur_db -= db_sim_step

        if geom_type == "Point":
            for _, row in group_gdf.iterrows():
                point_from = row.geometry
                point_buffer = point_from.buffer(max_dist, resolution=16)
                local_obstacles = obstacles[obstacles.intersects(point_buffer)]
                vis_poly = get_visibility_accurate(point_from, obstacles=local_obstacles, view_distance=max_dist)
                noise_from_feature = _eval_donuts_gdf(point_from, dist_db, local_crs, vis_poly)
                frame_result.append(noise_from_feature)

        elif geom_type == "LineString":
            layer_points = distribute_points_on_linestrings(group_gdf, radius=linestring_point_radius, lloyd_relax_n=1)
            noise_from_feature = _process_lines_or_polygons(
                group_gdf, max_dist, obstacles, layer_points, dist_db, local_crs
            )
            frame_result.append(noise_from_feature)
        elif geom_type == "Polygon":
            group_gdf.geometry = group_gdf.buffer(0.1, resolution=1)
            layer_points = distribute_points_on_polygons(
                group_gdf, only_exterior=False, radius=polygon_point_radius, lloyd_relax_n=1
            )
            noise_from_feature = _process_lines_or_polygons(
                group_gdf, max_dist, obstacles, layer_points, dist_db, local_crs
            )
            frame_result.append(noise_from_feature)
        else:
            pass

    noise_gdf = gpd.GeoDataFrame(pd.concat(frame_result, ignore_index=True), crs=local_crs)
    polygons = gpd.GeoDataFrame(
        geometry=list(polygonize(noise_gdf.geometry.apply(polygons_to_multilinestring).union_all())), crs=local_crs
    )
    polygons_points = polygons.copy()
    polygons_points.geometry = polygons.representative_point()
    sim_result = polygons_points.sjoin(noise_gdf, predicate="within").reset_index()
    sim_result = sim_result.groupby("index").agg({"noise_level": "max"})
    sim_result["geometry"] = polygons
    sim_result = (
        gpd.GeoDataFrame(sim_result, geometry="geometry", crs=local_crs).dissolve(by="noise_level").reset_index()
    )

    return sim_result.to_crs(original_crs)


def _process_lines_or_polygons(group_gdf, max_dist, obstacles, layer_points, dist_db, local_crs) -> gpd.GeoDataFrame:
    features_vision_polys = []
    layer_buffer = group_gdf.buffer(max_dist, resolution=16).union_all()
    local_obstacles = obstacles[obstacles.intersects(layer_buffer)]
    for _, row in layer_points.iterrows():
        point_from = row.geometry
        vis_poly = get_visibility_accurate(point_from, obstacles=local_obstacles, view_distance=max_dist)
        features_vision_polys.append(vis_poly)
    features_vision_polys = unary_union(features_vision_polys)
    return _eval_donuts_gdf(group_gdf.union_all(), dist_db, local_crs, features_vision_polys)


def _eval_donuts_gdf(initial_geometry, dist_db, local_crs, clip_poly) -> gpd.GeoDataFrame:
    donuts = []
    don_values = []
    to_cut_off = initial_geometry
    for i in range(len(dist_db[:-1])):
        cur_buffer = initial_geometry.buffer(dist_db[i + 1][0])
        donuts.append(cur_buffer.difference(to_cut_off))
        don_values.append(dist_db[i][1])
        to_cut_off = cur_buffer
    noise_from_feature = (
        gpd.GeoDataFrame(geometry=donuts, data={"noise_level": don_values}, crs=local_crs)
        .clip(clip_poly, keep_geom_type=True)
        .explode(ignore_index=True)
    )
    return noise_from_feature
