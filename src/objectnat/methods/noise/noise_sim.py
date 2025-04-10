import concurrent.futures
import math
import multiprocessing
import time

import geopandas as gpd
import pandas as pd
from shapely import GEOSException
from shapely.geometry import GeometryCollection, MultiPolygon, Point, Polygon
from shapely.ops import polygonize, unary_union
from tqdm import tqdm

from objectnat import config
from objectnat.methods.noise.noise_exceptions import InvalidStepError
from objectnat.methods.noise.noise_reduce import dist_to_target_db, green_noise_reduce_db
from objectnat.methods.utils.geom_utils import (
    gdf_to_circle_zones_from_point,
    get_point_from_a_thorough_b,
    polygons_to_multilinestring,
)
from objectnat.methods.visibility.visibility_analysis import get_visibility_accurate

logger = config.logger


def simulate_noise(
    source_points: gpd.GeoDataFrame, obstacles: gpd.GeoDataFrame, source_noise_db, geometric_mean_freq_hz, **kwargs
):
    """
    Simulates noise propagation from a set of source points considering obstacles, trees, and environmental factors.

    Args:
        source_points (gpd.GeoDataFrame): A GeoDataFrame containing one or more points representing the noise sources.
            A separate simulation will be run for each point.
        obstacles (gpd.GeoDataFrame): A GeoDataFrame representing obstacles in the environment. If a column with
            sound absorption coefficients is present, its name should be provided in the `absorb_ratio_column` argument.
            Missing values will be filled with the `standart_absorb_ratio`.
        source_noise_db (float): The noise level of the point source in decibels (dB). Decibels are logarithmic units
            used to measure sound intensity. A value of 20 dB represents a barely audible whisper, while 140 dB
            is comparable to the noise of jet engines.
        geometric_mean_freq_hz (float): The geometric mean frequency of the sound (in Hz). This parameter influences
            the sound wave's propagation and scattering in the presence of trees. Lower frequencies travel longer
            distances than higher frequencies. It's recommended to use values between 63 Hz and 8000 Hz; values outside
            this range will be clamped to the nearest boundary for the sound absorption coefficient calculation.

    Optional kwargs:
        absorb_ratio_column (str, optional): The name of the column in the `obstacles` GeoDataFrame that contains the
            sound absorption coefficients for each obstacle. Default is None. If not specified, all obstacles will have
            the `standart_absorb_ratio`.
        standart_absorb_ratio (float, optional): The default sound absorption coefficient to use for obstacles without
            specified values in the `absorb_ratio_column`. Default is 0.05, which is a typical value for concrete walls.
        trees (gpd.GeoDataFrame, optional): A GeoDataFrame containing trees or dense vegetation along the sound wave's
            path. Trees will scatter and absorb sound waves.
        tree_resolution (int, optional): A resolution parameter for simulating tree interactions with sound waves.
            Recommended values are between 2 and 16, with higher values providing more accurate simulation results.
        air_temperature (float, optional): The air temperature in degrees Celsius. The recommended range is from 0 to
            30 degrees Celsius, as temperatures outside this range will be clipped. Temperature affects the sound
            propagation in the air.
        target_noise_db (float, optional): The target noise level (in dB) for the simulation. Default is 40 dB.
            Lower values may not be relevant for further analysis, as they are near the threshold of human hearing.
        db_sim_step (float, optional): The step size in decibels for the noise simulation. Default is 1. For more
            precise analysis, this can be adjusted. If the difference between `source_noise_db` and `target_noise_db`
            is not divisible by the step size, the function will raise an error.
        reflection_n (int, optional): The maximum number of reflections (bounces) to simulate for each sound wave.
            Recommended values are between 1 and 3. Larger values will result in longer simulation times.
        dead_area_r (float, optional): A debugging parameter that defines the radius of the "dead zone" for reflections.
            Points within this area will not generate reflections. This is useful to prevent the algorithm from getting
            stuck in corners or along building walls.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the noise simulation results, including noise levels and geometries
        of the affected areas. Each point's simulation results will be merged into a single GeoDataFrame.
    """
    # Obstacles args
    absorb_ratio_column = kwargs.get("absorb_ratio_column", None)
    standart_absorb_ratio = kwargs.get("standart_absorb_ratio", 0.05)

    # Trees args
    trees = kwargs.get("trees", None)
    tree_res = kwargs.get("tree_resolution", 4)

    # Simulation conditions
    air_temperature = kwargs.get("air_temperature", 20)
    target_noise_db = kwargs.get("target_noise_db", 40)

    # Simulation params
    db_sim_step = kwargs.get("db_sim_step", 1)
    reflection_n = kwargs.get("reflection_n", 3)
    dead_area_r = kwargs.get("dead_area_r", 5)

    original_crs = source_points.crs

    div_ = (source_noise_db - target_noise_db) % db_sim_step
    if div_ != 0:
        raise InvalidStepError(source_noise_db, target_noise_db, db_sim_step, div_)
    # Choosing crs and simplifying obs if any
    source_points = source_points.copy()
    if len(obstacles) > 0:
        obstacles = obstacles.copy()
        obstacles.geometry = obstacles.geometry.simplify(tolerance=1)
        local_crs = obstacles.estimate_utm_crs()
        obstacles.to_crs(local_crs, inplace=True)
        source_points.to_crs(local_crs, inplace=True)
    else:
        local_crs = source_points.estimate_utm_crs()
        source_points.to_crs(local_crs, inplace=True)
        source_points.reset_index(drop=True)
        source_points.geometry = source_points.centroid

    # Simplifying trees
    if trees is not None:
        trees = trees.copy()
        trees.to_crs(local_crs, inplace=True)
        trees.geometry = trees.geometry.simplify(tolerance=1)
    else:
        trees = gpd.GeoDataFrame()

    if absorb_ratio_column is None:
        obstacles["absorb_ratio"] = standart_absorb_ratio
    else:
        obstacles["absorb_ratio"] = obstacles[absorb_ratio_column]
        obstacles["absorb_ratio"] = obstacles["absorb_ratio"].fillna(standart_absorb_ratio)
    obstacles = obstacles[["absorb_ratio", "geometry"]]

    logger.info(
        dist_to_target_db(
            source_noise_db,
            target_noise_db,
            geometric_mean_freq_hz,
            air_temperature,
            return_desc=True,
            check_temp_freq=True,
        )
    )
    # calculating layer dist and db values
    dist_db = [(0, source_noise_db)]
    cur_db = source_noise_db - db_sim_step
    while cur_db != target_noise_db - db_sim_step:
        max_dist = dist_to_target_db(source_noise_db, cur_db, geometric_mean_freq_hz, air_temperature)
        dist_db.append((max_dist, cur_db))
        cur_db = cur_db - db_sim_step

    # creating initial task and simulating for each point
    all_p_res = []
    for ind, row in source_points.iterrows():
        logger.info(f"Started simulation for point {ind+1} / {len(source_points)}")
        source_point = row.geometry
        task_queue = multiprocessing.Queue()
        args = (source_point, obstacles, trees, 0, 0, dist_db)
        kwargs = {
            "reflection_n": reflection_n,
            "geometric_mean_freq_hz": geometric_mean_freq_hz,
            "tree_res": tree_res,
            "min_db": target_noise_db,
        }
        task_queue.put((_noise_from_point_task, args, kwargs))

        noise_gdf = _parallel_split_queue(
            task_queue, dead_area=source_point.buffer(dead_area_r, resolution=2), dead_area_r=dead_area_r
        )

        noise_gdf = gpd.GeoDataFrame(pd.concat(noise_gdf, ignore_index=True), crs=local_crs)
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
        sim_result["source_point_ind"] = ind
        all_p_res.append(sim_result)

    return gpd.GeoDataFrame(pd.concat(all_p_res, ignore_index=True), crs=local_crs).to_crs(original_crs)


def _noise_from_point_task(task, **kwargs) -> tuple[gpd.GeoDataFrame, list[tuple] | None]:  # pragma: no cover
    # Unpacking task
    point_from, obstacles, trees_orig, passed_dist, deep, dist_db = task

    def donuts_dist_values(dist_db, passed_dist, max_view_dist):
        new_dist_db = dist_db + [(passed_dist, None), (max_view_dist + passed_dist, None)]
        new_dist_db = sorted(new_dist_db, key=lambda x: x[0])
        start = None
        end = None
        for i, (dist, db) in enumerate(new_dist_db[:-1]):
            if db is None:
                if start is None:
                    new_dist_db[i] = (dist, new_dist_db[i - 1][1])
                    start = i
                else:
                    new_dist_db[i] = (dist, new_dist_db[i + 1][1])
                    end = i + 1
                    break
        return [(dist - passed_dist, db) for dist, db in new_dist_db[start:end]]

    max_dist = max(dist_db, key=lambda x: x[0])[0]
    min_db = kwargs.get("min_db")
    reflection_n = kwargs.get("reflection_n")
    geometric_mean_freq_hz = kwargs.get("geometric_mean_freq_hz")
    tree_res = kwargs.get("tree_res")
    local_crs = obstacles.crs
    dist = round(max_dist - passed_dist, 1)

    obstacles = obstacles[obstacles.intersects(point_from.buffer(dist, resolution=8))]

    if len(obstacles) == 0:
        obstacles_union = Polygon()
    else:
        obstacles_union = obstacles.union_all()

    vis_poly, max_view_dist = get_visibility_accurate(point_from, obstacles, dist, return_max_view_dist=True)

    donuts_dist_values = donuts_dist_values(dist_db, passed_dist, max_view_dist)

    allowed_geom_types = ["MultiPolygon", "Polygon"]

    # Trees noise reduce
    reduce_polygons = []
    if len(trees_orig) > 0:
        trees_orig = trees_orig[trees_orig.intersects(point_from.buffer(dist, resolution=8))]
        if len(trees_orig) > 0:
            try:
                trees = gdf_to_circle_zones_from_point(trees_orig, point_from, dist, resolution=tree_res)
                trees = trees.clip(vis_poly, keep_geom_type=True).explode(index_parts=False)
            except TypeError:
                trees = gpd.GeoDataFrame()

            for _, row in trees.iterrows():
                tree_geom = row.geometry
                if tree_geom.area < 1:
                    continue
                dist_to_centroid = tree_geom.centroid.distance(point_from)

                points_with_angle = [
                    (
                        Point(pt),
                        round(abs(math.atan2(pt[1] - point_from.y, pt[0] - point_from.x)), 5),
                        Point(pt).distance(point_from),
                    )
                    for pt in tree_geom.exterior.coords
                ]

                p0_1 = max(points_with_angle, key=lambda x: (x[1], x[2]))
                p0_2 = min(points_with_angle, key=lambda x: (x[1], -x[2]))
                delta_angle = 2 * math.pi + p0_1[1] - p0_2[1]
                if delta_angle > math.pi:
                    delta_angle = 2 * math.pi - delta_angle

                a = math.sqrt((dist**2) * (1 + (math.tan(delta_angle / 2) ** 2)))
                p1 = get_point_from_a_thorough_b(point_from, p0_1[0], a)
                p2 = get_point_from_a_thorough_b(point_from, p0_2[0], a)
                red_polygon = unary_union([Polygon([p0_1[0], p1, p2, p0_2[0]]).intersection(vis_poly), tree_geom])
                if isinstance(red_polygon, GeometryCollection):
                    red_polygon = max(
                        ((poly, poly.area) for poly in red_polygon.geoms if isinstance(poly, (MultiPolygon, Polygon))),
                        key=lambda x: x[1],
                    )[0]
                if isinstance(red_polygon, MultiPolygon):
                    red_polygon = red_polygon.buffer(0.1, resolution=1).buffer(-0.1, resolution=1)
                if isinstance(red_polygon, MultiPolygon):
                    red_polygon = max(((poly, poly.area) for poly in red_polygon.geoms), key=lambda x: x[1])[0]
                if isinstance(red_polygon, Polygon) and not red_polygon.is_empty:
                    red_polygon = Polygon(red_polygon.exterior)
                    r_tree_new = round(
                        tree_geom.area / (2 * dist_to_centroid * math.sin(abs(p0_1[1] - p0_2[1]) / 2)), 2
                    )

                    noise_reduce = int(round(green_noise_reduce_db(geometric_mean_freq_hz, r_tree_new)))
                    reduce_polygons.append((red_polygon, noise_reduce))

    # Generating donuts - db values
    donuts = []
    don_values = []
    to_cut_off = point_from
    for i in range(len(donuts_dist_values[:-1])):
        cur_buffer = point_from.buffer(donuts_dist_values[i + 1][0])
        donuts.append(cur_buffer.difference(to_cut_off))
        don_values.append(donuts_dist_values[i][1])
        to_cut_off = cur_buffer

    noise_from_point = (
        gpd.GeoDataFrame(geometry=donuts, data={"noise_level": don_values}, crs=local_crs)
        .clip(vis_poly, keep_geom_type=True)
        .explode(ignore_index=True)
    )

    # intersect noise poly with noise reduce
    if len(reduce_polygons) > 0:
        reduce_polygons = gpd.GeoDataFrame(
            reduce_polygons, columns=["geometry", "reduce"], geometry="geometry", crs=local_crs
        )

        all_lines = (
            reduce_polygons.geometry.apply(polygons_to_multilinestring).tolist()
            + noise_from_point.geometry.apply(polygons_to_multilinestring).tolist()
        )

        cutted_polygons = gpd.GeoDataFrame(geometry=list(polygonize(unary_union(all_lines))), crs=local_crs)

        cutted_polygons_points = cutted_polygons.copy()
        cutted_polygons_points.geometry = cutted_polygons.representative_point()

        joined = (
            cutted_polygons_points.sjoin(noise_from_point, predicate="within", how="left")
            .drop(columns="index_right")
            .sjoin(reduce_polygons, predicate="within", how="left")
            .drop(columns="index_right")
        )
        joined.geometry = cutted_polygons.geometry
        joined = (
            joined.reset_index().groupby("index").agg({"geometry": "first", "reduce": "sum", "noise_level": "first"})
        )
        joined = gpd.GeoDataFrame(joined, geometry="geometry", crs=local_crs)
        noise_from_point = joined.copy()

        noise_from_point = noise_from_point.dropna(subset=["noise_level"])

        noise_from_point["reduce"] = noise_from_point["reduce"].fillna(0)
        noise_from_point["noise_level"] = noise_from_point["noise_level"] - noise_from_point["reduce"]
    else:
        noise_from_point["reduce"] = 0
    noise_from_point = noise_from_point[noise_from_point.geom_type.isin(allowed_geom_types)]
    noise_from_point = noise_from_point[noise_from_point["noise_level"] >= min_db]
    if deep == reflection_n:
        return noise_from_point, None

    if isinstance(vis_poly, Polygon):
        vis_poly_points = [Point(coords) for coords in vis_poly.exterior.coords]
    else:
        vis_poly_points = [Point(coords) for geom in vis_poly.geoms for coords in geom.exterior.coords]
    vis_poly_points = gpd.GeoDataFrame(geometry=vis_poly_points, crs=local_crs)

    # Generating reflection points
    vis_poly_points["point"] = vis_poly_points["geometry"].copy()
    vis_poly_points.geometry = vis_poly_points.geometry.buffer(1, resolution=1)
    vis_poly_points = vis_poly_points.sjoin(obstacles, predicate="intersects").drop(columns="index_right")
    vis_poly_points = vis_poly_points[~vis_poly_points.index.duplicated(keep="first")]
    vis_poly_points.dropna(subset=["absorb_ratio"], inplace=True)
    nearby_poly = point_from.buffer(1.1, resolution=2)
    try:
        vis_poly_points.geometry = (
            vis_poly_points.difference(vis_poly).difference(obstacles_union).difference(nearby_poly)
        )
    except GEOSException:
        return noise_from_point, None
    vis_poly_points = vis_poly_points[~vis_poly_points.is_empty]
    vis_poly_points = vis_poly_points[vis_poly_points.area >= 0.01]
    vis_poly_points["geometry"] = vis_poly_points["point"]
    vis_poly_points["dist"] = vis_poly_points.distance(point_from)
    vis_poly_points = vis_poly_points[vis_poly_points["dist"] < max_dist - 5]
    vis_poly_points = vis_poly_points.sjoin(noise_from_point, predicate="intersects", how="left")

    if len(vis_poly_points) == 0:
        return noise_from_point, None

    new_obs = pd.concat([obstacles, gpd.GeoDataFrame(geometry=[vis_poly], crs=local_crs)], ignore_index=True)

    # Creating new reflection tasks
    new_tasks = []
    for _, loc in vis_poly_points.iterrows():
        if not isinstance(loc.geometry, Point):
            continue
        new_passed_dist = round(loc.dist + passed_dist, 2)
        dist_last = max_dist - new_passed_dist
        if dist_last > 1:
            db_change = loc["reduce"]
            dist_change = loc["absorb_ratio"] * dist_last
            new_dist_db = [(dist - dist_change, db - db_change) for dist, db in dist_db]
            task_obs = new_obs.copy()
            task_obs.geometry = task_obs.difference(loc.geometry.buffer(1, resolution=1))
            new_tasks.append(
                (
                    _noise_from_point_task,
                    (loc.geometry, task_obs, trees_orig, new_passed_dist, deep + 1, new_dist_db),
                    kwargs,
                )
            )

    return noise_from_point, new_tasks


def _parallel_split_queue(task_queue: multiprocessing.Queue, dead_area: Polygon, dead_area_r: int):
    results = []
    total_tasks = task_queue.qsize()

    with tqdm(total=total_tasks, desc="Simulating noise") as pbar:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            future_to_task = {}
            while True:
                while not task_queue.empty() and len(future_to_task) < executor._max_workers:
                    func, task, kwargs = task_queue.get_nowait()
                    future = executor.submit(func, task, **kwargs)
                    future_to_task[future] = task

                done, _ = concurrent.futures.wait(future_to_task.keys(), return_when=concurrent.futures.FIRST_COMPLETED)

                for future in done:
                    future_to_task.pop(future)
                    result, new_tasks = future.result()
                    if new_tasks:
                        new_tasks_n = 0
                        new_dead_area_points = [dead_area]
                        for func, new_task, kwargs in new_tasks:
                            if not dead_area.covers(new_task[0]):
                                new_tasks_n = new_tasks_n + 1
                                task_queue.put((func, new_task, kwargs))
                                new_dead_area_points.append(new_task[0].buffer(dead_area_r, resolution=2))

                        dead_area = unary_union(new_dead_area_points)
                        total_tasks += new_tasks_n
                        pbar.total = total_tasks
                        pbar.refresh()
                    results.append(result)
                    pbar.update(1)
                time.sleep(0.01)
                if not future_to_task and task_queue.empty():
                    break

    return results
