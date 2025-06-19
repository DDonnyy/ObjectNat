import math
from multiprocessing import cpu_count

import geopandas as gpd
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from shapely import LineString, MultiPolygon, Point, Polygon
from shapely.ops import unary_union
from tqdm.contrib.concurrent import process_map

from objectnat import config
from objectnat.methods.utils.geom_utils import (
    combine_geometry,
    explode_linestring,
    get_point_from_a_thorough_b,
    point_side_of_line,
    polygons_to_multilinestring,
)
from objectnat.methods.utils.math_utils import min_max_normalization

logger = config.logger


def get_visibility_accurate(
    point_from: Point | gpd.GeoDataFrame, obstacles: gpd.GeoDataFrame, view_distance, return_max_view_dist=False
) -> Polygon | gpd.GeoDataFrame | tuple[Polygon | gpd.GeoDataFrame, float]:
    """
    Function to get accurate visibility from a given point to buildings within a given distance.

    Parameters
    ----------
    point_from : Point | gpd.GeoDataFrame
        The point or GeoDataFrame with 1 point from which the line of sight is drawn.
        If Point is provided it should be in the same crs as obstacles
    obstacles : gpd.GeoDataFrame
        A GeoDataFrame containing the geometry of the obstacles.
    view_distance : float
        The distance of view from the point.
    return_max_view_dist
        If True, the max view distance is returned with view polygon in tuple.

    Returns
    -------
    Polygon | gpd.GeoDataFrame | tuple[Polygon | gpd.GeoDataFrame, float]
        A polygon representing the area of visibility from the given point or polygon with max view distance.
        if point_from was a GeoDataFrame, return GeoDataFrame with one feature, else Polygon.

    Notes
    -----
    If a quick result is important, consider using the `get_visibility()` function instead.
    However, please note that `get_visibility()` may provide less accurate results.

    Examples
    --------
    >>> from objectnat import get_visibility_accurate
    >>> obstacles = gpd.read_parquet('examples_data/buildings.parquet')
    >>> point_from = gpd.GeoDataFrame(geometry=[Point(30.2312112, 59.9482336)], crs=4326)
    >>> result = get_visibility_accurate(point_from, obstacles, 500)
    """

    def find_furthest_point(point_from, view_polygon):
        try:
            res = round(max(Point(coords).distance(point_from) for coords in view_polygon.exterior.coords), 1)
        except Exception as e:
            print(view_polygon)
            raise e
        return res

    local_crs = None
    original_crs = None
    return_gdf = False
    if isinstance(point_from, gpd.GeoDataFrame):
        original_crs = point_from.crs
        return_gdf = True
        if len(obstacles) > 0:
            local_crs = obstacles.estimate_utm_crs()
        else:
            local_crs = point_from.estimate_utm_crs()
        obstacles = obstacles.to_crs(local_crs)
        point_from = point_from.to_crs(local_crs)
        if len(point_from) > 1:
            logger.warning(
                f"This method processes only single point. The GeoDataFrame contains {len(point_from)} points - "
                "only the first geometry will be used for isochrone calculation. "
            )
        point_from = point_from.iloc[0].geometry
    else:
        obstacles = obstacles.copy()
    if obstacles.contains(point_from).any():
        return Polygon()
    obstacles.reset_index(inplace=True, drop=True)
    point_buffer = point_from.buffer(view_distance, resolution=32)
    allowed_geom_types = ["MultiPolygon", "Polygon", "LineString", "MultiLineString"]
    obstacles = obstacles[obstacles.geom_type.isin(allowed_geom_types)]
    s = obstacles.intersects(point_buffer)
    obstacles_in_buffer = obstacles.loc[s[s].index].geometry

    buildings_lines_in_buffer = gpd.GeoSeries(
        pd.Series(
            obstacles_in_buffer.apply(polygons_to_multilinestring).explode(index_parts=False).apply(explode_linestring)
        ).explode()
    )

    buildings_lines_in_buffer = buildings_lines_in_buffer.loc[buildings_lines_in_buffer.intersects(point_buffer)]

    buildings_in_buffer_points = gpd.GeoSeries(
        [Point(line.coords[0]) for line in buildings_lines_in_buffer.geometry]
        + [Point(line.coords[-1]) for line in buildings_lines_in_buffer.geometry]
    )

    max_dist = max(view_distance, buildings_in_buffer_points.distance(point_from).max())
    polygons = []
    buildings_lines_in_buffer = gpd.GeoDataFrame(geometry=buildings_lines_in_buffer, crs=obstacles.crs).reset_index()
    logger.debug("Calculation vis polygon")
    while not buildings_lines_in_buffer.empty:
        gdf_sindex = buildings_lines_in_buffer.sindex
        # TODO check if 2 walls are nearest and use the widest angle between points
        nearest_wall_sind = gdf_sindex.nearest(point_from, return_all=False, max_distance=max_dist)
        nearest_wall = buildings_lines_in_buffer.loc[nearest_wall_sind[1]].iloc[0]
        wall_points = [Point(coords) for coords in nearest_wall.geometry.coords]

        # Calculate angles and sort by angle
        points_with_angle = sorted(
            [(pt, math.atan2(pt.y - point_from.y, pt.x - point_from.x)) for pt in wall_points], key=lambda x: x[1]
        )
        delta_angle = 2 * math.pi + points_with_angle[0][1] - points_with_angle[-1][1]
        if round(delta_angle, 10) == round(math.pi, 10):
            wall_b_centroid = obstacles_in_buffer.loc[nearest_wall["index"]].centroid
            p1 = get_point_from_a_thorough_b(point_from, points_with_angle[0][0], max_dist)
            p2 = get_point_from_a_thorough_b(point_from, points_with_angle[1][0], max_dist)
            polygon = LineString([p1, p2])
            polygon = polygon.buffer(
                distance=max_dist * point_side_of_line(polygon, wall_b_centroid), single_sided=True
            )
        else:
            if delta_angle > math.pi:
                delta_angle = 2 * math.pi - delta_angle
            a = math.sqrt((max_dist**2) * (1 + (math.tan(delta_angle / 2) ** 2)))
            p1 = get_point_from_a_thorough_b(point_from, points_with_angle[0][0], a)
            p2 = get_point_from_a_thorough_b(point_from, points_with_angle[-1][0], a)
            polygon = Polygon([points_with_angle[0][0], p1, p2, points_with_angle[1][0]])

        polygons.append(polygon)
        buildings_lines_in_buffer.drop(nearest_wall_sind[1], inplace=True)

        if not polygon.is_valid or polygon.area < 1:
            buildings_lines_in_buffer.reset_index(drop=True, inplace=True)
            continue

        lines_to_kick = buildings_lines_in_buffer.within(polygon)
        buildings_lines_in_buffer = buildings_lines_in_buffer.loc[~lines_to_kick]
        buildings_lines_in_buffer.reset_index(drop=True, inplace=True)
    logger.debug("Done calculating!")
    res = point_buffer.difference(unary_union(polygons + obstacles_in_buffer.to_list()))

    if isinstance(res, MultiPolygon):
        res = list(res.geoms)
        for polygon in res:
            if polygon.intersects(point_from):
                res = polygon
                break

    if return_gdf:
        res = gpd.GeoDataFrame(geometry=[res], crs=local_crs).to_crs(original_crs)

    if return_max_view_dist:
        return res, find_furthest_point(point_from, res)
    return res


def get_visibility(
    point_from: Point | gpd.GeoDataFrame, obstacles: gpd.GeoDataFrame, view_distance: float, resolution: int = 32
) -> Polygon | gpd.GeoDataFrame:
    """
    Function to get a quick estimate of visibility from a given point to buildings within a given distance.

    Parameters
    ----------
    point_from : Point | gpd.GeoDataFrame
        The point or GeoDataFrame with 1 point from which the line of sight is drawn.
        If Point is provided it should be in the same crs as obstacles
    obstacles : gpd.GeoDataFrame
        A GeoDataFrame containing the geometry of the buildings.
    view_distance : float
        The distance of view from the point.
    resolution: int
        Buffer resolution for more accuracy (may give result slower)

    Returns
    -------
    Polygon | gpd.GeoDataFrame
        A polygon representing the area of visibility from the given point.
        if point_from was a GeoDataFrame, return GeoDataFrame with one feature, else Polygon.

    Notes
    -----
    This function provides a quicker but less accurate result compared to `get_visibility_accurate()`.
    If accuracy is important, consider using `get_visibility_accurate()` instead.

    Examples
    --------
    >>> from objectnat import get_visibility
    >>> obstacles = gpd.read_parquet('examples_data/buildings.parquet')
    >>> point_from = gpd.GeoDataFrame(geometry=[Point(30.2312112, 59.9482336)], crs=4326)
    >>> result = get_visibility(point_from, obstacles, 500)
    """
    return_gdf = False
    if isinstance(point_from, gpd.GeoDataFrame):
        original_crs = point_from.crs
        return_gdf = True
        if len(obstacles) > 0:
            local_crs = obstacles.estimate_utm_crs()
        else:
            local_crs = point_from.estimate_utm_crs()
        obstacles = obstacles.to_crs(local_crs)
        point_from = point_from.to_crs(local_crs)
        if len(point_from) > 1:
            logger.warning(
                f"This method processes only single point. The GeoDataFrame contains {len(point_from)} points - "
                "only the first geometry will be used for isochrone calculation. "
            )
        point_from = point_from.iloc[0].geometry
    else:
        obstacles = obstacles.copy()
    point_buffer = point_from.buffer(view_distance, resolution=resolution)
    s = obstacles.intersects(point_buffer)
    buildings_in_buffer = obstacles.loc[s[s].index].reset_index(drop=True)
    buffer_exterior_ = list(point_buffer.exterior.coords)
    line_geometry = [LineString([point_from, ext]) for ext in buffer_exterior_]
    buffer_lines_gdf = gpd.GeoDataFrame(geometry=line_geometry)
    united_buildings = buildings_in_buffer.union_all()
    if united_buildings:
        splited_lines = buffer_lines_gdf["geometry"].apply(lambda x: x.difference(united_buildings))
    else:
        splited_lines = buffer_lines_gdf["geometry"]

    splited_lines_gdf = gpd.GeoDataFrame(geometry=splited_lines).explode(index_parts=True)
    splited_lines_list = []

    for _, v in splited_lines_gdf.groupby(level=0):
        splited_lines_list.append(v.iloc[0]["geometry"].coords[-1])
    circuit = Polygon(splited_lines_list)
    if united_buildings:
        circuit = circuit.difference(united_buildings)

    if return_gdf:
        circuit = gpd.GeoDataFrame(geometry=[circuit], crs=local_crs).to_crs(original_crs)
    return circuit


def get_visibilities_from_points(
    points: gpd.GeoDataFrame,
    obstacles: gpd.GeoDataFrame,
    view_distance: int,
    sectors_n=None,
    max_workers: int = cpu_count(),
) -> list[Polygon]:
    """
    Calculate visibility polygons from a set of points considering obstacles within a specified view distance.

    Parameters
    ----------
    points : gpd.GeoDataFrame
        GeoDataFrame containing the points from which visibility is calculated.
    obstacles : gpd.GeoDataFrame
        GeoDataFrame containing the obstacles that block visibility.
    view_distance : int
        The maximum distance from each point within which visibility is calculated.
    sectors_n : int, optional
        Number of sectors to divide the view into for more detailed visibility calculations. Defaults to None.
    max_workers: int, optional
        Maximum workers in multiproccesing, multipocessing.cpu_count() by default.

    Returns
    -------
    list[Polygon]
        A list of visibility polygons for each input point.

    Notes
    -----
    This function uses `get_visibility_accurate()` in multiprocessing way.

    """
    if points.crs != obstacles.crs:
        raise ValueError(f"CRS mismatch, points crs:{points.crs} != obstacles crs:{obstacles.crs}")
    if points.crs.is_geographic:
        logger.warning("Points crs is geographic, it may produce invalid results")
    # remove points inside polygons
    joined = gpd.sjoin(points, obstacles, how="left", predicate="intersects")
    points = joined[joined.index_right.isnull()]

    # remove unused obstacles
    points_view = points.geometry.buffer(view_distance).union_all()
    s = obstacles.intersects(points_view)
    buildings_in_buffer = obstacles.loc[s[s].index].reset_index(drop=True)

    buildings_in_buffer.geometry = buildings_in_buffer.geometry.apply(
        lambda geom: MultiPolygon([geom]) if isinstance(geom, Polygon) else geom
    )
    args = [(point, buildings_in_buffer, view_distance, sectors_n) for point in points.geometry]
    all_visions = process_map(
        _multiprocess_get_vis,
        args,
        chunksize=5,
        desc="Calculating Visibility Catchment Area from each Point, it might take a while for a "
        "big amount of points",
        max_workers=max_workers,
    )

    # could return sectorized visions if sectors_n is set
    return all_visions


def calculate_visibility_catchment_area(
    points: gpd.GeoDataFrame, obstacles: gpd.GeoDataFrame, view_distance: int | float, max_workers: int = cpu_count()
) -> gpd.GeoDataFrame:
    """
    Calculate visibility catchment areas for a large urban area based on given points and obstacles.
    This function is designed to work with at least 1000 points spaced 10-20 meters apart for optimal results.
    Points can be generated using a road graph.

    Parameters
    ----------
    points : gpd.GeoDataFrame
        GeoDataFrame containing the points from which visibility is calculated.
    obstacles : gpd.GeoDataFrame
        GeoDataFrame containing the obstacles that block visibility.
    view_distance : int
        The maximum distance from each point within which visibility is calculated.
    max_workers: int
        Maximum workers in multiproccesing, multipocessing.cpu_count() by default.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing the calculated visibility catchment areas.

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point, Polygon
    >>> points = gpd.read_file('points.shp')
    >>> obstacles = gpd.read_file('obstacles.shp')
    >>> view_distance = 1000

    >>> visibility_areas = calculate_visibility_catchment_area(points, obstacles, view_distance)
    >>> visibility_areas
    """

    def filter_geoms(x):
        if x.geom_type == "GeometryCollection":
            return MultiPolygon([y for y in x.geoms if y.geom_type in ["Polygon", "MultiPolygon"]])
        return x

    def calc_group_factor(x):
        # pylint: disable-next=redefined-outer-name,reimported,import-outside-toplevel
        import numpy as np

        return np.mean(x.new_ratio) * x.count_n

    def unary_union_groups(x):
        # pylint: disable-next=redefined-outer-name,reimported,import-outside-toplevel
        from shapely import MultiPolygon

        # pylint: disable-next=redefined-outer-name,reimported,import-outside-toplevel
        from shapely.ops import unary_union

        return unary_union(MultiPolygon(list(x["geometry"])).buffer(0))

    pandarallel.initialize(progress_bar=True, verbose=0)

    local_crs = obstacles.estimate_utm_crs()
    obstacles = obstacles.to_crs(local_crs)
    points = points.to_crs(local_crs)

    sectors_n = 12
    logger.info("Calculating Visibility Catchment Area from each point")
    all_visions_sectorized = get_visibilities_from_points(points, obstacles, view_distance, sectors_n, max_workers)
    all_visions_sectorized = gpd.GeoDataFrame(
        geometry=[item for sublist in all_visions_sectorized for item in sublist], crs=local_crs
    )
    logger.info("Calculating non-vision part...")
    all_visions_unary = all_visions_sectorized.union_all()
    convex = all_visions_unary.convex_hull
    dif = convex.difference(all_visions_unary)

    del convex, all_visions_unary

    buf_area = (math.pi * view_distance**2) / sectors_n
    all_visions_sectorized["ratio"] = all_visions_sectorized.area / buf_area
    all_visions_sectorized["ratio"] = min_max_normalization(
        all_visions_sectorized["ratio"].values, new_min=1, new_max=10
    )
    groups = all_visions_sectorized.sample(frac=1).groupby(all_visions_sectorized.index // 6000)
    groups = [group for _, group in groups]

    del all_visions_sectorized

    groups_result = process_map(
        _process_group,
        groups,
        desc="Counting intersections in each group...",
        max_workers=max_workers,
    )
    logger.info("Calculating all groups intersection...")
    all_in = combine_geometry(gpd.GeoDataFrame(data=pd.concat(groups_result), geometry="geometry", crs=local_crs))

    del groups_result

    all_in["count_n"] = all_in["index_right"].apply(len)

    logger.info("Calculating intersection's parameters")
    all_in["factor"] = all_in.parallel_apply(calc_group_factor, axis=1)
    threshold = all_in["factor"].quantile(0.3)
    all_in = all_in[all_in["factor"] > threshold]

    all_in["factor_normalized"] = np.round(
        min_max_normalization(np.sqrt(all_in["factor"].values), new_min=1, new_max=5)
    ).astype(int)
    logger.info("Calculating normalized groups geometry...")
    all_in = all_in.groupby("factor_normalized").parallel_apply(unary_union_groups).reset_index()
    all_in = gpd.GeoDataFrame(data=all_in.rename(columns={0: "geometry"}), geometry="geometry", crs=32636)

    all_in = all_in.explode(index_parts=True).reset_index(drop=True)
    all_in["area"] = all_in.area
    threshold = all_in["area"].quantile(0.9)
    all_in = all_in[all_in["area"] > threshold]
    all_in = all_in.groupby("factor_normalized").apply(unary_union_groups).reset_index()
    all_in = gpd.GeoDataFrame(data=all_in.rename(columns={0: "geometry"}), geometry="geometry", crs=32636)

    all_in.geometry = all_in.geometry.buffer(20).buffer(-20).difference(dif)

    all_in.sort_values(by="factor_normalized", ascending=False, inplace=True)
    all_in.reset_index(drop=True, inplace=True)
    logger.info("Smoothing normalized groups geometry...")
    for ind, row in all_in.iloc[:-1].iterrows():
        for ind2 in range(ind + 1, len(all_in)):
            current_geometry = all_in.at[ind2, "geometry"]
            all_in.at[ind2, "geometry"] = current_geometry.difference(row.geometry)
            all_in["geometry"] = all_in["geometry"].apply(filter_geoms)

    all_in = all_in.explode(index_parts=True)
    logger.info("Done!")
    return all_in


def _multiprocess_get_vis(args):  # pragma: no cover
    point, buildings, view_distance, sectors_n = args
    result = get_visibility_accurate(point, buildings, view_distance)

    if sectors_n is not None:
        sectors = []

        cx, cy = point.x, point.y

        angle_increment = 2 * math.pi / sectors_n
        view_distance = math.sqrt((view_distance**2) * (1 + (math.tan(angle_increment / 2) ** 2)))
        for i in range(sectors_n):
            angle1 = i * angle_increment
            angle2 = (i + 1) * angle_increment

            x1, y1 = cx + view_distance * math.cos(angle1), cy + view_distance * math.sin(angle1)
            x2, y2 = cx + view_distance * math.cos(angle2), cy + view_distance * math.sin(angle2)

            sector_triangle = Polygon([point, (x1, y1), (x2, y2)])
            sector = result.intersection(sector_triangle)

            if not sector.is_empty:
                sectors.append(sector)
        result = sectors
    return result


def _process_group(group):  # pragma: no cover
    geom = group
    combined_geometry = combine_geometry(geom)
    combined_geometry.drop(columns=["index", "index_right"], inplace=True)
    combined_geometry["count_n"] = combined_geometry["ratio"].apply(len)
    combined_geometry["new_ratio"] = combined_geometry.apply(
        lambda x: np.power(np.prod(x.ratio), 1 / x.count_n) * x.count_n, axis=1
    )

    threshold = combined_geometry["new_ratio"].quantile(0.25)
    combined_geometry = combined_geometry[combined_geometry["new_ratio"] > threshold]

    combined_geometry["new_ratio_normalized"] = min_max_normalization(
        combined_geometry["new_ratio"].values, new_min=1, new_max=10
    )

    combined_geometry["new_ratio_normalized"] = np.round(combined_geometry["new_ratio_normalized"]).astype(int)

    result_union = (
        combined_geometry.groupby("new_ratio_normalized")
        .agg({"geometry": lambda x: unary_union(MultiPolygon(list(x)).buffer(0))})
        .reset_index(drop=True)
    )
    result_union.set_geometry("geometry", inplace=True)
    result_union.set_crs(geom.crs, inplace=True)

    result_union = result_union.explode("geometry", index_parts=False).reset_index(drop=True)

    representative_points = combined_geometry.copy()
    representative_points["geometry"] = representative_points["geometry"].representative_point()

    joined = gpd.sjoin(result_union, representative_points, how="inner", predicate="contains").reset_index()
    joined = joined.groupby("index").agg({"geometry": "first", "new_ratio": lambda x: np.mean(list(x))})

    joined.set_geometry("geometry", inplace=True)
    joined.set_crs(geom.crs, inplace=True)
    return joined
