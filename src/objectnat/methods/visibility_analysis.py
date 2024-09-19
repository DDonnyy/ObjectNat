import math
from multiprocessing import cpu_count

import geopandas as gpd
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from shapely import LineString, MultiPolygon, Point, Polygon
from shapely.ops import polygonize, unary_union
from tqdm.contrib.concurrent import process_map

from objectnat import config

logger = config.logger


def get_visibility_accurate(point_from: Point, obstacles: gpd.GeoDataFrame, view_distance) -> Polygon:
    """
    Function to get accurate visibility from a given point to buildings within a given distance.

    Parameters
    ----------
    point_from : Point
        The point from which the line of sight is drawn.
    obstacles : gpd.GeoDataFrame
        A GeoDataFrame containing the geometry of the obstacles.
    view_distance : float
        The distance of view from the point.

    Returns
    -------
    Polygon
        A polygon representing the area of visibility from the given point.

    Notes
    -----
    If a quick result is important, consider using the `get_visibility_result()` function instead.
    However, please note that `get_visibility_result()` may provide less accurate results.

    Examples
    --------
    >>> point_from = Point(1, 1)
    >>> buildings = gpd.read_file('buildings.shp')
    >>> view_distance = 1000
    >>> visibility = get_visibility_accurate(point_from, obstacles, view_distance)
    """

    def get_point_from_a_thorough_b(a: Point, b: Point, dist):
        """
        Func to get Point from point a thorough point b on dist
        """
        direction = math.atan2(b.y - a.y, b.x - a.x)
        c_x = a.x + dist * math.cos(direction)
        c_y = a.y + dist * math.sin(direction)
        return Point(c_x, c_y)

    def polygon_to_linestring(geometry: Polygon):
        """A function to return all segments of a polygon as a list of linestrings"""
        coords_ext = geometry.exterior.coords  # Create a list of all line node coordinates
        polygons_inter = [Polygon(x) for x in geometry.interiors]
        result = [LineString(part) for part in zip(coords_ext, coords_ext[1:])]
        for poly in polygons_inter:
            poly_coords = poly.exterior.coords
            result.extend([LineString(part) for part in zip(poly_coords, poly_coords[1:])])
        return result

    point_buffer = point_from.buffer(view_distance, resolution=32)
    s = obstacles.intersects(point_buffer)
    buildings_in_buffer = obstacles.loc[s[s].index]
    # TODO kick all geoms except Polygons/MultiPolygons
    buildings_in_buffer = buildings_in_buffer.geometry.apply(
        lambda x: list(x.geoms) if isinstance(x, MultiPolygon) else x
    ).explode()

    buildings_lines_in_buffer = gpd.GeoSeries(buildings_in_buffer.apply(polygon_to_linestring).explode())
    buildings_lines_in_buffer = buildings_lines_in_buffer.loc[buildings_lines_in_buffer.intersects(point_buffer)]

    buildings_in_buffer_points = gpd.GeoSeries(
        [Point(line.coords[0]) for line in buildings_lines_in_buffer.geometry]
        + [Point(line.coords[-1]) for line in buildings_lines_in_buffer.geometry]
    )

    max_dist = max(view_distance, buildings_in_buffer_points.distance(point_from).max())
    polygons = []
    buildings_lines_in_buffer = gpd.GeoDataFrame(geometry=buildings_lines_in_buffer, crs=obstacles.crs).reset_index(
        drop=True
    )
    iteration = 0
    while not buildings_lines_in_buffer.empty:
        iteration += 1
        gdf_sindex = buildings_lines_in_buffer.sindex
        # TODO check if 2 walls are nearest and use the widest angle between points
        nearest_wall_sind = gdf_sindex.nearest(point_from, return_all=False)
        nearest_wall = buildings_lines_in_buffer.loc[nearest_wall_sind[1]].iloc[0]
        wall_points = [Point(coords) for coords in nearest_wall.geometry.coords]

        # Calculate angles and sort by angle
        points_with_angle = sorted(
            [(pt, math.atan2(pt.y - point_from.y, pt.x - point_from.x)) for pt in wall_points], key=lambda x: x[1]
        )

        delta_angle = 2 * math.pi + points_with_angle[0][1] - points_with_angle[-1][1]
        if delta_angle > math.pi:
            delta_angle = 2 * math.pi - delta_angle
        a = math.sqrt((max_dist**2) * (1 + (math.tan(delta_angle / 2) ** 2)))
        p1 = get_point_from_a_thorough_b(point_from, points_with_angle[0][0], a)
        p2 = get_point_from_a_thorough_b(point_from, points_with_angle[1][0], a)
        polygon = Polygon([points_with_angle[0][0], p1, p2, points_with_angle[1][0]])

        polygons.append(polygon)

        buildings_lines_in_buffer.drop(nearest_wall_sind[1], inplace=True)

        lines_to_kick = buildings_lines_in_buffer.within(polygon)
        buildings_lines_in_buffer = buildings_lines_in_buffer.loc[~lines_to_kick]
        buildings_lines_in_buffer.reset_index(drop=True, inplace=True)
    res = point_buffer.difference(unary_union(polygons))
    if isinstance(res, Polygon):
        return res
    res = list(res.geoms)
    polygon_containing_point = None
    for polygon in res:
        if polygon.contains(point_from):
            polygon_containing_point = polygon
            break
    return polygon_containing_point


def get_visibility(point: Point, obstacles: gpd.GeoDataFrame, view_distance: float, resolution: int = 32) -> Polygon:
    """
    Function to get a quick estimate of visibility from a given point to buildings within a given distance.

    Parameters
    ----------
    point : Point
        The point from which the line of sight is drawn.
    obstacles : gpd.GeoDataFrame
        A GeoDataFrame containing the geometry of the buildings.
    view_distance : float
        The distance of view from the point.
    resolution: int
        Buffer resolution for more accuracy (may give result slower)

    Returns
    -------
    Polygon
        A polygon representing the estimated area of visibility from the given point.

    Notes
    -----
    This function provides a quicker but less accurate result compared to `get_visibility_accurate()`.
    If accuracy is important, consider using `get_visibility_accurate()` instead.

    Examples
    --------
    >>> point = Point(1, 1)
    >>> buildings = gpd.read_file('buildings.shp')
    >>> view_distance = 1000
    >>> visibility = get_visibility(point, obstacles, view_distance)
    """

    point_buffer = point.buffer(view_distance, resolution=resolution)
    s = obstacles.within(point_buffer)
    buildings_in_buffer = obstacles.loc[s[s].index].reset_index(drop=True)
    buffer_exterior_ = list(point_buffer.exterior.coords)
    line_geometry = [LineString([point, ext]) for ext in buffer_exterior_]
    buffer_lines_gdf = gpd.GeoDataFrame(geometry=line_geometry)
    united_buildings = buildings_in_buffer.unary_union
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
    return circuit


def _multiprocess_get_vis(args):
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


def _polygons_to_linestring(geom):
    # pylint: disable-next=redefined-outer-name,reimported,import-outside-toplevel
    from shapely import LineString, MultiLineString, MultiPolygon

    def convert_polygon(polygon: Polygon):
        lines = []
        exterior = LineString(polygon.exterior.coords)
        lines.append(exterior)
        interior = [LineString(p.coords) for p in polygon.interiors]
        lines = lines + interior
        return lines

    def convert_multipolygon(polygon: MultiPolygon):
        return MultiLineString(sum([convert_polygon(p) for p in polygon.geoms], []))

    if geom.geom_type == "Polygon":
        return MultiLineString(convert_polygon(geom))
    if geom.geom_type == "MultiPolygon":
        return convert_multipolygon(geom)
    return geom


def _combine_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Combine geometry of intersecting layers into a single GeoDataFrame.
    Parameters
    ----------
    gdf: gpd.GeoDataFrame
        A GeoPandas GeoDataFrame

    Returns
    -------
    gpd.GeoDataFrame
        The combined GeoDataFrame with aggregated in lists columns.

    Examples
    --------
    >>> gdf = gpd.read_file('path_to_your_file.geojson')
    >>> result = _combine_geometry(gdf)
    """

    crs = gdf.crs
    polygons = polygonize(gdf["geometry"].apply(_polygons_to_linestring).unary_union)
    enclosures = gpd.GeoSeries(list(polygons), crs=crs)
    enclosures_points = gpd.GeoDataFrame(enclosures.representative_point(), columns=["geometry"], crs=crs)
    joined = gpd.sjoin(enclosures_points, gdf, how="inner", predicate="within").reset_index()
    cols = joined.columns.tolist()
    cols.remove("geometry")
    joined = joined.groupby("index").agg({column: list for column in cols})
    joined["geometry"] = enclosures
    joined = gpd.GeoDataFrame(joined, geometry="geometry", crs=crs)
    return joined


def _min_max_normalization(data, new_min=0, new_max=1):
    """
    Min-max normalization for a given array of data.

    Parameters
    ----------
    data: numpy.ndarray
        Input data to be normalized.
    new_min: float, optional
        New minimum value for normalization. Defaults to 0.
    new_max: float, optional
        New maximum value for normalization. Defaults to 1.

    Returns
    -------
    numpy.ndarray
        Normalized data.

    Examples
    --------
    >>> import numpy as np
    >>> data = np.array([1, 2, 3, 4, 5])
    >>> normalized_data = min_max_normalization(data, new_min=0, new_max=1)
    """

    min_value = np.min(data)
    max_value = np.max(data)
    normalized_data = (data - min_value) / (max_value - min_value) * (new_max - new_min) + new_min
    return normalized_data


def _process_group(group):
    geom = group
    combined_geometry = _combine_geometry(geom)
    combined_geometry.drop(columns=["index", "index_right"], inplace=True)
    combined_geometry["count_n"] = combined_geometry["ratio"].apply(len)
    combined_geometry["new_ratio"] = combined_geometry.apply(
        lambda x: np.power(np.prod(x.ratio), 1 / x.count_n) * x.count_n, axis=1
    )

    threshold = combined_geometry["new_ratio"].quantile(0.25)
    combined_geometry = combined_geometry[combined_geometry["new_ratio"] > threshold]

    combined_geometry["new_ratio_normalized"] = _min_max_normalization(
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

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Point, Polygon
    >>> points = gpd.GeoDataFrame({'geometry': [Point(0, 0), Point(1, 1)]}, crs='epsg:4326')
    >>> obstacles = gpd.GeoDataFrame({'geometry': [Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])]}, crs='epsg:4326')
    >>> view_distance = 100

    >>> visibilities = get_visibilities_from_points(points, obstacles, view_distance)
    >>> visibilities
    """
    # remove points inside polygons
    joined = gpd.sjoin(points, obstacles, how="left", predicate="intersects")
    points = joined[joined.index_right.isnull()]

    # remove unused obstacles
    points_view = points.geometry.buffer(view_distance).unary_union
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

    assert points.crs == obstacles.crs
    crs = obstacles.crs
    sectors_n = 12
    logger.info("Calculating Visibility Catchment Area from each point")
    all_visions_sectorized = get_visibilities_from_points(points, obstacles, view_distance, sectors_n, max_workers)
    all_visions_sectorized = gpd.GeoDataFrame(
        geometry=[item for sublist in all_visions_sectorized for item in sublist], crs=crs
    )
    logger.info("Calculating non-vision part...")
    all_visions_unary = all_visions_sectorized.unary_union
    convex = all_visions_unary.convex_hull
    dif = convex.difference(all_visions_unary)

    del convex, all_visions_unary

    buf_area = (math.pi * view_distance**2) / sectors_n
    all_visions_sectorized["ratio"] = all_visions_sectorized.area / buf_area
    all_visions_sectorized["ratio"] = _min_max_normalization(
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
    all_in = _combine_geometry(gpd.GeoDataFrame(data=pd.concat(groups_result), geometry="geometry", crs=crs))

    del groups_result

    all_in["count_n"] = all_in["index_right"].apply(len)

    logger.info("Calculating intersection's parameters")
    all_in["factor"] = all_in.parallel_apply(calc_group_factor, axis=1)
    threshold = all_in["factor"].quantile(0.3)
    all_in = all_in[all_in["factor"] > threshold]

    all_in["factor_normalized"] = np.round(
        _min_max_normalization(np.sqrt(all_in["factor"].values), new_min=1, new_max=5)
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
