import numpy as np
import pandas as pd
from shapely import Point, LineString, Polygon, MultiPolygon, MultiLineString
from shapely.ops import unary_union, polygonize
import geopandas as gpd
import math
from pandarallel import pandarallel
import tqdm.auto as tqdm
import time

from tqdm.contrib.concurrent import process_map


def get_visibility_accurate(point_from: Point, buildings: gpd.GeoDataFrame, view_distance) -> Polygon:
    """
    Function to get accurate visibility from a given point to buildings within a given distance.

    Parameters
    ----------
    point_from : Point
        The point from which the line of sight is drawn.
    buildings : gpd.GeoDataFrame
        A GeoDataFrame containing the geometry of the buildings.
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
    >>> visibility = get_visibility_accurate(point_from, buildings, view_distance)
    """

    def get_point_from_a_thorough_b(A: Point, B: Point, dist):
        """
        Func to get line from point A thorough point B on dist
        """
        direction = math.atan2(B.y - A.y, B.x - A.x)
        C_x = A.x + dist * math.cos(direction)
        C_y = A.y + dist * math.sin(direction)
        return Point(C_x, C_y)

    def polygon_to_linestring(geometry):
        """A function to return all segments of a line as a list of linestrings"""
        coords = geometry.exterior.coords  # Create a list of all line node coordinates
        return [LineString(part) for part in zip(coords, coords[1:])]

    point_buffer = point_from.buffer(view_distance, resolution=32)
    s = buildings.intersects(point_buffer)
    buildings_in_buffer = buildings.loc[s[s].index]
    buildings_in_buffer = buildings_in_buffer.geometry.apply(
        lambda x: list(x.geoms) if isinstance(x, MultiPolygon) else x
    ).explode()

    buildings_lines_in_buffer = gpd.GeoSeries(buildings_in_buffer.apply(lambda x: polygon_to_linestring(x)).explode())
    buildings_lines_in_buffer = buildings_lines_in_buffer.loc[buildings_lines_in_buffer.intersects(point_buffer)]

    buildings_in_buffer_points = gpd.GeoSeries(
        [Point(line.coords[0]) for line in buildings_lines_in_buffer.geometry]
        + [Point(line.coords[-1]) for line in buildings_lines_in_buffer.geometry]
    )

    max_dist = buildings_in_buffer_points.distance(point_from).max()
    polygons = []
    buildings_lines_in_buffer = gpd.GeoDataFrame(geometry=buildings_lines_in_buffer, crs=buildings.crs).reset_index(
        drop=True
    )
    while buildings_lines_in_buffer.shape[0] > 0:
        gdf_sindex = buildings_lines_in_buffer.sindex
        nearest_wall_sind = gdf_sindex.nearest(point_from, return_all=False)
        nearest_wall = buildings_lines_in_buffer.loc[nearest_wall_sind[1]].iloc[0]
        wall_points = [Point(coords) for coords in nearest_wall.geometry.coords]

        points_with_angle = []
        for point_to in wall_points:
            angle = math.atan2(point_to.y - point_from.y, point_to.x - point_from.x)
            points_with_angle.append((point_to, angle))

        points_with_angle.sort(key=lambda x: x[1])

        delta_angle = 2 * math.pi + points_with_angle[0][1] - points_with_angle[-1][1]
        if delta_angle > math.pi:
            delta_angle = 2 * math.pi - delta_angle
        a = math.sqrt((max_dist ** 2) * (1 + (math.tan(delta_angle / 2) ** 2)))
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


def get_visibility(
        point: Point, buildings: gpd.GeoDataFrame, view_distance: float, resolution: int = 32
) -> Polygon:
    """
    Function to get a quick estimate of visibility from a given point to buildings within a given distance.

    Parameters
    ----------
    point : Point
        The point from which the line of sight is drawn.
    buildings : gpd.GeoDataFrame
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
    >>> visibility = get_visibility(point, buildings, view_distance)
    """

    point_buffer = point.buffer(view_distance, resolution=resolution)
    s = buildings.within(point_buffer)
    buildings_in_buffer = buildings.loc[s[s].index].reset_index(drop=True)
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

    for u, v in splited_lines_gdf.groupby(level=0):
        splited_lines_list.append(v.iloc[0]["geometry"].coords[-1])
    circuit = Polygon(splited_lines_list)
    if united_buildings:
        circuit = circuit.difference(united_buildings)
    return circuit


def _multiprocess_get_vis(args):
    point, buildings, view_distance = args
    return get_visibility_accurate(point, buildings, view_distance)


def _polygons_to_linestring(geom):
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
    cols.remove('geometry')
    joined = joined.groupby("index").agg({column: list for column in cols})
    joined["geometry"] = enclosures
    joined = gpd.GeoDataFrame(joined, geometry='geometry', crs=crs)
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


def process_group(group):
    geom = group[0]
    result = _combine_geometry(geom)
    threshold = group[1]
    result['count_n'] = result['index_right'].apply(lambda x: len(x))
    result['area_normalized'] = (_min_max_normalization(result.geometry.area, new_min=1, new_max=10))
    result['new_ratio'] = result.apply(
        lambda x: np.power(np.prod(x.ratio), 1 / len(x.ratio)) * x.area_normalized * x.count_n, axis=1)
    result['new_ratio_normalized'] = (_min_max_normalization(result['new_ratio'].values, new_min=1, new_max=10))
    result = result[result['new_ratio_normalized'] > threshold]

    result['new_ratio_normalized'] = np.round(result['new_ratio_normalized']).astype(int)
    result = result.groupby("new_ratio_normalized").agg(
        {'geometry': lambda x: unary_union(MultiPolygon(list(x)).buffer(0))}).reset_index()
    result.set_geometry("geometry", inplace=True)
    result.set_crs(geom.crs, inplace=True)
    result = result.explode('geometry', index_parts=False).reset_index(drop=True)
    first = result.area.quantile(.75)
    result = result[result.area >= first]

    return result


def calculate_visibility_catchment_area(points, buildings: gpd.GeoDataFrame, view_distance, chunk_size, group_size):
    pandarallel.initialize(progress_bar=True)

    crs = buildings.crs
    points_view = points.geometry.buffer(view_distance).unary_union
    s = buildings.intersects(points_view)
    buildings_in_buffer = buildings.loc[s[s].index].reset_index(drop=True)
    buildings_in_buffer.geometry = buildings_in_buffer.geometry.apply(
        lambda geom: MultiPolygon([geom]) if isinstance(geom, Polygon) else geom)

    args = [(point, buildings_in_buffer, view_distance) for point in points.geometry]
    all_visions = process_map(_multiprocess_get_vis, args, chunksize=chunk_size,
                              desc='Calculating Visibility Catchment Area from each Point, it might take a while for a '
                                   'big amount of points')

    # groups = all_visions.sample(frac=1).groupby(all_visions.index // group_size)
    # groups = [group for _, group in groups]
    # groups_result = process_map(process_group, groups)
    # all_in = _combine_geometry(gpd.GeoDataFrame(data=pd.concat(groups_result), geometry='geometry', crs=crs))
    # all_in['factor'] = all_in.apply(lambda x: np.prod(x.normalized), axis=1)
    # all_in['factor_normalized'] = np.round(
    #     _min_max_normalization(all_in['factor'].values, new_min=1, new_max=10)).astype(int)
    # result = all_in.groupby("factor_normalized").agg(
    #     {'geometry': lambda x: unary_union(MultiPolygon(list(x)).buffer(0))}).reset_index()
    # result = gpd.GeoDataFrame(data=result, geometry='geometry', crs=crs)
    return all_visions


def start(groups):
    return process_map(process_group, groups)
