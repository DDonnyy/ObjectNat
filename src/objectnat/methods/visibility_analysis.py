import pandas as pd
from shapely import Point, LineString, Polygon, MultiPolygon
from shapely.ops import unary_union
import geopandas as gpd
import math
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
    if isinstance(res,Polygon):
        return res
    res = list(res.geoms)
    polygon_containing_point = None
    for polygon in res:
        if polygon.contains(point_from):
            polygon_containing_point = polygon
            break
    return polygon_containing_point


def get_visibility_result(
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
    >>> visibility = get_visibility_result(point, buildings, view_distance)
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

def multiprocess(args):
    point, buildings, view_distance = args
    return get_visibility_accurate(point,buildings,view_distance)

def calculate_visibility_catchment_area(points, buildings: gpd.GeoDataFrame, view_distance, chunksize):
    points_view = points.geometry.buffer(view_distance).unary_union
    s = buildings.intersects(points_view)
    buildings_in_buffer = buildings.loc[s[s].index].reset_index(drop=True)
    buildings_in_buffer.geometry = buildings_in_buffer.geometry.apply(lambda geom: MultiPolygon([geom]) if isinstance(geom, Polygon) else geom)
    args =[(point, buildings_in_buffer,view_distance) for point in points.geometry]
    results = process_map(multiprocess,args,chunksize=chunksize)
    return results