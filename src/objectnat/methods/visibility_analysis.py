import pandas as pd
from shapely import Point, LineString, Polygon, intersection, MultiPolygon
import geopandas as gpd
import math
import time


def get_visibility_from_point(mid_point: Point, buildings: gpd.GeoDataFrame, view_distance: int):

    def get_line_from_a_thorough_b(A: Point, B: Point, dist):
        """
        Func to get line from point A thorough point B on dist
        """
        direction = math.atan2(B.y - A.y, B.x - A.x)
        C_x = A.x + dist * math.cos(direction)
        C_y = A.y + dist * math.sin(direction)
        return LineString([Point(A.x, A.y), Point(C_x, C_y)])

    crs = buildings.crs
    point_buffer = Point(mid_point).buffer(view_distance)

    s = buildings.intersects(point_buffer)

    # getting all building's points in buffer

    buildings_in_buffer = buildings.loc[s[s].index].reset_index(drop=True)
    builds_points_inside = [
        Point(coords) for geom in buildings_in_buffer.geometry.unary_union.geoms for coords in geom.exterior.coords
    ]
    united_buildings = buildings_in_buffer.unary_union

    # raytracing lines from mid point to building's corners
    lines_to_corners = [LineString([mid_point, ext]) for ext in builds_points_inside]
    lines_to_corners = gpd.GeoDataFrame(geometry=lines_to_corners, crs=crs)
    lines_to_corners = lines_to_corners["geometry"].apply(lambda x: x.difference(united_buildings))
    lines_to_corners = gpd.GeoDataFrame(geometry=lines_to_corners, crs=crs).explode(index_parts=True)
    # filter lines, removing all intersecting parts
    points_to_build = pd.DataFrame()
    for u, v in lines_to_corners.groupby(level=0):
        points_to_build = pd.concat(
            [points_to_build, pd.DataFrame({"geometry": [Point(v.iloc[0].geometry.coords[-1])]})]
        )

    # raytracing lines from mid point thorough building's corners to buffer limit
    lines_to_buf = [get_line_from_a_thorough_b(mid_point, p_obn_b, view_distance) for p_obn_b in builds_points_inside]

    buffer_exterior_ = list(point_buffer.exterior.coords)
    lines_to_buf2 = [LineString([mid_point, ext]) for ext in buffer_exterior_]

    lines_to_buf = gpd.GeoDataFrame(geometry=lines_to_buf + lines_to_buf2, crs=crs)

    # decreasing building size, due to geometry issues
    united_buildings_buf = united_buildings.buffer(-0.005)
    lines_to_buf = lines_to_buf["geometry"].apply(lambda x: x.difference(united_buildings_buf))

    lines_to_buf = gpd.GeoDataFrame(geometry=lines_to_buf, crs=crs).explode(index_parts=True)
    # filter lines, removing all intersecting parts
    for u, v in lines_to_buf.groupby(level=0):
        points_to_build = pd.concat(
            [points_to_build, pd.DataFrame({"geometry": [Point(v.iloc[0].geometry.coords[-1])]})]
        )

    points_to_build.drop_duplicates(subset=["geometry"], inplace=True, keep="first")

    # calculating the angle of each line to determine the connection order
    points_with_angle = []
    for point_to in points_to_build.geometry:
        angle = math.atan2(point_to.y - mid_point.y, point_to.x - mid_point.x)
        points_with_angle.append((point_to, angle))

    points_with_angle = [(p, round(angle, 3), p.distance(mid_point)) for p, angle in points_with_angle]
    points_with_angle.sort(key=lambda x: -x[1])

    # for the equal angles, it is necessary to determine the connection order.
    # Connecting not just nearest points, but points that are closer in distance from the center
    new_points_with_angle = []
    point_pre_last = points_with_angle[-1]
    point_last = points_with_angle[0]

    new_tuple = (
        points_with_angle[0][0],
        points_with_angle[0][1],
        points_with_angle[0][2],
        abs(point_pre_last[2] - points_with_angle[0][2]),
    )
    new_points_with_angle.append(new_tuple)

    for i in range(1, len(points_with_angle)):
        current_point = points_with_angle[i]
        cur_angle = current_point[1]

        if cur_angle == point_last[1]:
            diff = abs(point_pre_last[2] - current_point[2])
            new_tuple = (current_point[0], current_point[1], current_point[2], diff)
            new_points_with_angle.append(new_tuple)
            point_last = current_point

        if cur_angle != point_last[1]:
            diff = abs(point_last[2] - current_point[2])
            new_tuple = (current_point[0], current_point[1], current_point[2], diff)
            new_points_with_angle.append(new_tuple)
            point_pre_last = point_last
            point_last = current_point

    new_points_with_angle.sort(key=lambda x: (-x[1], x[3]))
    # creating polygon with properly ordered points
    result = Polygon([t[0] for t in new_points_with_angle])
    result = intersection(point_buffer, result)
    result = (
        MultiPolygon([y for y in result.geoms if (y.geom_type == "Polygon" or y.geom_type == "MultiPolygon")])
        if result.geom_type == "GeometryCollection"
        else result
    )

    result = gpd.GeoDataFrame(geometry=[result], crs=crs)
    return result
