import math

import geopandas as gpd
from shapely import LineString, MultiPolygon, Point, Polygon
from shapely.ops import polygonize, unary_union

from objectnat import config

logger = config.logger


def polygons_to_multilinestring(geom: Polygon | MultiPolygon):
    # pylint: disable-next=redefined-outer-name,reimported,import-outside-toplevel
    from shapely import LineString, MultiLineString, MultiPolygon

    def convert_polygon(polygon: Polygon):
        lines = []
        exterior = LineString(polygon.exterior)
        lines.append(exterior)
        interior = [LineString(p) for p in polygon.interiors]
        lines = lines + interior
        return lines

    def convert_multipolygon(polygon: MultiPolygon):
        return MultiLineString(sum([convert_polygon(p) for p in polygon.geoms], []))

    if geom.geom_type == "Polygon":
        return MultiLineString(convert_polygon(geom))
    return convert_multipolygon(geom)


def explode_linestring(geometry: LineString) -> list[LineString]:
    """A function to return all segments of a linestring as a list of linestrings"""
    coords_ext = geometry.coords  # Create a list of all line node coordinates
    result = [LineString(part) for part in zip(coords_ext, coords_ext[1:])]
    return result


def point_side_of_line(line: LineString, point: Point) -> int:
    """A positive indicates the left-hand side a negative indicates the right-hand side"""
    x1, y1 = line.coords[0]
    x2, y2 = line.coords[-1]
    x, y = point.coords[0]
    cross_product = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
    if cross_product > 0:
        return 1
    return -1


def get_point_from_a_thorough_b(a: Point, b: Point, dist):
    """
    Func to get Point from point a thorough point b on dist
    """
    direction = math.atan2(b.y - a.y, b.x - a.x)
    c_x = a.x + dist * math.cos(direction)
    c_y = a.y + dist * math.sin(direction)
    return Point(c_x, c_y)


def gdf_to_circle_zones_from_point(
    gdf: gpd.GeoDataFrame, point_from: Point, zone_radius, resolution=4, explode_multigeom=True
) -> gpd.GeoDataFrame:
    """n_segments = 4*resolution,e.g. if resolution = 4 that means there will be 16 segments"""
    crs = gdf.crs
    buffer = point_from.buffer(zone_radius, resolution=resolution)
    gdf_unary = gdf.clip(buffer, keep_geom_type=True).union_all()
    gdf_geometry = (
        gpd.GeoDataFrame(geometry=[gdf_unary], crs=crs)
        .explode(index_parts=True)
        .geometry.apply(polygons_to_multilinestring)
        .union_all()
    )
    zones_lines = [LineString([Point(coords1), Point(point_from)]) for coords1 in buffer.exterior.coords[:-1]]
    if explode_multigeom:
        return (
            gpd.GeoDataFrame(geometry=list(polygonize(unary_union([gdf_geometry] + zones_lines))), crs=crs)
            .clip(gdf_unary, keep_geom_type=True)
            .explode(index_parts=False)
        )
    return gpd.GeoDataFrame(geometry=list(polygonize(unary_union([gdf_geometry] + zones_lines))), crs=crs).clip(
        gdf_unary, keep_geom_type=True
    )


def remove_inner_geom(polygon: Polygon | MultiPolygon):
    """function to get rid of inner polygons"""
    if isinstance(polygon, Polygon):
        return Polygon(polygon.exterior.coords)
    if isinstance(polygon, MultiPolygon):
        polys = []
        for poly in polygon.geoms:
            polys.append(Polygon(poly.exterior.coords))
        return MultiPolygon(polys)
    else:
        return Polygon()


def combine_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
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
    >>> result = combine_geometry(gdf)
    """

    crs = gdf.crs

    enclosures = gpd.GeoDataFrame(
        geometry=list(polygonize(gdf["geometry"].apply(polygons_to_multilinestring).union_all())), crs=crs
    )
    enclosures_points = enclosures.copy()
    enclosures_points.geometry = enclosures.representative_point()
    joined = gpd.sjoin(enclosures_points, gdf, how="inner", predicate="within").reset_index()
    cols = joined.columns.tolist()
    cols.remove("geometry")
    joined = joined.groupby("index").agg({column: list for column in cols})
    joined["geometry"] = enclosures
    joined = gpd.GeoDataFrame(joined, geometry="geometry", crs=crs)
    return joined


def distribute_points_on_linestrings(lines: gpd.GeoDataFrame, radius, lloyd_relax_n=2) -> gpd.GeoDataFrame:
    lines = lines.copy()
    lines = lines.explode(ignore_index=True)
    lines = lines[lines.geom_type == "LineString"]
    original_crs = lines.crs
    lines = lines.to_crs(crs=lines.estimate_utm_crs())
    lines = lines.reset_index(drop=True)
    lines = lines[["geometry"]]
    radius = radius * 1.1
    segmentized = lines.geometry.apply(lambda x: x.simplify(radius).segmentize(radius))
    points = [Point(pt) for line in segmentized for pt in line.coords]

    points = gpd.GeoDataFrame(geometry=points, crs=lines.crs)
    lines["lines"] = lines.geometry
    geom_concave = lines.buffer(5, resolution=1).union_all()

    for i in range(lloyd_relax_n):
        points.geometry = points.voronoi_polygons().clip(geom_concave).centroid
        points = points.sjoin_nearest(lines, how="left")
        points = points[~points.index.duplicated(keep="first")]
        points["geometry"] = points["lines"].interpolate(points["lines"].project(points.geometry))
        points.drop(columns=["lines", "index_right"], inplace=True)

    return points.dropna().to_crs(original_crs)


def distribute_points_on_polygons(
    polygons: gpd.GeoDataFrame, radius, only_exterior=True, lloyd_relax_n=2
) -> gpd.GeoDataFrame:
    polygons = polygons.copy()
    polygons = polygons.explode(ignore_index=True)
    polygons = polygons[polygons.geom_type == "Polygon"]

    if only_exterior:
        polygons.geometry = polygons.geometry.apply(lambda x: LineString(x.exterior))
    else:
        polygons = gpd.GeoDataFrame(
            geometry=list(polygons.geometry.apply(polygons_to_multilinestring)), crs=polygons.crs
        )

    return distribute_points_on_linestrings(polygons, radius, lloyd_relax_n=lloyd_relax_n)
