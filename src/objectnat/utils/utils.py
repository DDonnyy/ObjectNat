import geopandas as gpd
from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info


def get_utm_crs_for_4326_gdf(gdf: gpd.GeoDataFrame) -> CRS:
    assert gdf.crs == CRS.from_epsg(4326), "provided GeoDataFrame is not in EPSG 4326"
    minx, miny, maxx, maxy = gdf.total_bounds
    utm_crs_list = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
            west_lon_degree=minx,
            south_lat_degree=miny,
            east_lon_degree=maxx,
            north_lat_degree=maxy,
        ),
    )
    return CRS.from_epsg(utm_crs_list[0].code)
