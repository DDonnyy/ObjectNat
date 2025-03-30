import geopandas as gpd
import numpy as np


def get_radius_coverage(gdf_from:gpd.GeoDataFrame,radius:float,resolution:int=32):
    local_crs =  gdf_from.estimate_utm_crs()
    gdf_from = gdf_from.to_crs(local_crs)
    bounds = gdf_from.buffer(radius).union_all()
    coverage_polys = gpd.GeoDataFrame(geometry=gdf_from.voronoi_polygons().clip(bounds, keep_geom_type=True))
    coverage_polys = coverage_polys.sjoin(gdf_from)
    coverage_polys['area'] = coverage_polys.area
    coverage_polys['buffer'] = np.pow(coverage_polys['area'], 1 / 3)
    coverage_polys.geometry = coverage_polys.buffer(-coverage_polys['buffer'], resolution=1, join_style='mitre').buffer(
        coverage_polys['buffer'] * 0.9, resolution=resolution)
    coverage_polys.drop(columns=['buffer','area'], inplace=True)
    return coverage_polys