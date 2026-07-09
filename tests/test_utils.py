import geopandas as gpd


def test_edges_gdf(intermodal_osm_1114252_edges_gdf):
    assert isinstance(intermodal_osm_1114252_edges_gdf, gpd.GeoDataFrame)
    assert not intermodal_osm_1114252_edges_gdf.empty
    assert "geometry" in intermodal_osm_1114252_edges_gdf.columns
