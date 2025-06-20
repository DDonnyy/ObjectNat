import networkx as nx

from objectnat import gdf_to_graph
import geopandas as gpd


def test_graph_to_gdf(intermodal_osm_1114252_edges_gdf):
    assert isinstance(intermodal_osm_1114252_edges_gdf, gpd.GeoDataFrame)
    assert not intermodal_osm_1114252_edges_gdf.empty
    assert "geometry" in intermodal_osm_1114252_edges_gdf.columns


def test_gdf_to_graph(intermodal_osm_1114252_edges_gdf):
    walk_edges = intermodal_osm_1114252_edges_gdf[intermodal_osm_1114252_edges_gdf['type'] == 'walk']
    graph = gdf_to_graph(walk_edges, project_gdf_attr=False, check_intersections=True)
    assert isinstance(graph, nx.Graph)
