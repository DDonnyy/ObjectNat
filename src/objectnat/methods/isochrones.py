import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from shapely import Point
from loguru import logger
from shapely.ops import unary_union


def get_accessibility_isochrones(
    points: Point | gpd.GeoSeries,
    weight_value: int,
    weight_type: str,
    graph_nx: nx.Graph,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """
    Calculate accessibility isochrones based on the provided city graph for the selected graph_type from point(s).

    Parameters
    ----------

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]
        A tuple containing the accessibility isochrones, and optionally the routes and stops.

    Examples
    --------

    """

    nodes_with_data = list(graph_nx.nodes(data=True))
    logger.info("Calculating isochrones distances...")
    coordinates = np.array([(data["x"], data["y"]) for node, data in nodes_with_data])
    tree = KDTree(coordinates)

    if isinstance(points,Point):
        points = [points]

    target_coord = [(p.x, p.y) for p in points]
    distance, indices = tree.query(target_coord)

    nearest_nodes = [(index, nodes_with_data[idx][0]) for index, idx in enumerate(indices)] # TODO use get_closest_nodes from IduEdu

    dist_nearest = pd.DataFrame(data=distance, index=pd.MultiIndex.from_tuples(nearest_nodes), columns=["dist"])
    walk_speed = 5 * 1000 / 60 # TODO Скорость из графа взять, если есть
    dist_nearest = dist_nearest / walk_speed if weight_type == "time_min" else dist_nearest

    if (dist_nearest > weight_value).all().all():
        raise RuntimeError(
            "The point(s) lie further from the graph than weight_value, it's impossible to "
            "construct isochrones. Check the coordinates of the point(s)/their projection"
        )

    data = []
    for point_index, source in nearest_nodes:
        dist, path = nx.single_source_dijkstra(
            graph_nx, source, weight=weight_type, cutoff=weight_value
        )
        for target_node, way in path.items():
            source = way[0]
            distance = dist.get(target_node, np.nan)
            data.append((point_index, source, target_node, distance))

    dist_matrix = pd.DataFrame(data, columns=["point_index", "source", "destination", "distance"])
    dist_matrix = dist_matrix.pivot_table(index=["point_index", "source"], columns="destination", values="distance")

    dist_matrix = dist_matrix.add(dist_nearest.dist, axis=0).transpose()
    dist_matrix = dist_matrix.mask(dist_matrix >= weight_value, np.nan)
    dist_matrix.dropna(how="all", inplace=True)

    results = []
    point_num = 0
    logger.info("Building isochrones geometry...")

    for unique_ind, node_ind in dist_matrix.columns:
        geometry = []
        for ind in dist_matrix.index:
            value = dist_matrix.loc[ind, (unique_ind, node_ind)]
            if not pd.isna(value):
                node = graph_nx.nodes[ind]
                point = Point(node["x"], node["y"])
                geometry.append(
                    point.buffer(round(weight_value - value, 2) * walk_speed * 0.8)
                    if weight_type == "time_min"
                    else point.buffer(round(weight_value - value, 2) * 0.8)
                )
        geometry = unary_union(geometry)
        results.append({"geometry": geometry, "point": str(points.iloc[point_num]), "point_number": point_num})
        point_num += 1

    isochrones = gpd.GeoDataFrame(data=results, geometry="geometry", crs=city_crs)
    # isochrones["travel_type"] = str([d.russian_name for d in graph_type])
    isochrones["weight_type"] = weight_type
    isochrones["weight_value"] = weight_value
    stops, routes = (None, None)
    if GraphType.PUBLIC_TRANSPORT in self.graph_type and self.weight_type == "time_min":
        node_data = {
            node: mobility_graph.nodes[node]
            for node in dist_matrix.index
            if mobility_graph.nodes[node]["stop"] == "True"
        }
        if len(node_data) > 0:
            logger.info("Building public transport geometry...")
            stops, routes = self._get_routes(node_data, mobility_graph)
        else:
            logger.info("No public transport nodes in accessibility")
    return isochrones, routes, stops


def _get_routes(self, stops, mobility_graph):
    stops = pd.DataFrame(stops).T
    stops["geometry"] = stops.apply(lambda x: Point(x.x, x.y), axis=1)
    stops.drop(columns=["x", "y"], inplace=True)

    subgraph = mobility_graph.subgraph(stops.index)
    routes = pd.DataFrame.from_records([e[-1] for e in subgraph.edges(data=True)])
    if routes.empty:
        return None, None
    routes["geometry"] = routes["geometry"].apply(lambda x: from_wkt(str(x)))
    routes_result = gpd.GeoDataFrame(data=routes, geometry="geometry", crs=self.city_crs)
    stops_result = gpd.GeoDataFrame(data=stops, geometry="geometry", crs=self.city_crs)
    return stops_result, routes_result
