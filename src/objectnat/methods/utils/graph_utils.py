import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from loguru import logger
from scipy.spatial import KDTree
from shapely import LineString, MultiLineString, line_merge, node
from shapely.geometry.point import Point


def _edges_to_gdf(graph: nx.Graph, crs) -> gpd.GeoDataFrame:
    """
    Converts nx graph to gpd.GeoDataFrame as edges.
    """
    graph_df = pd.DataFrame(list(graph.edges(data=True)), columns=["u", "v", "data"])
    edge_data_expanded = pd.json_normalize(graph_df["data"])
    graph_df = pd.concat([graph_df.drop(columns=["data"]), edge_data_expanded], axis=1)
    graph_df = gpd.GeoDataFrame(graph_df, geometry="geometry", crs=crs).set_index(["u", "v"])
    graph_df["geometry"] = graph_df["geometry"].fillna(LineString())
    return graph_df


def _nodes_to_gdf(graph: nx.Graph, crs: int) -> gpd.GeoDataFrame:
    """
    Converts nx graph to gpd.GeoDataFrame as nodes.
    """

    ind, data = zip(*graph.nodes(data=True))
    node_geoms = (Point(d["x"], d["y"]) for d in data)
    gdf_nodes = gpd.GeoDataFrame(data, index=ind, crs=crs, geometry=list(node_geoms))

    return gdf_nodes


def _restore_edges_geom(nodes_gdf, edges_gdf) -> gpd.GeoDataFrame:
    edges_wout_geom = edges_gdf[edges_gdf["geometry"].is_empty].reset_index()
    edges_wout_geom["geometry"] = [
        LineString((s, e))
        for s, e in zip(
            nodes_gdf.loc[edges_wout_geom["u"], "geometry"], nodes_gdf.loc[edges_wout_geom["v"], "geometry"]
        )
    ]
    edges_wout_geom.set_index(["u", "v"], inplace=True)
    edges_gdf.update(edges_wout_geom)
    return edges_gdf


def graph_to_gdf(
    graph: nx.MultiDiGraph | nx.Graph | nx.DiGraph, edges: bool = True, nodes: bool = True, restore_edge_geom=False
) -> gpd.GeoDataFrame | tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Converts nx graph to gpd.GeoDataFrame as edges.

    Parameters
    ----------
    graph : nx.MultiDiGraph
        The graph to convert.
    edges: bool, default to True
        Keep edges in GoeDataFrame.
    nodes: bool, default to True
        Keep nodes in GoeDataFrame.
    restore_edge_geom: bool, default to False
        if True, will try to restore edge geometry from nodes.
    Returns
    -------
    gpd.GeoDataFrame | tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        Graph representation in GeoDataFrame format, either nodes or nodes,or tuple of them nodes,edges.
    """
    try:
        crs = graph.graph["crs"]
    except KeyError as exc:
        raise ValueError("Graph does not have crs attribute") from exc
    if not edges and not nodes:
        raise AttributeError("Neither edges or nodes were selected")
    if nodes and not edges:
        nodes_gdf = _nodes_to_gdf(graph, crs)
        return nodes_gdf
    if not nodes and edges:
        edges_gdf = _edges_to_gdf(graph, crs)
        if restore_edge_geom:
            nodes_gdf = _nodes_to_gdf(graph, crs)
            edges_gdf = _restore_edges_geom(nodes_gdf, edges_gdf)
        return edges_gdf

    nodes_gdf = _nodes_to_gdf(graph, crs)
    edges_gdf = _edges_to_gdf(graph, crs)
    if restore_edge_geom:
        edges_gdf = _restore_edges_geom(nodes_gdf, edges_gdf)
    return nodes_gdf, edges_gdf


def gdf_to_graph(
    gdf: gpd.GeoDataFrame, project_gdf_attr=True, reproject_to_utm_crs=True, speed=5, check_intersections=True
) -> nx.DiGraph:
    """
    Converts a GeoDataFrame of LineStrings into a directed graph (nx.DiGraph).

    This function transforms a set of linear geometries (which may or may not form a planar graph)
    into a directed graph where each edge corresponds to a LineString (or its segment) from the GeoDataFrame.
    Intersections are optionally checked and merged. Attributes from the original GeoDataFrame
    can be projected onto the graph edges using spatial matching.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        A GeoDataFrame containing at least one LineString geometry.
    project_gdf_attr : bool, default=True
        If True, attributes from the input GeoDataFrame will be spatially projected
        onto the resulting graph edges. This can be an expensive operation for large datasets.
    reproject_to_utm_crs : bool, default=True
        If True, reprojects the GeoDataFrame to the estimated local UTM CRS
        to ensure accurate edge length calculations in meters.
        If False, edge lengths are still computed in UTM CRS, but the final graph
        will remain in the original CRS of the input GeoDataFrame.
    speed : float, default=5
        Assumed travel speed in km/h used to compute edge traversal time (in minutes).
    check_intersections : bool, default=True
        If True, merges geometries to ensure topological correctness.
        Can be disabled if the input geometries already form a proper planar graph
        with no unintended intersections.

    Returns
    -------
    nx.DiGraph
        A directed graph where each edge corresponds to a line segment from the input GeoDataFrame.
        Edge attributes include geometry, length in meters, travel time (in minutes), and any
        additional projected attributes from the original GeoDataFrame.

    Raises
    ------
    ValueError
        If the input GeoDataFrame contains no valid LineStrings.
    """

    def unique_list(agg_vals):
        agg_vals = list(set(agg_vals.dropna()))
        if len(agg_vals) == 1:
            return agg_vals[0]
        return agg_vals

    original_crs = gdf.crs
    gdf = gdf.to_crs(gdf.estimate_utm_crs())

    gdf = gdf.explode(ignore_index=True)
    gdf = gdf[gdf.geom_type == "LineString"]

    if len(gdf) == 0:
        raise ValueError("Provided GeoDataFrame contains no valid LineStrings")

    if check_intersections:
        lines = line_merge(node(MultiLineString(gdf.geometry.to_list())))
    else:
        lines = line_merge(MultiLineString(gdf.geometry.to_list()))

    lines = gpd.GeoDataFrame(geometry=list(lines.geoms), crs=gdf.crs)

    if len(gdf.columns) > 1 and project_gdf_attr:
        lines_centroids = lines.copy()
        lines_centroids.geometry = lines_centroids.apply(
            lambda row: row.geometry.line_interpolate_point(row.geometry.length / 2), axis=1
        ).buffer(0.05, resolution=2)
        lines_with_attrs = gpd.sjoin(lines_centroids, gdf, how="left", predicate="intersects")
        aggregated_attrs = (
            lines_with_attrs.drop(columns=["geometry", "index_right"])  # убрать геометрию буфера
            .groupby(lines_with_attrs.index)
            .agg(unique_list)
        )
        lines = pd.concat([lines, aggregated_attrs], axis=1)

    lines["length_meter"] = np.round(lines.length, 2)
    if not reproject_to_utm_crs:
        lines = lines.to_crs(original_crs)

    coords = lines.geometry.get_coordinates()
    coords_grouped_by_index = coords.reset_index(names="old_index").groupby("old_index")
    start_coords = coords_grouped_by_index.head(1).apply(lambda a: (a.x, a.y), axis=1).rename("start")
    end_coords = coords_grouped_by_index.tail(1).apply(lambda a: (a.x, a.y), axis=1).rename("end")
    coords = pd.concat([start_coords.reset_index(), end_coords.reset_index()], axis=1)[["start", "end"]]
    lines = pd.concat([lines, coords], axis=1)
    unique_coords = pd.concat([coords["start"], coords["end"]], ignore_index=True).unique()
    coord_to_index = {coord: idx for idx, coord in enumerate(unique_coords)}

    lines["u"] = lines["start"].map(coord_to_index)
    lines["v"] = lines["end"].map(coord_to_index)

    speed = speed * 1000 / 60
    lines["time_min"] = np.round(lines["length_meter"] / speed, 2)

    graph = nx.Graph()
    for coords, node_id in coord_to_index.items():
        x, y = coords
        graph.add_node(node_id, x=float(x), y=float(y))

    columns_to_attr = lines.columns.difference(["start", "end", "u", "v"])
    for _, row in lines.iterrows():
        edge_attrs = {}
        for col in columns_to_attr:
            edge_attrs[col] = row[col]
        graph.add_edge(row.u, row.v, **edge_attrs)

    graph.graph["crs"] = lines.crs
    graph.graph["speed m/min"] = speed
    return nx.DiGraph(graph)


def get_closest_nodes_from_gdf(gdf: gpd.GeoDataFrame, nx_graph: nx.Graph) -> tuple:
    """
    Finds the closest graph nodes to the geometries in a GeoDataFrame.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame with geometries for which the nearest graph nodes will be found.
    nx_graph : nx.Graph
        A NetworkX graph where nodes have 'x' and 'y' attributes (coordinates).

    Returns
    -------
    tuple
        A tuple of (distances, nearest_nodes), where:
        - distances: List of distances from each geometry to the nearest node.
        - nearest_nodes: List of node IDs closest to each geometry in the input GeoDataFrame.

    Raises
    ------
    ValueError
        If any node in the graph is missing 'x' or 'y' attributes.
    """
    nodes_with_data = list(nx_graph.nodes(data=True))
    try:
        coordinates = np.array([(data["x"], data["y"]) for node, data in nodes_with_data])
    except KeyError as e:
        raise ValueError("Graph does not have coordinates attribute") from e
    tree = KDTree(coordinates)
    target_coord = [(p.x, p.y) for p in gdf.representative_point()]
    distances, indices = tree.query(target_coord)
    nearest_nodes = [nodes_with_data[idx][0] for idx in indices]
    return distances, nearest_nodes


def remove_weakly_connected_nodes(graph: nx.DiGraph) -> nx.DiGraph:
    """
    Removes all nodes that are not part of the largest strongly connected component in the graph.

    Parameters
    ----------
    graph : nx.DiGraph
        A directed NetworkX graph.

    Returns
    -------
    nx.DiGraph
        A new graph with only the largest strongly connected component retained.

    Notes
    -----
    - Also logs a warning if multiple weakly connected components are detected.
    - Logs the number of nodes removed and size of the remaining component.
    """
    graph = graph.copy()

    weakly_connected_components = list(nx.weakly_connected_components(graph))
    if len(weakly_connected_components) > 1:
        logger.warning(
            f"Found {len(weakly_connected_components)} disconnected subgraphs in the network. "
            f"These are isolated groups of nodes with no connections between them. "
            f"Size of components: {[len(c) for c in weakly_connected_components]}"
        )

    all_scc = sorted(nx.strongly_connected_components(graph), key=len)
    nodes_to_del = set().union(*all_scc[:-1])

    if nodes_to_del:
        logger.warning(
            f"Removing {len(nodes_to_del)} nodes that form {len(all_scc) - 1} trap components. "
            f"These are groups where you can enter but can't exit (or vice versa). "
            f"Keeping the largest strongly connected component ({len(all_scc[-1])} nodes)."
        )
        graph.remove_nodes_from(nodes_to_del)

    return graph


def reverse_graph(nx_graph: nx.Graph, weight: str) -> tuple[nx.Graph, nx.DiGraph]:
    """
    Generate a reversed version of a directed or weighted graph.

    If the input graph is undirected, the original graph is returned as-is.
    For directed graphs, the function returns a new graph with all edge directions reversed,
    preserving the specified edge weight.

    Parameters
    ----------
    nx_graph : nx.Graph
        Input NetworkX graph (can be directed or undirected).
    weight : str
        Name of the edge attribute to use as weight in graph conversion.

    Returns
    -------
    tuple[nx.Graph, nx.DiGraph]
        A tuple containing:
        - normalized_graph: Original graph with relabeled nodes (if needed)
        - reversed_graph: Directed graph with reversed edges and preserved weights
    """

    if nx_graph.is_multigraph():
        nx_graph = nx.DiGraph(nx_graph) if nx_graph.is_directed() else nx.Graph(nx_graph)
    if not nx_graph.is_multigraph() and not nx_graph.is_directed():
        return nx_graph, nx_graph

    nx_graph = remove_weakly_connected_nodes(nx_graph)

    mapping = {old_label: new_label for new_label, old_label in enumerate(nx_graph.nodes())}
    nx_graph = nx.relabel_nodes(nx_graph, mapping)

    sparse_matrix = nx.to_scipy_sparse_array(nx_graph, weight=weight)
    transposed_matrix = sparse_matrix.transpose()
    reversed_graph = nx.from_scipy_sparse_array(transposed_matrix, edge_attribute=weight, create_using=type(nx_graph))
    return nx_graph, reversed_graph
