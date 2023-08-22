# -*- coding: utf-8 -*-

import logging
import pickle
from pathlib import Path
from typing import Tuple, Union

import geopandas as gpd
import momepy
import networkx as nx
from pyproj import CRS, Transformer
from shapely.geometry import LineString, Point, box
from shapely.ops import transform

logger = logging.getLogger(__name__)


__all__ = [
    # creat
    "preprocess_graph",
    "validate_graph",
    # convert
    "gpd_to_digraph",
    "network_to_graph",
    "graph_to_network",
    # io
    "read_graph",
    "write_graph",
    # property
    "graph_region",
    "graph_utils.graph_edges",
    "graph_nodes",
]


"""
graph: nx.Graph
    graph property ["crs"]
    edge property ["id", "geometry", "node_start", "node_end"]
    node property ["id", "geometry"]
    
graph_edges: gpd.GeoDataFrame
    columns: ["id", "geometry", "node_start", "node_end"]

graph_nodes: gpd.GeoDataFrame
    columns: ["id", "geometry"]
    
graph_region: gpd.GeoDataFrame
"""


def gpd_to_digraph(data: gpd.GeoDataFrame) -> nx.DiGraph():
    """Convert a `gpd.GeoDataFrame` to a `nx.DiGraph`.

    This is done by taking the first and last coordinate in a row as source and target.

    Parameters
    ----------
    data : gpd.GeoDataFrame
        The data to convert.

    Returns
    -------
    nx.DiGraph
        The converted directed graph.
    """
    DeprecationWarning("will be replaced by network_to_graph")

    _ = data.copy()

    _["from_node"] = [row.geometry.coords[0] for index, row in _.iterrows()]
    _["to_node"] = [row.geometry.coords[-1] for index, row in _.iterrows()]
    G = nx.from_pandas_edgelist(
        _,
        source="from_node",
        target="to_node",
        create_using=nx.DiGraph,
        edge_attr=True,
    )

    if _.crs is not None:
        G.graph["crs"] = _.crs.to_epsg()

    return G


def preprocess_graph(graph: nx.Graph, to_crs: CRS = None) -> nx.Graph:
    """
    Preprocess graph geometry and indexes.

    Parameters
    ----------
    graph: nx.Graph
        The input graph. Must provide crs in global properties
    to_crs: CRS
        The desired CRS for the graph's geometry (if reprojecting is desired).

    Returns
    -------
    nx.Graph
        The prepared graph.
        Includes global properties ['crs']
        edge properties ['id', 'geometry','node_start', and 'node_end']
        node properties ['id', 'geometry']
    """
    crs = graph.graph.get("crs")
    if not crs:
        raise ValueError("must provide crs in graph.")
    crs = CRS.from_user_input(crs)

    # assign graph geometry
    graph = assign_graph_geometry(graph)

    # If CRS is provided, reproject the graph
    if to_crs and to_crs != crs:
        graph = reproject_graph_geometry(graph, crs, to_crs)
        graph.graph["crs"] = graph.graph["crs"].to_epsg()

    # Assign node and edge indexes
    graph = assign_graph_index(graph)

    validate_graph(graph)

    return graph


# adapted from ra2ce funcs
def _add_missing_geoms_one_by_one(
    graph: nx.Graph, geom_name: str = "geometry"
) -> nx.Graph:
    # Not all nodes have geometry attributed (some only x and y coordinates) so add a geometry columns
    nodes_without_geom = [n[0] for n in graph.nodes(data=geom_name) if n[1] is None]
    for nd in nodes_without_geom:
        graph.nodes[nd][geom_name] = Point(graph.nodes[nd]["x"], graph.nodes[nd]["y"])

    edges_without_geom = [
        e for e in graph.edges.data(data=True) if geom_name not in e[-1]
    ]
    for ed in edges_without_geom:
        graph[ed[0]][ed[1]][0][geom_name] = LineString(
            [graph.nodes[ed[0]][geom_name], graph.nodes[ed[1]][geom_name]]
        )

    return graph


def assign_graph_geometry(graph: nx.Graph) -> nx.Graph:
    """Add geometry to graph nodes and edges.

    Note that this function does not add missing geometries one by one.
    """
    _exist_node_geometry = any(nx.get_node_attributes(graph, "geometry"))
    _exist_edge_geometry = any(nx.get_edge_attributes(graph, "geometry"))

    if _exist_node_geometry and _exist_edge_geometry:
        graph = _add_missing_geoms_one_by_one(graph)
        return graph

    if (not _exist_node_geometry) and _exist_edge_geometry:
        if any(nx.get_node_attributes(graph, "x")):
            # use xy from nodes
            geometry = {n: Point(a["x"], a["y"]) for n, a in graph.nodes(data=True)}
            nx.set_node_attributes(graph, geometry, name="geometry")
        else:
            # add node from egde geometry
            assert not graph.is_multigraph()  # FIXME support multigraph
            geometry = {}
            for s, t, e in graph.edges(data="geometry"):
                geometry.update({s: Point(e.coords[0]), t: Point(e.coords[-1])})
            nx.set_node_attributes(graph, geometry, name="geometry")
        return graph

    if _exist_node_geometry and (not _exist_edge_geometry):
        # add edge from node geometry
        assert not graph.is_multigraph()  # FIXME support multigraph
        geometry = {}
        for s, t, e in graph.edges(data="geometry"):
            geometry.update(
                {
                    (s, t): LineString(
                        [graph.nodes[s]["geometry"], graph.nodes[t]["geometry"]]
                    )
                }
            )
        nx.set_edge_attributes(graph, geometry, name="geometry")
        return graph

    if not (_exist_node_geometry or _exist_edge_geometry):
        # no geometry, check if the geometry can be derived from node tuple
        assert not graph.is_multigraph()  # FIXME support multigraph
        assert len(list(graph.nodes)[0]) == 2  # must have tuple (x,y) as nodes
        geometry = {n: Point(n[0], n[1]) for n in graph.nodes}
        nx.set_node_attributes(graph, geometry, name="geometry")
        return graph


def reproject_graph_geometry(graph: nx.Graph, crs_from: CRS, crs_to: CRS) -> nx.Graph:
    """
    Reprojects the geometry attributes of the nodes and edges of the graph.

    Parameters
    ----------
    graph: nx.Graph
        The input graph with geometry attributes.
    crs_from: CRS
        The current CRS of the graph's geometry.
    crs_to: CRS
        The desired CRS for the graph's geometry.

    Returns
    -------
    nx.Graph
        A new graph with reprojected geometry attributes.
    """
    # Create a new graph instance to avoid modifying the original graph
    reprojected_graph = graph.copy()

    # Initialize the transformer
    transformer = Transformer.from_crs(crs_from, crs_to, always_xy=True)

    # Reproject the nodes
    for node, attributes in graph.nodes(data=True):
        if "geometry" in attributes:
            reprojected_geometry = transform(
                transformer.transform, attributes["geometry"]
            )
            reprojected_graph.nodes[node].update({"geometry": reprojected_geometry})

    # Reproject the edges
    if graph.is_multigraph():
        for u, v, key, attributes in graph.edges(keys=True, data=True):
            reprojected_geometry = transform(
                transformer.transform, attributes["geometry"]
            )
            reprojected_graph[u][v][key].update({"geometry": reprojected_geometry})
    else:
        for u, v, attributes in graph.edges(data=True):
            if "geometry" in attributes:
                reprojected_geometry = transform(
                    transformer.transform, attributes["geometry"]
                )
                reprojected_graph[u][v].update({"geometry": reprojected_geometry})

    reprojected_graph.graph["crs"] = crs_to

    return reprojected_graph


def assign_graph_index(graph: nx.Graph) -> nx.Graph:
    """
    Assign unique indexes to the nodes and edges of the graph.

    Parameters
    ----------
    graph: nx.Graph
        The input graph.

    Returns
    -------
    nx.Graph
        The prepared graph.
        edge properties ['id', 'geometry','node_start', and 'node_end']
        node properties ['id', 'geometry']
    """
    # Create a new graph instance to avoid modifying the original graph
    indexed_graph = graph.copy()

    # Assign ids
    for index, node in enumerate(graph.nodes()):
        indexed_graph.nodes[node]["id"] = index

    # Assign ids
    if graph.is_multigraph():
        for u, v, key in graph.edges(keys=True):
            node_start = str(indexed_graph.nodes[u]["id"])
            node_end = str(indexed_graph.nodes[v]["id"])
            edge_id = node_start + "_" + node_end + "_" + str(key)
            indexed_graph[u][v][key]["id"] = edge_id
            indexed_graph[u][v][key]["node_start"] = node_start
            indexed_graph[u][v][key]["node_end"] = node_end
    else:
        for u, v in graph.edges():
            node_start = str(indexed_graph.nodes[u]["id"])
            node_end = str(indexed_graph.nodes[v]["id"])
            edge_id = node_start + "_" + node_end
            indexed_graph[u][v]["id"] = edge_id
            indexed_graph[u][v]["node_start"] = node_start
            indexed_graph[u][v]["node_end"] = node_end
    return indexed_graph


def network_to_graph(
    edges: gpd.GeoDataFrame,
    nodes: Union[gpd.GeoDataFrame, None] = None,
    create_using=nx.DiGraph,
) -> nx.Graph:
    """
    Convert the network edges and/or nodes to a graph.

    Parameters
    ----------
    edges: gpd.GeoDataFrame
        The network lines.
    nodes: gpd.GeoDataFrame or None
        The network nodes, optional.
    create_using : NetworkX graph constructor, optional (default=nx.Graph)
        Graph type to create. If graph instance, then cleared before populated.

    Returns
    -------
    nx.Graph
        The prepared graph.
        global properties ['crs']
        edge properties ['id', 'geometry','node_start', and 'node_end']
        node properties ['id', 'geometry']
    """
    # use memopy convention but not momepy.gdf_to_nx
    # create graph from edges
    graph = nx.from_pandas_edgelist(
        edges,
        source="node_start",
        target="node_end",
        create_using=create_using,
        edge_attr=True,
    )

    # add node attribtues
    if nodes is not None:
        nx.set_node_attributes(graph, {nid: {"id": nid} for nid in nodes["id"]})
        nx.set_node_attributes(graph, nodes.set_index("id").to_dict("index"))

    graph.graph["crs"] = edges.crs
    return graph


def graph_to_network(
    graph: nx.Graph, crs: CRS = None
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Convert the graph to network edges and nodes.

    Parameters
    ----------
    graph : nx.Graph
        The input graph.
        Required variables: ["geometry"] for graph.edges
        Optional variables: ["geometry"] for graph.nodes

        If "geometry" is not present, graph nodes should have "x" and "y" as attributes.
        Optional global property crs (pyproj.CRS).
    crs: CRS
        The crs of the output network.
        Optional if graph.graph["crs"] exsit.

    Returns
    -------
    Tuple
        gpd.GeoDataFrame
            edges of the graph.
            Contains attributes ['id', 'geometry', node_start', 'node_end', '_graph_edge_index']
        gpd.GeoDataFrame
            nodes of the graph.
            Contains attributes ['id', 'geometry', '_graph_node_index']

    See Also
    --------
    momepy.nx_to_gdf
    """
    # valida graph
    validate_graph(graph)

    # get crs
    crs = crs if crs is not None else graph.graph.get("crs")
    if not crs:
        raise ValueError("must provide crs.")

    # assign geometry
    graph = assign_graph_geometry(graph)

    # convert to geodataframe
    edges = gpd.GeoDataFrame(nx.to_pandas_edgelist(graph))
    nodes = gpd.GeoDataFrame(data=[attrs for _, attrs in graph.nodes(data=True)])

    # assign crs
    crs = CRS.from_user_input(crs)
    if not edges.crs:
        edges = edges.set_crs(crs)
        nodes = nodes.set_crs(crs)
    else:
        edges = edges.to_crs(crs)
        nodes = nodes.to_crs(crs)

    return edges, nodes


def validate_graph(graph: nx.Graph) -> None:
    """
    Validate the graph.

    The method checks that the graph has the required attribtues
    """
    _exist_attribute = "crs" in graph.graph
    assert _exist_attribute, "CRS does not exist in graph."
    _correct_attribute = isinstance(graph.graph["crs"], int)
    assert _correct_attribute, "crs must be epsg int."

    _exist_attribute = any(nx.get_node_attributes(graph, "id"))
    assert _exist_attribute, "Missing id in nodes"
    _exist_attribute = any(nx.get_node_attributes(graph, "geometry"))
    assert _exist_attribute, "Missing geometries in nodes"

    _exist_attribute = any(nx.get_edge_attributes(graph, "id"))
    assert _exist_attribute, "Missing id in edges"
    _exist_attribute = any(nx.get_edge_attributes(graph, "node_start"))
    assert _exist_attribute, "Missing node_start in edges"
    _exist_attribute = any(nx.get_edge_attributes(graph, "node_end"))
    assert _exist_attribute, "Missing node_end in edges"
    _exist_attribute = any(nx.get_edge_attributes(graph, "geometry"))
    assert _exist_attribute, "Missing geometries in edges"


def read_graph(graph_fn: str) -> nx.Graph:
    """
    Read the graph from a file.

    Parameters
    ----------
    graph_fn: str or Path
        The path to the file where the graph is stored.
        Supported extentions: gml, graphml, gpickle, geojson

    Returns
    -------
    nx.Graph
    """
    filepath = Path(graph_fn)
    extension = filepath.suffix.lower()

    if extension == ".gml":
        return nx.read_gml(filepath)
    elif extension == ".graphml":
        return nx.read_graphml(filepath)
    elif extension == ".gpickle":
        with open(filepath, "rb") as f:
            return pickle.load(f)
    elif extension == ".geojson":
        edges = gpd.read_file(filepath, driver="GeoJSON")
        return network_to_graph(edges)
    else:
        raise ValueError(f"Unsupported file format: {extension}")


def _pop_geometry_from_graph(graph: nx.Graph) -> nx.Graph:
    graph_copy = graph.copy()

    x = {}
    y = {}

    for node, attrs in graph_copy.nodes(data=True):
        geometry = attrs.pop("geometry", None)
        if geometry:
            x[node] = geometry.x
            y[node] = geometry.y

    for _, _, attrs in graph_copy.edges(data=True):
        attrs.pop("geometry", None)

    nx.set_node_attributes(graph_copy, x, name="x")
    nx.set_node_attributes(graph_copy, y, name="y")

    return graph_copy


def _pop_list_from_graph(graph: nx.Graph) -> nx.Graph:
    graph_copy = graph.copy()
    for u, v, attrs in graph_copy.edges(data=True):
        for key, value in list(
            attrs.items()
        ):  # Convert to list to avoid "dictionary size changed during iteration" error
            if isinstance(value, list):
                attrs.pop(key)
    for node, attrs in graph_copy.nodes(data=True):
        for key, value in list(
            attrs.items()
        ):  # Convert to list to avoid "dictionary size changed during iteration" error
            if isinstance(value, list):
                attrs.pop(key)
    return graph_copy


def write_graph(graph: nx.Graph, graph_fn: str) -> None:
    """
    Write the graph to a file.

    Parameters
    ----------
    graph: nx.Graph
        The graph to be written.
    graph_fn: str or Path
        The path to the file where the graph will be written.
        Supported extentions: gml, graphml, gpickle, geojson
        To write complete information, use gpickle.
        To write graph without geometry info, use gml.
        To write graph without geometry and any lists, use graphml.
        To write graph with geometry and index only, use geojson.

    Returns
    -------
    None
    """
    filepath = Path(graph_fn)
    extension = filepath.suffix.lower()

    if extension == ".gml":
        nx.write_gml(_pop_geometry_from_graph(graph), filepath)
    elif extension == ".graphml":
        nx.write_graphml(
            _pop_list_from_graph(_pop_geometry_from_graph(graph)), filepath
        )
    elif extension == ".gpickle":
        with open(filepath, "wb") as f:
            pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)
    elif extension == ".geojson":
        edges, _ = graph_to_network(graph)
        edges = edges[["id", "node_start", "node_end", "geometry"]]
        edges.to_file(filepath, driver="GeoJSON")
    else:
        raise ValueError(f"Unsupported file format: {extension}")


def graph_edges(graph: nx.Graph) -> gpd.GeoDataFrame:
    """Get graph edges as geodataframe"""
    return graph_to_network(graph)[0]


def graph_nodes(graph: nx.Graph) -> gpd.GeoDataFrame:
    """Get graph nodes as geodataframe"""
    return graph_to_network(graph)[1]


def graph_region(graph: nx.Graph) -> gpd.GeoDataFrame:
    """Get graph region as geodataframe"""
    edges = graph_to_network(graph)[0]
    region = gpd.GeoDataFrame(geometry=[box(*edges.total_bounds)], crs=edges.crs)
    return region
