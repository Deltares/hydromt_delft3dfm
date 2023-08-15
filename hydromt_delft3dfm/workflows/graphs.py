# -*- coding: utf-8 -*-

import logging
from typing import Tuple, Union
import geopandas as gpd
import networkx as nx
import momepy
from pyproj.crs import CRS
from shapely.geometry import (
    LineString,
    Point,
)

logger = logging.getLogger(__name__)


__all__ = [
    "gpd_to_digraph",
    "network_to_graph",
    "graph_to_network",
    "validate_graph",
    "read_graph",
    "write_graph",
]


def gpd_to_digraph(data: gpd.GeoDataFrame) -> nx.DiGraph():
    """Convert a `gpd.GeoDataFrame` to a `nx.DiGraph` by taking the first and last coordinate in a row as source and target, respectively.

    Parameters
    ----------
    data : gpd.GeoDataFrame
        The data to convert.

    Returns
    -------
    nx.DiGraph
        The converted directed graph.
    """
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
    return G


# TODO: func from ra2ce utils due to installation issue
def _add_missing_geoms_graph(graph: nx.Graph, geom_name: str = "geometry") -> nx.Graph:
    # Not all nodes have geometry attributed (some only x and y coordinates) so add a geometry columns
    nodes_without_geom = [
        n[0] for n in graph.nodes(data=True) if geom_name not in n[-1]
    ]
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


def network_to_graph(
    edges: gpd.GeoDataFrame,
    nodes: Union[gpd.GeoDataFrame, None] = None,
    create_using=nx.DiGraph,
) -> nx.Graph:
    """
    Converts the network (edges and nodes) to a graph.

    The resulting graph includes global properties ('crs'),  edge properties (such as 'edgeid', 'geometry',
    'node_start', and 'node_end') and node properties ('nodeid', 'geometry')

    Parameters:
    -----------
    edges: gpd.GeoDataFrame
        The network lines.
    nodes: gpd.GeoDataFrame or None
        The network nodes, optional.
    create_using : NetworkX graph constructor, optional (default=nx.Graph)
        Graph type to create. If graph instance, then cleared before populated.


    Returns:
    --------
    nx.Graph
    """

    # use memopy convention but not momepy.gdf_to_nx because it does not create nodes using nodeID
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
        nx.set_node_attributes(graph, {nid: {"nodeid": nid} for nid in nodes["nodeid"]})
        nx.set_node_attributes(graph, nodes.set_index("nodeid").to_dict("index"))

    graph.graph["crs"] = edges.crs
    return graph


def graph_to_network(
    graph: nx.Graph, crs: CRS = None
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Converts the graph to a network representation.

    The resulting network includes global properties (such as 'crs') and edge properties (such as 'edgeid', 'geometry',
    'node_start', and 'node_end').

    Parameters
    ----------
    graph : nx.Graph
        The graph with "geometry" as edge and node attribute.
        If "geometry" is not present, graph nodes should have "x" and "y" as attributes.
        Optional global property crs (pyproj.CRS).
    crs: CRS
        The crs of the output network.
        Optional if graph.graph["crs"] exsit.

    Returns
    -------
    Tuple
        gpd.GeoDataFrame: edges of the graph. Contains attributes ['edgeid', 'geometry', node_start', 'node_end']
        gpd.GeoDataFrame: nodes of the graph. Contains attributes ['nodeid', 'geometry']

    See Also
    --------
    momepy.nx_to_gdf
    """
    if not crs:
        crs = graph.graph.get("crs")

    crs = CRS.from_user_input(crs)
    if not crs:
        raise ValueError("must provide crs.")

    # obtain edges and nodes from the simplfied graph
    graph = _add_missing_geoms_graph(graph)
    nodes, edges = momepy.nx_to_gdf(graph, nodeID="nodeid")
    edges["edgeid"] = (
        edges["node_start"].astype(str) + "_" + edges["node_end"].astype(str)
    )
    # convert network back to region crs
    edges = edges.to_crs(crs)
    nodes = nodes.to_crs(crs)
    return edges, nodes


def validate_graph(graph: nx.Graph) -> bool:
    """
    Validates the graph.

    The method checks that the graph has the required global property ('crs') and node and edge properties
    ('nodeid', 'geometry', 'edgeid').

    Returns:
    --------
    bool:
        True if the graph is valid, False otherwise.
    """
    # method implementation goes here
    pass


def read_graph(graph_fn) -> nx.Graph:
    """
    Reads a graph from a file.

    Parameters:
    -----------
    filepath: str or Path
        The path to the file from which to read the graph.

    Returns:
    --------
    None
    """
    # method implementation goes here
    pass


def write_graph(graph: nx.graph, graph_fn) -> None:
    """
    Writes the graph to a file.

    Parameters:
    -----------
    filepath: str or Path
        The path to the file where the graph will be written.

    Returns:
    --------
    nx.graph
    """
    # method implementation goes here
    pass
