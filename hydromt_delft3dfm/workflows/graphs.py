# -*- coding: utf-8 -*-

import logging
from typing import Tuple, Union
import geopandas as gpd
import networkx as nx

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


def network_to_graph(
    edges: gpd.GeoDataFrame,
    nodes: Union[gpd.GeoDataFrame, None] = None,
) -> nx.Graph:
    """
    Converts the network to a graph.

    The resulting graph includes global properties (such as 'crs') and edge properties (such as 'edgeid', 'geometry',
    'node_start', and 'node_end').

    Parameters:
    -----------
    edges: gpd.GeoDataFrame
        The network lines.
    nodes: gpd.GeoDataFrame or None
        The network nodes, optional.

    Returns:
    --------
    nx.Graph
    """
    # method implementation goes here
    pass


def graph_to_network(graph: nx.Graph) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Converts the graph to a network representation.

    The resulting network includes global properties (such as 'crs') and edge properties (such as 'edgeid', 'geometry',
    'node_start', and 'node_end').

    Returns:
    --------
    Tuple:
        gpd.GeoDataFrame: edges
        gpd.GeoDataFrame: nodes
    """
    # method implementation goes here
    pass


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
