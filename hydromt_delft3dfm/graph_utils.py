# -*- coding: utf-8 -*-

import logging

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import Point

logger = logging.getLogger(__name__)


__all__ = ["gpd_to_digraph", "get_endnodes_from_lines"]


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


def get_endnodes_from_lines(
    lines: gpd.GeoDataFrame, where: str = "both"
) -> gpd.GeoDataFrame:
    """Get the possible boundary locations from the branches with id.

    Parameters
    ----------
    where : {'both', 'upstream', 'downstream'}
        Where at the branches should the boundaries be derived.
        An upstream end node is defined as a node which has 0 incoming branches and 1 outgoing branch.
        A downstream end node is defined as a node which has 1 incoming branch and 0 outgoing branches.

    Returns
    -------
    gpd.GeoDataFrame
        A data frame containing all the upstream and downstream end nodes of the branches
    """
    # convert branches to graph
    G = gpd_to_digraph(lines)

    # get boundary locations at where
    if where == "downstream":
        endnodes = {
            dn: {**d, **{"where": "downstream"}}
            for up, dn, d in G.edges(data=True)
            if G.out_degree[dn] == 0 and G.degree[dn] == 1
        }
    elif where == "upstream":
        endnodes = {
            up: {**d, **{"where": "upstream"}}
            for up, dn, d in G.edges(data=True)
            if G.in_degree[up] == 0 and G.degree[up] == 1
        }
    elif where == "both":
        endnodes = {
            dn: {**d, **{"where": "downstream"}}
            for up, dn, d in G.edges(data=True)
            if G.out_degree[dn] == 0 and G.degree[dn] == 1
        }
        endnodes.update(
            {
                up: {**d, **{"where": "upstream"}}
                for up, dn, d in G.edges(data=True)
                if G.in_degree[up] == 0 and G.degree[up] == 1
            }
        )
    else:
        pass

    if len(endnodes) == 0:
        logger.error(f"cannot generate endnodes for given condition {where}")

    endnodes_pd = (
        pd.DataFrame().from_dict(endnodes, orient="index").drop(columns=["geometry"])
    )
    endnodes_gpd = gpd.GeoDataFrame(
        data=endnodes_pd,
        geometry=[Point(endnode) for endnode in endnodes],
        crs=lines.crs,
    )
    endnodes_gpd.reset_index(inplace=True)
    return endnodes_gpd
