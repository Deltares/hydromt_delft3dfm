# -*- coding: utf-8 -*-

import logging
from typing import List, Union

import numpy as np
from shapely.geometry import (
    LineString,
    MultiLineString,
)
from shapely.wkt import loads, dumps

from hydrolib.core.dflowfm import Branch, Network

logger = logging.getLogger(__name__)


__all__ = [
    "mesh1d_add_branch",
]


def mesh1d_add_branch(
    network: Network,
    branches: Union[
        LineString, MultiLineString, List[Union[LineString, MultiLineString]]
    ],
    node_distance: Union[float, int],
    branch_names: Union[str, List[str]] = None,
    branch_orders: Union[float, int, List[Union[float, int]]] = -1,
) -> List[str]:
    """Add branch to 1d mesh, from a (list of) (Multi)LineString geometry.
    The branch is discretized with the given node distance.
    if node distance is given as infinity, no discretization will be performed at mid point of the branch, #TODO, if minimum node distance is given, no discretization will be performed at mid point of the branch
    i.e. branch is treated as a pipe
    Args:
        network (Network): Network to which the branch is added
        branches (Union[ LineString, MultiLineString, List[Union[LineString, MultiLineString]] ]): Geometry object(s) for which the branch is created
        node_distance (Union[float, int]): Preferred node distance between branch nodes
        branch_names (Union[str, list[str]]): Branch names to be used in the mesh1d object
        branch_orders (Union[float, int, list[Union[float, int]]]): Branch orders to be used in the mesh1d object

    Returns:
        List[str]: List of names of added branches
    """
    if node_distance == np.inf:
        force_midpoint = False
    else:
        force_midpoint = True

    branchids = []

    if branch_names is None or isinstance(branch_names, str):
        branch_names = np.repeat(branch_names, len(branches))
    if isinstance(branch_orders, int) or isinstance(branch_orders, float):
        branch_orders = np.repeat(branch_orders, len(branches))

    for line, branch_name, branch_order in zip(
        branches, branch_names, branch_orders
    ):
        branch = Branch(
            geometry=np.array(round_geometry(line).coords[:])
        )  # avoid error caused by rounding precision
        branch.generate_nodes(node_distance)
        branchid = network.mesh1d_add_branch(
            branch,
            name=branch_name,
            branch_order=int(branch_order),
            force_midpoint=force_midpoint,
        )
        branchids.append(branchid)
    
    return branchids


def round_geometry(geometry, rounding_precision: int = 6):
    """
    Round the coordinates of the geometry object to the provided precision.
    Parameters
    ----------
    geometry
        The geometry object.
    rounding_preicision: int, optional
        Round coordinates to the specified number of digits.
        Defaults to 6.
    Returns
    -------
    A shapely geometry object.
    """
    return loads(dumps(geometry, rounding_precision=rounding_precision))
