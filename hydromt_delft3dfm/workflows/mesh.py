# -*- coding: utf-8 -*-

import logging
from typing import List, Union

import numpy as np
import xugrid as xu
from hydrolib.core.dflowfm import Branch, Network, Mesh1d
from meshkernel import GeometryList
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon, box
from shapely.wkt import dumps, loads

logger = logging.getLogger(__name__)


__all__ = [
    "mesh1d_add_branch",
    "links1d2d_add_links_1d_to_2d",
    "links1d2d_add_links_2d_to_1d_embedded",
    "links1d2d_add_links_2d_to_1d_lateral",
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
        branch_orders (Union[float, int, list[Union[float, int]]]): Branch orders to be used in the mesh1d object.

    Returns
    -------
        Network: Network with added branches
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

    for line, branch_name, branch_order in zip(branches, branch_names, branch_orders):
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

    return network, branchids


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


# below function is a copy of functions from hydrolib, to avoid installation error.


def polygon_to_geometrylist(polygon):
    # Create a list of coordinate lists
    cls = GeometryList
    # Add exterior
    x_ext, y_ext = np.array(polygon.exterior.coords[:]).T
    x_crds = [x_ext]
    y_crds = [y_ext]
    # Add interiors, seperated by inner_outer_separator
    for interior in polygon.interiors:
        x_int, y_int = np.array(interior.coords[:]).T
        x_crds.append([cls.inner_outer_separator])
        x_crds.append(x_int)
        y_crds.append([cls.inner_outer_separator])
        y_crds.append(y_int)
    gl = cls(x_coordinates=np.concatenate(x_crds), y_coordinates=np.concatenate(y_crds))
    return gl


def links1d2d_add_links_1d_to_2d(
    network: Network,
    branchids: List[str] = None,
    within: Union[Polygon, MultiPolygon] = None,
    max_length: float = np.inf,
) -> None:
    """Function to add 1d2d links to network, by generating them from 1d to 2d.
    Branchids can be specified for 1d branches that need to be linked.
    A (Multi)Polygon can be provided were links should be made.

    Note: The boundary nodes of Mesh1d (those sharing only one Mesh1d edge) are not connected to any Mesh2d face.

    Args:
        network (Network): Network in which the connections are made
        branchids (List[str], optional): List of branchid's to connect. If None, all branches are connected. Defaults to None.
        within (Union[Polygon, MultiPolygon], optional): Area within which connections are made. Defaults to None.
        max_length (float, optional): Max edge length. Defaults to None.
    """
    # Load 1d and 2d in meshkernel
    network._mesh1d._set_mesh1d()
    network._mesh2d._set_mesh2d()

    if within is None:
        # If not provided, create a box from the maximum bounds
        xmin = min(
            network._mesh1d.mesh1d_node_x.min(), network._mesh2d.mesh2d_node_x.min()
        )
        xmax = max(
            network._mesh1d.mesh1d_node_x.max(), network._mesh2d.mesh2d_node_x.max()
        )
        ymin = min(
            network._mesh1d.mesh1d_node_y.min(), network._mesh2d.mesh2d_node_y.min()
        )
        ymax = max(
            network._mesh1d.mesh1d_node_y.max(), network._mesh2d.mesh2d_node_y.max()
        )

        within = box(xmin, ymin, xmax, ymax)

    # If a 'within' polygon was provided, convert it to a geometrylist
    geometrylist = polygon_to_geometrylist(within)

    # Get the nodes for the specific branch ids
    node_mask = network._mesh1d.get_node_mask(branchids)

    # Get the already present links. These are not filtered on length
    npresent = len(network._link1d2d.link1d2d)

    # Generate links
    network._link1d2d._link_from_1d_to_2d(node_mask, polygon=geometrylist)

    # Filter the links that are longer than the max distance
    id1d = network._link1d2d.link1d2d[npresent:, 0]
    id2d = network._link1d2d.link1d2d[npresent:, 1]
    nodes1d = np.stack(
        [network._mesh1d.mesh1d_node_x[id1d], network._mesh1d.mesh1d_node_y[id1d]],
        axis=1,
    )
    faces2d = np.stack(
        [network._mesh2d.mesh2d_face_x[id2d], network._mesh2d.mesh2d_face_y[id2d]],
        axis=1,
    )
    lengths = np.hypot(nodes1d[:, 0] - faces2d[:, 0], nodes1d[:, 1] - faces2d[:, 1])
    keep = np.concatenate(
        [np.arange(npresent), np.where(lengths < max_length)[0] + npresent]
    )
    _filter_links_on_idx(network, keep)


def links1d2d_add_links_1d_to_2d_include_boundary(
    network: Network,
    branchids: List[str] = None,
    within: Union[Polygon, MultiPolygon] = None,
    max_length: float = np.inf,
) -> None:
    """Function to add 1d2d links to network, by generating them from 1d to 2d.
    Branchids can be specified for 1d branches that need to be linked.
    A (Multi)Polygon can be provided were links should be made.
    Modified from links1d2d_add_links_1d_to_2d to include also boundary locations.

    Note: The boundary nodes of Mesh1d (those sharing only one Mesh1d edge) are also connected to Mesh2d face.

    Args:
        network (Network): Network in which the connections are made
        branchids (List[str], optional): List of branchid's to connect. If None, all branches are connected. Defaults to None.
        within (Union[Polygon, MultiPolygon], optional): Area within which connections are made. Defaults to None.
        max_length (float, optional): Max edge length. Defaults to None.

    See Also
    --------
        links1d2d_add_links_1d_to_2d
    """
    # Load 1d and 2d in meshkernel
    network._mesh1d._set_mesh1d()
    network._mesh2d._set_mesh2d()

    if within is None:
        # If not provided, create a box from the maximum bounds
        xmin = min(
            network._mesh1d.mesh1d_node_x.min(), network._mesh2d.mesh2d_node_x.min()
        )
        xmax = max(
            network._mesh1d.mesh1d_node_x.max(), network._mesh2d.mesh2d_node_x.max()
        )
        ymin = min(
            network._mesh1d.mesh1d_node_y.min(), network._mesh2d.mesh2d_node_y.min()
        )
        ymax = max(
            network._mesh1d.mesh1d_node_y.max(), network._mesh2d.mesh2d_node_y.max()
        )

        within = box(xmin, ymin, xmax, ymax)

    # If a 'within' polygon was provided, convert it to a geometrylist
    geometrylist = polygon_to_geometrylist(within)

    # Get the nodes for the specific branch ids
    node_mask = network._mesh1d.get_node_mask(branchids)

    # Get the already present links. These are not filtered on length
    npresent = len(network._link1d2d.link1d2d)

    # Generate links
    network._link1d2d._link_from_1d_to_2d(node_mask, polygon=geometrylist)

    # generate 1d2d links #FIXME does not work yet
    network._link1d2d.meshkernel.contacts_compute_boundary(
        node_mask=node_mask, polygons=geometrylist, search_radius=max_length * 10
    )
    network._link1d2d._process()

    # Filter the links that are longer than the max distance
    id1d = network._link1d2d.link1d2d[npresent:, 0]
    id2d = network._link1d2d.link1d2d[npresent:, 1]
    nodes1d = np.stack(
        [network._mesh1d.mesh1d_node_x[id1d], network._mesh1d.mesh1d_node_y[id1d]],
        axis=1,
    )
    faces2d = np.stack(
        [network._mesh2d.mesh2d_face_x[id2d], network._mesh2d.mesh2d_face_y[id2d]],
        axis=1,
    )
    lengths = np.hypot(nodes1d[:, 0] - faces2d[:, 0], nodes1d[:, 1] - faces2d[:, 1])
    keep = np.concatenate(
        [np.arange(npresent), np.where(lengths < max_length)[0] + npresent]
    )
    _filter_links_on_idx(network, keep)


def _filter_links_on_idx(network: Network, keep: np.ndarray) -> None:
    # Select the remaining links
    network._link1d2d.link1d2d = network._link1d2d.link1d2d[keep]
    network._link1d2d.link1d2d_contact_type = network._link1d2d.link1d2d_contact_type[
        keep
    ]
    network._link1d2d.link1d2d_id = network._link1d2d.link1d2d_id[keep]
    network._link1d2d.link1d2d_long_name = network._link1d2d.link1d2d_long_name[keep]


def links1d2d_add_links_2d_to_1d_embedded(
    network: Network,
    branchids: List[str] = None,
    within: Union[Polygon, MultiPolygon] = None,
) -> None:
    """Generates links from 2d to 1d, where the 2d mesh intersects the 1d mesh: the 'embedded' links.

    To find the intersecting cells in an efficient way, we follow we the next steps. 1) Get the
    maximum length of a face edge. 2) Buffer the branches with this length. 3) Find all face nodes
    within this buffered geometry. 4) Check for each of the corresponding faces if it crossed the
    branches.

    Args:
        network (Network): Network in which the links are made. Should contain a 1d and 2d mesh
        branchids (List[str], optional): List is branch id's for which the connections are made. Defaults to None.
        within (Union[Polygon, MultiPolygon], optional): Clipping polygon for 2d mesh that is. Defaults to None.

    """
    # Load 1d and 2d in meshkernel
    network._mesh1d._set_mesh1d()
    network._mesh2d._set_mesh2d()

    # Get the max edge distance
    nodes2d = np.stack(
        [network._mesh2d.mesh2d_node_x, network._mesh2d.mesh2d_node_y], axis=1
    )
    edge_node_crds = nodes2d[network._mesh2d.mesh2d_edge_nodes]

    diff = edge_node_crds[:, 0, :] - edge_node_crds[:, 1, :]
    maxdiff = np.hypot(diff[:, 0], diff[:, 1]).max()

    # Create multilinestring from branches
    # branchnrs = np.unique(network._mesh1d.mesh1d_node_branch_id)
    nodes1d = np.stack(
        [network._mesh1d.mesh1d_node_x, network._mesh1d.mesh1d_node_y], axis=1
    )

    # Create a prepared multilinestring of the 1d network, to check for intersections
    mls = MultiLineString(nodes1d[network._mesh1d.mesh1d_edge_nodes].tolist())
    mls_prep = prep(mls)

    # Buffer the branches with the max cell distances
    area = mls.buffer(maxdiff)

    # If a within polygon is provided, clip the buffered area with this polygon.
    if within is not None:
        area = area.intersection(within)

    # Create an array with 2d facecenters and check which intersect the (clipped) area
    faces2d = np.stack(
        [network._mesh2d.mesh2d_face_x, network._mesh2d.mesh2d_face_y], axis=1
    )
    mpgl = GeometryList(*faces2d.T.copy())
    idx = np.zeros(len(faces2d), dtype=bool)
    if isinstance(area, MultiPolygon):
        area = [a for a in area]
    else:
        area = [area]
    for subarea in area:
        subarea = polygon_to_geometrylist(subarea)
        idx |= (
            network.meshkernel.polygon_get_included_points(subarea, mpgl).values == 1.0
        )

    # Check for each of the remaining faces, if it actually crosses the branches
    nodes2d = np.stack(
        [network._mesh2d.mesh2d_node_x, network._mesh2d.mesh2d_node_y], axis=1
    )
    where = np.where(idx)[0]
    for i, face_crds in enumerate(nodes2d[network._mesh2d.mesh2d_face_nodes[idx]]):
        if not mls_prep.intersects(LineString(face_crds)):
            idx[where[i]] = False

    # Use the remaining points to create the links
    multipoint = GeometryList(
        x_coordinates=faces2d[idx, 0], y_coordinates=faces2d[idx, 1]
    )

    # Get the nodes for the specific branch ids
    node_mask = network._mesh1d.get_node_mask(branchids)

    # Generate links
    network._link1d2d._link_from_2d_to_1d_embedded(node_mask, points=multipoint)


def links1d2d_add_links_2d_to_1d_lateral(
    network: Network,
    dist_factor: Union[float, None] = 2.0,
    branchids: List[str] = None,
    within: Union[Polygon, MultiPolygon] = None,
    max_length: float = np.inf,
) -> None:
    """Generate 1d2d links from the 2d mesh to the 1d mesh, with a lateral connection.
    If a link is kept, is determined based on the distance between the face center and
    the intersection with the 2d mesh exterior. By default, links with an intersection
    distance larger than 2 times the center to edge distance of the cell, are removed.
    Note that for a square cell with a direct link out of the cell (without passing any
    other cells) this max distance is sqrt(2) = 1.414. The default value of 2 provides
    some flexibility. Note that a link with more than 1 intersection with the 2d mesh
    boundary is removed anyway.

    Furthermore:
    - Branch ids can be specified to connect only specific branches.
    - A 'within' polygon can be given to only connect 2d cells within this polygon.
    - A max link length can be given to limit the link length.

    Args:
        network (Network): Network in which the links are made. Should contain a 1d and 2d mesh
        dist_factor (Union[float, None], optional): Factor to determine which links are kept (see description above). Defaults to 2.0.
        branchids (List[str], optional): List is branch id's for which the conncetions are made. Defaults to None.
        within (Union[Polygon, MultiPolygon], optional): Clipping polygon for 2d mesh that is. Defaults to None.
        max_length (float, optional): Max edge length. Defaults to None.
    """
    # Load 1d and 2d in meshkernel
    network._mesh1d._set_mesh1d()
    network._mesh2d._set_mesh2d()

    geometrylist = network.meshkernel.mesh2d_get_mesh_boundaries_as_polygons()
    mpboundaries = GeometryList(**geometrylist.__dict__).to_geometry()
    if within is not None:
        # If a 'within' polygon was provided, get the intersection with the meshboundaries
        # and convert it to a geometrylist
        # Note that the provided meshboundaries is a (list of) polygon(s). Holes are provided
        # as polygons as well, which dont make it a valid MultiPolygon
        geometrylist = GeometryList.from_geometry(
            MultiPolygon([geom.intersection(within) for geom in mpboundaries.geoms])
        )

    # Get the nodes for the specific branch ids
    node_mask = network._mesh1d.get_node_mask(branchids)

    # Get the already present links. These are not filtered subsequently
    npresent = len(network._link1d2d.link1d2d)

    # Generate links
    network._link1d2d._link_from_2d_to_1d_lateral(
        node_mask, polygon=geometrylist, search_radius=max_length
    )

    # If the provided distance factor was None, no further selection is needed, all links are kept.
    if dist_factor is None:
        return

    # Create multilinestring
    multilinestring = MultiLineString([poly.exterior for poly in mpboundaries.geoms])

    # Find the links that intersect the boundary close to the origin
    id1d = network._link1d2d.link1d2d[npresent:, 0]
    id2d = network._link1d2d.link1d2d[npresent:, 1]

    nodes1d = np.stack(
        [network._mesh1d.mesh1d_node_x[id1d], network._mesh1d.mesh1d_node_y[id1d]],
        axis=1,
    )
    faces2d = np.stack(
        [network._mesh2d.mesh2d_face_x[id2d], network._mesh2d.mesh2d_face_y[id2d]],
        axis=1,
    )
    nodes2d = np.stack(
        [network._mesh2d.mesh2d_node_x, network._mesh2d.mesh2d_node_y], axis=1
    )
    face_node_crds = nodes2d[network._mesh2d.mesh2d_face_nodes[id2d]]

    # Calculate distance between face edge and face center
    x1 = np.take(face_node_crds, 0, axis=2)
    y1 = np.take(face_node_crds, 1, axis=2)
    face_node_crds[:] = np.roll(face_node_crds, 1, axis=1)
    x2 = np.take(face_node_crds, 0, axis=2)
    y2 = np.take(face_node_crds, 1, axis=2)
    x0, y0 = faces2d[:, 0], faces2d[:, 1]
    distance = (
        np.absolute((x2 - x1) * (y1 - y0[:, None]) - (x1 - x0[:, None]) * (y2 - y1))
        / np.hypot(x2 - x1, y2 - y1)
    ).mean(axis=1)

    # Check which links to keep
    keep = list(range(npresent))
    for i, (node1d, face2d, comp_dist) in enumerate(
        zip(nodes1d, faces2d, distance * dist_factor)
    ):
        isect = multilinestring.intersection(LineString([face2d, node1d]))

        # If the intersection is for some reason not a Point of Multipoint, skip it.
        if not isinstance(isect, (Point, MultiPoint)):
            continue

        # Skip the link if it has more than one intersection with the boundary
        if isinstance(isect, MultiPoint):
            isect_list = [a for a in isect]
        else:
            isect_list = [isect]
        if len(isect_list) != 1:
            continue

        # If the distance to the mesh 2d exterior intersection is smaller than
        # the compared distance, keep it.
        dist = np.hypot(*(face2d - isect_list[0]))
        if dist < comp_dist:
            keep.append(i)

    # Select the remaining links
    network._link1d2d.link1d2d = network._link1d2d.link1d2d[keep]
    network._link1d2d.link1d2d_contact_type = network._link1d2d.link1d2d_contact_type[
        keep
    ]
    network._link1d2d.link1d2d_id = network._link1d2d.link1d2d_id[keep]
    network._link1d2d.link1d2d_long_name = network._link1d2d.link1d2d_long_name[keep]
