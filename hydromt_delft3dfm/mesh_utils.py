# -*- coding: utf-8 -*-

import logging
from typing import List, Union, Tuple

import numpy as np
import xarray as xr
import geopandas as gpd
import xugrid as xu
from hydrolib.core.dflowfm import Branch, Network
import meshkernel as mk
from pyproj import CRS

logger = logging.getLogger(__name__)


__all__ = [
    "hydrolib_network_from_mesh",
    "mesh1d_network1d_from_hydrolib_network",
    "links1d2d_from_hydrolib_network",
    "mesh_from_hydrolib_network",
    "mesh1d_nodes_geodataframe",
    "network1d_nodes_geodataframe",
]


def hydrolib_network_from_mesh(
    mesh: xu.UgridDataset,
) -> Network:
    """
    Converts from xugrid mesh to hydrolib-core network object.

    Parameters
    ----------
    mesh : xu.UgridDataset
        Mesh UgridDataset.

    Returns
    -------
    network : Network
        Network object.
    """
    # split grids
    grids = dict()
    for grid in mesh.ugrid.grids:
        grids[grid.name] = grid

    # create network
    dfm_network = Network(is_geographic=mesh.ugrid.grids[0].crs.is_geographic)

    # add mesh2d
    if "mesh2d" in grids:
        dfm_network._mesh2d._process(grids["mesh2d"].mesh)

    # add mesh1d and networkd1d
    if "mesh1d" in grids:
        # set hydrolib Mesh1d atribute one by one and then add to network
        # mesh1d and network1d variables
        mesh1d_dict = {
            "mesh1d_node_id": "mesh1d_node_id",
            "mesh1d_node_long_name": "mesh1d_node_long_name",
            "mesh1d_node_x": "mesh1d_node_x",
            "mesh1d_node_y": "mesh1d_node_y",
            "mesh1d_node_branch_id": "mesh1d_node_branch",
            "mesh1d_node_branch_offset": "mesh1d_node_offset",
            "mesh1d_edge_nodes": "mesh1d_edge_nodes",
            "mesh1d_edge_x": "mesh1d_edge_x",
            "mesh1d_edge_y": "mesh1d_edge_y",
            "mesh1d_edge_branch_id": "mesh1d_edge_branch",
            "mesh1d_edge_branch_offset": "mesh1d_edge_offset",
            "network1d_node_id": "network1d_node_id",
            "network1d_node_long_name": "network1d_node_long_name",
            "network1d_node_x": "network1d_node_x",
            "network1d_node_y": "network1d_node_y",
            "network1d_branch_id": "network1d_branch_id",
            "network1d_branch_long_name": "network1d_branch_long_name",
            "network1d_branch_length": "network1d_edge_length",
            "network1d_branch_order": "network1d_branch_order",
            "network1d_edge_nodes": "network1d_edge_nodes",
            "network1d_part_node_count": "network1d_geom_node_count",
            "network1d_geom_x": "network1d_geom_x",
            "network1d_geom_y": "network1d_geom_y",
        }
        for hydrolibkey, meshkey in mesh1d_dict.items():
            # Key in the UgridDataset
            if meshkey in mesh:
                setattr(dfm_network._mesh1d, hydrolibkey, mesh[meshkey].values)
            # Key in Ugrid1D
            elif meshkey in grids["mesh1d"].to_dataset():
                setattr(
                    dfm_network._mesh1d,
                    hydrolibkey,
                    grids["mesh1d"].to_dataset()[meshkey].values,
                )
            elif meshkey in grids["network1d"].to_dataset():
                setattr(
                    dfm_network._mesh1d,
                    hydrolibkey,
                    grids["network1d"].to_dataset()[meshkey].values,
                )

        # process
        dfm_network._mesh1d._process_network1d()

    # add 1d2dlinks
    if "link1d2d" in mesh:
        link1d2d_dict = {
            "link1d2d": "link1d2d",
            "link1d2d_id": "link1d2d_ids",
            "link1d2d_long_name": "link1d2d_long_names",
            "link1d2d_contact_type": "link1d2d_contact_type",
        }
        for hydrolibkey, meshkey in link1d2d_dict.items():
            if meshkey in mesh:
                setattr(dfm_network._link1d2d, hydrolibkey, mesh[meshkey].values)

    return dfm_network


def mesh1d_network1d_from_hydrolib_network(
    network: Network,
    crs: CRS,
) -> Tuple[xu.UgridDataset, xu.UgridDataset]:
    """
    Creates xugrid mesh1d and network1d UgridDataset from hydrolib-core network object.

    Parameters
    ----------
    network : Network
        Network hydrolib-core object.
    crs : pyproj.CRS
        Coordinate reference system of the network.

    Returns
    -------
    uds_mesh1d : xu.UgridDataset
        Mesh1d UgridDataset.
    uds_network1d : xu.UgridDataset
        Network1d UgridDataset.
    """
    # Mesh1d to mesh1d and network1d xugrid
    mesh1d = network._mesh1d

    if not mesh1d.is_empty():
        uds_mesh1d = xu.Ugrid1d(
            node_x=mesh1d.mesh1d_node_x,
            node_y=mesh1d.mesh1d_node_y,
            fill_value=-1,
            edge_node_connectivity=mesh1d.mesh1d_edge_nodes,
            name="mesh1d",
            projected=crs.is_projected,
            crs=crs,
        )
        # Convert to UgridDataset
        uds_mesh1d = xu.UgridDataset(uds_mesh1d.to_dataset()).ugrid.assign_edge_coords()
        uds_mesh1d.ugrid.set_crs(crs)
        # Add extra properties
        edge_dim = uds_mesh1d.ugrid.grid.edge_dimension
        node_dim = uds_mesh1d.ugrid.grid.node_dimension
        uds_mesh1d["mesh1d_node_id"] = (node_dim, mesh1d.mesh1d_node_id)
        uds_mesh1d["mesh1d_node_long_name"] = (node_dim, mesh1d.mesh1d_node_long_name)
        uds_mesh1d["mesh1d_node_branch"] = (node_dim, mesh1d.mesh1d_node_branch_id)
        uds_mesh1d["mesh1d_node_offset"] = (node_dim, mesh1d.mesh1d_node_branch_offset)
        uds_mesh1d["mesh1d_edge_branch"] = (edge_dim, mesh1d.mesh1d_edge_branch_id)
        uds_mesh1d["mesh1d_edge_offset"] = (edge_dim, mesh1d.mesh1d_edge_branch_offset)

        # derive network1d
        uds_network1d = xu.Ugrid1d(
            node_x=mesh1d.network1d_node_x,
            node_y=mesh1d.network1d_node_y,
            fill_value=-1,
            edge_node_connectivity=mesh1d.network1d_edge_nodes,
            name="network1d",
            projected=crs.is_projected,
            crs=crs,
        )
        # Convert to UgridDataset
        uds_network1d = xu.UgridDataset(
            uds_network1d.to_dataset()
        ).ugrid.assign_edge_coords()
        uds_network1d.ugrid.set_crs(crs)
        # Add extra properties
        edge_dim = uds_network1d.ugrid.grid.edge_dimension
        node_dim = uds_network1d.ugrid.grid.node_dimension
        uds_network1d["network1d_node_id"] = (node_dim, mesh1d.network1d_node_id)
        uds_network1d["network1d_node_long_name"] = (
            node_dim,
            mesh1d.network1d_node_long_name,
        )
        uds_network1d["network1d_branch_id"] = (edge_dim, mesh1d.network1d_branch_id)
        uds_network1d["network1d_branch_long_name"] = (
            edge_dim,
            mesh1d.network1d_branch_long_name,
        )
        uds_network1d["network1d_edge_length"] = (
            edge_dim,
            mesh1d.network1d_branch_length,
        )
        uds_network1d["network1d_branch_order"] = (
            edge_dim,
            mesh1d.network1d_branch_order,
        )
        uds_network1d["network1d_geom_node_count"] = (
            edge_dim,
            mesh1d.network1d_part_node_count,
        )
        uds_network1d["network1d_geom_x"] = (
            "network1d_nGeometryNodes",
            mesh1d.network1d_geom_x,
        )
        uds_network1d["network1d_geom_y"] = (
            "network1d_nGeometryNodes",
            mesh1d.network1d_geom_y,
        )

    else:
        uds_mesh1d = None
        uds_network1d = None

    return uds_mesh1d, uds_network1d


def links1d2d_from_hydrolib_network(
    network: Network,
) -> xr.Dataset():
    """
    Extract link1d2d from hydrolib-core network object.

    Parameters
    ----------
    network : Network
        Network hydrolib-core object.

    Returns
    -------
    link1d2d : xr.Dataset
        Link1d2d Dataset.
    """
    # extract links from network object
    link1d2d = xr.Dataset()
    link1d2d["link1d2d"] = (["nLink1D2D_edge", "Two"], network._link1d2d.link1d2d)
    # extra variables
    link1d2d["link1d2d_ids"] = ("nLink1D2D_edge", network._link1d2d.link1d2d_id)
    link1d2d["link1d2d_long_names"] = (
        "nLink1D2D_edge",
        network._link1d2d.link1d2d_long_name,
    )
    link1d2d["link1d2d_contact_type"] = (
        "nLink1D2D_edge",
        network._link1d2d.link1d2d_contact_type,
    )

    return link1d2d


def mesh_from_hydrolib_network(
    network: Network,
    crs: CRS,
) -> xu.UgridDataset:
    """
    Creates xugrid mesh from hydrolib-core network object.

    Parameters
    ----------
    network : Network
        Network hydrolib-core object.
    crs : pyproj.CRS
        Coordinate reference system of the network.

    Returns
    -------
    mesh : xu.UgridDataset
        Mesh UgridDataset.
    """
    # initialise mesh
    mesh = None

    # Mesh1d to mesh1d and network1d xugrid
    mesh1d = network._mesh1d
    if not mesh1d.is_empty():
        uds_mesh1d, uds_network1d = mesh1d_network1d_from_hydrolib_network(network, crs)

        # add to mesh
        mesh = xu.UgridDataset(
            xr.merge([uds_mesh1d.ugrid.to_dataset(), uds_network1d.ugrid.to_dataset()])
        )

    # Mesh2d
    if not network._mesh2d.is_empty():
        # network._mesh2d._set_mesh2d()
        mesh2d = network._mesh2d

        # meshkernel to xugrid Ugrid2D
        uds_mesh2d = xu.Ugrid2d(
            node_x=mesh2d.mesh2d_node_x,
            node_y=mesh2d.mesh2d_node_y,
            fill_value=-1,
            face_node_connectivity=mesh2d.mesh2d_face_nodes,
            edge_node_connectivity=mesh2d.mesh2d_edge_nodes,
            name="mesh2d",
            projected=crs.is_projected,
            crs=crs,
        )
        # Convert to UgridDataset
        uds_mesh2d = xu.UgridDataset(uds_mesh2d.to_dataset())
        uds_mesh2d = uds_mesh2d.ugrid.assign_face_coords()
        uds_mesh2d.ugrid.set_crs(crs)

        if mesh is None:
            mesh = uds_mesh2d
        else:
            mesh = xu.UgridDataset(
                xr.merge([mesh.ugrid.to_dataset(), uds_mesh2d.ugrid.to_dataset()])
            )

    # 1d2dlinks
    if not network._link1d2d.is_empty():
        link1d2d = links1d2d_from_hydrolib_network(network)
        # Add to mesh (links should only exist if mesh1d and mesh2d exist)
        for v in link1d2d.data_vars:
            mesh[v] = link1d2d[v]

    # Set crs
    for grid in mesh.ugrid.grids:
        grid.set_crs(crs)

    return mesh


def mesh1d_nodes_geodataframe(
    uds_mesh1d: xu.UgridDataset,
    branches: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Returns the nodes of mesh 1D as geodataframe.

    Parameters
    ----------
    uds_mesh1d : xu.UgridDataset
        Mesh1d UgridDataset including mesh1d data_vars.
    branches : gpd.GeoDataFrame
        Branches GeoDataFrame.

    Returns
    -------
    mesh1d_nodes : gpd.GeoDataFrame
        Mesh1d nodes GeoDataFrame.
    """
    # create node points
    mesh1d_nodes = gpd.points_from_xy(
        x=uds_mesh1d.ugrid.grid.node_x,
        y=uds_mesh1d.ugrid.grid.node_y,
        crs=uds_mesh1d.ugrid.grid.crs,
    )
    # Creates gdf with some extra data
    mesh1d_nodes = gpd.GeoDataFrame(
        data={
            "branch_id": uds_mesh1d["mesh1d_node_branch"],
            "branch_name": [
                branches.branchid[i] for i in uds_mesh1d["mesh1d_node_branch"].values
            ],
            "branch_chainage": uds_mesh1d["mesh1d_node_offset"],
            "geometry": mesh1d_nodes,
        }
    )

    return mesh1d_nodes


def network1d_nodes_geodataframe(
    uds_network1d: xu.UgridDataset,
) -> gpd.GeoDataFrame:
    """
    Get network1d nodes as gdp.

    Parameters
    ----------
    uds_network1d : xu.UgridDataset
        Network1d UgridDataset including network1d data_vars.

    Returns
    -------
    network1d_nodes : gpd.GeoDataFrame
        Network1d nodes GeoDataFrame.
    """
    # get networkids to complete the boundaries
    network1d_nodes = gpd.points_from_xy(
        x=uds_network1d.ugrid.grid.node_x,
        y=uds_network1d.ugrid.grid.node_y,
        crs=uds_network1d.ugrid.grid.crs,
    )
    network1d_nodes = gpd.GeoDataFrame(
        data={
            "nodeid": uds_network1d["network1d_node_id"],
            "geometry": network1d_nodes,
        }
    )

    return network1d_nodes
