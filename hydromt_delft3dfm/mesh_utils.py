# -*- coding: utf-8 -*-

import logging
from typing import Tuple

import geopandas as gpd
import xarray as xr
import xugrid as xu
from hydrolib.core.dflowfm import Network
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
    _mesh1d_attrs = dfm_network._mesh1d.__dict__.keys()

    # add mesh2d
    if "mesh2d" in grids:
        dfm_network._mesh2d._process(
            grids["mesh2d"].mesh
        )  # FIXME: test what if the mesh had bedlevel variable

    # add mesh1d (including mesh1d and networkd1d)
    if "mesh1d" in grids:
        mesh1d = mesh.ugrid.to_dataset(
            optional_attributes=True
        )  # optional_attributes includes edge_x and edge_y
        for var, val in mesh1d.variables.items():
            if var in _mesh1d_attrs:
                # use hydrolib-core conventions as it does harmonization when reading.
                # check conventions at hydrolib.core.dflowfm.net.ugrid_conventions.json
                setattr(dfm_network._mesh1d, var, val.values)
        # process
        dfm_network._mesh1d._process_network1d()
        dfm_network._mesh1d._set_mesh1d()

    # add 1d2dlinks
    _link1d2d_attrs = dfm_network._link1d2d.__dict__.keys()
    if "link1d2d" in mesh:
        for var, val in mesh.variables.items():
            if var in _link1d2d_attrs:
                # use hydrolib-core conventions as it does harmonization when reading.
                setattr(dfm_network._link1d2d, var, val.values)

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
        # get grid
        grid_mesh1d = xu.Ugrid1d(
            node_x=mesh1d.mesh1d_node_x,
            node_y=mesh1d.mesh1d_node_y,
            fill_value=-1,
            edge_node_connectivity=mesh1d.mesh1d_edge_nodes,
            name="mesh1d",
            projected=crs.is_projected,
            crs=crs,
        )
        grid_mesh1d.set_crs(crs)
        edge_dim = grid_mesh1d.edge_dimension
        node_dim = grid_mesh1d.node_dimension

        # get data
        ds = xr.Dataset()
        ds["mesh1d_node_id"] = (node_dim, mesh1d.mesh1d_node_id)
        ds["mesh1d_node_long_name"] = (node_dim, mesh1d.mesh1d_node_long_name)
        ds["mesh1d_node_branch_id"] = (node_dim, mesh1d.mesh1d_node_branch_id)
        ds["mesh1d_node_branch_offset"] = (node_dim, mesh1d.mesh1d_node_branch_offset)
        ds["mesh1d_edge_branch_id"] = (edge_dim, mesh1d.mesh1d_edge_branch_id)
        ds["mesh1d_edge_branch_offset"] = (edge_dim, mesh1d.mesh1d_edge_branch_offset)

        # get ugrid dataset
        uds_mesh1d = xu.UgridDataset(ds, grids=grid_mesh1d)

        # derive network1d
        # The 1D network topology serves as the coordinate space in which a 1D mesh discretization
        # will later be defined. The network is largely based on the UGRID conventions for its topology
        # (i.e., nodes and edges) and additionally uses an optional edge_geometry to define the
        # precise network branch geometries (more about this in the next Section).

        grid_network1d = xu.Ugrid1d(
            node_x=mesh1d.network1d_node_x,
            node_y=mesh1d.network1d_node_y,
            fill_value=-1,
            edge_node_connectivity=mesh1d.network1d_edge_nodes,
            name="network1d",
            projected=crs.is_projected,
            crs=crs,
        )
        grid_network1d.set_crs(crs)
        network_edge_dim = grid_network1d.edge_dimension
        network_node_dim = grid_network1d.node_dimension

        # get data
        ds_network1d = xr.Dataset()
        ds_network1d["network1d_node_id"] = (network_node_dim, mesh1d.network1d_node_id)
        ds_network1d["network1d_node_long_name"] = (
            network_node_dim,
            mesh1d.network1d_node_long_name,
        )
        ds_network1d["network1d_branch_id"] = (
            network_edge_dim,
            mesh1d.network1d_branch_id,
        )
        ds_network1d["network1d_branch_long_name"] = (
            network_edge_dim,
            mesh1d.network1d_branch_long_name,
        )
        ds_network1d["network1d_branch_length"] = (
            network_edge_dim,
            mesh1d.network1d_branch_length,  # real length of the branches
        )
        # network1d_geometry related
        ds_network1d["network1d_part_node_count"] = (
            network_edge_dim,
            mesh1d.network1d_part_node_count,
        )
        ds_network1d["network1d_geom_x"] = (
            "network1d_nGeometryNodes",
            mesh1d.network1d_geom_x,
        )
        ds_network1d["network1d_geom_y"] = (
            "network1d_nGeometryNodes",
            mesh1d.network1d_geom_y,
        )
        ds_network1d["network1d_branch_order"] = (
            network_edge_dim,
            mesh1d.network1d_branch_order,
        )
        # might be supported in the future https://github.com/Deltares/HYDROLIB-core/issues/561
        # ds["network1d_branch_type"] = (
        #     edge_dim,
        #     mesh1d.network1d_branch_type,
        # )

        # get ugrid dataset
        uds_network1d = xu.UgridDataset(ds_network1d, grids=grid_network1d)

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
    link1d2d["link1d2d_id"] = ("nLink1D2D_edge", network._link1d2d.link1d2d_id)
    link1d2d["link1d2d_long_name"] = (
        "nLink1D2D_edge",
        network._link1d2d.link1d2d_long_name,
    )
    link1d2d["link1d2d_contact_type"] = (
        "nLink1D2D_edge",
        network._link1d2d.link1d2d_contact_type,
    )

    return link1d2d


def mesh2d_from_hydrolib_network(
    network: Network,
    crs: CRS,
) -> xu.UgridDataset:
    """
    Creates xugrid mesh2d UgridDataset from hydrolib-core network object.

    Parameters
    ----------
    network : Network
        Network hydrolib-core object.
    crs : pyproj.CRS
        Coordinate reference system of the network.

    Returns
    -------
    uds_mesh2d : xu.UgridDataset
        Mesh2d UgridDataset.
    """
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

    # set crs
    uds_mesh2d.ugrid.set_crs(crs)
    return uds_mesh2d


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
    if not network._mesh1d.is_empty():
        uds_mesh1d, uds_network1d = mesh1d_network1d_from_hydrolib_network(network, crs)

        # add to mesh
        mesh = xu.UgridDataset(
            xr.merge([uds_mesh1d.ugrid.to_dataset(), uds_network1d.ugrid.to_dataset()])
        )

    # Mesh2d
    if not network._mesh2d.is_empty():
        # network._mesh2d._set_mesh2d()
        uds_mesh2d = mesh2d_from_hydrolib_network(network, crs)

        if mesh is None:
            mesh = uds_mesh2d
        else:
            mesh = xu.UgridDataset(
                xr.merge(
                    [
                        mesh.ugrid.to_dataset(optional_attributes=True),
                        uds_mesh2d.ugrid.to_dataset(optional_attributes=True),
                    ]
                )
            )

    # 1d2dlinks
    if not network._link1d2d.is_empty():
        link1d2d = links1d2d_from_hydrolib_network(network)
        # Add to mesh (links should only exist if mesh1d and mesh2d exist)
        if mesh is not None:
            for v in link1d2d.data_vars:
                mesh[v] = link1d2d[v]

    # Set crs
    if mesh is not None:
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
    # Creates gdf with some extra data (using hydrolib-core convention)
    mesh1d_nodes = gpd.GeoDataFrame(
        data={
            "branch_id": uds_mesh1d["mesh1d_node_branch_id"],
            "branch_name": [
                branches.branchid[i] for i in uds_mesh1d["mesh1d_node_branch_id"].values
            ],
            "branch_chainage": uds_mesh1d["mesh1d_node_branch_offset"],
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
