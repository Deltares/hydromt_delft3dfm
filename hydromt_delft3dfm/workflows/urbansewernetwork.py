# -*- coding: utf-8 -*-

import logging
import tempfile
from pathlib import Path
from typing import Union, Optional, List, Dict

import pandas as pd
import geopandas as gpd
import networkx as nx
import numpy as np
import osmnx
import pyproj
import xarray as xr

# hydromt
from hydromt import flw
from hydromt import DataCatalog
from hydromt.gis_utils import nearest_merge

from hydromt_delft3dfm import graph_utils
from hydromt_delft3dfm import workflows


logger = logging.getLogger(__name__)

__all__ = [
    "setup_urban_sewer_network_topology_from_osm",
    # graph additions
    # "setup_graph_nodes",
    # "setup_graph_egdes",
    # graph workflows
    "optimize_graph",
    "update_graph_from_dem",
    # TODO to be categorized
    "setup_network_connections_based_on_flowdirections",
    "setup_network_parameters_from_rasters",
    "setup_network_topology_optimization",
    "setup_network_dimentions_from_rainfallstats",
]


def setup_urban_sewer_network_topology_from_osm(
    region,
    highway_types: list[str] = None,
    waterway_types: list[str] = None,
    logger: logging.Logger = logger,
):
    """
    Set up a complete urban sewer network topology OpenStreetMap.

    Include data for waterways and highways for open and closed systems.
    The two systems are connected and outlets are added at intersection points.

    Adds/Updates model layers:

    * **rivers** geom: 1D rivers vector
    * **pipes** geom: 1D pipes vector
    * **branches** geom: 1D branches vector
    * **outlets** geom: 1D outlets vector

    Parameters:
    - region : object
        The geographical region for which the sewer network is to be constructed.
        This could be a bounding box, polygon, or any other representation of a region.
    - highway_types : list
        List of highway types (e.g., ["motorway", "primary"]) to include from the OpenStreetMap data.
    - waterway_types : list
        List of waterway types (e.g., ["river", "stream"]) to include from the OpenStreetMap data.
    - logger : Logger object
        An instance of a logger to capture logs and messages during the process.
        Useful for debugging and understanding the workflow steps.

    Returns:
    - graph : networkx.Graph
        The final preprocessed graph representing the urban sewer network.
        It integrates both the open and closed drainage systems and includes outlets at
        the intersection points.
    """
    # input
    if waterway_types is None:
        waterway_types = ["river", "stream", "brook", "canal", "ditch"]
    if highway_types is None:
        highway_types = [
            "motorway",
            "motorway_link",
            "primary",
            "primary_link",
            "secondary",
            "secondary_link",
            "tertiary",
            "tertiary_link",
            "residential",
        ]
    # output
    _required_columns = None

    # 1. Build the graph for waterways (open system)
    graph_osm_open = workflows.create_graph_from_openstreetmap(
        region=region,
        buffer=1000,
        osm_key="waterway",
        osm_values=waterway_types,
        logger=logger,
    )
    logger.info(f"Download waterway {waterway_types} form osm")

    # 2. Build the graph for highways (closed system)
    graph_osm_closed = workflows.create_graph_from_openstreetmap(
        region=region,
        buffer=500,
        osm_key="highway",
        osm_values=highway_types,
        logger=logger,
    )
    logger.info(f"Download highway {highway_types} form osm")

    # 3. obtain open and closed systems and outlets by intersecting
    open_system, closed_system, outlets = workflows.intersect_lines(
        graph_utils.graph_to_network(graph_osm_open)[0],
        graph_utils.graph_to_network(graph_osm_closed)[0],
    )
    # add attributes
    open_system["branchid"] = open_system["id"]
    open_system["branchtype"] = "river"
    closed_system["branchid"] = closed_system["id"]
    closed_system["branchtype"] = "pipe"
    outlets["nodeid"] = outlets["id"]
    outlets["nodetype"] = "outlet"

    # TODO select required columns

    # 4. Recreate the graph from the complete system topology
    graph = workflows.create_graph_from_geodataframe(
        pd.concat([open_system, closed_system])
    )
    branches = graph_utils.graph_to_network(graph)[0]

    # TODO set branches
    # TODO set graph

    # 5. Add outlets to graph
    graph = workflows.setup_graph_from_geodataframe(
        graph,
        data_catalog=None,  # TODO replace by self.data_catalog
        vector_fn=outlets,
        variables=["nodeid", "nodetype"],
        max_dist=1.0,
        rename=dict(nodeid="nodeid", nodetype="nodetype"),
        graph_component="nodes",
        logger=logger,
    )

    return graph


def optimize_graph(graph, logger=logger):
    """Optimize graph based on various methods

    if only elevtn is available at nodes, then the cost along path is used
    need to determin which nodes are sources, which nodes are targets


    """

    logger.info("compute gradient from elevtn")
    graph = _compute_gradient_from_elevtn(graph)

    logger.info("update graph direction")
    graph = reverse_edges_on_negative_weight(graph, weight="gradient")
    return graph


def _compute_gradient_from_elevtn(graph):
    """
    Compute the gradient for each edge in the graph based on node elevations and edge length.

    Gradient is calculated as:
    gradient = (elevation of source node - elevation of target node) / edge length

    Modifies the input graph in place.
    """
    for e in graph.edges:
        graph.edges[e].update(
            {
                "gradient": (
                    (graph.nodes[e[0]]["elevtn"] - graph.nodes[e[1]]["elevtn"])
                    / graph.edges[e]["length"]
                )
            }
        )
    return graph


def update_graph_from_dem(
    graph: nx.Graph,  # TODO replace by self.graphs
    data_catalog: DataCatalog,  # TODO replace by self.data_catalog
    dem_fn: Union[str, Path, xr.DataArray, xr.Dataset],
    fill_method: Optional[str] = None,
    logger: logging.Logger = logger,
):
    """
    Update the graph with elevation data from a DEM, compute gradient for edges,
    and reverse the direction of edges with negative gradient.

    The function samples elevation data from the provided DEM to the nodes of the graph.
    It then computes the gradient of edges based on the elevation difference between the nodes
    and the edge length. If the computed gradient is negative, the direction of the edge is reversed.

    Parameters
    ----------
    graph : nx.Graph
        Input graph to be updated.
    data_catalog : DataCatalog
        Data catalog to read and process raster data.
    dem_fn : Union[str, Path, xr.DataArray, xr.Dataset]
        Data catalog key, path to DEM file or DEM xarray data object.
    fill_method : str, optional
        If specified, fills no data values using fill_nodata method.
        Available methods are {'linear', 'nearest', 'cubic', 'rio_idw'}.
    logger : logging.Logger, optional
        Logger object to log messages. Default is the global logger.

    Returns
    -------
    nx.Graph
        Graph with updated node elevations, edge gradients, and potentially reversed edge directions.
    """
    logger.info(
        "Update the graph with elevation data from a DEM and compute gradient for edges"
    )

    # Sample DEM
    # TODO replaced by super method
    graph = setup_graph_from_rasterdataset(
        graph=graph,
        data_catalog=data_catalog,
        raster_fn=dem_fn,
        fill_method=fill_method,
        graph_component="nodes",
        logger=logger,
    )

    logger.info("compute gradient from elevtn")
    graph = _compute_gradient_from_elevtn(graph)

    logger.info("update graph direction")
    graph = reverse_edges_on_negative_weight(graph, weight="gradient")
    return graph


def setup_network_connections_based_on_flowdirections(
    user_input_graph: nx.graph,
    dem_fn: Union[str, Path],
    logger: logging.Logger = logger,
    **kwargs,
) -> nx.graph:
    """
    This method sets up connections in the urban drainage network graph based on flow directions derived from a Digital Elevation Model (DEM).

    The method performs the following steps:
    1. Retrieves a flow direction raster from a DEM.
    2. Converts the flow direction raster into a graph and a network, using the network to represent flow directions.
    3. Snaps nodes from the `user_input_graph` to the flow direction network to ensure that the network aligns with natural flow paths.
    4. Identifies and fills in any missing links in the `user_input_graph`, ensuring a fully connected network.
    5. Fills in any missing geometry in the `user_input_graph` and converts it into a network.

    Parameters
    ----------
    user_input_graph: nx.graph
        The initial graph provided by the user, which will be updated to align with the flow direction network.

    dem_fn: Union[str, Path]
        The file path to the digital elevation model (DEM) used to derive the flow direction raster. This can be provided as a string or a Path object.

    **kwargs:
        Other keyword arguments that may be required for specific implementations of the function.

    Returns
    -------
    nx.graph:
        The processed `user_input_graph`, now an instance of nx.graph class, with any missing links and geometry filled in based on flow directions. The graph represents the urban drainage network.

    Notes
    -----
    The function reprojects the DEM-derived flow network to the local CRS (Coordinate Reference System) of the `user_input_graph` to ensure that all spatial data aligns correctly.
    """
    # method implementation goes here
    pass


def setup_network_parameters_from_rasters(
    graph: nx.graph,
    dem_fn,
    landuse_fn,
    water_demand_fn,
    population_fn,
    building_footprint_fn,
    logger: logging.Logger = logger,
    **kwargs,
) -> nx.graph:
    """
    This method sets up physical parameters of the urban drainage network represented by the graph, using various raster data sources.

    The method performs the following steps:
    1. Adds upstream and downstream street levels to the graph using a digital elevation model (DEM).
    2. Adds upstream areas to the graph using land use data.
    3. Adds water demand, population, and building footprint data to the graph using corresponding raster data sources.
    4. Approximates the size of the network based on the street size/type using data from the graph itself. --> use as weight for the next step.
    5. Approximates the slope of the network based on upstream and downstream street levels, the depths, and the length of the edge geometry in the graph. --> use as weight for the next step.

    Parameters
    ----------
    graph: nx.graph
        The network graph to be filled with physical parameters. This is an instance of the nx.graph class, which wraps around the NetworkX Graph class.

    dem_fn: Union[str, Path]
        The file path to the digital elevation model (DEM) used to derive the upstream and downstream street levels.

    landuse_fn: Union[str, Path]
        The file path to the land use data used to derive the upstream area.

    water_demand_fn: Union[str, Path]
        The file path to the water demand data used to add water demand information to the graph.

    population_fn: Union[str, Path]
        The file path to the population data used to add population information to the graph.

    building_footprint_fn: Union[str, Path]
        The file path to the building footprint data used to add building footprint information to the graph.

    **kwargs:
         Additional keyword arguments for more specific implementations of the function.

    Returns
    -------
    nx.graph:
        The updated graph, an instance of nx.graph class, representing the urban drainage network with physical parameters filled in.
    """
    # method implementation goes here
    pass


def setup_network_topology_optimization(
    graph: nx.graph, method: str, logger: logging.Logger = logger, **kwargs
) -> nx.graph:
    """
    This method optimizes the topology of the urban drainage network represented by the graph.

    The method performs the following steps:
    1. Computes the selected graph metric (defined by 'method') for all nodes in the graph.
    2. Removes unnecessary links based on the computed graph metric.
    3. Iterates the above steps as needed until the network is optimized.

    Parameters
    ----------
    graph: nx.graph
        The network graph to be optimized. This is an instance of the nx.graph class, which wraps around the NetworkX Graph class.

    method: str
        The method to use for optimization. This determines the graph metric that will be computed and used for the optimization process.

    **kwargs:
        Additional keyword arguments for more specific implementations of the function.

    Returns
    -------
    nx.graph:
        The optimized graph, an instance of nx.graph class, representing the urban drainage network.

    Notes
    -----
    The betweenness centrality of a node in a graph is a measure of how often it appears on the shortest paths between nodes. Nodes with high betweenness centrality can be considered 'hubs' in the network. By removing links with low betweenness centrality, we can potentially simplify the network topology without significantly impacting the overall connectivity.

    The optimal way to optimize network topology will highly depend on the specific characteristics of the network, which are often determined by previous steps in the network analysis process.

    Other candidate graph metrics that might be useful for network topology optimization include:
    - Degree centrality: Measures the number of edges connected to a node. This can be useful for identifying nodes that are most connected.
    - Closeness centrality: Measures how close a node is to all other nodes in the network. This can be useful for identifying nodes that are centrally located.
    - Eigenvector centrality: Measures a node's influence based on the number of links it has to other influential nodes. This can be useful for identifying influential nodes in the network.

    If the network's weights have physical meanings, shortest path or maximum flow methods can also be considered for optimization. These methods, such as the Ford-Fulkerson or Edmonds-Karp algorithm, can be used to identify critical paths in the network. This can be particularly useful for identifying vulnerabilities in the network or for identifying potential areas for network enhancement.

    See Also
    --------
    https://github.com/xldeltares/hydrolib-nowcasting/tree/master

    """
    # method implementation goes here
    pass


def setup_network_dimentions_from_rainfallstats(
    graph: nx.graph,
    rainfall_fn: Union[str, Path],
    rainfall_assumption,
    capacity_assumption=None,
    diameter_range=None,
    velocity_range=None,
    slope_range=None,
    logger: logging.Logger = logger,
) -> nx.graph:
    """
    This method updates the dimensions of the pipes in the urban drainage network represented by the graph, using historical rainfall data.

    After the graph topology is optimized, this method updates the physical parameters of the network, such as the diameter of the pipes. Rainfall statistics are obtained from historical rainfall datasets, and assumptions are made about the capacity of the sewer network (for example, based on the rainfall return period). Alternatively, the problem can be simplified by directly assuming a capacity (in mm/hr) for the network, bypassing the need for historical rainfall data.

    The dimensions of the pipe (such as diameter) are then approximated based on the edge weights used for network optimization, subject to reasonable ranges for diameter, velocity, and slope.

    Parameters
    ----------
    graph: nx.graph
        The network graph to be updated with new pipe dimensions. This is an instance of the nx.graph class, which wraps around the NetworkX Graph class.

    rainfall_fn: Union[str, Path]
        The file path to the historical rainfall data used to derive rainfall statistics.

    rainfall_assumption: float
        An assumption about the capacity of the sewer system based on the rainfall return period (in years). This is used to derive statistics from the historical rainfall data.

    capacity_assumption: float, optional
        An assumption about the capacity of the sewer network (in mm/hr). If provided, this value is used instead of deriving capacity from historical rainfall data.

    diameter_range: tuple, optional
        A tuple specifying the minimum and maximum feasible pipe diameters. Used to constrain the approximation of pipe dimensions.

    velocity_range: tuple, optional
        A tuple specifying the minimum and maximum feasible pipe velocities. Used to constrain the approximation of pipe dimensions.

    slope_range: tuple, optional
        A tuple specifying the minimum and maximum feasible pipe slopes. Used to constrain the approximation of pipe dimensions.

    Returns
    -------
    nx.graph:
        The updated graph, an instance of nx.graph class, representing the urban drainage network with new pipe dimensions.

    References
    ----------
    - hydromt.stats.extremes: Function used to derive statistics from historical rainfall data.
    """
    # method implementation goes here
    pass


# graph workflows
def setup_graph_from_rasterdataset(
    graph: nx.Graph,  # TODO replace by self.graphs
    data_catalog: DataCatalog,  # TODO replace by self.data_catalog
    raster_fn: Union[str, Path, xr.DataArray, xr.Dataset],
    variables: Optional[List[str]] = None,
    fill_method: Optional[str] = None,
    resampling_method: Optional[str] = "mean",
    all_touched: Optional[bool] = True,
    rename: Optional[Dict[str, str]] = dict(),
    graph_component: Optional[str] = "both",
    logger: logging.Logger = logger,
) -> nx.graph:
    """Add data variable(s) from ``raster_fn`` to attribute(s) in graph object.

    Raster data is sampled to the graph edges and nodes using the ``resampling_method``.
    If raster is a dataset, all variables will be added unless ``variables`` list
    is specified.

    Parameters
    ----------
    raster_fn: str, Path, xr.DataArray, xr.Dataset
        Data catalog key, path to raster file or raster xarray data object.
    variables: list, optional
        List of variables to add to mesh from raster_fn. By default all.
    fill_method : str, optional
        If specified, fills no data values using fill_nodata method.
        Available methods are {'linear', 'nearest', 'cubic', 'rio_idw'}.
    rename: dict, optional
        Dictionary to rename variable names in raster_fn before adding to mesh
        {'name_in_raster_fn': 'name_in_mesh'}. By default empty.
    graph_component: str, optional
        Specifies which component of the graph to process. Can be one of the following:
        * "edges" - Only processes and updates the edges of the graph.
        * "nodes" - Only processes and updates the nodes of the graph.
        * "both" - Processes and updates both nodes and edges of the graph.
        By default, it processes both nodes and edges ("both").
    resampling_method: str, optional
        Method to sample from raster data to mesh. By default mean. Options include
        {'count', 'min', 'max', 'sum', 'mean', 'std', 'median', 'q##'}.
        Only used when ``graph_component`` is "edges" or "both".
    all_touched : bool, optional
        If True, all pixels touched by geometries will used to define the sample.
        If False, only pixels whose center is within the geometry or that are
        selected by Bresenham's line algorithm will be used. By default True.
        Only used when ``graph_component`` is "edges" or "both".
    Returns
    -------
    list
        List of variables added to mesh.
    """  # noqa: E501

    assert graph_component in [
        "edges",
        "nodes",
        "both",
    ], "Invalid graph_component value."

    logger.info(f"Preparing graph data from raster source {raster_fn}")
    region = graph_utils.graph_region(graph)

    # Read raster data, select variables, and interpolate na if needed
    ds = data_catalog.get_rasterdataset(
        raster_fn, region=region, buffer=1000, variables=variables
    )
    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset()
    if fill_method is not None:
        ds = ds.raster.interpolate_na(method=fill_method)

    # Sample raster data
    # Rename variables
    rm_dict = {f"{var}_{resampling_method}": var for var in ds.data_vars}
    # get sample at edges
    if graph_component in ["edges", "both"]:
        edges = graph_utils.graph_edges(graph)
        ds_sample = ds.raster.zonal_stats(
            gdf=edges,
            stats=resampling_method,
            all_touched=all_touched,
        )
        ds_sample = ds_sample.rename(rm_dict).rename(rename)
        ds_sample_df = ds_sample.to_dataframe().set_index(edges["id"])
        graph = update_edges_attributes(graph, ds_sample_df, id_col="id")
    # get sample at nodes
    if graph_component in ["nodes", "both"]:
        nodes = graph_utils.graph_nodes(graph)
        ds_sample = ds.raster.sample(
            nodes
        )  # FIXME dem for sample data too small, need to add nodedata value and correct crs
        ds_sample = ds_sample.rename(rename)
        ds_sample_df = ds_sample.to_dataframe().set_index(nodes["id"])
        graph = update_nodes_attributes(graph, ds_sample_df)

    # TODO Convert to UgridDataset
    # uds_sample = xu.UgridDataset(ds_sample, grids=self.mesh_grids[grid_name])

    return graph


def setup_graph_from_geodataframe(
    graph: nx.Graph,  # TODO replace by self.graphs
    data_catalog: DataCatalog,  # TODO replace by self.data_catalog
    vector_fn: Union[str, Path, xr.DataArray, xr.Dataset],
    variables: Optional[List[str]] = None,
    max_dist: float = np.inf,
    rename: Optional[Dict[str, str]] = dict(),
    graph_component: Optional[str] = "both",
    logger: logging.Logger = logger,
) -> nx.graph:
    """Add data variable(s) from ``vector_fn`` to attribute(s) in graph object using nearest_join.

    Raster data is sampled to the graph edges and nodes using the ``resampling_method``.
    If raster is a dataset, all variables will be added unless ``variables`` list
    is specified.

    Parameters
    ----------
    raster_fn: str, Path, xr.DataArray, xr.Dataset
        Data catalog key, path to raster file or raster xarray data object.
    variables: list, optional
        List of variables to add to mesh from raster_fn. By default all.
    fill_method : str, optional
        If specified, fills no data values using fill_nodata method.
        Available methods are {'linear', 'nearest', 'cubic', 'rio_idw'}.
    rename: dict, optional
        Dictionary to rename variable names in raster_fn before adding to mesh
        {'name_in_raster_fn': 'name_in_mesh'}. By default empty.
    graph_component: str, optional
        Specifies which component of the graph to process. Can be one of the following:
        * "edges" - Only processes and updates the edges of the graph.
        * "nodes" - Only processes and updates the nodes of the graph.
        * "both" - Processes and updates both nodes and edges of the graph.
        By default, it processes both nodes and edges ("both").
    resampling_method: str, optional
        Method to sample from raster data to mesh. By default mean. Options include
        {'count', 'min', 'max', 'sum', 'mean', 'std', 'median', 'q##'}.
        Only used when ``graph_component`` is "edges" or "both".
    all_touched : bool, optional
        If True, all pixels touched by geometries will used to define the sample.
        If False, only pixels whose center is within the geometry or that are
        selected by Bresenham's line algorithm will be used. By default True.
        Only used when ``graph_component`` is "edges" or "both".
    Returns
    -------
    list
        List of variables added to mesh.
    """  # noqa: E501

    assert graph_component in [
        "edges",
        "nodes",
        "both",
    ], "Invalid graph_component value."

    logger.info(f"Preparing graph data from raster source {vector_fn}")
    region = graph_utils.graph_region(graph)

    # Read vector data, select variables
    gdf = data_catalog.get_geodataframe(
        vector_fn, region=region, buffer=1000, variables=variables
    )

    # relate the geodataframe to graph using nearet join
    # get gdf attributes at edges
    if graph_component in ["edges", "both"]:
        edges = graph_utils.graph_edges(graph)
        gdf_sample = find_edge_ids_by_snapping(graph, gdf, snap_offset=max_dist)
        gdf_sample = gdf_sample.rename(rename)
        gdf_sample = gdf_sample.set_index(edges["id"])
        graph = update_edges_attributes(graph, gdf_sample, id_col="id")
    # get sample at nodes
    if graph_component in ["nodes", "both"]:
        nodes = graph_utils.graph_nodes(graph)
        gdf_sample = find_node_ids_by_snapping(graph, gdf, snap_offset=max_dist)
        gdf_sample = gdf_sample.rename(rename)
        gdf_sample = gdf_sample.set_index(nodes["id"])
        graph = update_nodes_attributes(graph, gdf_sample)

    # TODO Convert to UgridDataset
    # uds_sample = xu.UgridDataset(ds_sample, grids=self.mesh_grids[grid_name])

    return graph


# func from hybridurb


def find_edge_ids_by_snapping(
    G: nx.Graph,
    edges: gpd.GeoDataFrame,
    snap_offset: float = 1,
    snap_method: str = "overall",
) -> gpd.GeoDataFrame:
    """This function adds "id" to edges GeoDataFrame"""

    # graph
    _ = gpd.GeoDataFrame(nx.to_pandas_edgelist(G).set_index("id"))

    # wrapper to use delft3dfmpy function to find "branch_id"
    _ = _.rename({"id": "branch_id"}).assign(branchType=None)
    find_nearest_branch(
        _,
        edges,
        method=snap_method,
        maxdist=snap_offset,
        move_geometries=True,
    )

    # rename "branch_id" to "edge_id"
    edges_with_ids = edges.rename(columns={"branch_id": "_id"})

    return edges_with_ids


def find_node_ids_by_snapping(
    G: nx.Graph,
    nodes: gpd.GeoDataFrame,
    snap_offset: float = 1,
    snap_method: str = "overall",
) -> gpd.GeoDataFrame:
    """This function adds "id" to nodes GeoDataFrame"""

    # graph
    G_nodes = gpd.GeoDataFrame(
        {
            "geometry": [Point(p) for p in G.nodes],
            "node_id": [f"{p[0]:.6f}_{p[1]:.6f}" for p in G.nodes],
            "_id": G.nodes(data="id"),
        }
    ).set_index("node_id")

    # nodes
    nodes.loc[:, "node_id"] = [
        f"{x:.6f}_{y:.6f}" for x, y in zip(nodes.geometry.x, nodes.geometry.y)
    ]
    nodes = nodes.set_index("node_id")

    # map user nodes to derived nodes
    if set(nodes.index).issubset(set(G_nodes.index)):
        # check if 1-to-1 match
        G_nodes_new = G_nodes.join(nodes.drop(columns="geometry"))
    else:
        # use snap_nodes_to_nodes function to find "node_id"
        G_nodes_new = snap_nodes_to_nodes(nodes, G_nodes, snap_offset)
        logger.debug("performing snap nodes to graph nodes")

    # assign id from graph to nodes
    nodes = nodes.join(G_nodes_new["_id"])

    return nodes


def update_edges_attributes(
    graph: nx.Graph,
    edges: gpd.GeoDataFrame,
    id_col: str = "id",
) -> nx.Graph():
    """This function updates the graph by adding new edges attributes specified in edges.

    Only edges with matching ids specified in "id" will be updated."""

    # graph df
    _graph_df = nx.to_pandas_edgelist(graph).set_index("id")
    _graph_df["_graph_edge_tuple"] = list(graph.edges)

    # check if edges id in attribute df
    if edges.index.name == id_col:
        edges.index.name = "id"
    elif id_col in edges.columns:
        edges = edges.set_index(id_col)
        edges.index.name = "id"
    else:
        raise ValueError(
            "attributes could not be updated to graph: could not perform join"
        )

    # last item that isnt NA
    graph_df = _graph_df.reindex(
        columns=_graph_df.columns.union(edges.columns, sort=False)
    )
    graph_df.update(edges)

    # add each attribute
    for c in edges.columns:
        dict = {row._graph_edge_tuple: row[c] for i, row in graph_df.iterrows()}
        nx.set_node_attributes(graph, dict, c)

    return graph


# func from hybridurb
def update_nodes_attributes(
    graph: nx.Graph,
    nodes: gpd.GeoDataFrame,
    id_col: str = "id",
) -> nx.Graph():
    """This function updates the graph by adding new edges attributes specified in edges.

    Only edges with matching ids specified in "id" will be updated."""

    # graph df
    _graph_df = pd.DataFrame(graph.nodes(data="id"), columns=["tuple", "id"]).set_index(
        "id"
    )

    # check if edges id in attribute df
    if nodes.index.name == id_col:
        nodes.index.name = "id"
    elif id_col in nodes.columns:
        nodes = nodes.set_index(id_col)
        nodes.index.name = "id"
    else:
        raise ValueError(
            "attributes could not be updated to graph: could not perform join"
        )

    # last item that isnt NA
    _graph_df = _graph_df.reindex(
        columns=_graph_df.columns.union(nodes.columns, sort=False)
    )
    graph_df = pd.concat([_graph_df, nodes]).groupby(level=0).last()
    graph_df = graph_df.loc[_graph_df.index]

    # add each attribute
    for c in nodes.columns:
        dict = {row.tuple: row[c] for i, row in graph_df.iterrows()}
        nx.set_node_attributes(graph, dict, c)

    return graph


def reverse_edges_on_negative_weight(
    graph: Union[nx.MultiDiGraph, nx.DiGraph], weight: str = "gradient"
) -> Union[nx.MultiDiGraph, nx.DiGraph]:
    """Reverse graph edges based on a negative weight attribute.

    Parameters:
    -----------
    graph : Union[nx.DiGraph,nx.MultiDiGraph]
        Input graph.
    weight : str
        Name of the edge attribute to check. Default is 'gradient'.

    Returns:
    --------
    Union[nx.DiGraph,nx.MultiDiGraph]
        The function modifies the graph in-place and returns it.
    """

    if isinstance(graph, nx.MultiDiGraph):
        # Collect edges that need to be reversed for MultiDiGraph
        edges_to_reverse = [
            (u, v, key, data)
            for u, v, key, data in graph.edges(keys=True, data=True)
            if data.get(weight, 0) < 0
        ]
        for u, v, key, data in edges_to_reverse:
            # Remove original edge
            graph.remove_edge(u, v, key=key)
            # Add reversed edge with preserved attributes
            graph.add_edge(v, u, key=key, **data)

    elif isinstance(graph, nx.DiGraph):
        # Collect edges that need to be reversed for DiGraph
        edges_to_reverse = [
            (u, v, data)
            for u, v, data in graph.edges(data=True)
            if data.get(weight, 0) < 0
        ]
        for u, v, data in edges_to_reverse:
            # Remove original edge
            graph.remove_edge(u, v)
            # Add reversed edge with preserved attributes
            graph.add_edge(v, u, **data)

    else:
        raise ValueError("The graph should be either a DiGraph or a MultiDiGraph.")

    return graph
