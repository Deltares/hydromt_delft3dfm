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

# TODO #65: ra2ce installation issue
from ra2ce.graph.network_config_data.network_config_data import NetworkConfigData
from ra2ce.graph.network_wrappers.osm_network_wrapper.osm_network_wrapper import (
    OsmNetworkWrapper,
)

from hydromt_delft3dfm import graph_utils

logger = logging.getLogger(__name__)


__all__ = [
    # graph creations
    "setup_graph_from_openstreetmap",
    "setup_graph_from_hydrography",
    # graph update
    "update_graph_from_dem",
    # TODO to be categorized
    "setup_network_connections_based_on_flowdirections",
    "setup_network_parameters_from_rasters",
    "setup_network_topology_optimization",
    "setup_network_dimentions_from_rainfallstats",
    # graph manipulations
    "setup_graph_from_rasterdataset",
]


# FIXME maybe rename to create*
# TODO remove the non-main component outside of region polygon.
def setup_graph_from_openstreetmap(
    region: gpd.GeoDataFrame,
    network_type: str = None,
    road_types: list[str] = None,
    logger: logging.Logger = logger,
) -> nx.MultiDiGraph:
    """
    Download and process Open Street Map (OSM) data within region into a graph network.

    The steps involved are:
    1. Fetch vector lines from OSM data based on a provided query and perform cleaning.
    2. simplify the graph's topology by removing interstitial nodes
    3. extract network edges and nodes from the graph (extract geometries from graph)
    4. recreate graph using network edges and nodes (add geometries to graph)

    Parameters
    ----------
    region: GeoDataFrame
        The region polygon to consider for fetching OSM data.
        Must contain crs.

    network_type: str {"all_private", "all", "bike", "drive", "drive_service", "walk"})
        The type of street network to consider. This helps filter the OSM data.
        If None, use road_types to filter features.
        By default None.

    road_types: list[str], optional
        A list of road types to consider during the creation of the graph.
        This is a refined filter to get data that is included from the OSM dataset.
        A complete list can be found in: https://wiki.openstreetmap.org/wiki/Key:highway.
        by default None.

    Returns
    -------
    nx.MultiDiGraph:
        An instance of the graph as a multi-digraph, with geometry information for edges
        and nodes.
        - Global property: 'crs'.
        - Edge property: ['edgeid', 'geometry', 'node_start', 'node_end','osmid',
                          'highway', 'oneway', 'reversed', 'length', 'rfid_c', , 'rfid']
        - Node property: ['nodeid', 'geometry', 'y', 'x', 'street_count']

    See Also
    --------
    - ra2ce.OsmNetworkWrapper.get_clean_graph_from_osm
    - osmnx.simplfy_graph
    - ra2ce.OsmNetworkWrapper.get_clean_graph
    - graph_utils.preprocess_graph
    """
    # method implementation goes here

    # this function use OsmNetworkWrapper from race to download OSM data into network

    # get crs
    crs = region.crs

    def _get_osm_wrapper() -> OsmNetworkWrapper:
        # creates and configures a osm network wrapper
        _network_config_data = NetworkConfigData()
        _osm_network_wrapper = OsmNetworkWrapper(_network_config_data)
        # configure wrapper properties
        _osm_network_wrapper.is_directed = True
        _osm_network_wrapper.network_type = network_type
        _osm_network_wrapper.road_types = road_types
        _osm_network_wrapper.polygon_path = _get_temp_polygon_path()
        _osm_network_wrapper.output_graph_dir = _get_temp_output_graph_dir()
        return _osm_network_wrapper

    # TODO #65: replace funcs to get temp paths after ra2ce update
    def _get_temp_polygon_path() -> Path:
        # create a temporary file for region polygon in EPSG:4326
        with tempfile.TemporaryFile(suffix=".geojson") as temp:
            # dump the data to the file
            region.to_crs(pyproj.CRS.from_user_input(4326)).geometry.reset_index(
                drop=True
            ).to_file(temp.name, driver="GeoJSON")
        return Path(temp.name)

    def _get_temp_output_graph_dir() -> Path:
        # create a temporary directory
        temp = tempfile.TemporaryDirectory()
        return Path(temp.name)

    # download from osm and perform cleaning (drop_duplicates, add_missing_geoms_graph, snap_nodes_to_nodes)
    _osm_network_wrapper = _get_osm_wrapper()
    _osm_graph = _osm_network_wrapper.get_clean_graph_from_osm()
    logger.info("Get clean graph from osm.")

    # simplify the graph's topology by removing interstitial nodes
    _osm_simplified_graph = osmnx.simplify_graph(
        _osm_graph
    )  # TODO #64: needs testing, too longe/too short are not beneficial for dem based direction.
    _osm_simplified_graph = _osm_network_wrapper.get_clean_graph(
        _osm_simplified_graph
    )  # clean again for the simplified graph
    logger.info("Get simplified graph from osm.")

    # preprocess to desired graph format
    graph = graph_utils.preprocess_graph(_osm_simplified_graph, to_crs=crs)

    # get edges and nodes from graph (momepy convention)
    # edges, nodes = graph_utils.graph_to_network(_osm_simplified_graph, crs=crs)

    # get new graph with correct crs
    # graph = graph_utils.network_to_graph(
    #     edges=edges, nodes=nodes, create_using=nx.MultiDiGraph
    # )

    # unit test
    # _edges, _nodes = graph_utils.graph_to_network(graph)
    # _graph = graph_utils.network_to_graph(
    #     edges=_edges, nodes=_nodes, create_using=nx.MultiDiGraph
    # )
    # nx.is_isomorphic(_graph, graph) -->  True

    return graph


# FIXME maybe rename to create*
def setup_graph_from_hydrography(
    region: gpd.GeoDataFrame,
    ds_hydro: xr.Dataset,
    min_sto: int = 1,
    logger: logging.Logger = logger,
) -> nx.Graph:
    """Preprocess DEM to graph representing flow directions.

    The steps involved are:
    1. Obtain or derive D8 flow directions from a dem.
    2. Get flow direction vectors as geodataframe based on a minumum stream order.
    3. convert the geodataframe into graph.

    Parameters
    ----------
    region: gpd.GeoDataframe
        Region geodataframe that provide model extent and crs.
    ds_hydro : xr.Dataset
        Hydrography data (or dem) to derive river shape and characteristics from.
        * Required variables: ['elevtn']
        * Optional variables: ['flwdir', 'uparea']
    min_sto : int, optional
        minimum stream order of subbasins, by default the stream order is set to
        one under the global maximum stream order (slow).

    Returns
    -------
    nx.DiGraph:
        The processed flow direction as a graph.

    """
    crs = region.crs

    # extra fix of data if user provided tiff
    elevtn = ds_hydro["elevtn"]
    if elevtn.raster.res[1] > 0:
        elevtn = elevtn.raster.flipud()
    if np.isnan(elevtn.raster.nodata):
        elevtn.raster.set_nodata(-9999.0)
        logger.debug("Missing nodata value, setting to -9999.0")

    # check if flwdir and uparea in ds_hydro
    if "flwdir" not in ds_hydro.data_vars:
        da_flw = flw.d8_from_dem(
            elevtn,
            max_depth=-1,  # no local pits
            outlets="edge",
            idxs_pit=None,
        )
        logger.info("flwdir computed from elevtn")
    else:
        da_flw = ds_hydro["flwdir"]

    flwdir = flw.flwdir_from_da(da_flw, ftype="d8")
    if min_sto > 1:
        # will add stream order key "strord" to feat
        feat = flwdir.streams(min_sto=min_sto)
    else:
        # add stream order key "strord" manually
        feat = flwdir.streams(strord=flwdir.stream_order())
    logger.info("Obtain stream vectors")

    # get stream vector geodataframe
    gdf = gpd.GeoDataFrame.from_features(feat, crs=ds_hydro.raster.crs)
    # convert to local crs, because flwdir is performed on the raster crs
    gdf = gdf.to_crs(crs)

    # convert to graph
    graph = graph_utils.gpd_to_digraph(gdf)
    graph = graph_utils.preprocess_graph(graph, to_crs=crs)

    # TODO #63: test a few different scales and compare with real cases

    # plot
    # import matplotlib.pyplot as plt
    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot()
    # ds_hydro["elevtn"].plot(
    #     ax=ax, zorder=1, cbar_kwargs=dict(aspect=30, shrink=0.5), alpha=0.5, **kwargs
    # )
    # gdf.to_crs(ds_hydro.raster.crs).plot(
    #     ax=ax,
    #     color="blue",
    #     linewidth=gdf["strord"] / 3,
    #     label="Original flow directions",
    # )
    # f = Path().joinpath("temp.geojson")
    # gdf.to_file(f)

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


# func from hybridurb
def update_edges_attributes(
    graph: nx.Graph,
    edges: gpd.GeoDataFrame,
    id_col: str = "id",
) -> nx.Graph():
    """This function updates the graph by adding new edges attributes specified in edges"""

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
    """This function updates the graph by adding new edges attributes specified in edges"""

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
