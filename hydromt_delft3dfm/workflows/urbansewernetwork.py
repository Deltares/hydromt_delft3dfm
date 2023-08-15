# -*- coding: utf-8 -*-

import logging
from typing import List, Union, Tuple
from pathlib import Path
import tempfile

import numpy as np
import geopandas as gpd
import osmnx
import networkx as nx
import pyproj

from shapely.wkt import dumps, loads
from hydromt_delft3dfm import workflows

# TODO: ra2ce installation issue
from ra2ce.graph.network_config_data.network_config_data import NetworkConfigData
from ra2ce.graph.network_wrappers.osm_network_wrapper.osm_network_wrapper import (
    OsmNetworkWrapper,
)

logger = logging.getLogger(__name__)


__all__ = [
    "setup_graph_from_openstreetmap",
    "setup_network_connections_based_on_flowdirections",
    "setup_network_parameters_from_rasters",
    "setup_network_topology_optimization",
    "setup_network_dimentions_from_rainfallstats",
]


def setup_graph_from_openstreetmap(
    region: gpd.GeoDataFrame,
    network_type: str = "drive",
    road_types: list[str] = None,
) -> nx.MultiDiGraph:
    """
    This method transforms Open Street Map (OSM) data into a graph network representing the road network within region.

    The steps involved are:
    1. Fetch vector lines from OSM data based on a provided query and perform cleaning to remove any irrelevant or erroneous data.
    2. simplify the graph's topology by removing interstitial nodes
    3. extract network edges and nodes from the graph (extract geometries from graph)
    4. recreate graph using network edges and nodes (add geometries to graph)

    Parameters:
    -----------
    region: GeoDataFrame
        The region polygon to consider for fetching OSM data. This defines the geographical area of interest.
        Must contain crs.

    network_type: str {"all_private", "all", "bike", "drive", "drive_service", "walk"})
        The type of street network to consider. This helps filter the OSM data to include only relevant road types.
        By default "drive"

    road_types: list[str], optional
        A list of road types to consider during the creation of the graph. This further refines the data that is included from the OSM dataset.
        A complete list can be found in: https://wiki.openstreetmap.org/wiki/Key:highway.
        by default None.

    Returns:
    --------
    nx.MultiDiGraph:
        An instance of the graph as a multi-digraph, with geometry information for edges and nodes.
        - Global property: 'crs'.
        - Edge property: ['edgeid', 'geometry', 'node_start', 'node_end','osmid', 'highway', 'oneway', 'reversed', 'length', 'rfid_c', , 'rfid']
        - Node property: ['nodeid', 'geometry', 'y', 'x', 'street_count']

    See Also:
    -----------
    - ra2ce.OsmNetworkWrapper.get_clean_graph_from_osm
    - osmnx.simplfy_graph
    - workflows.graph_to_network
    - workflows.network_to_graph
    """
    # method implementation goes here

    # this function use OsmNetworkWrapper from race to download OSM data into network

    # get crs
    crs = region.crs

    # funcs to get temp paths
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

    # get graph from osm data
    # create an emtpy wrapper
    _network_config_data = NetworkConfigData()
    _osm_network_wrapper = OsmNetworkWrapper(_network_config_data)
    # configure wrapper properties
    _osm_network_wrapper.is_directed = True
    _osm_network_wrapper.network_type = network_type
    _osm_network_wrapper.road_types = road_types
    _osm_network_wrapper.polygon_path = _get_temp_polygon_path()
    _osm_network_wrapper.output_graph_dir = _get_temp_output_graph_dir()
    # download from osm and perform cleaning (drop_duplicates, add_missing_geoms_graph, snap_nodes_to_nodes)
    _osm_graph = _osm_network_wrapper.get_clean_graph_from_osm()

    # simplify the graph's topology by removing interstitial nodes
    _osm_simplified_graph = osmnx.simplify_graph(_osm_graph)

    # get edges and nodes from graph (momepy convention)
    edges, nodes = workflows.graph_to_network(_osm_simplified_graph, crs=crs)

    # get new graph with correct crs
    graph = workflows.network_to_graph(
        edges=edges, nodes=nodes, create_using=nx.MultiDiGraph
    )

    # unit test
    # _edges, _nodes = workflows.graph_to_network(graph)
    # _graph = workflows.network_to_graph(
    #     edges=_edges, nodes=_nodes, create_using=nx.MultiDiGraph
    # )
    # nx.is_isomorphic(_graph, graph) -->  True

    return graph


def setup_network_connections_based_on_flowdirections(
    user_input_graph: nx.graph, dem_fn: Union[str, Path], **kwargs
) -> nx.graph:
    """
    This method sets up connections in the urban drainage network graph based on flow directions derived from a Digital Elevation Model (DEM).

    The method performs the following steps:
    1. Retrieves a flow direction raster from a DEM.
    2. Converts the flow direction raster into a graph and a network, using the network to represent flow directions.
    3. Snaps nodes from the `user_input_graph` to the flow direction network to ensure that the network aligns with natural flow paths.
    4. Identifies and fills in any missing links in the `user_input_graph`, ensuring a fully connected network.
    5. Fills in any missing geometry in the `user_input_graph` and converts it into a network.

    Parameters:
    -----------
    user_input_graph: nx.graph
        The initial graph provided by the user, which will be updated to align with the flow direction network.

    dem_fn: Union[str, Path]
        The file path to the digital elevation model (DEM) used to derive the flow direction raster. This can be provided as a string or a Path object.

    **kwargs:
        Other keyword arguments that may be required for specific implementations of the function.

    Returns:
    --------
    nx.graph:
        The processed `user_input_graph`, now an instance of nx.graph class, with any missing links and geometry filled in based on flow directions. The graph represents the urban drainage network.

    Notes:
    ------
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

    Parameters:
    -----------
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

    Returns:
    --------
    nx.graph:
        The updated graph, an instance of nx.graph class, representing the urban drainage network with physical parameters filled in.
    """
    # method implementation goes here
    pass


def setup_network_topology_optimization(
    graph: nx.graph, method: str, **kwargs
) -> nx.graph:
    """
    This method optimizes the topology of the urban drainage network represented by the graph.

    The method performs the following steps:
    1. Computes the selected graph metric (defined by 'method') for all nodes in the graph.
    2. Removes unnecessary links based on the computed graph metric.
    3. Iterates the above steps as needed until the network is optimized.

    Parameters:
    -----------
    graph: nx.graph
        The network graph to be optimized. This is an instance of the nx.graph class, which wraps around the NetworkX Graph class.

    method: str
        The method to use for optimization. This determines the graph metric that will be computed and used for the optimization process.

    **kwargs:
        Additional keyword arguments for more specific implementations of the function.

    Returns:
    --------
    nx.graph:
        The optimized graph, an instance of nx.graph class, representing the urban drainage network.

    Notes:
    ------
    The betweenness centrality of a node in a graph is a measure of how often it appears on the shortest paths between nodes. Nodes with high betweenness centrality can be considered 'hubs' in the network. By removing links with low betweenness centrality, we can potentially simplify the network topology without significantly impacting the overall connectivity.

    The optimal way to optimize network topology will highly depend on the specific characteristics of the network, which are often determined by previous steps in the network analysis process.

    Other candidate graph metrics that might be useful for network topology optimization include:
    - Degree centrality: Measures the number of edges connected to a node. This can be useful for identifying nodes that are most connected.
    - Closeness centrality: Measures how close a node is to all other nodes in the network. This can be useful for identifying nodes that are centrally located.
    - Eigenvector centrality: Measures a node's influence based on the number of links it has to other influential nodes. This can be useful for identifying influential nodes in the network.

    If the network's weights have physical meanings, shortest path or maximum flow methods can also be considered for optimization. These methods, such as the Ford-Fulkerson or Edmonds-Karp algorithm, can be used to identify critical paths in the network. This can be particularly useful for identifying vulnerabilities in the network or for identifying potential areas for network enhancement.

    See Also:
    ----------
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
) -> nx.graph:
    """
    This method updates the dimensions of the pipes in the urban drainage network represented by the graph, using historical rainfall data.

    After the graph topology is optimized, this method updates the physical parameters of the network, such as the diameter of the pipes. Rainfall statistics are obtained from historical rainfall datasets, and assumptions are made about the capacity of the sewer network (for example, based on the rainfall return period). Alternatively, the problem can be simplified by directly assuming a capacity (in mm/hr) for the network, bypassing the need for historical rainfall data.

    The dimensions of the pipe (such as diameter) are then approximated based on the edge weights used for network optimization, subject to reasonable ranges for diameter, velocity, and slope.

    Parameters:
    -----------
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

    Returns:
    --------
    nx.graph:
        The updated graph, an instance of nx.graph class, representing the urban drainage network with new pipe dimensions.

    References:
    -----------
    - hydromt.stats.extremes: Function used to derive statistics from historical rainfall data.
    """
    # method implementation goes here
    pass
