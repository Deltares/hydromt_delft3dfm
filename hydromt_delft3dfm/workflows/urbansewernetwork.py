# -*- coding: utf-8 -*-

import logging
from typing import List, Union
from pathlib import Path

import numpy as np
import geopandas as gpd
import networkx as nx
import pyproj
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon, box
from shapely.wkt import dumps, loads

logger = logging.getLogger(__name__)


__all__ = [
    "setup_network_from_openstreetmap",
    "setup_network_connections_based_on_flowdirections",
    "setup_network_parameters_from_rasters",
    "setup_network_topology_optimization",
    "setup_network_dimentions_from_rainfallstats",
]


def setup_network_from_openstreetmap(
    polygon: gpd.GeoDataFrame, network_type: str, road_types: list[str], crs: pyproj.CRS
) -> nx.graph:
    """
    This method transforms Open Street Map (OSM) data into a graph network representing an urban drainage system.

    The graph is an object that wraps around the NetworkX Graph (nx.Graph) class, providing additional functionality such as conversion methods and reading/writing capabilities.

    The steps involved are:
    1. Fetch vector lines from OSM data based on a provided query related to the specified network type and road types.
    2. Clean up the retrieved lines to remove any irrelevant or erroneous data.
    3. Compose the cleaned lines into a graph.
    4. Extract a network from the graph, where the network corresponds to the graph edges.

    Parameters:
    -----------
    polygon: GeoDataFrame
        The region polygon to consider for fetching OSM data. This defines the geographical area of interest.

    network_type: str
        The type of street network to consider. This helps filter the OSM data to include only relevant road types.

    road_types: list[str]
        A list of road types to consider during the creation of the graph. This further refines the data that is included from the OSM dataset.

    crs: pyproj.CRS
        The Coordinate Reference System to use for the geospatial data.

    Returns:
    --------
    nx.graph:
        An instance of the nx.graph class, which is a wrapper around nx.Graph. The graph represents the urban drainage network with the following properties:
        - Global property: 'crs'
        - Edge property: 'edgeid', 'geometry'
        - Node property: 'nodeid', 'geometry'

    References:
    -----------
    - ra2ce.OsmNetworkWrapper.get_network: https://github.com/Deltares/ra2ce/blob/master/ra2ce/graph/network_wrappers/osm_network_wrapper/osm_network_wrapper.py
    - osmnx.graph.graph_from_polygon: https://osmnx.readthedocs.io/en/latest/internals-reference.html
    """
    # method implementation goes here
    pass


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
