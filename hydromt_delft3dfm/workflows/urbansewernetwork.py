# -*- coding: utf-8 -*-

import logging
from typing import List, Union, Tuple
from pathlib import Path
import tempfile

import xarray as xr
import numpy as np
import geopandas as gpd
import osmnx
import networkx as nx
import pyproj
from shapely.geometry import Point, LineString

from shapely.wkt import dumps, loads
from hydromt_delft3dfm import workflows

# TODO: ra2ce installation issue
from ra2ce.graph.network_config_data.network_config_data import NetworkConfigData
from ra2ce.graph.network_wrappers.osm_network_wrapper.osm_network_wrapper import (
    OsmNetworkWrapper,
)

# hydromt
from hydromt import DataCatalog, flw


logger = logging.getLogger(__name__)


__all__ = [
    "setup_graph_from_openstreetmap",
    "setup_graph_from_hydrography",
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

    # TODO: replace funcs to get temp paths after ra2ce update
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

    # TODO: func from ra2ce utils due to installation issue
    def _add_missing_geoms_osmgraph(
        graph: nx.Graph, geom_name: str = "geometry"
    ) -> nx.Graph:
        # Not all nodes have geometry attributed (some only x and y coordinates) so add a geometry columns
        nodes_without_geom = [
            n[0] for n in graph.nodes(data=True) if geom_name not in n[-1]
        ]
        for nd in nodes_without_geom:
            graph.nodes[nd][geom_name] = Point(
                graph.nodes[nd]["x"], graph.nodes[nd]["y"]
            )

        edges_without_geom = [
            e for e in graph.edges.data(data=True) if geom_name not in e[-1]
        ]
        for ed in edges_without_geom:
            graph[ed[0]][ed[1]][0][geom_name] = LineString(
                [graph.nodes[ed[0]][geom_name], graph.nodes[ed[1]][geom_name]]
            )

        return graph

    # download from osm and perform cleaning (drop_duplicates, add_missing_geoms_graph, snap_nodes_to_nodes)
    _osm_network_wrapper = _get_osm_wrapper()
    _osm_graph = _osm_network_wrapper.get_clean_graph_from_osm()

    # simplify the graph's topology by removing interstitial nodes
    _osm_simplified_graph = osmnx.simplify_graph(_osm_graph)

    # preprocess to desired graph format
    # add missing geoms using func customised for osm data
    _graph = _add_missing_geoms_osmgraph(_osm_simplified_graph)
    graph = workflows.preprocess_graph(_graph, to_crs=crs)

    # get edges and nodes from graph (momepy convention)
    # edges, nodes = workflows.graph_to_network(_osm_simplified_graph, crs=crs)

    # get new graph with correct crs
    # graph = workflows.network_to_graph(
    #     edges=edges, nodes=nodes, create_using=nx.MultiDiGraph
    # )

    # unit test
    # _edges, _nodes = workflows.graph_to_network(graph)
    # _graph = workflows.network_to_graph(
    #     edges=_edges, nodes=_nodes, create_using=nx.MultiDiGraph
    # )
    # nx.is_isomorphic(_graph, graph) -->  True

    return graph


def setup_graph_from_hydrography(
    region: gpd.GeoDataFrame,
    ds_hydro: xr.Dataset,
    min_sto: int = 1,
    **kwargs,
) -> nx.Graph:
    """Preprocess DEM to graph representing flow directions.

    The steps involved are:
    1. Obtain or derive D8 flow directions from a dem.
    2. Get flow direction vectors as geodataframe based on a minumum stream order (`min_sto`)
    3. convert the geodataframe into graph
    4. recreate graph using network edges and nodes (add geometries to graph)

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
    **kwargs:
        Other keyword arguments that may be required for specific implementations of the function.

    Returns:
    --------
    nx.graph:
        The processed `user_input_graph`, now an instance of nx.graph class, with any missing links and geometry filled in based on flow directions. The graph represents the urban drainage network.

    """

    crs = region.crs

    # FIXME: discuss with Helene
    # question regarding hydromt.flw.d8_from_dem:
    # I am trying to use the function on a raster dataset resembling a "local dem"
    # - raster data containing elevtn variable that is read from a tiff in datacatalog.
    # The resulting raster has positive resolutions on both x,y direction: (0.0008, 0.0008).
    # However, the `hydromt.flw.d8_from_dem`  requires the y resolution to be negative.
    # hence the flipud()
    # also the function does not allow no data as nan
    # hence the set_nodata(-9999)
    # But I also wonder if I should even support local data at all,
    # considering the function is to derive that from global data

    # extra fix of data if user provided tiff
    elevtn = ds_hydro["elevtn"]
    if elevtn.raster.res[1] > 0:
        elevtn = elevtn.raster.flipud()
    if np.isnan(elevtn.raster.nodata):
        elevtn.raster.set_nodata(-9999.0)

    # check if flwdir and uparea in ds_hydro
    if "flwdir" not in ds_hydro.data_vars:
        da_flw = flw.d8_from_dem(
            elevtn,
            max_depth=-1,  # no local pits
            outlets="edge",
            idxs_pit=None,
        )
    else:
        da_flw = ds_hydro["flwdir"]

    flwdir = flw.flwdir_from_da(da_flw, ftype="d8")
    if min_sto > 1:
        # will add stream order key "strord" to feat
        feat = flwdir.streams(min_sto=min_sto)
    else:
        # add stream order key "strord" manually
        feat = flwdir.streams(strord=flwdir.stream_order())

    # get stream vector geodataframe
    gdf = gpd.GeoDataFrame.from_features(feat, crs=ds_hydro.raster.crs)
    # convert to local crs, because flwdir is performed on the raster crs
    gdf = gdf.to_crs(crs)

    # convert to graph
    graph = workflows.gpd_to_digraph(gdf)
    graph = workflows.preprocess_graph(graph, to_crs=crs)

    # plot
    # FIXME discuss the results with helene
    # it seems like the flow is going upstream
    # but nice basin delineation --> check how it compares with the bangkok model
    # might indicate the challenges for getting flow directions for flat areas.

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

    # FIXME do I need the below?
    # if "uparea" not in ds_hydro.data_vars:
    #     da_upa = xr.DataArray(
    #         dims=ds_hydro["elevtn"].raster.dims,
    #         coords=ds_hydro["elevtn"].raster.coords,
    #         data=flwdir.upstream_area(unit="km2"),
    #         name="uparea",
    #     )
    #     da_upa.raster.set_nodata(-9999)
    #     ds_hydro["uparea"] = da_upa

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
