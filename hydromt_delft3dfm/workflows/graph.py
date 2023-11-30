"""Workflows to prepare graphs from different sources (hydromt networkmixins + other workflows)"""

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

from hydromt_delft3dfm import graph_utils
from hydromt_delft3dfm.workflows import (
    find_nearest_branch,
    explode_and_deduplicate_geometries,
)
from hydromt.gis_utils import nearest_merge

logger = logging.getLogger(__name__)


__all__ = [
    # create graph topology
    "create_graph_from_geodataframe",
    "create_graph_from_openstreetmap",
    "create_graph_from_hydrography",
    # setup graph attributes
    "setup_graph_from_rasterdataset",
    "setup_graph_from_geodataframe",
]


def create_graph_from_openstreetmap(
    region: gpd.GeoDataFrame,
    osm_key: str,
    osm_values: list[str],
    buffer: float = 0.0,
    logger: logging.Logger = logger,
) -> nx.MultiDiGraph:
    """
    Download and process Open Street Map (OSM) data within region into a graph network.

    Only map features selected by the `osm_key` and 'osm_values' are used.
    A complete describtion of `osm_key` and 'osm_values' can be found in:
    https://wiki.openstreetmap.org/wiki/.


    Parameters
    ----------
    region: GeoDataFrame
        The region polygon to consider for fetching OSM data.
        Must contain crs.

    osm_key: str
        OSM "key" to categorize or describe a particular characteristic of map features.
        The concept of "key" in OSM is part of the key-value pairing system,
        which is central to how OSM data is structured.
        Some common keys include:
            highway: Describes roads and paths.
            waterway: Describes rivers, streams, canals, etc.

    osm_values: list[str]
        OSM "values" to specify attributes assigned to a key to further detail or describe a map feature.
        Together, a key-value pair provides descriptive information about OSM entities.
        Some common values associated to the common keys include:
            highway: motorway,motorway_link,primary,primary_link,secondary,secondary_link,tertiary,tertiary_link,residential
            waterway:river,stream,brook,canal,ditch

    buffer: float, optional
        A buffer applied to the region polygon.
        By default 0.0.

    Returns
    -------
    nx.MultiDiGraph:
        An instance of the graph as a multi-digraph, with geometry information for edges
        and nodes.
        - Global property: 'crs'.
        - Edge property: ['edgeid', 'geometry', 'node_start', 'node_end','osmid']
        - Node property: ['nodeid', 'geometry']

    See Also
    --------
    - osmnx.graph_from_polygon
    - graph_utils.preprocess_graph
    """

    # get crs
    crs = region.crs

    # download osm data as graph
    logger.info("Download from osm.")
    _osm_graph = osmnx.graph_from_polygon(
        polygon=region.buffer(buffer)
        .to_crs(pyproj.CRS.from_epsg(4326))
        .geometry.unary_union,
        custom_filter=f'["{osm_key}"~"{"|".join(osm_values)}"]',
        simplify=False,
        retain_all=True,
    )

    # preprocess to desired graph format
    graph = graph_utils.preprocess_graph(_osm_graph, to_crs=crs)

    # # simplify the graph's topology by removing interstitial nodes
    # _osm_simplified_graph = osmnx.simplify_graph(
    #     _osm_graph
    # )  # TODO #64: needs testing, too longe/too short are not beneficial for dem based direction.
    # logger.info("Get simplified graph from osm.")

    # # preprocess to desired graph format
    # graph = graph_utils.preprocess_graph(_osm_simplified_graph, to_crs=crs)

    return graph


def create_graph_from_hydrography(
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


# TODO imcomplete
def create_graph_from_geodataframe(
    gdf: gpd.GeoDataFrame, is_directed: bool = True, logger=logger
) -> nx.Graph:
    """
    Convert a GeoDataFrame of line geometries into a graph.

    This function first converts a GeoDataFrame into a DiGraph using the `gpd_to_digraph` function
    which treats the first and last coordinates of each row's geometry as source and target, respectively.
    Then, it preprocesses the graph using the `preprocess_graph` function which may involve geometry assignments,
    reprojections based on CRS, and graph validations.
    Finally, a digraph or a undirected graph will be returned based on `is_directed`.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        The input GeoDataFrame that is to be converted into a DiGraph.
    is_directed: bool, optional
        Whether the output graph has directions.
        By default, True

    Returns
    -------
    nx.Graph
        The converted and preprocessed graph.
        Directed or undirected given `is_directed`

    Raises
    ------
    ValueError
        If the graph does not have a 'crs' attribute during preprocessing.

    Examples
    --------
    >>> import geopandas as gpd
    >>> import shapely.geometry as geom
    >>> gdf = gpd.GeoDataFrame({'geometry': [geom.LineString([(0,0), (1,1)]), geom.LineString([(1,1), (2,2)])]})
    >>> graph = create_graph_from_geodataframe(gdf)
    >>> print(type(graph))
    <class 'networkx.classes.digraph.DiGraph'>
    """
    logger.info("Preprocessing geodataframes")
    gdf = explode_and_deduplicate_geometries(
        gdf
    )  # TODO more cleanup funcs might be needed

    logger.info("Creating digraph")
    graph = graph_utils.gpd_to_digraph(gdf)  # TODO allow different graph types

    logger.info("Processing graph")
    graph = graph_utils.preprocess_graph(graph)

    if is_directed:
        return graph
    else:
        return graph.to_undirected()


def setup_graph_from_rasterdataset(
    graph: nx.Graph,  # TODO replace by self.graphs
    raster_fn: Union[str, Path, xr.DataArray, xr.Dataset],
    data_catalog: DataCatalog = None,  # TODO replace by self.data_catalog
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
        List of variables to add to graph from raster_fn. By default all.
    fill_method : str, optional
        If specified, fills no data values using fill_nodata method.
        Available methods are {'linear', 'nearest', 'cubic', 'rio_idw'}.
    rename: dict, optional
        Dictionary to rename variable names in raster_fn before adding to graph
        {'name_in_raster_fn': 'name_in_graph'}. By default empty.
    graph_component: str, optional
        Specifies which component of the graph to process. Can be one of the following:
        * "edges" - Only processes and updates the edges of the graph.
        * "nodes" - Only processes and updates the nodes of the graph.
        * "both" - Processes and updates both nodes and edges of the graph.
        By default, it processes both nodes and edges ("both").
    resampling_method: str, optional
        Method to sample from raster data to graph. By default mean. Options include
        {'count', 'min', 'max', 'sum', 'mean', 'std', 'median', 'q##'}.
        Only used when ``graph_component`` is "edges" or "both".
    all_touched : bool, optional
        If True, all pixels touched by geometries will used to define the sample.
        If False, only pixels whose center is within the geometry or that are
        selected by Bresenham's line algorithm will be used. By default True.
        Only used when ``graph_component`` is "edges" or "both".
    Returns
    -------
    graph
        the new graph
    """  # noqa: E501

    assert graph_component in [
        "edges",
        "nodes",
        "both",
    ], "Invalid graph_component value."

    logger.info(f"Preparing graph data from raster source {raster_fn}")
    region = graph_utils.graph_region(graph)

    # Read raster data, select variables, and interpolate na if needed
    if isinstance(raster_fn, str) or isinstance(raster_fn, Path):
        ds = data_catalog.get_rasterdataset(
            raster_fn, region=region, buffer=1000, variables=variables
        )
    else:
        ds = raster_fn
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

    # TODO set graph
    # TODO Convert to UgridDataset
    # uds_sample = xu.UgridDataset(ds_sample, grids=self.mesh_grids[grid_name])

    return graph


def setup_graph_from_geodataframe(
    graph: nx.Graph,  # TODO replace by self.graphs
    vector_fn: Union[str, Path, gpd.GeoDataFrame],
    data_catalog: DataCatalog = None,  # TODO replace by self.data_catalog
    variables: Optional[List[str]] = None,
    max_dist: float = np.inf,
    rename: Optional[Dict[str, str]] = dict(),
    graph_component: Optional[str] = "both",
    logger: logging.Logger = logger,
) -> nx.Graph:
    """Add data variable(s) from ``vector_fn`` to attribute(s) in graph object using nearest_merge.

    Raster data is sampled to the graph edges and nodes using the ``resampling_method``.
    If raster is a dataset, all variables will be added unless ``variables`` list
    is specified.

    Parameters
    ----------
    vector_fn: str, Path, gpd.GeoDataFrame
        Data catalog key, path to vector file or gpd.GeoDataFrame data object.
    variables: list, optional
        List of variables to add to graph from vector_fn. By default all.
    max_dist: float
        Maixmum distance within which the geodataframe would be used to setup graph.
        By default, infinite.
    rename: dict, optional
        Dictionary to rename variable names in vector_fn before adding to graph
        {'name_in_vector_fn': 'name_in_graph'}. By default empty.
    graph_component: str, optional
        Specifies which component of the graph to process. Can be one of the following:
        * "edges" - Only processes and updates the edges of the graph.
        * "nodes" - Only processes and updates the nodes of the graph.
        * "both" - Processes and updates both nodes and edges of the graph.
        By default, it processes both nodes and edges ("both").
    Returns
    -------
    graph
        the new graph
    """  # noqa: E501

    assert graph_component in [
        "edges",
        "nodes",
        "both",
    ], "Invalid graph_component value."

    logger.info(f"Preparing graph data from raster source {vector_fn}")
    region = graph_utils.graph_region(graph)

    # Read vector data, select variables
    if isinstance(vector_fn, gpd.GeoDataFrame):
        gdf = vector_fn
    else:
        gdf = data_catalog.get_geodataframe(
            vector_fn, region=region, buffer=1000, variables=variables
        )

    # relate the geodataframe to graph using nearet join
    # get gdf attributes at edges
    if graph_component in ["edges", "both"]:
        edges = nearest_merge(
            graph_utils.graph_edges(graph), gdf, columns=variables, max_dist=max_dist
        )
        edges = edges.rename(columns=rename)
        graph = graph_utils.network_to_graph(edges, graph_utils.graph_nodes(graph))
        # set graph
    # get sample at nodes
    if graph_component in ["nodes", "both"]:
        nodes = nearest_merge(
            graph_utils.graph_nodes(graph), gdf, columns=variables, max_dist=max_dist
        )
        graph = graph_utils.network_to_graph(graph_utils.graph_edges(graph), nodes)

    # TODO set graph
    # TODO Convert to UgridDataset
    # uds_sample = xu.UgridDataset(ds_sample, grids=self.mesh_grids[grid_name])

    return graph


# func from hybridurb
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


# func from hybridurb network.py
def _find_target_nodes(self, G: nx.DiGraph, node_query=None, edges_query=None):
    """helper to find the target nodes"""

    targets = None

    if all(v is None for v in [node_query, edges_query]):
        targets = [n for n in G.nodes if G.out_degree[n] == 0]

    if isinstance(edges_query, str):
        try:
            _target_G = graph.query_graph_edges_attributes(
                G,
                edge_query=edges_query,
            )
            targets = [v for u, v in _target_G.edges()]
            self.logger.debug("Find targets in edges")
        except Exception as e:
            self.logger.debug(e)
            pass

    if isinstance(node_query, str):
        try:
            _target_G = graph.query_graph_nodes_attributes(
                G,
                node_query=node_query,
            )
            targets = [n for n in _target_G.nodes()]
            self.logger.debug("Find targets in nodes")
        except Exception as e:
            self.logger.debug(e)
            pass

    if targets is None:
        raise ValueError("could not find targets.")

    return targets


def _find_target_graph(self, G: nx.DiGraph, node_query=None, edges_query=None):
    """helper to find the target nodes"""

    targets = None

    if all(v is None for v in [node_query, edges_query]):
        targets = [n for n in G.nodes if G.out_degree[n] == 0]

    if isinstance(edges_query, str):
        try:
            self.logger.debug("Finding targets in edges")
            _target_G = graph.query_graph_edges_attributes(
                G,
                edge_query=edges_query,
            )
            targets = [v for u, v in _target_G.edges()]
        except Exception as e:
            self.logger.debug(e)
            pass

    if isinstance(node_query, str):
        try:
            self.logger.debug("Finding targets in nodes")
            _target_G = graph.query_graph_nodes_attributes(
                G,
                node_query=node_query,
            )
            targets = [n for n in _target_G.nodes()]
        except Exception as e:
            self.logger.debug(e)
            pass

    if targets is None:
        raise ValueError("could not find targets.")

    return _target_G


def setup_dag(
    self,
    G: nx.Graph = None,
    targets=None,
    target_query: str = None,
    weight: str = None,
    loads: list = [],
    report: str = None,
    algorithm: str = "simple",
    **kwargs,
):
    """This component prepares subgraph as Directed Acyclic Graphs (dag) using shortest path
    step 1: add a supernode to subgraph (representing graph - subgraph)
    step 2: use shortest path to prepare dag edges (source: node in subgraph; target: super node)

    in progress

    Parameters
    ----------
    G : nx.Graph
    targets : None or String or list, optional (default = None)
        DAG super targets.
        If None, a target node will be any node with out_degree == 0.
        If a string, use this to query part of graph as a super target node.
        If a list, use targets as a list of target nodes.
    weight : None or string, optional (default = None)
        Weight used for shortest path.
        If None, every edge has weight/distance/cost 1.
        If a string, use this edge attribute as the edge weight.
        Any edge attribute not present defaults to 1.
    loads : None or list of strings, optional (default - None)
        Load used from edges/nodes attributes.
        If None, every node/edge has a load equal to the total number of nodes upstream (nnodes), and number of edges upstream (nedges).
        If a list, use the attributes in the list as load.
        The attribute can be either from edges or from nodes.
        Any attributes that are not present defaults to 0.
    algorithm : string, optional (default = 'simple')
        The algorithm to use to compute the dag.
        Supported options: 'simple', 'flowpath'.
        Other inputs produce a ValueError.
        if 'simple', the input graph is treated as a undirected graph, the resulting graph might alter the original edges direction
        If 'flowpath', the input graph is treated as a directed graph, the resulting graph do not alter the original edges direction.
        Both option will have less edges. the missing edges indicates their insignificance in direction, i.e. both way are possible. # FIXME

    Arguments
    ---------
    use_super_target : bool
        whether to add a super target at the ends of all targets.
        True if the weight of DAG exist for all edges.
        False if the weight of DAG also need to consider attribute specified for targets.

    See Also
    --------
    self._setup_dag

    """
    # convert Digraph to Graph
    if G is None:
        G = self._graphmodel.copy()
        self.logger.debug("Apply dag on graphmodel.")

    # check targets of the dag
    if isinstance(target_query, str):
        targets = self._find_target_nodes(G, target_query, target_query)
    self.logger.debug(f"{len(targets)} targets are selected")

    # remove connected components without targets
    _remain_nodes = set()
    for comp in nx.connected_components(G.to_undirected()):
        if not any(comp.intersection(targets)):
            self.logger.warning(f"Removing nodes disconnected from targets: {comp}.")
        else:
            _remain_nodes.update(comp)
    G = G.subgraph(_remain_nodes)

    # for each connected component with a target
    for comp in nx.connected_components(G.to_undirected()):
        _graph = G.subgraph(comp).copy()
        _targets = comp.intersection(targets)
        self._setup_dag(
            G=_graph,
            targets=_targets,
            weight=weight,
            loads=loads,
            report=report,
            algorithm=algorithm,
        )


def _setup_dag(
    self,
    G: nx.Graph = None,
    targets=None,
    weight: str = None,
    loads: list = [],
    report: str = None,
    algorithm: str = "simple",
    **kwargs,
):
    """This component prepares subgraph as Directed Acyclic Graphs (dag) using shortest path
    step 1: add a supernode to subgraph (representing graph - subgraph)
    step 2: use shortest path to prepare dag edges (source: node in subgraph; target: super node)

    in progress

    Parameters
    ----------
    G : nx.Graph
    targets : None or String or list, optional (default = None)
        DAG super targets.
        If None, a target node will be any node with out_degree == 0.
        If a string, use this to query part of graph as a super target node.
        If a list, use targets as a list of target nodes.
    weight : None or string, optional (default = None)
        Weight used for shortest path.
        If None, every edge has weight/distance/cost 1.
        If a string, use this edge attribute as the edge weight.
        Any edge attribute not present defaults to 1.
    loads : None or list of strings, optional (default - None)
        Load used from edges/nodes attributes.
        If None, every node/edge has a load equal to the total number of nodes upstream (nnodes), and number of edges upstream (nedges).
        If a list, use the attributes in the list as load.
        The attribute can be either from edges or from nodes.
        Any attributes that are not present defaults to 0.
    algorithm : string, optional (default = 'simple')
        The algorithm to use to compute the dag.
        Supported options: 'simple', 'flowpath'.
        Other inputs produce a ValueError.
        if 'simple', the input graph is treated as a undirected graph, the resulting graph might alter the original edges direction
        If 'flowpath', the input graph is treated as a directed graph, the resulting graph do not alter the original edges direction.
        Both option will have less edges. the missing edges indicates their insignificance in direction, i.e. both way are possible. # FIXME

    Arguments
    ----------
    use_super_target : bool
        whether to add a super target at the ends of all targets.
        True if the weight of DAG exist for all edges.
        False if the weight of DAG also need to consider attribute specified for targets.

    """

    # check algorithm of the setup
    if algorithm not in ("simple", "flowpath"):
        raise ValueError(f"algorithm not supported: {algorithm}")
    self.logger.debug(f"Performing {algorithm}")

    # check method of the dag
    method = "dijkstra"
    if method not in ("dijkstra", "bellman-ford"):
        raise ValueError(f"method not supported: {method}")
    self.logger.debug(f"Performing {method} dag")

    # started making dag
    DAG = nx.DiGraph()

    # try adding super nodes
    G.add_edges_from([(n, -1) for n in targets])
    self.logger.debug(f"connecting targets to supernode")
    if len([_ for _ in nx.connected_components(G.to_undirected())]) > 1:
        raise TypeError("Cannot apply dag on disconnected graph.")

    # 1. add path
    # FIXME: if don't do this step: networkx.exception.NetworkXNoPath: No path to **.
    if algorithm == "simple":
        self.logger.info("dag based on method simple. undirected graph will be used. ")
        G = G.to_undirected()
        DAG = graph.make_dag(G, targets=targets, weight=weight, drop_unreachable=False)
    elif algorithm == "flowpath":
        self.logger.info(
            "dag based on method flowpath. unreachable path will be dropped. "
        )
        DAG = graph.make_dag(G, targets=targets, weight=weight, drop_unreachable=True)

    # 2. add back weights
    for u, v, new_d in DAG.edges(data=True):
        # get the old data from X
        old_d = G[u].get(v)
        non_shared = set(old_d) - set(new_d)
        if non_shared:
            # add old weights to new weights if not in new data
            new_d.update(dict((key, old_d[key]) for key in non_shared))

    # 3. add auxiliary calculations
    _ = DAG.reverse()
    nodes_attributes = [k for n in G.nodes for k in G.nodes[n].keys()]
    edges_attribute = [k for e in G.edges for k in G.edges[e].keys()]
    nodes_loads = [l for l in loads if l in nodes_attributes]
    edegs_loads = [l for l in loads if l in edges_attribute]

    for s, t in _.edges():
        # fixed loads
        upstream_nodes = list(nx.dfs_postorder_nodes(_, t))  # exclusive
        upstream_edges = list(G.subgraph(upstream_nodes).edges())
        DAG[t][s].update(
            {"upstream_nodes": upstream_nodes, "upstream_edges": upstream_edges}
        )
        DAG.nodes[s].update(
            {"upstream_nodes": upstream_nodes, "upstream_edges": upstream_edges}
        )
        nnodes = len(upstream_nodes)
        nedges = len(upstream_edges)
        DAG[t][s].update({"nnodes": nnodes, "nedges": nedges})
        DAG.nodes[s].update({"nnodes": nnodes, "nedges": nedges})
        # customised nodes
        sumload_from_nodes = 0
        sumload_from_edges = 0
        for l in loads:
            if l in nodes_loads:
                sumload_from_nodes = np.nansum([G.nodes[n][l] for n in upstream_nodes])
            elif l in edegs_loads:
                sumload_from_edges = np.nansum(
                    [G[e[0]][e[1]][l] for e in upstream_edges]
                )
            else:
                raise KeyError(f"Load {l} does exist in nodes/edges attributes.")
            sumload = sumload_from_nodes + sumload_from_edges
            DAG[t][s].update({l: sumload})
            DAG.nodes[s].update({l: sumload})

    # validate DAG
    if nx.is_directed_acyclic_graph(DAG):
        self.logger.debug("dag is directed acyclic graph")
    else:
        self.logger.error("dag is NOT directed acyclic graph")

    # visualise DAG
    if report:
        # plot graphviz
        graph.make_graphplot_for_targetnodes(DAG, targets, layout="graphviz")
        plt.title(report, wrap=True, loc="left")

        # plot xy
        plt.figure(figsize=(8, 8))
        plt.title(report, wrap=True, loc="left")
        plt.axis("off")
        pos = {xy: xy for xy in G.nodes()}

        # base
        nx.draw(
            DAG,
            pos=pos,
            node_size=0,
            with_labels=False,
            arrows=False,
            node_color="gray",
            edge_color="silver",
            width=0.5,
        )
        nx.draw_networkx_nodes(
            DAG,
            pos=pos,
            nodelist=targets,
            node_size=200,
            node_shape="*",
            node_color="r",
        )
        edge_width = [d[2] / 100 for d in DAG.edges(data="nnodes")]
        nx.draw_networkx_edges(
            DAG,
            pos=pos,
            edgelist=DAG.edges(),
            arrows=False,
            width=[float(i) / (max(edge_width) + 0.1) * 20 + 0.5 for i in edge_width],
        )

    self._graphmodel = DAG
    return DAG


def setup_partition(
    self,
    G: nx.Graph = None,
    subgraph_fn: str = None,
    algorithm: str = "simple",
    report: str = None,
    contracted: bool = False,
    **kwargs,
):
    """This component prepares the partition based on the connectivity of graph

    Parameters
    ----------
    subgraph_fn : str
        Name of the subgraph instance.
        If None, self._graphmodel will be used for partition; the function will update the self._graphmodel
        if String and new, self._graphmodel will be used for partition; the function will create the instance in self._subgraphmodels
        if String and old, self._subgraphmodels[subgraph_fn] will be used for partition; the function will update the instance in self._subgraphmodels[subgraph_fn]
    algorithm : str
        Algorithm to derive partitions from subgraph. Options: ['simple', 'louvain' ]
        testing methods:
        "simple" : based on connected components, every connected components will be considered as a partition (undirected)
        "flowpath" : based on direction of the edges, the graph is devided into a few partitions, each of which represents a target node with all of its sources.
        "louvain": based on louvain algorithm (work in progress)  (undirected). "It measures the relative density of edges inside communities with respect to edges outside communities. Optimizing this value theoretically results in the best possible grouping of the nodes of a given network."(from wikipedia)
    contracted : bool
        Specify whether to build contracted graph from the parititons. So a new contracted graph is created, with a node representing a part of the graph; edges represernting the connectivity between the parts
        If True, each partition will be contracted to a node
        If False, no contraction is performed
    """

    # get graph model
    if G is None:
        G = self._graphmodel.copy()
    G_targets = [n for n in G.nodes if G.out_degree[n] == 0]

    # get subgraph if applicable
    if subgraph_fn in self._subgraphmodels:
        self.logger.warning(
            f"subgraph instance {subgraph_fn} will be used for partition."
        )
        SG = self._subgraphmodels[subgraph_fn].copy()
    elif subgraph_fn is not None:
        self.logger.warning(f"graph will be used for partition.")
        SG = self._graphmodel.copy()
    else:
        self.logger.debug(f"graph will be used for partition.")
        SG = self._graphmodel.copy()
    SG_targets = [n for n in SG.nodes if SG.out_degree[n] == 0]

    # start partition
    partition = {n: -1 for n in G.nodes}
    partition_edges = {e: -1 for e in G.edges}

    if algorithm == "simple":  # based on connected subgraphs
        UG = SG.to_undirected()  # convert SG to UG for further partition
        for i, ig in enumerate(nx.connected_components(UG)):
            ig = UG.subgraph(ig)
            partition.update({n: i for n in ig.nodes})
        partition_edges.update(
            {
                (s, t): partition[s]
                for s, t in partition_edges
                if partition[s] == partition[t]
            }
        )
        logger.info(
            f"algorithm {algorithm} is applied. Note that different partitions are disconnected."
        )
    elif algorithm == "flowpath":
        assert isinstance(
            SG, nx.DiGraph
        ), f"algorithm {algorithm} can only be applied on directional graph"
        SG = graph.sort_direction(SG)
        endnodes = [n[0] for n in SG.nodes(data="_type") if n[-1] == "endnode"]
        partition.update(
            {
                nn: i
                for i, n in enumerate(endnodes)
                for nn in graph.get_predecessors(SG, n)
            }
        )
        partition_edges.update(
            {
                (s, t): partition[s]
                for s, t in partition_edges
                if partition[s] == partition[t]
            }
        )
        logger.info(
            f"algorithm {algorithm} is applied. Note that flowpath might be duplicated."
        )
    elif algorithm == "louvain":  # based on louvain algorithm
        UG = SG.to_undirected()  # convert SG to UG for further partition
        partition.update(graph.louvain_partition(UG))
        partition_edges.update(
            {
                (s, t): partition[s]
                for s, t in partition_edges
                if partition[s] == partition[t]
            }
        )
    else:
        raise ValueError(
            f"{algorithm} is not recognised. allowed algorithms: simple, louvain"
        )

    n_partitions = max(partition.values())
    self.logger.debug(
        f"{n_partitions} partitions is derived from subgraph using {algorithm} algorithm"
    )

    # update partition to graph
    nx.set_node_attributes(SG, partition, "part")
    nx.set_edge_attributes(SG, partition_edges, "part")
    if subgraph_fn in self._subgraphmodels:
        self.logger.warning(
            f"subgraph instance {subgraph_fn} will be updated with partition information (part)."
        )
        self._subgraphmodels[subgraph_fn] = SG
    elif subgraph_fn is not None:
        self.logger.warning(
            f"subgraph instance {subgraph_fn} will be created with partition information (part)."
        )
        self._subgraphmodels[subgraph_fn] = SG
    else:
        self.logger.warning(f"graph will be updated with partition information (part).")
        self._graphmodel = SG

    # contracted graph
    # induced graph from the partitions - faster but a bit confusing results
    # ind = community.induced_graph(partition, G)
    # ind.remove_edges_from(nx.selfloop_edges(ind))
    # induced by contracting - slower but neater
    if contracted == True:
        ind = G.copy()
        nx.set_node_attributes(ind, {n: {"ind_size": 1} for n in ind.nodes})
        for part in np.unique(list(partition.values())):
            part_nodes = [n for n in partition if partition[n] == part]
            if part == -1:
                # do not contract
                pass
            else:
                for to_node in [n for n in part_nodes if n in SG_targets]:
                    ind, targets = graph.contract_graph_nodes(ind, part_nodes, to_node)
                    ind.nodes[to_node]["ind_size"] = len(part_nodes)

    # visualise
    if report:
        # Cartesian coordinates centroid
        pos_G = {xy: xy for xy in G.nodes()}
        pos_SG = {xy: xy for xy in SG.nodes()}

        # partitions
        plt.figure(figsize=(8, 8))
        plt.title(report, wrap=True, loc="left")
        plt.axis("off")
        nx.draw(
            G,
            pos=pos_G,
            node_size=0,
            with_labels=False,
            arrows=False,
            node_color="k",
            edge_color="k",
            width=0.2,
        )
        nx.draw_networkx_nodes(
            G,
            pos=pos_G,
            node_size=30,
            cmap=plt.cm.RdYlBu,
            node_color=list(partition.values()),
        )

        # induced partition graph
        if contracted == True:
            pos_ind = {xy: xy for xy in ind.nodes()}

            graph.make_graphplot_for_targetnodes(ind, None, None, layout="graphviz")
            plt.title(report, wrap=True, loc="left")

            plt.figure(figsize=(8, 8))
            plt.title(report, wrap=True, loc="left")
            plt.axis("off")
            # base
            nx.draw(
                G,
                pos=pos_G,
                node_size=0,
                with_labels=False,
                arrows=False,
                node_color="gray",
                edge_color="silver",
                width=0.2,
            )
            nx.draw_networkx_nodes(
                ind,
                pos=pos_ind,
                node_size=list(dict(ind.nodes(data="ind_size")).values()),
                cmap=plt.cm.RdYlBu,
                node_color=range(len(ind)),
            )
            nx.draw_networkx_edges(ind, pos_ind, alpha=0.3)

    return SG


def setup_pruning(
    self,
    G: nx.Graph = None,
    subgraph_fn: str = None,
    algorithm: str = "simple",
    edge_prune_query: str = None,
    node_prune_query: str = None,
    weight: str = None,
    loads: list = [],
    max_loads: float = [],
    report: str = None,
    **kwargs,
):
    """This component prepares the pruning the 1D flow network based on aggregation

    Parameters
    ----------
    G: nx.Graph, optional (default = None)
        Networkx graph to allow the function being called by another function.
        If None, the network graph from self._graphmodel or self._subgraphmodels[subgraph_fn] will be used
    subgraph_fn : str, optional (default - None)
        Name of the subgraph instance.
        If None, self._graphmodel will be used for partition; the function will update the self._graphmodel
        if String and new, self._graphmodel will be used for partition; the function will create the instance in self._subgraphmodels
        if String and old, self._subgraphmodels[subgraph_fn] will be used for partition; the function will update the instance in self._subgraphmodels[subgraph_fn]
    algorithm : str, optional (default - simple)
        # FIXME: after creating dag, some edges are missing, also the edge attributes. Therefore, the loads are imcomplete. will not have influence if the loads are applied on nodes.
        Algorithm to derive partitions from subgraph. Options: ['simple', 'auto']
        testing methods:
        "simple" : based on user specifie.  The algorithm will prune the edges/nodes based on edge_prune_query and node_prune_query.
        "flowpath": based on direction defined for edges.
        "auto" : based on analysis of network. The algrithm will prune the arborescence of the tree. A threshold can be set to determine whether an arborescence is small enough to be pruned.
    weight : None or string, optional (default = None)
        Weight used for shortest path. Used in setup_dag.
        If None, every edge has weight/distance/cost 1.
        If a string, use this edge attribute as the edge weight.
        Any edge attribute not present defaults to 1.
    loads : None or list of strings, optional (default - None)
        Load used from edges/nodes attributes. Used in setup_dag.
        If None, every node/edge has a load equal to the total number of nodes upstream (nnodes), and number of edges upstream (nedges).
        If a list, use the attributes in the list as load.
        The attribute can be either from edges or from nodes.
        Any attributes that are not present defaults to 0.
    """

    # create the initial graph
    G = self._io_subgraph(subgraph_fn, G, "r")

    # check algorithm
    if algorithm is not None:
        if algorithm not in ("simple", "auto", "flowpath"):
            raise ValueError(f"algorithm not supported: {algorithm}")
    else:
        algorithm = "simple"
        self.logger.info(f"default algorithm {algorithm} is applied")

    # check prune query
    if algorithm == "simple":
        # check queries
        if all(v is None for v in [edge_prune_query, node_prune_query]):
            raise ValueError(
                "must specify one of the following: edge_prune_query OR node_prune_query"
            )
        if max_loads is not None:
            self.logger.debug(f"will ignore: max_loads")
    elif algorithm == "auto":
        if any(v is not None for v in [edge_prune_query, node_prune_query]):
            self.logger.debug(f"will ignore: edge_prune_query, node_prune_query")

    # check loads
    if loads is not None:
        if isinstance(loads, str):
            loads = list(loads)
    else:
        loads = []

    # check max_loads
    if max_loads is not None:
        if isinstance(max_loads, str):
            max_loads = list(max_loads)
            if len(loads) != len(max_loads):
                raise ValueError(f"max_loads must have the same length as loads")
    else:
        max_loads = []

    # get pruned graph
    self.logger.debug(f"getting pruned graph based on algorithm = {algorithm}")
    if algorithm == "auto":
        _PG = graph.get_arborescence(G)
        PG = self.setup_subgraph(
            _PG,
            edge_query="_arborescence == True"
            # report ='plot sg for pruning'
        )

    else:
        PG = self.setup_subgraph(
            G,
            edge_query=edge_prune_query,
            node_query=node_prune_query
            # report ='plot sg for pruning'
        )

    # remained graph
    self.logger.debug(f"getting remained graph based on difference")
    RG = graph.find_difference(G, PG)

    # adding loads to remained graph
    self.logger.debug(f"adding loads to remained graph using dag search")

    # graph connections pruned graph VS remained graph
    tree_roots = [n for n in RG if n in PG]
    tree_roots_edges = {}
    for n in tree_roots:
        keep = []
        dn = [e[2] for e in RG.edges(data="id") if e[0] == n]
        if len(dn) == 0:
            # try upstream
            up = [e[2] for e in RG.edges(data="id") if e[1] == n]
            if len(up) == 0:
                pass
            elif len(up) > 1:
                keep = [up[0]]  # just pick one
            else:
                keep = up
        elif len(dn) == 1:
            keep = dn
        elif len(dn) > 1:
            keep = [dn[0]]  # just pick one
        tree_roots_edges[n] = keep

    if algorithm == "flowpath":
        # no dag is performed
        PG_dag = self.setup_dag(
            PG,
            targets=tree_roots,
            loads=loads,
            weight=weight,
            report="plot dag for pruning",
            algorithm="flowpath",
        )
    else:
        # apply DAG to pruned graph -  get summed information for root nodes
        PG_dag = self.setup_dag(
            PG,
            targets=tree_roots,
            loads=loads,
            weight=weight,
            report="plot dag for pruning",
            algorithm="simple",
        )

    # add new id to PG nodes
    self.logger.debug(f"adding new id to nodes")
    new_id = {}
    for o in tree_roots:
        for n in PG.nodes():
            try:
                if nx.has_path(PG_dag, n, o):
                    new_id[n] = tree_roots_edges[o]
            except nx.exception.NodeNotFound:
                pass

    new_id = {k: v[0] if len(v) > 0 else None for k, v in new_id.items()}
    nx.set_node_attributes(PG, new_id, "new_id")

    # add new id to PG edges
    self.logger.debug(f"adding new id to edges")
    new_id = {}
    for o in tree_roots:
        for e in PG.edges():
            try:
                if nx.has_path(PG_dag, e[0], o):
                    new_id[e] = tree_roots_edges[o]
            except nx.exception.NodeNotFound:
                pass

    new_id = {k: v[0] if len(v) > 0 else None for k, v in new_id.items()}
    nx.set_edge_attributes(PG, new_id, "new_id")

    # add PG_dag loads back to tree roots in RG
    # add back missing edges (adding missing edges for weight, might alter flow path indeed)
    self.logger.debug(f"adding missing edges")
    for e in PG.edges:
        if e not in PG_dag.edges():
            if e not in PG_dag.reverse().edges():
                PG_dag.add_edges_from(e, **PG.edges[e[0], e[1]])

    # sum loads of the arborescence
    self.logger.debug(f"adding loads")
    for load in loads:
        for o in tree_roots:
            arborescence_graph = PG_dag.subgraph(
                graph.get_predecessors(PG_dag, o, inclusive=True)
            )
            load_from_nodes = sum(
                [
                    e[-1]
                    for e in arborescence_graph.edges(data=load)
                    if e[-1] is not None
                ]
            )
            loads_from_edges = sum(
                [
                    n[-1]
                    for n in arborescence_graph.nodes(data=load)
                    if n[-1] is not None
                ]
            )
            nx.set_node_attributes(RG, {o: load_from_nodes + loads_from_edges}, load)
            nx.set_node_attributes(RG, {o: len(arborescence_graph.nodes())}, "nnodes")
            nx.set_node_attributes(RG, {o: len(arborescence_graph.edges())}, "nedges")
    loads = set(loads + ["nnodes", "nedges"])

    # write into graph
    self._io_subgraph(subgraph_fn + "_RG", RG, "w")
    self._io_subgraph(subgraph_fn + "_PG", PG, "w")

    # # map arborescence to a node
    # self.logger.debug(f"adding mapping to remained graph for missing nodes")
    # for n in tree_roots:
    #     uns = graph.get_predecessors(G, n, inclusive=False) # upstream nodes
    #     ues = list(G.subgraph(uns + [n]).edges())
    #     nx.set_node_attributes(G, {un:G.nodes[n]['id'] for un in uns}, 'mapping_id')
    #     nx.set_edge_attributes(G, {ue:G.nodes[n]['id'] for ue in ues}, 'mapping_id')
    #
    # self._io_subgraph(subgraph_fn + '_mapping', G, 'w')

    # draw to see
    if report:
        # plot xy
        plt.figure(figsize=(8, 8))
        plt.title(report, wrap=True, loc="left")
        plt.axis("off")
        pos = {xy: xy for xy in G.nodes()}

        # base
        nx.draw(
            G,
            pos=pos,
            node_size=0,
            with_labels=False,
            arrows=False,
            node_color="gray",
            edge_color="silver",
            width=0.2,
        )
        nx.draw(
            RG,
            pos=pos,
            node_size=0,
            with_labels=False,
            arrows=False,
            node_color="gray",
            edge_color="gray",
            width=1,
        )

        # plot size of arborescence
        size = list(dict(RG.nodes(data=list(loads)[0])).values())
        size = np.array([0 if v is None else v for v in size])
        scaled_size = np.interp(size, (size.min(), size.max()), (0, 100))

        nx.draw_networkx_nodes(RG, pos, node_size=scaled_size, node_color="r")

    return None


def setup_diagram(
    self,
    G: nx.DiGraph = None,
    subgraph_fn: str = None,
    target_query: str = None,
    **kwargs,
):
    """function to derive flow diagram at certain target nodes/edges (in progress)

    The target nodes/edges include:
        nodes/edges identified using target_query
        nodes that are end nodes (out degree is 0)

    Arguments
    ---------
    G: nx.DiGraph
        Directional Graph
    subgraph_fn:
        Directional Graph node that are used as target to find predecessors
    target_query: bool
        Whether to include the input node in the results
    **kwargs:
        weight = None, for setup_dag # FIXME
        report = None, for plotting
    """

    # check arguments
    if G is not None:
        _G = G.copy()
        self.logger.info(f"Performing on given graph")
    else:
        _G = self._io_subgraph(subgraph_fn, G, "r")
        self.logger.info(f"Performing on subgraph instance {subgraph_fn}")

    if target_query is None:
        raise ValueError(
            "must specify target_query (diagram will be aggregated base don target_query)"
        )

    weight = kwargs.get("weight", None)

    # 1. find target nodes/edges and put them in a graph
    target_G = self._find_target_graph(_G, target_query, target_query)

    # 2. setup a subgraph instance without the target graph
    SG = _G.copy()
    SG.remove_edges_from(target_G.edges)

    # 3. setup a dag to delineate the upstream of the targets --> results is a dag split at target nodes
    targets = list(
        set(
            [n for n in target_G.nodes if target_G.in_degree[n] == 0]
            + [n for n in SG.nodes if SG.out_degree[n] == 0]
        )
    )
    G_dag = graph.make_dag(SG, targets=targets, weight=weight)
    self._graphmodel = G_dag.copy()
    # TODO: this method follows the flow direction specified in the DiGraph,
    #  there fore is incapable of making correction to any flow direcitons.
    #  to be implemented in setup_dag

    # 4. setup partition --> results is a graph with partition information as attributes
    G_part = self.setup_partition(G_dag, algorithm="simple")

    # 5. add target graph back --> result is the graph has reconstructed connecitivity
    G_part.add_edges_from(target_G.edges)

    # 6. contract graph based on partition to form tree diagram
    partition = pd.DataFrame.from_dict(G_part.nodes(data="part"))
    partition = partition.dropna().set_index(0)
    partition = partition.to_dict(orient="dict")[1]
    ind = graph.contract_graph(G_part, partition=partition, tonodes=targets)
    # TODO: improve this function in workflows, remove from setup partition

    # plot
    report = kwargs.get("report", "setup_diagram")

    if report is not None:
        G = G_part.copy()
        pos_G = {xy: xy for xy in G.nodes()}
        pos_ind = {xy: xy for xy in ind.nodes()}

        graph.make_graphplot_for_targetnodes(ind, None, None, layout="graphviz")
        plt.title(report, wrap=True, loc="left")

        plt.figure(figsize=(8, 8))
        plt.title(report, wrap=True, loc="left")
        plt.axis("off")
        # base
        nx.draw(
            G,
            pos=pos_G,
            node_size=0,
            with_labels=False,
            arrows=False,
            node_color="gray",
            edge_color="silver",
            width=0.2,
        )
        nx.draw_networkx_nodes(
            ind,
            pos=pos_ind,
            node_size=list(dict(ind.nodes(data="ind_size")).values()),
            cmap=plt.cm.RdYlBu,
            node_color=range(len(ind)),
        )
        nx.draw_networkx_edges(ind, pos_ind, alpha=0.3)

    return


# func from hybridurb gis_utils.py
# create graph
def add_edges_with_id(
    G: nx.Graph, edges: gpd.GeoDataFrame, id_col: str, snap_offset: float = 1e-6
) -> nx.Graph():
    """Return graph with edges and edges ids"""

    for index, row in edges.iterrows():
        from_node = row.geometry.coords[0]
        to_node = row.geometry.coords[-1]

        G.add_edge(from_node, to_node, id=row[id_col])

    return G


def add_nodes_with_id(
    G: nx.Graph, nodes: gpd.GeoDataFrame, id_col: str, snap_offset: float = 1e-6
) -> nx.Graph():
    """return graph with nodes and nodes ids"""

    # derived nodes and user nodes
    G_nodes = gpd.GeoDataFrame(
        {
            "geometry": [Point(p) for p in G.nodes],
            "id": [f"{p[0]:.6f}_{p[1]:.6f}" for p in G.nodes],
            "tuple": G.nodes,
        }
    ).set_index("id")
    if "id" in nodes.columns:
        logger.error("Abort: nodes contains id columns. Please remove the column.")
    F_nodes = nodes.rename(columns={id_col: "id"}).set_index("id")

    # map user nodes to derived nodes
    if set(F_nodes.index).issubset(set(G_nodes.index)):
        # check if 1-to-1 match
        G_nodes_new = G_nodes.join(F_nodes.drop(columns="geometry"))

    else:
        G_nodes_new = snap_nodes_to_nodes(F_nodes, G_nodes, snap_offset)
        logger.debug("performing snap nodes to graph nodes")

    # assign nodes id to graph
    dict = {row.tuple: i for i, row in G_nodes_new.iterrows()}
    nx.set_node_attributes(G, dict, "id")

    return G


def update_edges_attributes(
    G: nx.Graph,
    edges: gpd.GeoDataFrame,
    id_col: str,
) -> nx.Graph():
    """This function updates the graph by adding new edges attributes specified in edges"""

    # graph df
    _graph_df = nx.to_pandas_edgelist(G).set_index("id")

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
    _graph_df = _graph_df.reindex(
        columns=_graph_df.columns.union(edges.columns, sort=False)
    )
    graph_df = pd.concat([_graph_df, edges]).groupby(level=0).last()
    graph_df = graph_df.loc[_graph_df.index]

    G_updated = nx.from_pandas_edgelist(
        graph_df.reset_index(),
        source="source",
        target="target",
        edge_attr=True,
        create_using=type(G),
    )

    return G_updated


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


def update_nodes_attributes(
    G: nx.Graph,
    nodes: gpd.GeoDataFrame,
    id_col: str,
) -> nx.Graph():
    """This function updates the graph by adding new edges attributes specified in edges"""

    # graph df
    _graph_df = pd.DataFrame(G.nodes(data="id"), columns=["tuple", "id"]).set_index(
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
        nx.set_node_attributes(G, dict, c)

    return G


# process graph


def query_graph_edges_attributes(G, id_col: str = "id", edge_query: str = None):
    """This function queries the graph by selecting only the edges specified in edge_query"""

    if edge_query is None:
        G_query = G.copy()
    else:
        _graph_df = nx.to_pandas_edgelist(G).set_index(id_col)
        keep_df = _graph_df.query(edge_query)

        if len(keep_df) > 0:
            G_query = G.edge_subgraph(
                [(row.source, row.target) for row in keep_df.itertuples()]
            ).copy()

        else:
            raise ValueError("edges_query results in nothing left")

    return G_query


def query_graph_nodes_attributes(G, id_col: str = "id", node_query: str = None):
    """This function queries the graph by selecting only the nodes specified in node_query"""

    if node_query is None:
        G_query = G

    else:
        _graph_df = pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient="index")
        graph_df = _graph_df.query(node_query)

        if len(graph_df) != 0:
            G_query = G.subgraph(graph_df.index.tolist()).copy()
        else:
            raise ValueError("node_query results in nothing left")

    return G_query


def contract_graph_nodes(G, nodes, to_node=None):
    """This function contract the nodes into one node in G"""

    G_contracted = G.copy()
    node_contracted = []

    if len(nodes) > 1:
        nodes = sorted(nodes)
        if not to_node:  # use the first node
            to_node = nodes[0]
        nodes = [n for n in nodes if n != to_node]
        node_contracted.append(to_node)
        for node in nodes:
            G_contracted = nx.contracted_nodes(
                G_contracted, to_node, node, self_loops=False, copy=False
            )

    return G_contracted, node_contracted


def louvain_partition(G: nx.Graph):
    """This function is a wrapper around best partiton method in community
    See :py:meth:`~community.best_partition()` for more information.
    """
    return community.best_partition(G)


def sort_ids(G: nx.Graph):
    """Function to sort the ids of the graph.
    if there are no ids for the nodes, the ids will be generated based on the geometry
    """
    if set(dict(G.nodes(data="id")).values()) == {None}:
        nx.set_node_attributes(G, {p: f"{p[0]:.6f}_{p[1]:.6f}" for p in G.nodes}, "id")

    return G


def sort_ends(G: nx.Graph):
    """Function to sort the ends of the graph.

    Arguments
    ---------
    G: nx.Graph
        Networkx Graph

    Returns
    -------
    G: nx.Graph
        Networkx Graph with node attribute '_type'
    """
    if isinstance(G, nx.DiGraph):
        endnodes = {
            n: "endnode" for n in G.nodes if (G.degree[n] == 1 and G.out_degree[n] == 0)
        }
        startnodes = {
            n: "startnode"
            for n in G.nodes
            if (G.degree[n] == 1 and G.in_degree[n] == 0)
        }
        nx.set_node_attributes(G, endnodes, "_type")
        nx.set_node_attributes(G, startnodes, "_type")
    elif isinstance(G, nx.Graph):
        endnodes = {n: "endnode" for n in G.nodes if G.degree[n] == 1}
        nx.set_node_attributes(G, endnodes, "_type")
    return G


def sort_direction(G: nx.DiGraph) -> nx.DiGraph:
    """Function sort the start end direction of the graph and obtain start and end nodes.

    Arguments
    ---------
    G: nx.DiGraph
        Directional Graph

    Returns
    -------
    G: nx.DiGraph
        Directional Graph with node attributes endnodes and startnodes"""

    endnodes = {
        n: "endnode" for n in G.nodes if (G.degree[n] == 1 and G.out_degree[n] == 0)
    }
    startnodes = {
        n: "startnode" for n in G.nodes if (G.degree[n] == 1 and G.in_degree[n] == 0)
    }
    nx.set_node_attributes(G, endnodes, "_type")
    nx.set_node_attributes(G, startnodes, "_type")
    return G


def get_predecessors(G: nx.DiGraph, n, inclusive=True):
    """Function to find the predecessors of a node n
    See :py:meth:`~nx.bfs_predecessors()` for more information.

    Arguments
    ---------
    G: nx.DiGraph
        Directional Graph
    n:
        Directional Graph node that are used as target to find predecessors
    inclusive: bool
        Whether to include the input node in the results
    """

    RG = G.reverse()
    predecessors = list(dict(nx.bfs_predecessors(RG, n)).keys())
    if inclusive:
        predecessors = predecessors + [n]
    return predecessors


def find_difference(G, H):
    """function to find the difference between G and H (G-H) based on edges
    replace :py:meth:`~nx.difference()`"""
    c = G.copy()
    c.remove_edges_from(H.edges)
    c.remove_nodes_from(list(nx.isolates(c)))
    return c


def contract_graph(G: nx.Graph, partition, tonodes):
    """contract based on partition --> needs further improvements
    TODO: harmonize with setup partition
    """
    ind = G.copy()
    nx.set_node_attributes(ind, {n: {"ind_size": 1} for n in ind.nodes})
    for part in np.unique(list(partition.values())):
        part_nodes = [n for n in partition if partition[n] == part]
        if part == -1:
            # do not contract
            pass
        else:
            for to_node in [n for n in part_nodes if n in tonodes]:
                ind, targets = contract_graph_nodes(ind, part_nodes, to_node)
                ind.nodes[to_node]["ind_size"] = len(part_nodes)
    return ind


def make_dag(
    G: nx.DiGraph,
    targets: list,
    weight: str = None,
    algorithm="dijkstra",
    drop_unreachable=False,
    logger=logger,
):
    """dag making for digraph --> needs further improvements
    TODO: add to setup_dag
    # test
    # G = nx.DiGraph()
    # G.add_edge(0,1)
    # G.add_edge(1,2)
    # G.add_edge(1,3)
    # nx.draw(G)
    # targets = [3, 2]
    # G.add_edges_from ([(o, -1) for o in targets])
    # nx.draw(G)
    # X_new = nx.DiGraph()
    # for n in G.nodes:
    #     path = nx.shortest_path(G, n, -1,
    #              method = algorithm)
    #     nx.add_path(X_new, path)
    # X_new.remove_nodes_from([-1])
    """

    # copy
    X = G.copy()

    X_new = nx.DiGraph()
    X.add_edges_from([(o, -1) for o in targets])
    # for nodes that are reachable from target
    if drop_unreachable:
        XX = X.reverse()
        sub_network = X.subgraph(nx.dfs_tree(XX, -1).nodes)
        X = sub_network.copy()
        logger.debug("drop unreachable nodes")

    for n in X.nodes:
        path = nx.shortest_path(X, n, -1, weight=weight, method=algorithm)
        nx.add_path(X_new, path)

    X_new.remove_nodes_from([-1])

    # add back weights
    for u, v, weight in X_new.edges(data=True):
        # get the old data from X
        xdata = X[u].get(v)
        non_shared = set(xdata) - set(weight)
        if non_shared:
            # add old weights to new weights if not in new data
            weight.update(dict((key, xdata[key]) for key in non_shared))

    return X_new


# plot graph


def make_graphplot_for_targetnodes(
    G: nx.DiGraph,
    target_nodes: list,
    target_nodes_labeldict: dict = None,
    layout="xy",
    ax=None,
):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # layout graphviz
    if layout == "graphviz":
        # get position
        pos = graphviz_layout(G, prog="dot", args="")

        # draw network
        nx.draw(G, pos, node_size=20, alpha=0.5, node_color="blue", with_labels=False)

        if target_nodes_labeldict is not None:
            # draw labels
            nx.draw_networkx_labels(
                G, pos, target_nodes_labeldict, font_size=16, font_color="k"
            )

    # layout xy
    elif layout == "xy":
        # get position
        pos = {xy: xy for xy in G.nodes()}
        # make plot for each target node
        RG = G.reverse()

        for target in target_nodes:
            c = random_color()

            # make target upstream a graph
            target_G = G.subgraph(
                list(dict(nx.bfs_predecessors(RG, target)).keys()) + [target]
            )

            # draw graph
            nx.draw_networkx(
                target_G,
                pos,
                node_size=10,
                node_color=[c],
                width=2,
                edge_color=[c],
                with_labels=False,
                ax=ax,
            )

            # draw outlets
            nx.draw_networkx_nodes(
                target_G,
                pos,
                nodelist=[target],
                node_size=100,
                node_color="k",
                edgecolors=c,
                ax=ax,
            )

            # draw labels
            if target_nodes_labeldict is not None:
                nx.draw_networkx_labels(
                    target_G, pos, target_nodes_labeldict, font_size=16, font_color="k"
                )
    return fig, ax


def plot_xy(G: nx.DiGraph, plot_outfall=False):
    plt.figure(figsize=(8, 8))
    plt.axis("off")
    pos_G = {xy: xy for xy in G.nodes()}
    nx.draw_networkx_nodes(G, pos=pos_G, node_size=1, node_color="k")
    nx.draw_networkx_edges(G, pos=pos_G, edge_color="k", arrows=False)

    if plot_outfall:
        if isinstance(G, nx.DiGraph):
            endnodes = [n for n in G.nodes if G.out_degree[n] == 0 and G.degree[n] == 1]
            startnodes = [
                n for n in G.nodes if G.in_degree[n] == 0 and G.degree[n] == 1
            ]
            nx.draw_networkx_nodes(
                endnodes, pos=pos_G, node_size=10, node_color="r", node_shape="o"
            )
            nx.draw_networkx_nodes(
                startnodes, pos=pos_G, node_size=10, node_color="b", node_shape="v"
            )
        else:
            pass  # can only be used in case of digraph
    return


def plot_graphviz(G: nx.DiGraph):
    """This function makes plots for grahviz layout"""

    # convert to undirected graph
    UG = G.to_undirected()

    # find connected components in undirected graph
    outlets = []
    for i, SG in enumerate(nx.connected_components(UG)):
        # make components a subgraph
        SG = G.subgraph(SG)

        # find outlets of the subgraph
        outlets.append([n for n in SG.nodes() if G.out_degree(n) == 0])

    outlets = sum(outlets, [])
    outlet_ids = {}

    fig1, ax1 = make_graphplot_for_targetnodes(
        G, outlets, outlet_ids, layout="graphviz"
    )
    return (fig1, ax1)


def plot_graph(G: nx.DiGraph):
    """This function makes plots for two different layout"""

    # convert to undirected graph
    UG = G.to_undirected()

    # find connected components in undirected graph
    outlets = []
    for i, SG in enumerate(nx.connected_components(UG)):
        # make components a subgraph
        SG = G.subgraph(SG)

        # find outlets of the subgraph
        outlets.append([n for n in SG.nodes() if G.out_degree(n) == 0])

    outlets = sum(outlets, [])
    outlet_ids = {}

    fig1, ax1 = make_graphplot_for_targetnodes(
        G, outlets, outlet_ids, layout="graphviz"
    )
    fig2, ax2 = make_graphplot_for_targetnodes(G, outlets, outlet_ids, layout="xy")

    return (fig1, ax1), (fig2, ax2)


def random_color():
    return tuple(
        [
            random.randint(0, 255) / 255,
            random.randint(0, 255) / 255,
            random.randint(0, 255) / 255,
        ]
    )


# validate graph - old


def validate_1dnetwork_connectivity(
    branches: gpd.GeoDataFrame,
    plotit=False,
    ax=None,
    exportpath=os.getcwd(),
    logger=logging,
):
    """Function to validate the connectivity of provided branch"""

    # affirm datatype
    branches = gpd.GeoDataFrame(branches)

    # create digraph
    G = create_graph_from_branches(branches)
    pos = {xy: xy for xy in G.nodes()}

    # convert to undirected graph
    UG = G.to_undirected()

    # find connected components in undirected graph
    outlets = []
    for i, SG in enumerate(nx.connected_components(UG)):
        # make components a subgraph
        SG = G.subgraph(SG)

        # find outlets of the subgraph
        outlets.append([n for n in SG.nodes() if G.out_degree(n) == 0])

    outlets = sum(outlets, [])
    outlet_ids = {
        p: [li for li, l in branches.geometry.iteritems() if l.intersects(Point(p))]
        for p in outlets
    }

    # report
    if i == 0:
        logger.info(
            "Validation results: the 1D network are fully connected.  Supress plotit function."
        )
    else:
        logger.info(
            f"Validation results: the 1D network are disconnected have {i+1} connected components"
        )

    if plotit:
        ax = make_graphplot_for_targetnodes(G, outlets, outlet_ids, layout="graphviz")
        ax.set_title(
            "Connectivity of the 1d network, with outlets"
            + "(connectivity outlets, not neccessarily network outlets due to bi-directional flow, please check these)",
            wrap=True,
        )
        plt.savefig(exportpath.joinpath("validate_1dnetwork_connectivity"))

    return None


def validate_1dnetwork_flowpath(
    branches: gpd.GeoDataFrame,
    branchType_col="branchType",
    plotit=False,
    ax=None,
    exportpath=os.getcwd(),
    logger=logging,
):
    """function to validate flowpath (flowpath to outlet) for provided branch"""

    # affirm datatype
    branches = gpd.GeoDataFrame(branches)

    # create digraph
    G = gpd_to_digraph(branches)
    pos = {xy: xy for xy in G.nodes()}

    # create separate graphs for pipes and branches
    pipes = branches.query(f"{branchType_col} == 'Pipe'")
    channels = branches.query(f"{branchType_col} == 'Channel'")

    # validate 1d network based on pipes -> channel logic
    if len(pipes) > 0:
        # create graph
        PG = gpd_to_digraph(pipes)
        # pipes outlets
        pipes_outlets = [n for n in PG.nodes() if G.out_degree(n) == 0]
        pipes_outlet_ids = {
            p: [li for li, l in pipes.geometry.iteritems() if l.intersects(Point(p))]
            for p in pipes_outlets
        }
        logger.info(
            f"Validation result: the 1d network has {len(pipes_outlets)} pipe outlets."
        )

    if len(channels) > 0:
        # create graph
        CG = gpd_to_digraph(channels)
        # pipes outlets
        channels_outlets = [n for n in CG.nodes() if G.out_degree(n) == 0]
        channels_outlet_ids = {
            p: [li for li, l in channels.geometry.iteritems() if l.intersects(Point(p))]
            for p in channels_outlets
        }
        logger.info(
            f"Validation result: the 1d network has {len(channels_outlets)} channel outlets."
        )

    if (len(channels) > 0) and (len(pipes) > 0):
        isolated_outlets = [
            p
            for p in pipes_outlets
            if not any(Point(p).intersects(l) for _, l in channels.geometry.iteritems())
        ]
        isolated_outlet_ids = {}
        for p in isolated_outlets:
            isolated_outlet_id = [
                li for li, l in pipes.geometry.iteritems() if l.intersects(Point(p))
            ]
            isolated_outlet_ids[p] = isolated_outlet_id
            logger.warning(
                f"Validation result: downstream of {isolated_outlet_id} are not located on channels. Please double check. "
            )

    # plot
    if plotit:
        ax = make_graphplot_for_targetnodes(
            G,
            target_nodes={**isolated_outlet_ids, **channels_outlet_ids}.keys(),
            target_nodes_labeldict={**isolated_outlet_ids, **channels_outlet_ids},
        )
        ctx.add_basemap(
            ax=ax, url=ctx.providers.OpenStreetMap.Mapnik, crs=branches.crs.to_epsg()
        )
        ax.set_title(
            "Flow path of the 1d network, with outlets"
            + "(flowpath outlets, not neccessarily network outlets due to bi-directional flow , please check these)",
            wrap=True,
        )
        plt.savefig(exportpath.joinpath("validate_1dnetwork_flowpath"))

    return None


def get_arborescence(G: nx.DiGraph):
    """function to get arborescence from Digraph
    This function will loop through all bifurcation node and check if its predecessors forms a arborescence.
    If yes, _arborescence = True is assigned to nodes and edges.
    See :py:meth:`networkx.algorithms.tree.recognition.is_arborescence` for more.
    """

    if not isinstance(G, nx.DiGraph):
        raise TypeError("Must be a DiGraph")

    # prepare
    G_neutral = G.to_undirected()
    G_positive = G.copy()
    G_negative = G_positive.reverse()

    # get all bifurcation node
    bifurcations = [n for n in G.nodes if G.degree[n] > 2]

    # get edges that can be pruned and its load
    for n in bifurcations:
        # get its upstream as a subgraph (based on edges not nodes)
        n_predecessors = nx.dfs_predecessors(G_negative, n)
        if len(n_predecessors) > 1:
            subgraph_edges = []
            for nn in n_predecessors:
                _ = list(itertools.chain(*list(G_neutral.edges(nn))))
                subgraph_edges.append(_)
            subgraph_nodes = list(set(sum(subgraph_edges, [])))

            subgraph = G_negative.subgraph(subgraph_nodes)
            subgraph_edges = list(subgraph.edges())

            # check if its upstream subgraph is arborescence
            if nx.is_arborescence(subgraph):
                nx.set_node_attributes(
                    G, {k: True for k in subgraph_nodes}, "_arborescence"
                )
                nx.set_edge_attributes(
                    G, {e: True for e in subgraph_edges}, "_arborescence"
                )

    return G
