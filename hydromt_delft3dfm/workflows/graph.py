"""Workflows to prepare graphs from different sources (hydromt networkmixins + other workflows)."""  # noqa: E501

import itertools
import logging
import random
from pathlib import Path
from typing import Dict, List, Optional, Union

import geopandas as gpd
import matplotlib.pyplot as plt  # FIXME maybe remove
import networkx as nx
import numpy as np
import osmnx
import pandas as pd
import pyproj
import xarray as xr
from hydromt import DataCatalog, flw
from hydromt.gis_utils import nearest_merge

from hydromt_delft3dfm import graph_utils
from hydromt_delft3dfm.workflows import explode_and_deduplicate_geometries

logger = logging.getLogger(__name__)


__all__ = [
    # create graph topology
    "create_graph_from_geodataframe",
    "create_graph_from_openstreetmap",
    "create_graph_from_hydrography",
    # setup graph attributes
    "setup_graph_from_rasterdataset",
    "setup_graph_from_geodataframe",
    # workflows
    "query_graph_edges_attributes",
    "query_graph_nodes_attributes",
    "get_largest_component",
    "add_missing_edges_to_subgraph",
    "setup_dag",
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
        OSM "values" to specify attributes assigned to a key to further detail or
        describe a map feature.
        Together, a key-value pair provides descriptive information about OSM entities.
        Some common values associated to the common keys include:
            highway: motorway,motorway_link,primary,primary_link,secondary,
            secondary_link,tertiary,tertiary_link,residential
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
    # )  # TODO #64: needs testing,
    # too longe/too short are not beneficial for dem based direction.
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

    This function first converts a GeoDataFrame into a DiGraph using the
    `gpd_to_digraph` function, which treats the first and last coordinates
    of each row's geometry as source and target, respectively.
    Then, it preprocesses the graph using the `preprocess_graph` function,
    which may involve geometry assignments, reprojections based on CRS,
    and graph validations.
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
    >>> gdf = gpd.GeoDataFrame({'geometry': [geom.LineString([(0,0), (1,1)]),
    geom.LineString([(1,1), (2,2)])]})
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
        # FIXME dem for sample data too small, need to add nodata and correct crs
        ds_sample = ds.raster.sample(nodes)
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


# workflows
def get_largest_component(G):
    """
    Get the largest connected component from a network graph.

    Parameters
    ----------
        G (nx.Graph): The networkx graph from which to find the largest component.

    Returns
    -------
        nx.Graph: Subgraph representing the largest connected component.
    """
    UG = G.to_undirected()
    largest_cc = max(nx.connected_components(UG), key=len)
    return G.subgraph(largest_cc)


def add_missing_edges_to_subgraph(
    original_subgraph: nx.Graph, complete_graph: nx.Graph
) -> nx.Graph:
    """
    Add missing edges to a subgraph using a complete graph.

    This function identifies and adds the edges from the complete graph to the subgraph
    that are missing but share the same vertices as the subgraph. It preserves the
    original attributes of the edges and nodes.

    Parameters
    ----------
        original_subgraph (nx.Graph):
        The subgraph to which missing edges are to be added.
        complete_graph (nx.Graph):
        The complete graph used as a reference to find missing edges.

    Returns
    -------
        nx.Graph: The updated subgraph with missing edges added.
    """
    # Copy the original subgraph to avoid modifying it directly
    subgraph_copy = original_subgraph.copy()

    # Retrieve nodes from the original subgraph
    subgraph_nodes = subgraph_copy.nodes

    # Convert the complete graph to an undirected graph for processing
    undirected_complete_graph = complete_graph.to_undirected()

    # Create a new subgraph from the undirected complete graph
    # using the nodes of the original subgraph
    new_subgraph = undirected_complete_graph.__class__()
    new_subgraph.add_nodes_from(subgraph_copy.nodes(data=True))
    new_subgraph.add_edges_from(subgraph_copy.edges(data=True))

    # Add edges that has a note in the subgraph and neighbour not in subgraph
    if new_subgraph.is_multigraph():
        new_subgraph.add_edges_from(
            (node, neighbor, key, data)
            for node, neighbors in undirected_complete_graph.adj.items()
            if node in subgraph_nodes
            for neighbor, keydict in neighbors.items()
            if neighbor not in subgraph_nodes
            for key, data in keydict.items()
        )
    else:
        new_subgraph.add_edges_from(
            (node, neighbor, data)
            for node, neighbors in undirected_complete_graph.adj.items()
            if node in subgraph_nodes
            for neighbor, data in neighbors.items()
            if neighbor not in subgraph_nodes
        )
    # update new nodes properties
    nx.set_node_attributes(
        new_subgraph,
        {
            n[0]: n[1]
            for n in undirected_complete_graph.nodes(data=True)
            if n[0] in new_subgraph.nodes
        },
    )
    # Update the graph properties from the undirected complete graph
    new_subgraph.graph.update(undirected_complete_graph.graph)

    return new_subgraph


# func from hybridurb
def update_edges_attributes(
    graph: nx.Graph,
    edges: gpd.GeoDataFrame,
    id_col: str = "id",
) -> nx.Graph():
    """Update the graph by adding new edges attributes specified in edges.

    Only edges with matching ids specified in "id" will be updated.
    """
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
    """Update the graph by adding new edges attributes specified in edges.

    Only edges with matching ids specified in "id" will be updated.
    """
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
def _find_target_nodes(G: nx.DiGraph, node_query=None, edges_query=None):
    """Find the target nodes."""
    targets = None

    if all(v is None for v in [node_query, edges_query]):
        targets = [n for n in G.nodes if G.out_degree[n] == 0]

    if isinstance(edges_query, str):
        try:
            _target_G = query_graph_edges_attributes(
                G,
                edge_query=edges_query,
            )
            targets = [v for u, v in _target_G.edges()]
            logger.debug("Find targets in edges")
        except Exception as e:
            logger.debug(e)
            pass

    if isinstance(node_query, str):
        try:
            _target_G = query_graph_nodes_attributes(
                G,
                node_query=node_query,
            )
            targets = [n for n in _target_G.nodes()]
            logger.debug("Find targets in nodes")
        except Exception as e:
            logger.debug(e)
            pass

    if targets is None:
        raise ValueError("could not find targets.")

    return targets


def _find_target_graph(G: nx.DiGraph, node_query=None, edges_query=None):
    """Find the target nodes."""
    targets = None

    if all(v is None for v in [node_query, edges_query]):
        targets = [n for n in G.nodes if G.out_degree[n] == 0]

    if isinstance(edges_query, str):
        try:
            logger.debug("Finding targets in edges")
            _target_G = query_graph_edges_attributes(
                G,
                edge_query=edges_query,
            )
            targets = [v for u, v in _target_G.edges()]
        except Exception as e:
            logger.debug(e)
            pass

    if isinstance(node_query, str):
        try:
            logger.debug("Finding targets in nodes")
            _target_G = query_graph_nodes_attributes(
                G,
                node_query=node_query,
            )
            targets = [n for n in _target_G.nodes()]
        except Exception as e:
            logger.debug(e)
            pass

    if targets is None:
        raise ValueError("could not find targets.")

    return _target_G


def setup_dag(
    G: nx.Graph,
    targets=None,
    target_query: str = None,
    weight: str = None,
    loads: list = [],
    report: str = None,
    algorithm: str = "simple",
    **kwargs,
):
    """Prepare subgraph as Directed Acyclic Graphs (dag) using shortest path.

    step 1: add a supernode to subgraph (representing graph - subgraph)
    step 2: use shortest path to prepare dag edges
    (source: node in subgraph; target: super node).
    in progress.

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
        If None, every node/edge has a load equal to the total number of
        nodes upstream (nnodes), and number of edges upstream (nedges).
        If a list, use the attributes in the list as load.
        The attribute can be either from edges or from nodes.
        Any attributes that are not present defaults to 0.
    algorithm : string, optional (default = 'simple')
        The algorithm to use to compute the dag.
        Supported options: 'simple', 'flowpath'.
        Other inputs produce a ValueError.
        if 'simple', the input graph is treated as a undirected graph,
        the resulting graph might alter the original edges direction
        If 'flowpath', the input graph is treated as a directed graph,
        the resulting graph do not alter the original edges direction.
        Both option will have less edges. the missing edges indicates
        their insignificance in direction, i.e. both way are possible. # FIXME

    Arguments
    ---------
    use_super_target : bool
        whether to add a super target at the ends of all targets.
        True if the weight of DAG exist for all edges.
        False if the weight of DAG also need to consider attribute for targets.

    See Also
    --------
    _setup_dag

    """
    # check targets of the dag
    if isinstance(target_query, str):
        targets = _find_target_nodes(G, target_query, target_query)
    logger.debug(f"{len(targets)} targets are selected")

    # remove connected components without targets
    _remain_nodes = set()
    for comp in nx.connected_components(G.to_undirected()):
        if not any(comp.intersection(targets)):
            logger.warning(f"Removing nodes disconnected from targets: {comp}.")
        else:
            _remain_nodes.update(comp)
    G = G.subgraph(_remain_nodes)

    # for each connected component with a target
    for comp in nx.connected_components(G.to_undirected()):
        _graph = G.subgraph(comp).copy()
        _targets = comp.intersection(targets)
        dag = _setup_dag(
            G=_graph,
            targets=_targets,
            weight=weight,
            loads=loads,
            report=report,
            algorithm=algorithm,
        )

    return dag


def _setup_dag(
    G: nx.Graph,
    targets=None,
    weight: str = None,
    loads: list = [],
    report: str = None,
    algorithm: str = "simple",
    **kwargs,
) -> nx.Graph:
    """Prepare subgraph as Directed Acyclic Graphs (dag) using shortest path.

    step 1: add a supernode to subgraph (representing graph - subgraph)
    step 2: use shortest path to prepare dag edges
    (source: node in subgraph; target: super node).
    In progress.

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
        If None, every node/edge has a load equal to the total number of
        nodes upstream (nnodes), and number of edges upstream (nedges).
        If a list, use the attributes in the list as load.
        The attribute can be either from edges or from nodes.
        Any attributes that are not present defaults to 0.
    algorithm : string, optional (default = 'simple')
        The algorithm to use to compute the dag.
        Supported options: 'simple', 'flowpath'.
        Other inputs produce a ValueError.
        if 'simple', the input graph is treated as a undirected graph,
        the resulting graph might alter the original edges direction
        If 'flowpath', the input graph is treated as a directed graph,
        the resulting graph do not alter the original edges direction.
        Both option will have less edges. the missing edges indicates
        their insignificance in direction, i.e. both way are possible. # FIXME

    Arguments
    ----------
    use_super_target : bool
        whether to add a super target at the ends of all targets.
        True if the weight of DAG exist for all edges.
        False if the weight of DAG also need to consider attribute for targets.

    """
    # check algorithm of the setup
    if algorithm not in ("simple", "flowpath"):
        raise ValueError(f"algorithm not supported: {algorithm}")
    logger.debug(f"Performing {algorithm}")

    # check method of the dag
    method = "dijkstra"
    if method not in ("dijkstra", "bellman-ford"):
        raise ValueError(f"method not supported: {method}")
    logger.debug(f"Performing {method} dag")

    # started making dag
    DAG = nx.DiGraph()

    # try adding super nodes
    G.add_edges_from([(n, -1) for n in targets])
    logger.debug("connecting targets to supernode")
    if len([_ for _ in nx.connected_components(G.to_undirected())]) > 1:
        raise TypeError("Cannot apply dag on disconnected graph.")

    # 1. add path
    # FIXME: if don't do this step: networkx.exception.NetworkXNoPath: No path to **.
    if algorithm == "simple":
        logger.info("dag based on method simple. undirected graph will be used. ")
        G = G.to_undirected()
        DAG = make_dag(G, targets=targets, weight=weight, drop_unreachable=False)
    elif algorithm == "flowpath":
        logger.info("dag based on method flowpath. unreachable path will be dropped. ")
        DAG = make_dag(G, targets=targets, weight=weight, drop_unreachable=True)

    # 2. add back weights
    # TODO Xiaohan: for now the weight are not reversed
    for u, v, new_d in DAG.edges(data=True):
        # get the old data from X
        old_d = G[u].get(v)
        non_shared = set(old_d) - set(new_d)
        if non_shared:
            # add old weights to new weights if not in new data
            new_d.update(dict((key, old_d[key]) for key in non_shared))
    # add back node attributes
    attrs = {n[0]: n[1] for n in G.nodes(data=True) if n[0] in DAG.nodes}
    nx.set_node_attributes(DAG, attrs)

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
        logger.debug("dag is directed acyclic graph")
    else:
        logger.error("dag is NOT directed acyclic graph")

    # visualise DAG
    if report:
        # preprocess
        _DAG = DAG.copy()
        any_has_geometry = any("geometry" in G.nodes[node] for node in G)
        if any_has_geometry:
            pos = {
                xy[0]: xy[1].coords[0]
                for xy in G.nodes(data="geometry")
                if xy[0] in _DAG.nodes
            }
        else:
            pos = {xy: xy for xy in _DAG.nodes()}

        # plot graphviz
        make_graphplot_for_targetnodes(_DAG, targets, layout="graphviz")
        plt.title(report, wrap=True, loc="left")

        # plot xy
        plt.figure(figsize=(8, 8))
        plt.title(report, wrap=True, loc="left")
        plt.axis("off")

        # base
        nx.draw(
            _DAG,
            pos=pos,
            node_size=0,
            with_labels=False,
            arrows=False,
            node_color="gray",
            edge_color="silver",
            width=0.5,
        )
        nx.draw_networkx_nodes(
            _DAG,
            pos=pos,
            nodelist=targets,
            node_size=200,
            node_shape="*",
            node_color="r",
        )
        edge_width = [d[2] / 100 for d in _DAG.edges(data="nnodes")]
        nx.draw_networkx_edges(
            _DAG,
            pos=pos,
            edgelist=_DAG.edges(),
            arrows=False,
            width=[float(i) / (max(edge_width) + 0.1) * 20 + 0.5 for i in edge_width],
        )

    return DAG


# func from hybridurb gis_utils.py


# process graph


def query_graph_edges_attributes(G, id_col: str = "id", edge_query: str = None):
    """Query the graph by selecting only the edges specified in edge_query."""
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
    """Query the graph by selecting only the nodes specified in node_query."""
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
    """Contract the nodes into one node in G."""
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


def sort_ids(G: nx.Graph):
    """Sort the ids of the graph.

    if there are no ids for the nodes, the ids will be generated based on the geometry.
    """
    if set(dict(G.nodes(data="id")).values()) == {None}:
        nx.set_node_attributes(G, {p: f"{p[0]:.6f}_{p[1]:.6f}" for p in G.nodes}, "id")

    return G


def sort_ends(G: nx.Graph):
    """Sort the ends of the graph.

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
    """Sort the start end direction of the graph and obtain start and end nodes.

    Arguments
    ---------
    G: nx.DiGraph
        Directional Graph

    Returns
    -------
    G: nx.DiGraph
    Directional Graph with node attributes endnodes and startnodes
    """
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
    """Find the predecessors of a node n.

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
    """Find the difference between G and H (G-H) based on edges.

    replace :py:meth:`~nx.difference()`.
    """
    c = G.copy()
    c.remove_edges_from(H.edges)
    c.remove_nodes_from(list(nx.isolates(c)))
    return c


def contract_graph(G: nx.Graph, partition, tonodes):
    """Contract based on partition.

    Needs further improvements.

    TODO: harmonize with setup partition.
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
    """Make dag graph.

    Needs further improvements

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
    # X_new.remove_nodes_from([-1]).
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
    from networkx.drawing.nx_agraph import graphviz_layout

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
    try:
        pos_G = {xy: xy for xy in G.nodes()}
    except IndexError:
        # support also geometry
        pos_G = {xy[0]: xy[1].coords[0] for xy in G.nodes(data="geometry")}

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
    """Make plots for grahviz layout."""
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
    """Make plots for two different layout."""
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


def get_arborescence(G: nx.DiGraph):
    """Get arborescence from Digraph.

    This function will loop through all bifurcation node and check if its predecessors
    forms a arborescence.
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
