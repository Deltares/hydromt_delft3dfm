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

from hydromt_delft3dfm import graph_utils
from hydromt_delft3dfm.workflows import find_nearest_branch
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


# TODO imcomplete
def create_graph_from_geodataframe(gdf):
    graph = graph_utils.gpd_to_digraph(gdf)
    graph = graph_utils.preprocess_graph(graph)
    return graph


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
        polygon=region.buffer(1000)
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
    if raster_fn in [str, Path]:
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
) -> nx.graph:
    """Add data variable(s) from ``vector_fn`` to attribute(s) in graph object using nearest_join.

    Raster data is sampled to the graph edges and nodes using the ``resampling_method``.
    If raster is a dataset, all variables will be added unless ``variables`` list
    is specified.

    Parameters
    ----------
    vector_fn: str, Path, gpd.GeoDataFrame
        Data catalog key, path to vector file or gpd.GeoDataFrame data object.
    variables: list, optional
        List of variables to add to graph from vector_fn. By default all.
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
    if vector_fn in [str, Path]:
        gdf = data_catalog.get_geodataframe(
            vector_fn, region=region, buffer=1000, variables=variables
        )
    else:
        gdf = vector_fn

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
