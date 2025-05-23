"""Workflows to prepare branches for Delft3D FM model."""

import logging
from typing import List, Literal, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import shapely
from hydromt import gis_utils
from hydromt.gis_utils import nearest_merge
from scipy.spatial import distance
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point

from hydromt_delft3dfm import graph_utils, mesh_utils

from ..gis_utils import cut_pieces, split_lines

logger = logging.getLogger(__name__)


__all__ = [
    "prepare_branches",
    "process_branches",
    "validate_branches",
    "add_branches",
    "find_nearest_branch",
    "update_data_columns_attributes",
    "update_data_columns_attribute_from_query",
    "snap_newbranches_to_branches_at_snappednodes",
    "snap_geom_to_branches_and_drop_nonsnapped",
]


def prepare_branches(
    gdf_br: gpd.GeoDataFrame,
    params: pd.DataFrame,
    br_type: Literal["river", "channel", "pipe"],
    dst_crs: pyproj.CRS,
    filter: str = None,
    id_start: int = 1,
    spacing: pd.DataFrame = None,
    snap_offset: float = 0.0,
    allow_intersection_snapping: bool = False,
    allowed_columns: List[str] = [],
    logger: logging.Logger = logger,
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Set all common steps to add branches type of objects.

    Example for ie channels, rivers, pipes... branches.

    Default frictions and crossections will also be added.

    Parameters
    ----------
    gdf_br : gpd.GeoDataFrame
        gpd.GeoDataFrame of branches.
    params : pd.DataFrame
        pd.Dataframe of defaults values for gdf_br.
    br_type : str
        branches type. Either "river", "channel", "pipe".
    dst_crs : pyproj.CRS
        Destination crs for the branches and branch_nodes.
    id_start: int, optional
        Start index for branchid. By default 1.
    filter: str, optional
        Keyword in branchtype column of gdf_br used to filter lines. If None
        all lines in br_fn are used (default).
    spacing: float, optional
        Spacing value in meters to split the long pipelines lines into shorter pipes.
        By default inf - no splitting is applied.
    snap_offset: float, optional
        Snapping tolerance to automatically connecting branches. Tolerance
        must be smaller than the shortest pipe length.
        By default 0.0, no snapping is applied.
    allow_intersection_snapping: bool, optional
        Switch to choose whether snapping of multiple branch ends are
        allowed when ``snap_offset`` is used.
        By default True.
    allowed_columns: list, optional
        List of columns to filter in branches GeoDataFrame
    logger: logging.Logger, optional
        Logger.

    Returns
    -------
    branches: gpd.GeoDataFrame
        Prepared branches.
    branches_nodes: gpd.GeoDataFrame
        Nodes of the prepared branches.
    """
    # 1. Filter features based on filter
    if filter is not None and "branchtype" in gdf_br.columns:
        gdf_br = gdf_br[gdf_br["branchtype"].str.lower() == filter.lower()]
        logger.info(f"Set {filter} locations filtered from branches as {br_type} .")
    # Check if features are present
    if len(gdf_br) == 0:
        logger.warning(f"No 1D {type} locations found within domain")
        return (None, None)

    # 2. Add defaults
    # Add branchtype and branchid attributes if does not exist
    gdf_br["branchtype"] = br_type
    if "branchid" not in gdf_br.columns:
        data = [
            f"{br_type}_{i}" for i in np.arange(id_start, id_start + len(gdf_br))
        ]  # avoid duplicated ids being generated
        gdf_br["branchid"] = pd.Series(data, index=gdf_br.index, dtype=str)
    gdf_br.index = gdf_br["branchid"]
    gdf_br.index.name = "branchid"

    # filter for allowed columns
    allowed_columns = set(allowed_columns).intersection(gdf_br.columns)
    gdf_br = gpd.GeoDataFrame(gdf_br[list(allowed_columns)], crs=gdf_br.crs)

    # Add spacing to defaults
    if spacing is not None:
        params["spacing"] = spacing

    # add params
    gdf_br = update_data_columns_attributes(gdf_br, params, brtype=br_type)

    # 3. (geo-) process branches
    logger.info("Processing branches")
    if gdf_br.crs.is_geographic:  # needed for length and splitting
        gdf_br = gdf_br.to_crs(3857)
    branches, branches_nodes = process_branches(
        gdf_br,
        id_col="branchid",
        snap_offset=snap_offset,
        allow_intersection_snapping=allow_intersection_snapping,
        smooth_branches=br_type == "pipe",
        logger=logger,
    )
    logger.info("Validating branches")
    validate_branches(branches)

    # 4. convert to model crs
    branches = branches.to_crs(dst_crs)
    branches_nodes = branches_nodes.to_crs(dst_crs)

    return branches, branches_nodes


def add_branches(
    mesh1d,
    branches,
    new_branches,
    snap_newbranches_to_branches_at_snapnodes,
    snap_offset,
):
    """Add branches to exisitng open system branches at mesh1d node locations."""
    if not snap_newbranches_to_branches_at_snapnodes or mesh1d is None:
        # do not perform snap or no mesh nodes to perform snap
        branches = gpd.GeoDataFrame(
            pd.concat([branches, new_branches], ignore_index=True), crs=new_branches.crs
        )
        return snap_branch_ends(branches)

    # first snap nodes
    snappednodes = _snap_nodes(mesh1d, branches, new_branches, snap_offset)

    # snap branches
    (
        new_branches_snapped,
        branches_snapped,
    ) = snap_newbranches_to_branches_at_snappednodes(
        new_branches, branches, snappednodes
    )

    # update the branches
    branches = gpd.GeoDataFrame(
        pd.concat([branches_snapped, new_branches_snapped], ignore_index=True),
        crs=branches.crs,
    )
    return snap_branch_ends(branches)


def _snap_nodes(mesh1d, branches, newbranches, snap_offset):
    # snap nodes
    # get mesh1d nodes that allows to be snapped
    mesh1d_nodes = mesh_utils.mesh1d_nodes_geodataframe(mesh1d, branches)
    mesh1d_nodes_open = mesh1d_nodes.loc[
        mesh1d_nodes.branch_name.isin(
            branches[branches.branchtype.isin(["river", "channel"])].branchid.tolist()
        )
    ]

    # get locations that needs to be snapped
    unsnappednodes = _get_possible_unsnappednodes(newbranches)

    # snap the new endnodes to existing mesh1d_nodes_open
    snappednodes = _snap_unsnappednodes_to_nodes(
        unsnappednodes, mesh1d_nodes_open, snap_offset
    )
    return snappednodes


def _get_possible_unsnappednodes(newbranches):
    branchtype = newbranches["branchtype"].unique()[0]

    if branchtype in ["pipe", "tunnel"]:
        # add pipes - connect both ends to existing network
        endnodes = graph_utils.get_endnodes_from_lines(newbranches, where="downstream")
    else:
        # add rivers/channels - only connect downstream ends to existing network
        endnodes = graph_utils.get_endnodes_from_lines(newbranches, where="both")
    return endnodes


def _snap_unsnappednodes_to_nodes(
    unsnapped_nodes: gpd.GeoDataFrame, nodes: gpd.GeoDataFrame, snap_offset: float
) -> gpd.GeoDataFrame:
    snapped_nodes = gis_utils.nearest_merge(
        unsnapped_nodes, nodes, max_dist=snap_offset, overwrite=False
    )
    snapped_nodes = snapped_nodes[snapped_nodes.index_right != -1]  # drop not snapped
    snapped_nodes["geometry_left"] = snapped_nodes["geometry"]
    snapped_nodes["geometry_right"] = [
        nodes.at[i, "geometry"] for i in snapped_nodes["index_right"]
    ]
    return snapped_nodes


def update_data_columns_attributes(
    branches: gpd.GeoDataFrame,
    attributes: pd.DataFrame,
    brtype: str = None,
):
    """
    Add or update columns in the branches geodataframe.

    UPdate is based on column and values in attributes
    (excluding 1st col branchtype used for query).

    If brtype is set, only update the attributes of the specified type of branches.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        Branches.
    attribute : int or pd.DataFrame
        Values of the attribute. Either int of float for fixed value for
        all or a pd.DataFrame with values per "branchtype", "shape" and "width"
        (column names of the DataFrame).
    attribute_name : str
        Name of the new attribute column in branches.

    Returns
    -------
    branches : gpd.GeoDataFrame
        Branches with new attribute values.
    """
    # If brtype is specified, only update attributes for this brtype
    if brtype:
        attributes = attributes[attributes["branchtype"] == brtype]
    # Update attributes
    for i in range(len(attributes.index)):
        row = attributes.iloc[i, :]
        branch = row.loc["branchtype"]
        for colname in row.index[1:]:
            # If attribute is not at all in branches, add a new column
            if colname not in branches.columns:
                branches[colname] = pd.Series(dtype=attributes[colname].dtype)
            # Then fill in empty or NaN values with defaults
            branches.loc[
                np.logical_and(
                    branches[colname].isna(), branches["branchtype"] == branch
                ),
                colname,
            ] = row.loc[colname]

    return branches


def update_data_columns_attribute_from_query(
    branches: gpd.GeoDataFrame,
    attribute: pd.DataFrame,
    attribute_name: str,
    logger=logger,
):
    """
    Update an attribute column of branches.

    Update is based on query on "branchtype", "shape" and "width"/"diameter"
    values specified in attribute DataFrame.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        Branches.
    attribute : pd.DataFrame
        pd.DataFrame with specific attribute values per
        "branchtype", "shape" and "width"/"diameter" (column names of the DataFrame).
    attribute_name : str
        Name of the new attribute column in branches.

    Returns
    -------
    branches : gpd.GeoDataFrame
        Branches with new attribute values.
    """
    # Add a new empty attribute column to be filled in
    branches[attribute_name] = pd.Series(dtype=attribute[attribute_name].dtype)
    # Iterate over the attribute DataFrame lines
    for row in attribute.itertuples(index=False):
        # Width/diameter is not always mandatory
        if "width" in row._fields:
            if np.isnan(row.width):
                branches[attribute_name] = branches[attribute_name].where(
                    np.logical_and(
                        branches.branchtype != row.branchtype,
                        branches["shape"] != row.shape,  # shape is reserved
                    ),
                    getattr(row, attribute_name),
                )
            else:
                branches[attribute_name] = branches[attribute_name].where(
                    np.logical_and(
                        branches.branchtype != row.branchtype,
                        branches["shape"] != row.shape,  # shape is reserved
                        branches.width != row.width,
                    ),
                    getattr(row, attribute_name),
                )
        elif "diameter" in row._fields:
            if np.isnan(row.diameter):
                branches[attribute_name] = branches[attribute_name].where(
                    np.logical_and(
                        branches.branchtype != row.branchtype,
                        branches["shape"] != row.shape,  # shape is reserved
                    ),
                    getattr(row, attribute_name),
                )
            else:
                branches[attribute_name] = branches[attribute_name].where(
                    np.logical_and(
                        branches.branchtype != row.branchtype,
                        branches["shape"] != row.shape,  # shape is reserved
                        branches.diameter != row.diameter,
                    ),
                    getattr(row, attribute_name),
                )
        else:
            branches[attribute_name] = branches[attribute_name].where(
                np.logical_and(
                    branches.branchtype != row.branchtype,
                    branches["shape"] != row.shape,  # shape is reserved
                ),
                getattr(row, attribute_name),
            )

    return branches


def process_branches(
    branches: gpd.GeoDataFrame,
    id_col: str = "branchid",
    snap_offset: float = 0.01,
    allow_intersection_snapping: bool = True,
    smooth_branches: bool = False,
    logger=logger,
):
    """Process the branches.

    Process by cleaning up the branches, snapping them, splitting them
    and generating branchnodes.

    Parameters
    ----------
    branches: gpd.GeoDataFrame
        The branches to process.
    id_col: str, optional
        Defalt to branchid.
    snap_offset : float, optional
        Maximum distance in meters between end points. If the distance
        is larger, they are not snapped. Defaults to 0.01.
    allow_intersection_snapping : bool, optional
        Allow snapping at all branch ends, including intersections. Defaults to True.
    smooth_branches: bool, optional
        whether to return branches that are smoothed (straightend), needed for pipes
        Default to False.
    logger
        The logger to log messages with.

    Returns
    -------
    branches : gpd.GeoDataFrame
        Preprocessed branches.
    branches_nodes : gpd.GeoDataFrame
        Preprocessed branches' nodes.
    """
    logger.debug("Cleaning up branches")
    # TODO: maybe add arguments,use branch cross sections
    # global_controls = branches_ini.get("global", None)

    branches = cleanup_branches(
        branches,
        id_col=id_col,
        snap_offset=snap_offset,
        allow_intersection_snapping=allow_intersection_snapping,
        logger=logger,
    )

    logger.debug("Splitting branches based on spacing")
    # TODO: add check, if spacing is used,
    # then in branch cross section cannot be setup later
    branches = space_branches(branches, smooth_branches=smooth_branches, logger=logger)

    logger.debug("Generating branchnodes")
    branch_nodes = generate_branchnodes(branches, id_col, logger=logger)

    return branches, branch_nodes


def cleanup_branches(
    branches: gpd.GeoDataFrame,
    id_col: str = "branchid",
    snap_offset: float = 0.01,
    allow_intersection_snapping: bool = True,
    logger=logger,
):
    """Clean up the branches.

    Steps:
    * Removing null geomtry
    * Exploding branches with multiline strings
    * simply line geometry by removing Z coordinates
    * Removing branches with duplicated geometry
    * Removing branches that are shorter than 0.1 meters
    * Renaming branches with duplicate IDs
    * Reducing the precision of the branch geometry to 6 digits.
    * Snapping the branches.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        The branches to clean up.
    id_col : str, optional
        The branch id column name. Defaults to 'BRANCH_ID'.
    snap_offset : float, optional
        Maximum distance in meters between end points. If the distance
        is larger, they are not snapped. Defaults to 0.01.
    allow_intersection_snapping : bool, optional
        Allow snapping at all branch ends, including intersections.
        Defaults to True.
    logger
        The logger to log messages with.

    Returns
    -------
    gpd.GeoDataFrame
        The cleanup branches.
    """
    # remove null geometry
    branches = branches.loc[~branches.geometry.isna(), :]

    # explode multiline string
    _branches = branches.explode(index_parts=False)
    # 3 remove z coordinates
    _branches["geometry"] = _branches["geometry"].apply(
        lambda x: LineString([p[:2] for p in x.coords])
    )
    logger.debug("Exploding branches.")

    # remove duplicated geometry
    G = _branches["geometry"].apply(lambda geom: geom.wkb)
    n = len(G) - len(G.drop_duplicates().index)
    branches = _branches[_branches.index.isin(G.drop_duplicates().index)]
    logger.debug(f"Removing {n} branches which have duplicated geometry.")
    # remove branches that are too short
    if branches.crs.is_geographic:
        branches_length = branches.geometry.to_crs(3857).length
    else:
        branches_length = branches.geometry.length
    n = np.sum(list(branches_length <= 0.1))
    branches = branches[branches_length >= 0.1]
    logger.debug(f"Removing {n} branches that are shorter than 0.1 meter.")
    # remove branches with ring geometries
    branches = _remove_branches_with_ring_geometries(branches)

    # sort index
    if id_col in [
        "None",
        "NONE",
        "none",
        None,
        "",
    ]:  # TODO: id_column must be specified
        id_col = "BRANCH_ID"
        # regenerate ID based on ini
        # # NOTE BMA: this is the step to ensure unique id cross the network
        # TODO could not find anay example of this or default
        id_prefix = 100  # branches_ini["global"]["id_prefix"]
        id_suffix = 100  # branches_ini["global"]["id_suffix"]
        branches[id_col] = [
            f"{id_prefix}_{x}_{id_suffix}" for x in range(len(branches))
        ]
        logger.warning(
            "id_col is not specified. Branch id columns are"
            f"read/generated using default: {id_col}."
        )

    # check duplicated id
    _branches = branches.copy()
    _branches.reset_index(drop=True, inplace=True)
    _branches["count"] = (
        _branches.groupby(id_col)[id_col].transform("cumcount").astype(int)
    )
    for bi, b in _branches.iterrows():
        if b["count"] >= 1:
            _branches.rename(
                index={bi: b[id_col] + "-" + str(b["count"])}, inplace=True
            )
        else:
            _branches.rename(index={bi: b[id_col]}, inplace=True)
    _branches[id_col] = _branches.index
    _branches.index.name = id_col
    branches = _branches.copy()
    n = sum(_branches["count"])
    logger.debug(
        f"Renaming {n} id_col duplicates. Convention:"
        "BRANCH_1, BRANCH_1 --> BRANCH_1, BRANCH_1-2."
    )

    # precision correction --> check if needs to be changed by the user,
    # if so move that into a func argument
    branches = reduce_gdf_precision(
        branches, rounding_precision=6  # branches_ini["global"]["rounding_precision"]
    )  # recommned to be larger than e-8
    logger.debug("Reducing precision of the GeoDataFrame." "Rounding precision (e-6) .")

    # snap branches
    if allow_intersection_snapping is True:
        # snap points no matter it is at intersection or ends
        branches = snap_branch_ends(branches, offset=snap_offset)
        logger.debug(
            "Performing snapping at all branch ends, including intersections"
            "(To avoid messy results, please use a lower snap_offset)."
        )

    else:
        # snap points at ends only
        branches = snap_branch_ends(branches, offset=snap_offset, max_points=2)
        logger.debug(
            "Performing snapping at all branch ends, excluding intersections"
            "(To avoid messy results, please use a lower snap_offset).."
        )

    # Drop count column
    if "count" in branches.columns:
        branches = branches.drop(columns=["count"])

    return branches


def space_branches(
    branches: gpd.GeoDataFrame,
    spacing_col: str = "spacing",
    smooth_branches: bool = False,
    logger=logger,
):
    """
    Space the branches based on the spacing_col on the branch.

    Removes the spacing column from the branches afterwards.

    # TODO: seperate situation where interpolation
    is needed and interpolation is not needed

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        The branches to clean up.
    spacing_col : str, optional
        The branch id column name. Defaults to 'spacing'.
    logger
        The logger to log messages with.

    Returns
    -------
    gpd.GeoDataFrame
        The split branches.
    """
    # split branches based on spacing
    branches_ = split_branches(
        branches, spacing_col=spacing_col, smooth_branches=smooth_branches
    )
    logger.debug(f"clipping branches into {len(branches_)} segments")

    # remove spacing column
    branches_ = branches_.drop(columns=[spacing_col])

    return branches_


def generate_branchnodes(
    branches: gpd.GeoDataFrame,
    id_col: str = None,
    logger=logger,
):
    """Generate branch nodes at the branch ends.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        The branches to generate the end nodes for.
    id_col : str, optional
        The branch id column name. Defaults to None.
    logger
        The logger to log messages with.

    Returns
    -------
    gpd.GeoDataFrame
        The branch nodes.
    """
    # generate node up and downstream
    nodes = pd.DataFrame(
        [Point(l.coords[0]) for li, l in branches["geometry"].items()]
        + [Point(l.coords[-1]) for li, l in branches["geometry"].items()],
        columns=["geometry"],
    )

    if id_col is None:
        id_col = branches.index.name

    nodes = []
    for bi, b in branches.iterrows():
        nodes.append([Point(b.geometry.coords[0]), b[id_col]])  # start
        nodes.append([Point(b.geometry.coords[-1]), b[id_col]])  # end
    nodes = pd.DataFrame(nodes, columns=["geometry", "branch_id"])
    nodes = pd.merge(
        nodes,
        branches.reset_index(drop=True),
        left_on="branch_id",
        right_on=branches.index.name,
        suffixes=("", "_b"),
    )
    nodes = nodes.drop(columns="geometry_b")
    # remove duplicated geometry
    _nodes = nodes.copy()
    G = _nodes["geometry"].apply(lambda geom: geom.wkb)
    nodes = _nodes[_nodes.index.isin(G.drop_duplicates().index)]
    nodes = gpd.GeoDataFrame(nodes)
    nodes.crs = branches.crs
    # give index
    nodes.index = [f"branchnodes_{i}" for i in range(len(nodes))]
    nodes.index.name = "branchnodeid"
    # nodes[nodes.index.name] = nodes.index # creates a duplicate column
    return nodes


def validate_branches(
    branches: gpd.GeoDataFrame, logger=logger
):  # TODO: add more content and maybe make a seperate module
    """Validate the branches.

    Logs an error when one or more branches have a length of 0 meter.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        The branches to validate.
    logger
        The logger to log messages with.
    """
    # validate pipe geometry
    if sum(branches.geometry.length <= 0) == 0:
        logger.debug("Branches are valid.")
    else:
        logger.error(
            f"Branches {branches.index[branches.geometry.length <= 0]}"
            + "have length of 0 meter. "
            + "Issue might have been caused by using a snap_offset"
            + "that is too large. Please revise or modify the branches data layer. "
        )


def split_branches(
    branches: gpd.GeoDataFrame,
    spacing_const: float = float("inf"),
    spacing_col: str = None,
    smooth_branches: bool = False,
    logger=logger,
):
    """
    Split branches based on a given spacing.

    If spacing_col is used (default), apply spacing as a categorical
    variable-distance used to split branches.
    If spacing_const is used (overwrite), apply spacing as a
    constant -  distance used to split branches.
    Raise Error if neither exist.

    If ``smooth_branches``, split branches generated will be straight line.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
    spacing_const : float
        Constent spacing which will overwrite the spacing_col.
        Defaults to float("inf").
    spacing_col: str
        Name of the column in branchs that contains spacing information.
        Default to None.
    smooth_branches: bool, optional
        Switch to split branches into straight lines. By default False.
    logger
        The logger to log messages with.

    Returns
    -------
    split_branches : gpd.GeoDataFrame
        Branches after split, new ids will be overwritten for
        the branch index. Old ids are stored in "OLD_" + index.
    """
    id_col = branches.index.name
    if spacing_col is None:
        logger.info(f"Splitting branches with spacing of {spacing_const} [m]")
        split_branches = _split_branches_by_spacing_const(
            branches,
            spacing_const,
            id_col=id_col,
            smooth_branches=smooth_branches,
        )

    elif branches[spacing_col].astype(float).notna().any():
        logger.info(
            "Splitting branches with spacing specifed"
            f"in datamodel branches[{spacing_col}]"
        )
        split_branches = []
        for spacing_subset, branches_subset in branches.groupby(spacing_col):
            if spacing_subset:
                split_branches_subset = _split_branches_by_spacing_const(
                    branches_subset,
                    spacing_subset,
                    id_col=id_col,
                    smooth_branches=smooth_branches,
                )
            else:
                branches_subset.loc[:, f"ORIG_{id_col}"] = branches_subset[id_col]
                split_branches_subset = branches_subset
            split_branches.append(split_branches_subset)
        split_branches = pd.concat(split_branches)

    else:  # no spacing information specified anywhere, do not apply splitting
        branches.loc[:, f"ORIG_{id_col}"] = branches[id_col]
        split_branches = branches

    # reassign branch id for the generated
    split_branches.index = split_branches[id_col]
    return split_branches


def _split_branches_by_spacing_const(
    branches: gpd.GeoDataFrame,
    spacing_const: float,
    id_col: str = "BRANCH_ID",
    smooth_branches: bool = False,
):
    """
    Split branches based on a given spacing constant.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
    spacing_const : float
        Constant spacing which will overwrite the spacing_col.
    id_col: str
        Name of the column in branches that contains the id of the branches.
    smooth_branches: bool, optional
        Swith to split branches into straight lines. By default False.

    Returns
    -------
    split_branches : gpd.GeoDataFrame
        Branches after split, new ids will be stored in id_col.
        Original ids are stored in "ORIG_" + id_col.
    """
    if spacing_const == float("inf"):
        branches[f"ORIG_{id_col}"] = branches[id_col]
        branches.index = branches[id_col]
        return branches

    edge_geom = []
    edge_offset = []
    edge_invertup = []
    edge_invertdn = []
    edge_bedlevup = []
    edge_bedlevdn = []
    edge_index = []
    branch_index = []

    # Check for attributes
    interp_invlev = "invlev_up" and "invlev_dn" in branches.columns
    interp_bedlev = "bedlev_up" and "bedlev_dn" in branches.columns

    for bid, b in branches.iterrows():
        # prepare for splitting
        line = b.geometry
        num_new_lines = int(np.ceil(line.length / spacing_const))

        if num_new_lines <= 0:
            continue

        # interpolate geometry
        new_edges = split_lines(line, num_new_lines)
        if smooth_branches:
            for i in range(len(new_edges)):
                ed = new_edges[i]
                new_edges[i] = LineString([Point(ed.coords[0]), Point(ed.coords[-1])])
        offsets = np.linspace(0, line.length, num_new_lines + 1)

        # interpolate values
        edge_geom.extend(new_edges)
        edge_offset.extend(offsets[1:])
        if interp_invlev:
            edge_invertup.extend(
                np.interp(
                    offsets[:-1], [0, offsets[-1]], [b.invlev_up, b.invlev_dn]
                )  # TODO: renaming needed
            )
            edge_invertdn.extend(
                np.interp(offsets[1:], [0, offsets[-1]], [b.invlev_up, b.invlev_dn])
            )
        if interp_bedlev:
            edge_bedlevup.extend(
                np.interp(offsets[:-1], [0, offsets[-1]], [b.bedlev_up, b.bedlev_dn])
            )
            edge_bedlevdn.extend(
                np.interp(offsets[1:], [0, offsets[-1]], [b.bedlev_up, b.bedlev_dn])
            )
        edge_index.extend([bid + "_E" + str(i) for i in range(len(new_edges))])
        branch_index.extend([bid] * len(new_edges))

    edges = gpd.GeoDataFrame(
        {
            "EDGE_ID": edge_index,
            "geometry": edge_geom,
            id_col: branch_index,
            # "invlev_up": edge_invertup,
            # "invlev_dn": edge_invertdn,
            # "bedlev_up": edge_bedlevup,
            # "bedlev_dn": edge_bedlevdn,
        },
        crs=branches.crs,
    )
    if interp_invlev:
        edges["invlev_up"] = edge_invertup
        edges["invlev_dn"] = edge_invertdn
    if interp_bedlev:
        edges["bedlev_up"] = edge_bedlevup
        edges["bedlev_dn"] = edge_bedlevdn
    edges_attr = pd.concat(
        [branches.loc[idx, :] for idx in branch_index], axis=1
    ).transpose()
    edges = pd.concat(
        [
            edges,
            edges_attr.drop(
                columns=list(set(edges.columns) - set(["EDGE_ID"]))
            ).reset_index(),
        ],
        axis=1,
    )

    edges = edges.rename(columns={id_col: f"ORIG_{id_col}"})
    edges = edges.rename(columns={"EDGE_ID": id_col})
    edges.index = edges[id_col]
    split_branches = edges

    return split_branches


def reduce_gdf_precision(gdf: gpd.GeoDataFrame, rounding_precision: int = 8):
    """Reduce the geometry coordinate precision with the provided number of digits.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        The geo data frame to reduce the precision for.
    rounding_precision : int, optional
        The number of digits to round the coordinates with. Defaults to 8.

    Returns
    -------
    gpd.GeoDataFrame
        The geo data frame with the rounded geometry.

    Raises
    ------
    NotImplementedError
        If the geometry is not a LineString or Point.
    """
    if isinstance(gdf.geometry[0], LineString):
        branches = gdf.copy()
        for i_branch, branch in enumerate(branches.itertuples()):
            points = shapely.wkt.loads(
                shapely.wkt.dumps(
                    branch.geometry, rounding_precision=rounding_precision
                )
            ).coords[:]
            branches.at[i_branch, "geometry"] = LineString(points)

    elif isinstance(gdf.geometry[0], Point):
        points = gdf.copy()
        for i_point, point in enumerate(points.itertuples()):
            new_point = shapely.wkt.loads(
                shapely.wkt.dumps(point.geometry, rounding_precision=rounding_precision)
            ).coords[:]
            points.at[i_point, "geometry"] = Point(new_point)

    else:
        raise NotImplementedError

    return gdf


def snap_branch_ends(
    branches: gpd.GeoDataFrame,
    offset: float = 0.01,
    subsets=[],
    max_points: int = np.inf,
    id_col: str = "BRANCH_ID",
):
    """
    Snap branch ends to other branch ends within a given offset.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
    offset : float, optional
        Maximum distance in meters between end points. If the distance is larger,
        they are not snapped. Default to 0.01.
    subsets : list, optional
        A list of branch id subset to perform snapping (forced snapping).
        Default to an empty list.
    max_points : int, optional
        maximum points allowed in a group.
        if snapping branch ends only, use max_points = 2
        if not specified, branch intersections will also be snapped
        Defaults to np.inf.
    id_col : str, optional
        The branch id column name. Defaults to 'BRANCH_ID'.

    Returns
    -------
    branches : gpd.GeoDataFrame
        Branches updated with snapped geometry.
    """
    # Collect endpoints
    _endpoints = []
    for branch in branches.itertuples():
        _endpoints.append((branch.geometry.coords[0], branch.Index, 0))
        _endpoints.append((branch.geometry.coords[-1], branch.Index, -1))

    # determine which branches should be included
    if len(subsets) > 0:
        _endpoints = [[i for i in _endpoints if i[1] in subsets]]

    # # group branch ends based on off set
    groups = {}
    coords = [i[0] for i in _endpoints]
    dist = distance.squareform(distance.pdist(coords))
    bdist = dist <= offset
    for row_i, row in enumerate(bdist):
        groups[_endpoints[row_i]] = []
        for col_i, col in enumerate(row):
            if col:
                groups[_endpoints[row_i]].append(_endpoints[col_i])

    # remove duplicated group, group that does not satisfy max_points in groups.
    # Assign endpoints
    endpoints = {
        k: list(set(v))
        for k, v in groups.items()
        if (len(set(v)) >= 2) and (len(set(v)) <= max_points)
    }
    logger.debug(
        "Limit snapping to allow a max number of {max_points}"
        "contact points. If max number == 2, it means 1 to 1 snapping."
    )

    # Create a counter
    snapped = 0

    # snap each group (list) in endpoints together,
    # by using the coords from the first point
    for point_reference, points_to_snap in endpoints.items():
        # get the point_reference coords as reference point
        ref_crd = point_reference[0]
        # for each of the rest
        for j, (endpoint, branchid, side) in enumerate(points_to_snap):
            # Change coordinates of branch
            crds = branches.at[branchid, "geometry"].coords[:]
            if crds[side] != ref_crd:
                crds[side] = ref_crd
                branches.at[branchid, "geometry"] = LineString(crds)
                snapped += 1
    logger.debug(f"Snapped {snapped} points.")

    return branches


# TODO copied from dhydamo geometry.py, update when available in main
def possibly_intersecting(
    dataframebounds: np.ndarray, geometry: gpd.GeoDataFrame, buffer: int = 0
):
    """
    Find intersecting profiles.

    Finding intersecting profiles for each branch is a slow
    process in case of large datasets. To speed this up, we first determine
    which profile intersect a square box around the branch

    With the selection, the intersecting profiles can be determines much faster.

    Parameters
    ----------
    dataframebounds : numpy.array
    geometry : shapely.geometry.Polygon
    """
    geobounds = geometry.bounds
    idx = (
        (dataframebounds[0] - buffer < geobounds[2])
        & (dataframebounds[2] + buffer > geobounds[0])
        & (dataframebounds[1] - buffer < geobounds[3])
        & (dataframebounds[3] + buffer > geobounds[1])
    )
    # Get intersecting profiles
    return idx


def find_nearest_branch(
    branches: gpd.GeoDataFrame,
    geometries: gpd.GeoDataFrame,
    method: str = "overal",
    maxdist: float = 5.0,
) -> gpd.GeoDataFrame:
    """Determine the nearest branch for each geometry.

    The method of determination can vary.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        Geodataframe containing branch geometries.
    geometries : gpd.GeoDataFrame
        Geodataframe containing geometries for which the nearest branch needs to
        be found.
    method : str, optional
        Method to determine the nearest branch. Supports:
        - 'overal': Find the nearest branch based on the geometry's location.
        - 'intersecting': Convert the geometry to a centroid of its intersection
        points with branches.
        Default is 'overal'.
    maxdist : float, optional
        Maximum distance threshold for finding the nearest branch. Default is 5.0.

    Returns
    -------
    gpd.GeoDataFrame
        Geodataframe with additional columns:
        - 'branch_id': ID of the nearest branch.
        - 'branch_distance': Distance to the nearest branch.
        - 'branch_offset': Offset along the branch.

    Raises
    ------
    NotImplementedError
        If the specified method is not among the allowed methods.
    """
    # Check if method is in allowed methods
    allowed_methods = ["intersecting", "overal"]
    if method not in allowed_methods:
        raise NotImplementedError(f'Method "{method}" not implemented.')

    # Depending on method, modify geometries
    if method == "intersecting":
        # Get the intersection points directly
        for index, geom in geometries.iterrows():
            # Find branches that intersect with the current geometry
            intersected_branches = branches[branches.intersects(geom["geometry"])]

            if not intersected_branches.empty:
                # If there are multiple intersecting points, take the centroid
                intersection_points = [
                    geom["geometry"].intersection(branch["geometry"])
                    for _, branch in intersected_branches.iterrows()
                ]
                centroid = MultiPoint(intersection_points).centroid
                geometries.at[index, "geometry"] = centroid

    # Check for previous data and drop if exist
    geometries.drop(
        columns=["branch_id", "branch_distance", "branch_offset"],
        inplace=True,
        errors="ignore",
    )

    # Use nearest_merge to get the nearest branches
    result = nearest_merge(geometries, branches, max_dist=maxdist, columns=["geometry"])
    result.rename(
        columns={"index_right": "branch_id", "distance_right": "branch_distance"},
        inplace=True,
    )

    # Select ones that are merged
    valid_rows = result["branch_distance"] < maxdist

    # Interpolate the branch geometries based on index_right
    branchgeo = branches.loc[result.loc[valid_rows, "branch_id"], "geometry"].values
    maxdist = np.array([geo.length for geo in branchgeo])
    offset = np.array(
        [
            geo.project(result.loc[idx, "geometry"])
            for idx, geo in zip(valid_rows[valid_rows].index, branchgeo)
        ]
    )
    result.loc[valid_rows, "branch_offset"] = np.where(
        offset > maxdist, maxdist, offset
    )
    snapped_geometries = [geo.interpolate(o) for geo, o in zip(branchgeo, offset)]
    result.loc[valid_rows, "geometry"] = snapped_geometries

    # For rows where distance is greater than maxdist, set branch_id
    # to empty and branch_offset to NaN
    result.loc[~valid_rows, "branch_id"] = ""
    result.loc[~valid_rows, "branch_offset"] = np.nan
    result.loc[~valid_rows, "branch_distance"] = np.nan

    return result


def snap_newbranches_to_branches_at_snappednodes(
    new_branches: gpd.GeoDataFrame,
    branches: gpd.GeoDataFrame,
    snappednodes: gpd.GeoDataFrame,
):
    """Snap new_branches to branches at snappednodes.

    snapnodes are located at branches. new branches will be snapped,
    and branches will be splitted.

    # NOTE: no interpolation of crosssection is needed because inter branch
    interpolation is turned on using branchorder.

    Parameters
    ----------
    new_branches : geopandas.GeoDataFrame
        Geodataframe of new branches whose geometry will be modified:
        end nodes will be snapped to snapnodes
    branches : geopandas.GeoDataFrame
        Geodataframe who will be splitted at snapnodes to allow connection
        with the new_branches.
    snapnodes : geopandas.GeoDataFrame
        Geodataframe which contiains the spatial relation of the new_branches
        and branches.

    Returns
    -------
    new_branches_snapped : geopandas.GeoDataFrame
        Geodataframe of new branches with endnodes be snapped to snapnodes
        in branches_snapped.
    branches_snapped : geopandas.GeoDataFrame
        Geodataframe of branches splitted at snapnodes to allow connection
        with the new_branches_snapped.
    """
    new_branches.index = new_branches.branchid
    branches.index = branches.branchid

    # for each snapped endnodes
    new_branches_snapped = new_branches.copy()
    branches_snapped = branches.copy()

    # modify new branches
    for snapnode in snappednodes.itertuples():
        new_branch = new_branches.loc[snapnode.branchid]
        snapped_line = LineString(
            [
                (
                    snapnode.geometry_right
                    if Point(xy).equals(snapnode.geometry_left)
                    else Point(xy)
                )
                for xy in new_branch.geometry.coords[:]
            ]
        )
        new_branches_snapped.at[snapnode.branchid, "geometry"] = snapped_line

    # modify old branches
    _branch_order_base = max(max(branches_snapped.branchorder), 1)
    for branch_order, branch_name in enumerate(sorted(set(snappednodes.branch_name))):
        branch = branches.loc[branch_name]
        branch_order = _branch_order_base + branch_order
        distances = snappednodes[
            snappednodes.branch_name == branch_name
        ].branch_chainage.to_list()
        snapped_line = MultiLineString(cut_pieces(branch.geometry, distances))
        branches_snapped.at[branch_name, "geometry"] = snapped_line
        # allow interpolation on the snapped branch
        branches_snapped.at[branch_name, "branchorder"] = branch_order

    # explode multilinestring after snapping
    branches_snapped = branches_snapped.explode(index_parts=False)

    # reset the idex
    branches_snapped = cleanup_branches(branches_snapped)

    # precision correction
    branches_snapped = reduce_gdf_precision(branches_snapped, rounding_precision=6)

    return new_branches_snapped, branches_snapped


def _remove_branches_with_ring_geometries(
    branches: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    first_nodes = [l.coords[0] for l in branches.geometry]
    last_nodes = [l.coords[-1] for l in branches.geometry]
    duplicate_ids = np.isclose(first_nodes, last_nodes)
    duplicate_ids = [
        branches.index[i] for i in range(len(branches)) if np.all(duplicate_ids[i])
    ]
    branches = branches.drop(duplicate_ids, axis=0)
    logger.debug("Removing branches with ring geometries.")

    return branches


def snap_geom_to_branches_and_drop_nonsnapped(
    branches: gpd.GeoDataFrame, geoms: gpd.GeoDataFrame, snap_offset=0.0
):
    """
    Snap geoms to branches and drop the ones that are not snapped.

    Returns snapped geoms with branchid and chainage.
    Branches must have branchid.
    """
    geoms = find_nearest_branch(
        branches=branches,
        geometries=geoms,
        maxdist=snap_offset,
    )
    geoms = geoms.rename(columns={"branch_id": "branchid", "branch_offset": "chainage"})

    # drop ones non snapped
    _drop_geoms = geoms["chainage"].isna()
    if any(_drop_geoms):
        logger.debug(f"Unable to snap to branches: {geoms[_drop_geoms].index}")

    return geoms[~_drop_geoms]
