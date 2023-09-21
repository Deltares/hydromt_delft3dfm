# -*- coding: utf-8 -*-

import logging
from typing import Tuple, Union, Literal, List

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point

from .branches import find_nearest_branch, update_data_columns_attributes

# from delft3dfmpy.core import geometry
from ..gis_utils import check_gpd_attributes

logger = logging.getLogger(__name__)


__all__ = [
    "prepare_default_friction_and_crosssection",
    "init_crosssections_options",
    "set_branch_crosssections",
    "set_xyz_crosssections",
    "set_point_crosssections",
    "add_crosssections",
]  # , "process_crosssections", "validate_crosssections"]


def prepare_default_friction_and_crosssection(
    branches: gpd.GeoDataFrame,
    br_type: Literal["river", "channel", "pipe"],
    friction_type: str = "Manning",
    friction_value: float = 0.023,
    crosssections_shape: Literal["rectangle", "circle"] = None,
    crosssections_value: Union[List[float], float] = None,
    logger: logging.Logger = logger,
):
    """
    Prepare the default uniform friction type-value pairs and crosssection profiles for branches.

    Parameters
    ----------
    branches: gpd.GeoDataFrame
        Branches to add frictions and crosssections.
    br_type : str
        branches type. Either "river", "channel", "pipe".
    friction_type : str
        Type of friction to use. One of ["Manning", "Chezy", "wallLawNikuradse", "WhiteColebrook", "StricklerNikuradse", "Strickler", "deBosBijkerk"].
    friction_value : float
        Value corresponding to ''friction_type''. Units are ["Chézy C [m 1/2 /s]", "Manning n [s/m 1/3 ]", "Nikuradse k_n [m]", "Nikuradse k_n [m]", "Nikuradse k_n [m]", "Strickler k_s [m 1/3 /s]", "De Bos-Bijkerk γ [1/s]"]
    crosssections_shape : str, optional
        Shape of branch crosssections to overwrite defaults. Either "circle" or "rectangle".
    crosssections_value : float or list of float, optional
        Crosssections parameter value to overwrite defaults.
        If ``crosssections_shape`` = "circle", expects a diameter [m], used for br_type == "pipe"
        If ``crosssections_shape`` = "rectangle", expects a list with [width, height] (e.g. [1.0, 1.0]) [m]. used for br_type == "river" or "channel".
    logger: Logger, optional
        Logger.
    Return
    ------
    branches: gpd.GeoDataFrame
        Branches with frictions and crosssections added.

    """
    # intialise defaults with branch type
    defaults = pd.DataFrame({"branchtype": br_type}, index=[0])

    # Add friction to defaults
    defaults["frictiontype"] = friction_type
    defaults["frictionvalue"] = friction_value
    # Add crosssections to defaults
    if crosssections_shape == "circle":
        if isinstance(crosssections_value, float):
            defaults["shape"] = crosssections_shape
            defaults["diameter"] = crosssections_value
        else:
            logger.warning(
                "If crosssections_shape is circle, crosssections_value should be a single float for diameter. Keeping defaults"
            )
    elif crosssections_shape == "rectangle":
        if isinstance(crosssections_value, list) and len(crosssections_value) == 2:
            defaults["shape"] = crosssections_shape
            defaults["width"], defaults["height"] = crosssections_value
            defaults["closed"] = "no"
        else:
            logger.warning(
                "If crosssections_shape is rectangle, crosssections_value should be a list with [width, height] values. Keeping defaults"
            )

    logger.info("Adding/Filling branches attributes values")
    branches = update_data_columns_attributes(branches, defaults, brtype=br_type)

    # compose
    branches["frictionid"] = [
        f"{ftype}_{fvalue}"
        for ftype, fvalue in zip(branches["frictiontype"], branches["frictionvalue"])
    ]

    return branches


def init_crosssections_options(
    crosssections_type: Union[str, list], crosssections_fn: Union[str, list]
) -> Tuple[list, list]:
    """Initialise crosssection options from user input."""
    if crosssections_type is None:
        crosssections_type = []
        crosssections_fn = []
    elif isinstance(crosssections_type, str):
        crosssections_type = [crosssections_type]
        crosssections_fn = [crosssections_fn]
    elif isinstance(crosssections_type, list):
        assert len(crosssections_type) == len(crosssections_fn)
    # inser branch as default
    if "branch" not in crosssections_type:
        crosssections_type.insert(0, "branch")
        crosssections_fn.insert(0, None)
    return crosssections_type, crosssections_fn


def add_crosssections(
    crosssections: Union[gpd.GeoDataFrame, None], new_crosssections: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Adds new_crossections to crosssections.

    Besides adding, also removes duplicated crossections and overwrites generated with user specified.

    """
    # TODO: setup river crosssections, set contrains based on branch types
    if crosssections is not None:
        # add new crossections
        new_crosssections = gpd.GeoDataFrame(
            pd.concat([crosssections, new_crosssections]), crs=crosssections.crs
        )
        # seperate crossection locations (used by branches) and definitions (used by branches and structures)
        _crosssections_locations = new_crosssections[
            ~new_crosssections.crsloc_id.isna()
        ]
        _crosssections_definitions = new_crosssections[
            new_crosssections.crsloc_id.isna()
        ]
        # remove duplicated locations based on branchid_chainage (not on xy)
        _crosssections_locations["temp_id"] = _crosssections_locations.apply(
            lambda x: f'{x["crsloc_branchid"]}_{x["crsloc_chainage"]:.2f}', axis=1
        )
        if _crosssections_locations["temp_id"].duplicated().any():
            logger.warning(
                "Duplicate crosssections locations found, removing duplicates"
            )
            # Remove duplicates based on the branch_id, branch_offset column, keeping the first occurrence (with minimum branch_distance)
            _crosssections_locations = _crosssections_locations.drop_duplicates(
                subset=["temp_id"], keep="first"
            )
        # combine them
        new_crosssections = pd.concat(
            [_crosssections_locations, _crosssections_definitions]
        )
        # remove existing generated branch crosssections
        crossections_generated = new_crosssections[
            new_crosssections.crsdef_id.str.endswith("_branch")
        ]
        crossections_others = new_crosssections[
            ~new_crosssections.crsdef_id.str.endswith("_branch")
        ]
        mask_to_remove = crossections_generated["crsloc_branchid"].isin(
            crossections_others["crsloc_branchid"]
        )
        if mask_to_remove.sum() > 0:
            new_crosssections = gpd.GeoDataFrame(
                pd.concat(
                    [crossections_generated[~mask_to_remove], crossections_others]
                ),
                crs=crosssections.crs,
            )

    return new_crosssections


def set_branch_crosssections(
    branches: gpd.GeoDataFrame,
    midpoint: bool = True,
):
    """
    Function to set regular cross-sections for each branch.
    Only support rectangle, trapezoid (for rivers) and circle (for pipes).
    Crosssections are derived at branches mid points if ``midpoints`` is True. Use this to setup crossections for rivers.
    Crosssections are derived at up and downstream of branches if ``midpoints`` is False. Use this to setup crossections for pipes.
    else at both upstream and downstream extremities of branches if False.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        The branches.
        * Required variables: branchid, shape, shift, frictionid, frictiontype, frictionvalue
        * Optional variables:
            if shape = 'circle': 'diameter'
            if shape = 'rectangle': 'width', 'height', 'closed'
            if shape = 'trapezoid': 'width', 't_width', 'height', 'closed'
    midpoint : bool, Optional
        Whether the crossection should be derived at midpoint of the branches.
        if midpoint = True:
            support 'rectangle' and 'trapezoid' shapes.
        if midpoint = True:
            support 'circle' shape.
        True by default.

    Returns
    -------
    gpd.GeoDataFrame
        The cross sections.

    # FIXME: upstream and downstream crosssection type needs further improvement for channels
    """
    # Get the crs at the midpoint of branches if midpoint
    if midpoint:
        crosssections = branches.copy()
        crosssections["branch_id"] = crosssections["branchid"]
        crosssections["branch_offset"] = [
            l / 2 for l in crosssections["geometry"].length
        ]
        crosssections["shift"] = crosssections["bedlev"]
        crosssections["geometry"] = [
            l.interpolate(0.5, normalized=True) for l in crosssections.geometry
        ]

        crosssections_ = pd.DataFrame()

        # loop through the shapes
        all_shapes = crosssections["shape"].unique().tolist()
        for shape in all_shapes:
            if shape == "rectangle":
                rectangle_crs = crosssections.loc[crosssections["shape"] == shape, :]
                rectangle_crs["definitionid"] = rectangle_crs.apply(
                    lambda x: "rect_h{:,.3f}_w{:,.3f}_c{:s}_{:s}".format(
                        x["height"],
                        x["width"],
                        x["closed"],
                        "branch",
                    ),
                    axis=1,
                )
                check_gpd_attributes(
                    rectangle_crs,
                    required_columns=[
                        "branch_id",
                        "branch_offset",
                        "frictionid",
                        "frictiontype",
                        "frictionvalue",
                        "width",
                        "height",
                        "closed",
                    ],
                )
                crosssections_ = pd.concat(
                    [crosssections_, _set_rectangle_crs(rectangle_crs)]
                )
            if shape == "trapezoid":
                trapezoid_crs = crosssections.loc[crosssections["shape"] == shape, :]
                trapezoid_crs["definitionid"] = trapezoid_crs.apply(
                    lambda x: "trapz_h{:,.1f}_bw{:,.1f}_tw{:,.1f}_c{:s}_{:s}".format(
                        x["height"],
                        x["width"],
                        x["t_width"],
                        x["closed"],
                        "branch",
                    ),
                    axis=1,
                )
                check_gpd_attributes(
                    trapezoid_crs,
                    required_columns=[
                        "branch_id",
                        "branch_offset",
                        "frictionid",
                        "frictiontype",
                        "frictionvalue",
                        "width",
                        "height",
                        "t_width",
                        "closed",
                    ],
                )
                crosssections_ = pd.concat(
                    [crosssections_, _set_trapezoid_crs(trapezoid_crs)]
                )
    # Else prepares crosssections at both upstream and dowsntream extremities
    # for now only support circle profile for pipes
    else:
        # Upstream
        ids = [f"{i}_up" for i in branches.index]
        crosssections_up = gpd.GeoDataFrame(
            {"geometry": [Point(l.coords[0]) for l in branches.geometry]},
            index=ids,
            crs=branches.crs,
        )

        crosssections_up["crsloc_id"] = [
            f"crs_up_{bid}" for bid in branches["branchid"]
        ]
        crosssections_up["branch_id"] = branches["branchid"].values
        crosssections_up["branch_offset"] = [0.0 for l in branches.geometry]
        crosssections_up["shift"] = branches["invlev_up"].values
        crosssections_up["shape"] = branches["shape"].values
        crosssections_up["diameter"] = branches["diameter"].values
        crosssections_up["frictionid"] = branches["frictionid"].values
        crosssections_up["frictiontype"] = branches["frictiontype"].values
        crosssections_up["frictionvalue"] = branches["frictionvalue"].values
        # Downstream
        ids = [f"{i}_dn" for i in branches.index]
        crosssections_dn = gpd.GeoDataFrame(
            {"geometry": [Point(l.coords[-1]) for l in branches.geometry]},
            index=ids,
            crs=branches.crs,
        )
        crosssections_dn["crsloc_id"] = [
            f"crs_dn_{bid}" for bid in branches["branchid"]
        ]
        crosssections_dn["branch_id"] = branches["branchid"].values
        crosssections_dn["branch_offset"] = [l for l in branches["geometry"].length]
        crosssections_dn["shift"] = branches["invlev_dn"].values
        crosssections_dn["shape"] = branches["shape"].values
        crosssections_dn["diameter"] = branches["diameter"].values
        crosssections_dn["frictionid"] = branches["frictionid"].values
        crosssections_dn["frictiontype"] = branches["frictiontype"].values
        crosssections_dn["frictionvalue"] = branches["frictionvalue"].values
        # Merge
        crosssections = pd.concat([crosssections_up, crosssections_dn])

        crosssections_ = pd.DataFrame()

        # loop through the shapes
        all_shapes = crosssections["shape"].unique().tolist()
        for shape in all_shapes:
            if shape == "circle":
                circle_crs = crosssections.loc[crosssections["shape"] == shape, :]
                circle_crs["definitionid"] = circle_crs.apply(
                    lambda x: "circ_d{:,.3f}_{:s}".format(x["diameter"], "branch"),
                    axis=1,
                )  # note diameter is reserved keywords in geopandas
                check_gpd_attributes(
                    circle_crs,
                    required_columns=[
                        "branch_id",
                        "branch_offset",
                        "frictionid",
                        "frictiontype",
                        "frictionvalue",
                        "diameter",
                        "shift",
                    ],
                )
                crosssections_ = pd.concat(
                    [crosssections_, _set_circle_crs(circle_crs)]
                )

    # setup thalweg for GUI
    crosssections_["crsdef_thalweg"] = 0.0

    # support both string and boolean for closed column
    if "crsdef_closed" in crosssections_:
        crosssections_["crsdef_closed"].replace({"yes": 1, "no": 0}, inplace=True)

    crosssections_ = gpd.GeoDataFrame(crosssections_, crs=branches.crs)

    return crosssections_.drop_duplicates()


def set_xyz_crosssections(
    branches: gpd.GeoDataFrame,
    crosssections: gpd.GeoDataFrame,
):
    """Set up xyz crosssections.
    xyz crosssections should be points gpd, column z and column order.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        The branches.
    crosssections : gpd.GeoDataFrame
        The crosssections.

    Returns
    -------
    pd.DataFrame
        The xyz cross sections.
    """
    # check if require columns exist
    required_columns = ["geometry", "crsid", "order", "z"]
    if set(required_columns).issubset(crosssections.columns):
        crosssections = gpd.GeoDataFrame(crosssections[required_columns])
    else:
        logger.error(
            f"Cannot setup crosssections from branch. Require columns {required_columns}."
        )

    # apply data type
    crosssections.loc[:, "x"] = crosssections.geometry.x
    crosssections.loc[:, "y"] = crosssections.geometry.y
    crosssections.loc[:, "z"] = crosssections.z
    crosssections.loc[:, "order"] = crosssections.loc[:, "order"].astype("int")

    # convert xyz crosssection into yz profile
    crosssections = crosssections.groupby(level=0).apply(xyzp2xyzl, (["order"]))
    crosssections.crs = branches.crs

    # snap to branch
    # setup branch_id - snap bridges to branch (inplace of bridges, will add branch_id and branch_offset columns)
    find_nearest_branch(
        branches=branches, geometries=crosssections, method="intersecting"
    )
    logger.warning(
        "Snapping to branches using intersection: Please double check if the crossection is closely located to a bifurcation."
    )

    # setup failed - drop based on branch_offset that are not snapped to branch (inplace of yz_crosssections) and issue warning
    _old_ids = crosssections.index.to_list()
    crosssections.dropna(axis=0, inplace=True, subset=["branch_offset"])
    crosssections = crosssections.merge(
        branches[["frictionid", "frictiontype", "frictionvalue"]],
        left_on="branch_id",
        right_index=True,
    )
    _new_ids = crosssections.index.to_list()
    if len(_old_ids) != len(_new_ids):
        logger.warning(
            f"Crosssection with id: {list(set(_old_ids) - set(_new_ids))} are dropped: unable to find closest branch. "
        )

    # setup crsdef from xyz
    crsdefs = pd.DataFrame(
        {
            "crsdef_id": crosssections.index.to_list(),
            "crsdef_branchid": crosssections.branch_id.to_list(),
            "crsdef_type": "xyz",
            "crsdef_xyzcount": crosssections.x.map(len).to_list(),
            "crsdef_xcoordinates": [
                " ".join([f"{i}" for i in l]) for l in crosssections.x.to_list()
            ],
            "crsdef_ycoordinates": [
                " ".join([f"{i}" for i in l]) for l in crosssections.y.to_list()
            ],
            "crsdef_zcoordinates": [
                " ".join([f"{i}" for i in l]) for l in crosssections.z.to_list()
            ],
            # --- start of yz ---
            # "crsdef_type": "yz",
            # "crsdef_yzCount": crosssections.l.map(len).to_list(),
            # 'crsdef_ycoordinates': [
            #      " ".join(["{:.1f}".format(i) for i in l])
            #      for l in crosssections.l.to_list()
            #  ],
            # "crsdef_zcoordinates": [
            #     " ".join(["{:.1f}".format(i) for i in l])
            #     for l in crosssections.z.to_list()
            # ],
            # --- end of yz ---
            # lower case key means temp keys (not written to file)
            "crsdef_frictionids": branches.loc[
                crosssections.branch_id.to_list(), "frictionid"
            ],
            "frictiontype": branches.loc[
                crosssections.branch_id.to_list(), "frictiontype"
            ],
            "frictionvalue": branches.loc[
                crosssections.branch_id.to_list(), "frictionvalue"
            ],
            "crsdef_frictionpositions": [
                f"0 {l}" for l in crosssections.geometry.length.to_list()
            ],
        }
    )

    # setup crsloc from xyz
    # delete generated ones
    crslocs = pd.DataFrame(
        {
            "crsloc_id": [
                f"{bid}_{bc:.2f}"
                for bid, bc in zip(
                    crosssections.branch_id.to_list(),
                    crosssections.branch_offset.to_list(),
                )
            ],
            "crsloc_branchid": crosssections.branch_id.to_list(),
            "crsloc_chainage": crosssections.branch_offset.to_list(),
            "crsloc_shift": 0.0,
            "crsloc_definitionid": crosssections.index.to_list(),
            "geometry": crosssections.geometry.centroid.to_list(),  # line to centroid. because crossection geom has point feature.
        }
    )
    crosssections_ = pd.merge(
        crslocs,
        crsdefs,
        how="left",
        left_on=["crsloc_definitionid"],
        right_on=["crsdef_id"],
    )

    crosssections_["crsdef_thalweg"] = 0.0

    crosssections_ = gpd.GeoDataFrame(crosssections_, crs=branches.crs)
    return crosssections_.drop_duplicates()


def set_point_crosssections(
    branches: gpd.GeoDataFrame, crosssections: gpd.GeoDataFrame, maxdist: float = 1.0
):
    """
    Function to set regular cross-sections from point.
    only support rectangle, trapezoid, circle and yz.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        Require index to contain the branch id
        The branches.
    crosssections : gpd.GeoDataFrame
        Required columns: shape,shift
        The crosssections.
    maxdist : float, optional
        the maximum distant that a point crossection is snapped to the branch.
        By default 1.0
    Returns
    -------
    gpd.GeoDataFrame
        The cross sections.
    """
    # check if crs mismatch
    if crosssections.crs != branches.crs:
        logger.error("mismatch crs between cross-sections and branches")

    # remove duplicated geometries
    _nodes = crosssections.copy()
    G = _nodes["geometry"].apply(lambda geom: geom.wkb)
    # check for diff in numbers: n = len(G) - len(G.drop_duplicates().index)
    crosssections = _nodes[_nodes.index.isin(G.drop_duplicates().index)]

    # snap to branch
    # setup branch_id - snap bridges to branch (inplace of bridges, will add branch_id and branch_offset columns)
    find_nearest_branch(
        branches=branches, geometries=crosssections, method="overal", maxdist=maxdist
    )

    # setup failed - drop based on branch_offset that are not snapped to branch (inplace of yz_crosssections) and issue warning
    _old_ids = crosssections.index.to_list()
    crosssections.dropna(axis=0, inplace=True, subset=["branch_offset"])
    _new_ids = crosssections.index.to_list()
    if len(_old_ids) != len(_new_ids):
        logger.warning(
            f"Crosssection with id: {list(set(_old_ids) - set(_new_ids))} are dropped: unable to find closest branch. "
        )

    if len(crosssections.branch_offset.dropna()) == 0:
        logger.error("No crossections are set up.")
        return pd.DataFrame()

    # get branch friction
    crosssections = crosssections.merge(
        branches[["frictionid", "frictiontype", "frictionvalue"]],
        left_on="branch_id",
        right_index=True,
    )

    # NOTE: below is removed because in case of multiple structures at the same location, there can be multiple crossections
    # # get a temporary id based on the convention of branch_id and branch_offset(precision 2 decimals)
    # crosssections["temp_id"] = crosssections.apply(lambda x:f"{x.branch_id}_{x.branch_offset:.2f}", axis = 1)
    # # drop duplicated temp_id, keep the one with minimum branch_distance
    # if crosssections["temp_id"].duplicated().any():
    #     logger.warning(f"Duplicate crosssections found, removing duplicates")
    #     # Sort DataFrame by branch_distance in ascending order
    #     crosssections_sorted = crosssections.sort_values('branch_distance')
    #     # Remove duplicates based on the branch_id, branch_offset column, keeping the first occurrence (with minimum branch_distance)
    #     crosssections = crosssections_sorted.drop_duplicates(subset=["temp_id"], keep='first')

    crosssections_ = pd.DataFrame()
    # loop through the shapes
    all_shapes = crosssections["shape"].unique().tolist()
    for shape in all_shapes:
        if shape == "circle":
            circle_crs = crosssections.loc[crosssections["shape"] == shape, :]
            circle_crs["definitionid"] = circle_crs.apply(
                lambda x: "circ_d{:,.3f}_{:s}".format(x["diameter"], "point"),
                axis=1,
            )
            check_gpd_attributes(
                circle_crs,
                required_columns=[
                    "branch_id",
                    "branch_offset",
                    "frictionid",
                    "frictiontype",
                    "frictionvalue",
                    "diameter",
                    "shift",
                ],
            )
            crosssections_ = pd.concat([crosssections_, _set_circle_crs(circle_crs)])
        elif shape == "rectangle":
            rectangle_crs = crosssections.loc[crosssections["shape"] == shape, :]
            rectangle_crs["definitionid"] = rectangle_crs.apply(
                lambda x: "rect_h{:,.3f}_w{:,.3f}_c{:s}_{:s}".format(
                    x["height"],
                    x["width"],
                    x["closed"],
                    "point",
                ),
                axis=1,
            )
            check_gpd_attributes(
                rectangle_crs,
                required_columns=[
                    "branch_id",
                    "branch_offset",
                    "frictionid",
                    "frictiontype",
                    "frictionvalue",
                    "width",
                    "height",
                    "closed",
                ],
            )
            crosssections_ = pd.concat(
                [crosssections_, _set_rectangle_crs(rectangle_crs)]
            )
        elif shape == "trapezoid":
            trapezoid_crs = crosssections.loc[crosssections["shape"] == shape, :]
            trapezoid_crs["definitionid"] = trapezoid_crs.apply(
                lambda x: "trapz_h{:,.1f}_bw{:,.1f}_tw{:,.1f}_c{:s}_{:s}".format(
                    x["height"],
                    x["width"],
                    x["t_width"],
                    x["closed"],
                    "point",
                ),
                axis=1,
            )
            check_gpd_attributes(
                trapezoid_crs,
                required_columns=[
                    "branch_id",
                    "branch_offset",
                    "frictionid",
                    "frictiontype",
                    "frictionvalue",
                    "width",
                    "height",
                    "t_width",
                    "closed",
                ],
            )
            crosssections_ = pd.concat(
                [crosssections_, _set_trapezoid_crs(trapezoid_crs)]
            )
        elif shape == "zw":
            zw_crs = crosssections.loc[crosssections["shape"] == shape, :]
            check_gpd_attributes(
                trapezoid_crs,
                required_columns=[
                    "branch_id",
                    "branch_offset",
                    "frictionid",
                    "frictiontype",
                    "frictionvalue",
                    "numlevels",
                    "levels",
                    "flowwidths",
                    "totalwidths",
                    "closed",
                ],
            )
            crosssections_ = pd.concat([crosssections_, _set_zw_crs(zw_crs)])
        elif shape == "yz":
            yz_crs = crosssections.loc[crosssections["shape"] == shape, :]
            check_gpd_attributes(
                trapezoid_crs,
                required_columns=[
                    "branch_id",
                    "branch_offset",
                    "frictionid",
                    "frictiontype",
                    "frictionvalue",
                    "yzcount",
                    "ycoordinates",
                    "zcoordinates",
                    "closed",
                ],
            )
            crosssections_ = pd.concat([crosssections_, _set_yz_crs(yz_crs)])
        else:
            logger.error(
                "crossection shape not supported. For now only support rectangle, trapezoid, zw and yz"
            )

    # setup thaiweg for GUI
    crosssections_["crsdef_thalweg"] = 0.0

    # support both string and boolean for closed column
    crosssections_["crsdef_closed"].replace({"yes": 1, "no": 0}, inplace=True)

    crosssections_ = gpd.GeoDataFrame(crosssections_, crs=branches.crs)

    return crosssections_


def _set_circle_crs(crosssections: gpd.GeoDataFrame):
    """Circle crossection."""
    crsdefs = []
    crslocs = []
    for c in crosssections.itertuples():
        crsdefs.append(
            {
                "crsdef_id": c.definitionid,
                "crsdef_type": "circle",
                "crsdef_branchid": c.branch_id,
                "crsdef_diameter": c.diameter,
                "crsdef_closed": "yes",
                "crsdef_frictionid": c.frictionid,
                "frictiontype": c.frictiontype,
                "frictionvalue": c.frictionvalue,
            }
        )
        crslocs.append(
            {
                "crsloc_id": f"{c.branch_id}_{c.branch_offset:.2f}",
                "crsloc_branchid": c.branch_id,
                "crsloc_chainage": c.branch_offset,
                "crsloc_shift": c.shift,
                "crsloc_definitionid": c.definitionid,
                "geometry": c.geometry,
            }
        )

    crosssections_ = pd.merge(
        pd.DataFrame.from_records(crslocs),
        pd.DataFrame.from_records(crsdefs).drop_duplicates(),
        how="left",
        left_on=["crsloc_branchid", "crsloc_definitionid"],
        right_on=["crsdef_branchid", "crsdef_id"],
    )

    crosssections_.index = crosssections.index

    return crosssections_


def _set_rectangle_crs(crosssections: gpd.GeoDataFrame):
    """Rectangle crossection."""
    crsdefs = []
    crslocs = []
    for c in crosssections.itertuples():
        crsdefs.append(
            {
                "crsdef_id": c.definitionid,
                "crsdef_type": "rectangle",
                "crsdef_branchid": c.branch_id,
                "crsdef_height": c.height,
                "crsdef_width": c.width,
                "crsdef_closed": c.closed,
                "crsdef_frictionid": c.frictionid,
                "frictiontype": c.frictiontype,
                "frictionvalue": c.frictionvalue,
            }
        )
        crslocs.append(
            {
                "crsloc_id": f"{c.branch_id}_{c.branch_offset:.2f}",
                "crsloc_branchid": c.branch_id,
                "crsloc_chainage": c.branch_offset,
                "crsloc_shift": c.shift,
                "crsloc_definitionid": c.definitionid,
                "geometry": c.geometry,
            }
        )

    crosssections_ = pd.merge(
        pd.DataFrame.from_records(crslocs),
        pd.DataFrame.from_records(crsdefs).drop_duplicates(),
        how="left",
        left_on=["crsloc_branchid", "crsloc_definitionid"],
        right_on=["crsdef_branchid", "crsdef_id"],
    )

    crosssections_.index = crosssections.index

    return crosssections_


def _set_trapezoid_crs(crosssections: gpd.GeoDataFrame):
    """Trapezoid need to be converted into zw type."""
    # check for non-valid trapezoid crs
    if (
        (crosssections["width"] <= 0).any()
        or (crosssections["t_width"] <= 0).any()
        or (crosssections["height"] <= 0).any()
    ):
        logger.error(
            "Invalid DataFrame: Found non-positive values in the 'width', 't_width', or 'height' columns."
        )

    crsdefs = []
    crslocs = []
    for c in crosssections.itertuples():
        levels = f"0 {c.height:.6f}"
        flowwidths = f"{c.width:.6f} {c.t_width:.6f}"
        crsdefs.append(
            {
                "crsdef_id": c.definitionid,
                "crsdef_type": "zw",
                "crsdef_branchid": c.branch_id,
                "crsdef_numlevels": 2,
                "crsdef_levels": levels,
                "crsdef_flowwidths": flowwidths,
                "crsdef_totalwidths": flowwidths,
                "crsdef_frictionid": c.frictionid,
                "frictiontype": c.frictiontype,
                "frictionvalue": c.frictionvalue,
                "crsdef_closed": c.closed,
            }
        )
        crslocs.append(
            {
                "crsloc_id": f"{c.branch_id}_{c.branch_offset:.2f}",
                "crsloc_branchid": c.branch_id,
                "crsloc_chainage": c.branch_offset,
                "crsloc_shift": c.shift,
                "crsloc_definitionid": c.definitionid,
                "geometry": c.geometry,
            }
        )

    # merge
    crosssections_ = pd.merge(
        pd.DataFrame.from_records(crslocs),
        pd.DataFrame.from_records(crsdefs).drop_duplicates(),
        how="left",
        left_on=["crsloc_branchid", "crsloc_definitionid"],
        right_on=["crsdef_branchid", "crsdef_id"],
    )
    return crosssections_


def _set_zw_crs(crosssections: gpd.GeoDataFrame):
    """Set zw profile."""
    crsdefs = []
    crslocs = []
    for c in crosssections.itertuples():
        crsdefs.append(
            {
                "crsdef_id": c.Index,
                "crsdef_type": "zw",
                "crsdef_branchid": c.branch_id,
                "crsdef_numlevels": c.numlevels,
                "crsdef_levels": c.levels,
                "crsdef_flowwidths": c.flowwidths,
                "crsdef_totalwidths": c.totalwidths,
                "crsdef_frictionid": c.frictionid,
                "frictiontype": c.frictiontype,
                "frictionvalue": c.frictionvalue,
            }
        )
        crslocs.append(
            {
                "crsloc_id": f"{c.branch_id}_{c.branch_offset:.2f}",
                "crsloc_branchid": c.branch_id,
                "crsloc_chainage": c.branch_offset,
                "crsloc_shift": c.shift,
                "crsloc_definitionid": c.Index,
                "geometry": c.geometry,
            }
        )

    crosssections_ = pd.merge(
        pd.DataFrame.from_records(crslocs),
        pd.DataFrame.from_records(crsdefs).drop_duplicates(),
        how="left",
        left_on=["crsloc_branchid", "crsloc_definitionid"],
        right_on=["crsdef_branchid", "crsdef_id"],
    )

    crosssections_.index = crosssections.index

    return crosssections_


def _set_yz_crs(crosssections: gpd.GeoDataFrame):
    """Set yz profile."""
    crsdefs = []
    crslocs = []
    for c in crosssections.itertuples():
        crsdefs.append(
            {
                "crsdef_id": c.Index,
                "crsdef_type": "yz",
                "crsdef_branchid": c.branch_id,
                "crsdef_yzcount": c.yzcount,
                "crsdef_ycoordinates": c.ycoordinates,
                "crsdef_zcoordinates": c.zcoordinates,
                "crsdef_frictionid": c.frictionid,
                "frictiontype": c.frictiontype,
                "frictionvalue": c.frictionvalue,
            }
        )
        crslocs.append(
            {
                "crsloc_id": f"{c.branch_id}_{c.branch_offset:.2f}",
                "crsloc_branchid": c.branch_id,
                "crsloc_chainage": c.branch_offset,
                "crsloc_shift": c.shift,
                "crsloc_definitionid": c.Index,
                "geometry": c.geometry,
            }
        )

    crosssections_ = pd.merge(
        pd.DataFrame.from_records(crslocs),
        pd.DataFrame.from_records(crsdefs).drop_duplicates(),
        how="left",
        left_on=["crsloc_branchid", "crsloc_definitionid"],
        right_on=["crsdef_branchid", "crsdef_id"],
    )

    crosssections_.index = crosssections.index

    return crosssections_


# def parse_sobek_crs(filename, logger=logger):
#     """read sobek crosssection files as a dataframe. Include location and definition file.
#     #TODO: include parsing geometry as well

#     Parameters
#     ----------
#     filename : Path
#         Path to the sobek crosssection files. supported format: .DAT amd .DEF
#     logger : logger, Optional

#     Raise
#     -----
#     NotImplementedError
#         do not support other files than .dat and .def

#     Returns
#     -------
#     df.DataFrame
#         The data frame with each item as a row
#     """
#     import shlex
#     from pathlib import Path

#     import numpy as np
#     import pandas as pd

#     # check file
#     if Path(filename).name.lower().endswith(".def"):
#         logger.info("Parsing cross section definition")
#         prefix = "CRDS"
#         suffix = "crds"
#     elif Path(filename).name.lower().endswith(".dat"):
#         logger.info("Parsing cross section location")
#         prefix = "CRSN"
#         suffix = "crsn"
#     else:
#         raise NotImplementedError("do not support other files than .dat and .def")

#     with open(filename) as myFile:
#         text = myFile.read()
#         raw_lines = text.split(suffix + "\n")

#     lines = []
#     for l in raw_lines:
#         if l.startswith(prefix):  # new item
#             # preliminary handling
#             l = l.removeprefix(prefix)
#             t = None
#             table_dict = {}
#             # parse zw profile
#             if "lt lw\nTBLE" in l:
#                 # the table contains height, total width en flowing width.
#                 l, t = l.split("lt lw\nTBLE")
#                 levels, totalwidths, flowwidths = np.array(
#                     [shlex.split(r, posix=False) for r in t.split("<")][:-1]
#                 ).T  # last element is the suffix of tble
#                 table_dict["numlevels"] = len(levels)
#                 table_dict["levels"] = " ".join(str(n) for n in levels)
#                 table_dict["totalwidths"] = " ".join(str(n) for n in totalwidths)
#                 table_dict["flowwidths"] = " ".join(str(n) for n in flowwidths)
#             # parse yz profile
#             if "lt yz\nTBLE" in l:
#                 # Y horizontal distance increasing from the left to right,
#                 # Z vertical distance increasing from bottom to top in m.
#                 # In other words, use a coordinate system to define the Y-Z profile.
#                 l, t = l.split("lt yz\nTBLE")
#                 ycoordinates, zcoordinates = np.array(
#                     [shlex.split(r, posix=False) for r in t.split("<")][:-1]
#                 ).T
#                 table_dict["yzcount"] = len(ycoordinates)
#                 table_dict["ycoordinates"] = " ".join(str(n) for n in ycoordinates)
#                 table_dict["zcoordinates"] = " ".join(str(n) for n in zcoordinates)
#                 # storage width on surface in m
#                 if "lt sw 0" in l:
#                     l.replace("lt lw 0", "lt_lw_0")  # remove space
#                 else:
#                     logger.error(
#                         "storage width function is not supported. Check lt sw field"
#                     )
#             # parse line
#             line = shlex.split(l, posix=False)
#             line_dict = {line[i]: line[i + 1] for i in range(0, len(line), 2)}
#             # add table
#             if t is not None:
#                 line_dict.update(table_dict)
#             lines.append(line_dict)

#     df = pd.DataFrame.from_records(lines)
#     df["id"] = df["id"].str.strip("'")
#     df.set_index("id", inplace=True)

#     return df


def xyzp2xyzl(xyz: pd.DataFrame, sort_by: list = ["x", "y"]):
    """Convert xyz points to xyz lines.

    Parameters
    ----------
    xyz: pd.DataFrame
        The xyz points.
    sort_by: list, optional
        List of attributes to sort by. Defaults to ["x", "y"].

    Returns
    -------
    gpd.GeoSeries
        The xyz lines.
    """
    sort_by = [s.lower() for s in sort_by]

    if xyz is not None:
        yz_index = xyz.index.unique()
        xyz.columns = [c.lower() for c in xyz.columns]
        xyz.reset_index(drop=True, inplace=True)

        # sort
        xyz_sorted = xyz.sort_values(by=sort_by)

        new_z = xyz_sorted.z.to_list()
        # temporary
        # new_z[0] = 1.4
        # new_z[-1] = 1.4

        line = LineString([(px, py) for px, py in zip(xyz_sorted.x, xyz_sorted.y)])
        xyz_line = pd.Series(
            {
                "geometry": line,
                "l": list(
                    np.r_[
                        0.0,
                        np.cumsum(
                            np.hypot(
                                np.diff(line.coords, axis=0)[:, 0],
                                np.diff(line.coords, axis=0)[:, 1],
                            )
                        ),
                    ]
                ),
                "index": yz_index.to_list()[0],
                "x": xyz_sorted.x.to_list(),
                "y": xyz_sorted.y.to_list(),
                "z": new_z,
            }
        )
    return xyz_line
