# -*- coding: utf-8 -*-

import logging
from typing import List, Literal

import geopandas as gpd
import numpy as np
import pandas as pd

from .branches import find_nearest_branch
from .crosssections import set_point_crosssections
from ..gis_utils import update_data_columns_attributes_based_on_filter

logger = logging.getLogger(__name__)


__all__ = [
    "prepare_1dstructures",
]


def prepare_1dstructures(
    branches: gpd.GeoDataFrame,
    gdf_st: gpd.GeoDataFrame,
    params: pd.DataFrame,
    st_type: Literal["culvert", "bridge"],
    id_start: int = 1,
    filter: str = None,
    snap_offset: float = 0.0,
    logger: logging.Logger = logger,
) -> gpd.GeoDataFrame:
    """Prepare 1D structures from geodataframe.

    Include the universal weir (Not Implemented) , culvert and bridge (``str_type`` = 'bridge'), which can only be used in a single 1D channel.

    Structures are first filtered for value specified in ``filter`` on the column ``st_type``.
    They are then snapped to the existing network within a max distance defined in ``snap_offset`` and will be dropped if not snapped.
    Finally, crossections are read and set up for the remaining structures.

    Parameters
    ----------
    branches : gpd.GeoDataFrame
        gpd.GeoDataFrame of branches.
    gdf_st : gpd.GeoDataFrame
        gpd.GeoDataFrame of structures.
    params : pd.DataFrame
        pd.Dataframe of defaults values for gdf_st.
    st_type : str
        structure type. Either "culvert" or "bridge".
    id_start: int, optional
        Start index for structure id. By default 1.
    filter: str, optional
        Keyword in structure_type column of gdf_st used to filter features.
    snap_offset: float, optional
        Snapping tolerance to automatically snapping to branch.
        By default 0.0, no snapping is applied.
    logger: logging.Logger, optional
        Logger.

    Returns
    -------
    gpd.GeoDataFrame
        Prepared structures with structure_id, structure_type, branchid and chainage.
    """
    type_col = "structure_type"
    id_col = "structure_id"

    # 1. prepare branches and gdf_st
    branches = branches.set_index("branchid")
    gdf_st = gdf_st.to_crs(branches.crs)

    # 2. Filter features based on filter on type_col
    if filter is not None and type_col in gdf_st.columns:
        gdf_st = gdf_st[gdf_st[type_col].str.lower() == filter.lower()]
        logger.info(f"Set {filter} locations filtered from structurs as {st_type} .")
    # Check if features in region
    if len(gdf_st) == 0:
        logger.warning(f"No 1D {type} locations found within domain")
        return None

    # 3. Add defaults
    # overwrite type and add id attributes if does not exist
    gdf_st[type_col] = pd.Series(
        data=np.repeat(st_type, len(gdf_st)), index=gdf_st.index, dtype=str
    )
    if id_col not in gdf_st.columns:
        data = [
            f"{st_type}_{i}" for i in np.arange(id_start, id_start + len(gdf_st))
        ]  # avoid duplicated ids being generated
        gdf_st[id_col] = pd.Series(data, index=gdf_st.index, dtype=str)
    # assign id
    gdf_st.index = gdf_st[id_col]
    gdf_st.index.name = id_col
    # filter for allowed columns
    allowed_columns = set(params.columns).intersection(gdf_st.columns)
    allowed_columns.update({"geometry"})
    gdf_st = gpd.GeoDataFrame(gdf_st[list(allowed_columns)], crs=gdf_st.crs)
    logger.info("Adding/Filling default attributes values")
    gdf_st = update_data_columns_attributes_based_on_filter(
        gdf_st, params, type_col, st_type
    )

    # 4. snap structures to branches
    # setup branch_id - snap structures to branch (inplace of structures, will add branch_id and branch_offset columns)
    find_nearest_branch(branches=branches, geometries=gdf_st, maxdist=snap_offset)
    # setup failed - drop based on branch_offset that are not snapped to branch
    _old_ids = gdf_st.index.to_list()
    gdf_st.dropna(axis=0, inplace=True, subset=["branch_offset"])
    _new_ids = gdf_st.index.to_list()
    if len(_old_ids) != len(_new_ids):
        logger.warning(
            f"structure with id: {list(set(_old_ids) - set(_new_ids))} are dropped: unable to find closest branch. "
        )
    if len(_new_ids) == 0:
        logger.warning(
            f"No 1D {type} locations found within the proximity of the network"
        )
        return None
    else:
        # setup success, add branchid and chainage
        gdf_st["structure_branchid"] = gdf_st["branch_id"]
        gdf_st["structure_chainage"] = gdf_st["branch_offset"]

    # 5. add structure crossections
    # add a dummy "shift" for crossection locations if missing (e.g. culverts), because structures only needs crossection definitions.
    if "shift" not in gdf_st.columns:
        gdf_st["shift"] = np.nan
    # derive crosssections
    gdf_st_crossections = set_point_crosssections(branches, gdf_st, maxdist=snap_offset)
    # remove crossection locations and any friction from the setup
    gdf_st_crsdefs = gdf_st_crossections.drop(
        columns=[
            c
            for c in gdf_st_crossections.columns
            if c.startswith("crsloc") or "friction" in c
        ]
    )
    # add to structures
    gdf_st = gdf_st.merge(
        gdf_st_crsdefs.drop(columns="geometry"), left_index=True, right_index=True
    )

    # 6. replace np.nan as None
    gdf_st = gdf_st.replace(np.nan, None)

    # 7. remove index
    gdf_st = gdf_st.reset_index()

    return gdf_st
