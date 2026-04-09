"""Utilities GIS functions for Delft3D FM model."""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from hydromt.gis import vector_utils
from shapely.geometry import (
    LineString,
    Point,
)
from shapely.ops import snap, split

logger = logging.getLogger(f"hydromt.{__name__}")


__all__ = [
    "split_lines",
    "cut_pieces",
    "check_gpd_attributes",
    "update_data_columns_attributes_based_on_filter",
    "get_gdf_from_branches",
]


## geometry
def cut_pieces(line, distances):
    """Cut a line into pieces based on distances."""
    if distances[0] != 0:
        distances.insert(0, 0)
    if distances[-1] == line.length:
        distances.pop(-1)
    pieces = [line]
    for d in np.diff(np.sort(distances)):
        line = pieces.pop(-1)
        pieces.extend(cut(line, d))
    return pieces


def cut(line, distance):
    """Cut a line in two at a distance from its starting point.

    ref: https://shapely.readthedocs.io/en/stable/manual.html.
    """
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [LineString(coords[: i + 1]), LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:]),
            ]


def split_lines(line, num_new_lines):
    """Get a list of lines splitted from a line.

    Parameters
    ----------
    line
        The line to split
    num_new_lines : int
        The desired number of lines.

    Returns
    -------
    list
        The new lines.
    """
    _line = [line]
    points = [
        line.interpolate((i / num_new_lines), normalized=True)
        for i in range(0, num_new_lines + 1)
    ]

    new_lines = []
    for n, p in enumerate(points):
        split_line = split(snap(line, p, 1e-8), p)
        segments = [feature for feature in split_line.geoms]
        if n == 0:
            line = segments[0]
        elif n == num_new_lines:
            new_lines.append(segments[0])
        else:
            new_lines.append(segments[0])
            line = segments[-1]

    assert (
        len(new_lines) == num_new_lines
    ), "number of lines after splitting does not match input"
    assert np.isclose(
        sum([l.length for l in new_lines]), _line[0].length
    ), "length after splitting does not match input"

    return new_lines


def check_gpd_attributes(
    gdf: gpd.GeoDataFrame, required_columns: list, raise_error: bool = False
):
    """Check if the geodataframe contains all required columns.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame, required
        GeoDataFrame to be checked
    required_columns: list of strings, optional
        Check if the geodataframe contains all required columns
    raise_error: boolean, optional
        Raise error if the check failed
    """
    if not (set(required_columns).issubset(gdf.columns)):
        if raise_error:
            raise ValueError(
                "GeoDataFrame does not contain all required attributes:"
                f"{required_columns}."
            )
        else:
            logger.warning(
                "GeoDataFrame does not contain all required attributes:"
                f"{required_columns}."
            )
        return False
    return True


def update_data_columns_attributes_based_on_filter(
    gdf: gpd.GeoDataFrame,
    df: pd.DataFrame,
    filter_column: str,
    filter_value: str = None,
):
    """
    Add or update columns in gdf based on column and values in df.

    If filter_column and filter_value is set, only update the attributes of the filtered
    geodataframe.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        geodataframe containing user input
    df : attribute DataFrame
        a pd.DataFrame with attribute columns and values (e.g. width =  1) per
        filter_value in the filter_column (e.g. branch_type = pipe)
    filter_column : str
        Name of the column linking df to gdf.
    filter_value: str
        Value of the filter in the filter column.
        Must be used with filter_column

    Raises
    ------
    ValueError:
        if Name of the column linking df to gdf (`filter_column`) is not specified

    Returns
    -------
    gdf : gpd.GeoDataFrame
        geodataframe containing user input filled with new attribute values.

    """
    if filter_column is None:
        raise ValueError("Name of the column linking df to gdf must be specified")

    if filter_value is not None:
        attributes = df[df[filter_column] == filter_value]
    else:
        attributes = df

    # Update attributes
    for i in range(len(attributes.index)):
        row = attributes.iloc[i, :]
        filter_value = row.loc[filter_column]
        for colname in row.index[1:]:
            # If attribute is not at all in branches, add a new column
            if colname not in gdf.columns:
                gdf[colname] = pd.Series(dtype=attributes[colname].dtype)
            # Then fill in empty or NaN values with defaults
            gdf.loc[
                np.logical_and(gdf[colname].isna(), gdf[filter_column] == filter_value),
                colname,
            ] = row.loc[colname]

    return gdf


def get_gdf_from_branches(
    branches: gpd.GeoDataFrame, df: pd.DataFrame
) -> gpd.GeoDataFrame:
    """Get geodataframe from dataframe.

    Based on interpolation of branches, using columns ["branchid", "chainage" in df].

    Parameters
    ----------
    branches:gpd.GeoDataFrame
        line geometries of the branches
        Required varaibles: ["branchid"/"branchid", "geometry" ]
    df:pd.DataFrame
        dataframe containing the features located on branches
        Required varaibles: ["branchid"/"branchid", "chainage" ]

    Return
    ------
    gdf:gpd.GeoDataFrame
        dataframe cotaining the features located on branches, with point geometry

    """
    branches.columns = branches.columns.str.lower()
    branches = branches.set_index("branchid")

    df["geometry"] = None

    # Iterate over each point and interpolate a point along the corresponding line
    for i, row in df.iterrows():
        line_geometry = branches.loc[row.branchid, "geometry"]
        new_point_geometry = line_geometry.interpolate(row.chainage)
        df.loc[i, "geometry"] = new_point_geometry

    return gpd.GeoDataFrame(df, crs=branches.crs)


# temporary until hydromt fixes crs bug issue #1413
def nearest_merge(
    gdf1: gpd.GeoDataFrame,
    gdf2: gpd.GeoDataFrame,
    *,
    columns: list | None = None,
    max_dist: float | None = None,
    overwrite: bool = False,
    inplace: bool = False,
) -> gpd.GeoDataFrame:
    """Merge attributes of gdf2 with the nearest feature of gdf1.

    Output is optionally bounded by a maximum distance `max_dist`.
    Unless `overwrite = True`, gdf2 values are only
    merged where gdf1 has missing values.

    Parameters
    ----------
    gdf1, gdf2: geopandas.GeoDataFrame
        Source `gdf1` and destination `gdf2` geometries.
    columns : list of str, optional
        Names of columns in `gdf2` to merge, by default None
    max_dist : float, optional
        Maximum distance threshold for merge, by default None, i.e.: no threshold.
    overwrite : bool, optional
        If False (default) gdf2 values are only merged where gdf1 has missing values,
        i.e. NaN values for existing columns or missing columns.
    inplace : bool,
        If True, apply the merge to gdf1, otherwise return a new object.

    Returns
    -------
    gpd.GeoDataFrame
        Merged GeoDataFrames
    """
    idx_nn, dst = vector_utils.nearest(gdf1, gdf2)
    if not inplace:
        gdf1 = gdf1.copy()
    valid = dst < max_dist if max_dist is not None else np.ones_like(idx_nn, dtype=bool)
    columns = gdf2.columns if columns is None else columns
    gdf1["distance_right"] = dst
    gdf1["index_right"] = -1
    gdf1.loc[valid, "index_right"] = idx_nn[valid]
    skip = ["geometry"]
    for col in columns:
        if col in skip or col not in gdf2:
            if col not in gdf2:
                logger.warning(f"Column {col} not found in gdf2 and skipped.")
            continue
        new_vals = gdf2.loc[idx_nn[valid], col].values
        if col in gdf1 and not overwrite:
            old_vals = gdf1.loc[valid, col]
            replace = np.logical_or(old_vals.isnull(), old_vals.eq(""))
            new_vals = np.where(replace, new_vals, old_vals)
        gdf1.loc[valid, col] = new_vals

    return gdf1
