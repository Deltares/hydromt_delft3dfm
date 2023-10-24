"""Workflows to prepare boundaries for Delft3D-FM model."""

import logging
from pathlib import Path

import geopandas as gpd
import hydromt.io
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Point, Polygon

from hydromt_delft3dfm import graph_utils

logger = logging.getLogger(__name__)


__all__ = [
    "get_boundaries_with_nodeid",
    "generate_boundaries_from_branches",
    "select_boundary_type",
    "validate_boundaries",
    "compute_boundary_values",
    "compute_2dboundary_values",
    "compute_meteo_forcings",
    "compute_forcing_values_points",
    "compute_forcing_values_polygon",
    "compute_forcing_values_lines",
    "get_geometry_coords_for_linestrings",
    "get_geometry_coords_for_polygons",
]


def get_boundaries_with_nodeid(
    branches: gpd.GeoDataFrame, network1d_nodes: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Get boundary locations from branches and associate with node IDs.

    Parameters
    ----------
    branches: A GeoDataFrame containing the branches of the network.
    network1d_nodes: A GeoDataFrame containing the network1d nodes with node IDs.

    Returns
    -------
    A GeoDataFrame with boundary locations and their associated node IDs.
    """
    # generate all possible and allowed boundary locations
    _boundaries = generate_boundaries_from_branches(branches, where="both")

    boundaries = hydromt.gis_utils.nearest_merge(
        _boundaries, network1d_nodes, max_dist=0.1, overwrite=False
    )
    return boundaries


def generate_boundaries_from_branches(
    branches: gpd.GeoDataFrame, where: str = "both"
) -> gpd.GeoDataFrame:
    """Get the possible boundary locations from the branches with id.

    Parameters
    ----------
    where : {'both', 'upstream', 'downstream'}
        Where at the branches should the boundaries be derived.
        An upstream end node is defined as a node which has 0 incoming
        branches and 1 outgoing branch.
        A downstream end node is defined as a node which has 1 incoming
        branch and 0 outgoing branches.

    Returns
    -------
    gpd.GeoDataFrame
        A data frame containing all the upstream and downstream
        end nodes of the branches
    """
    # convert branches to graph
    G = graph_utils.gpd_to_digraph(branches)

    # get boundary locations at where
    if where == "downstream":
        endnodes = {
            dn: {**d, **{"where": "downstream"}}
            for up, dn, d in G.edges(data=True)
            if G.out_degree[dn] == 0 and G.degree[dn] == 1
        }
    elif where == "upstream":
        endnodes = {
            up: {**d, **{"where": "upstream"}}
            for up, dn, d in G.edges(data=True)
            if G.in_degree[up] == 0 and G.degree[up] == 1
        }
    elif where == "both":
        endnodes = {
            dn: {**d, **{"where": "downstream"}}
            for up, dn, d in G.edges(data=True)
            if G.out_degree[dn] == 0 and G.degree[dn] == 1
        }
        endnodes.update(
            {
                up: {**d, **{"where": "upstream"}}
                for up, dn, d in G.edges(data=True)
                if G.in_degree[up] == 0 and G.degree[up] == 1
            }
        )
    else:
        pass

    if len(endnodes) == 0:
        logger.error(f"cannot generate boundaries for given condition {where}")

    endnodes_pd = (
        pd.DataFrame().from_dict(endnodes, orient="index").drop(columns=["geometry"])
    )
    endnodes_gpd = gpd.GeoDataFrame(
        data=endnodes_pd,
        geometry=[Point(endnode) for endnode in endnodes],
        crs=branches.crs,
    )
    endnodes_gpd.reset_index(inplace=True)
    return endnodes_gpd


def select_boundary_type(
    boundaries: gpd.GeoDataFrame,
    branch_type: str,
    boundary_type: str,
    boundary_locs: str,
    logger=logger,
) -> pd.DataFrame:
    """Select boundary location per branch type and boundary type.

    Parameters
    ----------
    boundaries : gpd.GeoDataFrame
        The boundaries.
    branch_type : {'river', 'pipe'}
        The branch type.
    boundary_type : {'waterlevel', 'discharge'}
        For rivers 'waterlevel' and 'discharge' are supported.
        For pipes 'waterlevel' is supported.
    boundary_locs : {'both', 'upstream', 'downstream'}
        The boundary location to use.
    logger
        The logger to log messages with.

    Returns
    -------
    pd.DataFrame
        A data frame containing the boundary location per branch type
        and boundary type.
    """
    boundaries_branch_type = boundaries.loc[boundaries["branchtype"] == branch_type, :]
    if branch_type == "river":
        if boundary_type == "waterlevel":
            if boundary_locs != "both":
                boundaries_branch_type = boundaries_branch_type.loc[
                    boundaries_branch_type["where"] == boundary_locs, :
                ]
        elif boundary_type == "discharge":
            if boundary_locs != "upstream":
                logger.warning(
                    f"Applying boundary type {boundary_type} selected"
                    f"for {branch_type} boundaries might cause instabilities."
                )
            if boundary_locs != "both":
                boundaries_branch_type = boundaries_branch_type.loc[
                    boundaries_branch_type["where"] == boundary_locs, :
                ]
        else:
            logger.error(
                f"Wrong boundary type {boundary_type} selected"
                f"for {branch_type} boundaries."
            )
    # TODO: extend
    # for now only downstream boundaries not connected to openwater
    # later add connected to weir (upstream or downstream or intermediate)
    elif branch_type == "pipe":
        if boundary_type == "waterlevel":
            boundaries_branch_type = boundaries_branch_type.loc[
                boundaries_branch_type["where"] == "downstream", :
            ]
        else:
            logger.error(
                f"Wrong boundary type {boundary_type} selected"
                f"for {branch_type} boundaries."
            )

    return boundaries_branch_type


def validate_boundaries(boundaries: gpd.GeoDataFrame, branch_type: str = "river"):
    """Validate boundaries per branch type.

    Will log a warning if the validation fails.

    Parameters
    ----------
    boundaries : gpd.GeoDataFrame
        The boundaries.
    branch_type : {'river', 'pipe'}
        The branch type.

    """
    if branch_type == "river":  # TODO add other open system branch_type
        for _, bnd in boundaries.iterrows():
            # TODO extended
            if bnd["where"] == "downstream" and bnd["boundary_type"] == "discharge":
                logger.warning(
                    "Boundary type violates modeller suggestions: using"
                    f"downstream discharge boundary at branch {bnd['branchid']}"
                )

    if branch_type == "pipe":  # TODO add other close system branch_type
        for _, bnd in boundaries.iterrows():
            # TODO extended
            if bnd["where"] == "upstream":
                logger.warning(
                    "Boundary type violates modeller suggestions:"
                    f"using upstream boundary at branch {bnd['branchid']}"
                )


def compute_boundary_values(
    boundaries: gpd.GeoDataFrame,
    da_bnd: xr.DataArray = None,
    boundary_value: float = -2.5,
    boundary_type: str = "waterlevel",
    boundary_unit: str = "m",
    snap_offset: float = 0.1,
    logger=logger,
):
    """
    Compute 1d boundary values.

    Parameters
    ----------
    boundaries : gpd.GeoDataFrame
        Point locations of the 1D boundaries to which to add data.

        * Required variables: ['nodeid']
    da_bnd : xr.DataArray, optional
        xr.DataArray containing the boundary timeseries values.
        If None, uses a constant values for all boundaries.

        * Required variables if netcdf: [``boundary_type``]
    boundary_value : float, optional
        Constant value to use for all boundaries if ``da_bnd`` is
        None and to fill in missing data. By default -2.5 m.
    boundary_type : {'waterlevel', 'discharge'}
        Type of boundary to use. By default "waterlevel".
    boundary_unit : {'m', 'm3/s'}
        Unit corresponding to [boundary_type].
        If ``boundary_type`` = "waterlevel"
            Allowed unit is [m]
        if ''boundary_type`` = "discharge":
            Allowed unit is [m3/s]
        By default m.
    snap_offset : float, optional
        Snapping tolerance to automatically applying boundaries
        at the correct network nodes. By default 0.1,
        a small snapping is applied to avoid precision errors.
    logger
        Logger to log messages.
    """
    # Timeseries boundary values
    if da_bnd is not None:
        logger.info(f"Preparing 1D {boundary_type} boundaries from timeseries.")

        # snap user boundary to potential boundary locations to get nodeid
        gdf_bnd = da_bnd.vector.to_gdf()
        gdf_bnd.crs = (
            boundaries.crs
        )  # FIXME temp fix for hydromt reprojection issue #613
        gdf_bnd = hydromt.gis_utils.nearest_merge(
            gdf_bnd,
            boundaries,
            max_dist=snap_offset,
            overwrite=True,
        )
        gdf_bnd = gdf_bnd[~gdf_bnd["nodeid"].isna()]
        da_bnd = da_bnd.sel(index=gdf_bnd.index)

        # get forcing data time indes
        bd_times, freq_name = _standardize_forcing_timeindexes(da_bnd)

        # instantiate xr.DataArray for bnd data
        da_out = xr.DataArray(
            data=da_bnd.data,
            dims=["index", "time"],
            coords=dict(
                index=gdf_bnd["nodeid"],
                time=bd_times,
                x=("index", gdf_bnd.geometry.x.values),
                y=("index", gdf_bnd.geometry.y.values),
            ),
            attrs=dict(
                function="TimeSeries",
                timeInterpolation="Linear",
                quantity=f"{boundary_type}bnd",
                units=f"{boundary_unit}",
                time_unit=f"{freq_name} since {pd.to_datetime(da_bnd.time[0].values)}",
            ),
        )

        # fill in na using default
        da_out = da_out.fillna(boundary_value)

        # drop na in time
        da_out.dropna(dim="time")

        # add name
        da_out.name = f"{boundary_type}bnd"
    else:
        logger.info(
            f"Using constant value {boundary_value} {boundary_unit}"
            f"for all {boundary_type} boundaries."
        )
        # instantiate xr.DataArray for bnd data with boundary_value directly
        da_out = xr.DataArray(
            data=np.full((len(boundaries.index)), boundary_value, dtype=np.float32),
            dims=["index"],
            coords=dict(
                index=boundaries["nodeid"],
                x=("index", boundaries.geometry.x.values),
                y=("index", boundaries.geometry.y.values),
            ),
            attrs=dict(
                function="constant",
                offset=0.0,
                factor=1.0,
                quantity=f"{boundary_type}bnd",
                units=f"{boundary_unit}",
            ),
        )
        da_out.name = f"{boundary_type}bnd"

    return da_out.drop_duplicates(dim=...)


def compute_2dboundary_values(
    boundaries: gpd.GeoDataFrame = None,
    df_bnd: pd.DataFrame = None,
    boundary_value: float = 0.0,
    boundary_type: str = "waterlevel",
    boundary_unit: str = "m",
    logger=logger,
):
    """
    Compute 2d boundary timeseries.

    Line geometry will be converted into supporting points.
    Note that All quantities are specified per support point,
    except for discharges which are specified per polyline.

    Parameters
    ----------
    boundaries : gpd.GeoDataFrame, optional
        line geometry type of locations of the 2D boundaries to
        which to add data. Must be combined with ``df_bnd``.

        * Required variables: ["boundary_id"]
    df_bnd : pd.DataFrame, optional
        pd.DataFrame containing the boundary timeseries values.
        Must be combined with ``boundaries``. Columns must match the
        "boundary_id" in ``boundaries``.

        * Required variables: ["time"]
    boundary_value : float, optional
        Constant value to fill in missing data. By default 0 m.
    boundary_type : {'waterlevel', 'discharge'}
        Type of boundary to use. By default "waterlevel".
    boundary_unit : {'m', 'm3/s'}
        Unit corresponding to [boundary_type].
        If ``boundary_type`` = "waterlevel"
            Allowed unit is [m]
        if ''boundary_type`` = "discharge":
            Allowed unit is [m3/s]
        By default m.
    logger :
        Logger to log messages.

    Raises
    ------
    ValueError:
        if no boundary to compute.
    """
    # Timeseries boundary values
    if boundaries is None or len(boundaries) == 0:
        raise ValueError("No boundary to compute.")
    else:
        # prepare boundary data
        # get data freq in seconds
        _TIMESTR = {"D": "days", "H": "hours", "T": "minutes", "S": "seconds"}
        dt = df_bnd.time[1] - df_bnd.time[0]
        freq = dt.resolution_string
        multiplier = 1
        if freq == "D":
            logger.warning(
                "time unit days is not supported by the current GUI version: 2022.04"
            )  # converting to hours as temporary solution
            # FIXME: day is supported in version 2023.02,
            # general question: where to indicate gui version?
            multiplier = 24
        if len(
            pd.date_range(df_bnd.iloc[0, :].time, df_bnd.iloc[-1, :].time, freq=dt)
        ) != len(df_bnd.time):
            logger.error("does not support non-equidistant time-series.")
        freq_name = _TIMESTR[freq]
        freq_step = getattr(dt.components, freq_name)
        bnd_times = np.array([(i * freq_step) for i in range(len(df_bnd.time))])
        if multiplier == 24:
            bnd_times = np.array(
                [(i * freq_step * multiplier) for i in range(len(df_bnd.time))]
            )
            freq_name = "hours"

        # for each boundary apply boundary data
        da_out_dict = {}
        for _index, _bnd in boundaries.iterrows():
            bnd_id = _bnd["boundary_id"]

            # convert line to points
            support_points = pd.DataFrame(
                np.array([[x, y] for x, y in _bnd.geometry.coords[:]]),
                columns=["x", "y"],
            )
            support_points["_id"] = support_points.index + 1
            support_points["id"] = support_points["_id"].astype(str)
            support_points["id"] = support_points["id"].str.zfill(4)
            support_points["name"] = support_points.astype(str).apply(
                lambda x: f"{bnd_id}_{x.id}", axis=1
            )

            # instantiate xr.DataArray for bnd data with boundary_value directly
            da_out = xr.DataArray(
                data=np.full(
                    (len(support_points["name"]), len(bnd_times)),
                    np.tile(df_bnd[bnd_id].values, (len(support_points["name"]), 1)),
                    dtype=np.float32,
                ),
                dims=["index", "time"],
                coords=dict(
                    index=support_points["name"],
                    time=bnd_times,
                    x=("index", support_points.x.values),
                    y=("index", support_points.y.values),
                ),
                attrs=dict(
                    locationfile=bnd_id + ".pli",
                    function="TimeSeries",
                    timeInterpolation="Linear",
                    quantity=f"{boundary_type}bnd",
                    units=f"{boundary_unit}",
                    time_unit=f"{freq_name} since {pd.to_datetime(df_bnd.time[0])}",
                    # support only yyyy-mm-dd HH:MM:SS
                ),
            )
            # fill in na using default
            da_out = da_out.fillna(boundary_value)
            da_out.name = f"{bnd_id}"
            da_out_dict.update({f"{bnd_id}": da_out})

    return da_out_dict


def _standardize_forcing_timeindexes(da):
    """Standardize timeindexes frequency based on forcing DataArray"""
    _TIMESTR = {"D": "days", "H": "hours", "T": "minutes", "S": "seconds"}
    dt = pd.to_timedelta((da.time[1].values - da.time[0].values))
    freq = dt.resolution_string
    multiplier = 1
    if freq == "D":
        logger.warning(
            "time unit days is not supported by the current GUI version: 2022.04"
        )  # converting to hours as temporary solution # FIXME: day is converted to hours temporarily
        multiplier = 24
    if len(pd.date_range(da.time[0].values, da.time[-1].values, freq=dt)) != len(
        da.time
    ):
        logger.error("does not support non-equidistant time-series.")
    freq_name = _TIMESTR[freq]
    freq_step = getattr(dt.components, freq_name)
    bd_times = np.array([float(i * freq_step) for i in range(len(da.time))])
    if multiplier == 24:
        bd_times = np.array([(i * freq_step * multiplier) for i in range(len(da.time))])
        freq_name = "hours"
    return bd_times, freq_name


def get_geometry_coords_for_linestrings(gdf):
    """Gets xarray DataArray coordinates that describes linestring geometries.
    Inlcudes numcoordinates, xcoordinates and ycoordinates"""
    if gdf.geometry.type.iloc[0] == "LineString":
        # Get the maximum number of coordinates for any polygon
        max_coords = gdf["geometry"].apply(lambda x: len(x.coords[:])).max()

        def get_xcoords(geom):
            coords = [xy[0] for xy in geom.coords[:]]
            return np.pad(
                coords,
                (0, max_coords - len(coords)),
                "constant",
                constant_values=np.nan,
            )

        def get_ycoords(geom):
            coords = [xy[1] for xy in geom.coords[:]]
            return np.pad(
                coords,
                (0, max_coords - len(coords)),
                "constant",
                constant_values=np.nan,
            )

        # Create the 2D arrays
        x_2d = np.vstack(gdf["geometry"].apply(get_xcoords))
        y_2d = np.vstack(gdf["geometry"].apply(get_ycoords))

        return dict(
            index=gdf.index,
            numcoordinates=np.arange(max_coords),
            xcoordinates=(("index", "numcoordinates"), x_2d),
            ycoordinates=(("index", "numcoordinates"), y_2d),
        )


def compute_forcing_values_lines(
    gdf: gpd.GeoDataFrame,
    da: xr.DataArray = None,
    forcing_value: float = 0.0,
    forcing_type: str = "dischargebnd",
    forcing_unit: str = "m3/s",
    logger=logger,
):
    """
    Compute 2d forcing values.

    Used for 2d boundaries from line geometries.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame of lines to add 2D forcing data.

        * Required variables: ['geometry']
    da : xr.DataArray, optional
        xr.DataArray containing the forcing timeseries values.
        If None, uses a constant ``forcing_value`` for all forcings.

        * Required variables: ['forcing_type']

    forcing_value : float, optional
        Constant value to use for all forcings if ``da`` is None and to
        fill in missing data.
        By default 0.0 ``forcing_unit``
    forcing_type : {'dischargebnd', 'waterlevelbnd'}
        Type of forcing to use.
        For now only support 'dischargebnd' and 'waterlevelbnd'.
        By default 'discharge'
    forcing_unit : {'m3/s', 'm'}
        Unit corresponding to ``forcing_type``.
        For now only support 'm3/s' and 'm'.
        By default 'm3/s'
    logger
        Logger to log messages.
    """

    # first process data based on either timeseries or constant
    # then update data based on either nodes or branches
    # Timeseries forcing values
    if da is not None:
        logger.info(f"Preparing 1D forcing type {forcing_type} from timeseries.")

        # get forcing data freq in seconds
        bd_times, freq_name = _standardize_forcing_timeindexes(da)

        # instantiate xr.DataArray for forcing data
        coords_dict = get_geometry_coords_for_linestrings(gdf)

        # Prepare the data
        data_3d = np.tile(
            np.expand_dims(da.data, axis=-1), (1, 1, len(coords_dict["numcoordinates"]))
        )
        # instantiate xr.DataArray for forcing data
        # NOTE only support points on branches
        da_out = xr.DataArray(
            data=data_3d,
            dims=("index", "time", "numcoordinates"),
            coords=dict(
                index=coords_dict["index"],
                time=bd_times,
                numcoordinates=coords_dict["numcoordinates"],
                x=coords_dict["xcoordinates"],
                y=coords_dict["ycoordinates"],
            ),
            attrs=dict(
                function="TimeSeries",
                timeInterpolation="Linear",
                quantity=f"{forcing_type}",
                units=f"{forcing_unit}",
                time_unit=f"{freq_name} since {pd.to_datetime(da.time[0].values)}",  # support only yyyy-mm-dd HH:MM:SS
            ),
        )
        # fill in na using default
        da_out = da_out.fillna(forcing_value)

        # drop na in time
        da_out.dropna(dim="time")

        # add name
        da_out.name = f"{forcing_type}"
    else:
        logger.info(
            f"Using constant value {forcing_value} {forcing_unit} for all {forcing_type} forcings."
        )
        # instantiate xr.DataArray for forcing data with forcing_type directly
        coords_dict = get_geometry_coords_for_linestrings(gdf)
        data_3d = np.full(
            (len(coords_dict["index"]), len(coords_dict["numcoordinates"])),
            forcing_value,
            dtype=np.float32,
        )
        da_out = xr.DataArray(
            data=data_3d,
            dims=("index", "numcoordinates"),
            coords={
                "index": coords_dict["index"],
                "numcoordinates": coords_dict["numcoordinates"],
                "x": coords_dict["xcoordinates"],
                "y": coords_dict["ycoordinates"],
            },
            attrs=dict(
                function="constant",
                offset=0.0,
                factor=1.0,
                quantity=f"{forcing_type}",
                units=f"{forcing_unit}",
            ),
        )
        da_out.name = f"{forcing_type}"

    return da_out.drop_duplicates(dim=...)


def gpd_to_pli(gdf: gpd.GeoDataFrame, output_dir: Path):
    """Convert geopandas GeoDataFrame (gdf) into pli files.

    Pli files at 'output_dir' directory.

    the geodataframe must has index as stations and geometry
    of the stations.
    each row of the geodataframe will be converted into a single pli file.

    the file name and the station name will be the index of that row.
    """
    for _, g in gdf.iterrows():
        pli_name = g.index
        pli_coords = g.geometry.coords[:]
        with open(output_dir.joinpath(f"{pli_name}.pli"), "w") as f:
            f.write(f"{pli_name}\n")
            f.write(f"\t{len(pli_coords)} {2}\n")
            for p in pli_coords:
                f.write(f"\t{' '.join(str(pi) for pi in p)}\n")


def df_to_bc(
    df,
    output_dir,
    output_filename="boundary",
    quantity="discharge",
    unit="m3/s",
    freq="H",
):
    """Convert pandas timeseires 'df' into bc file.

    bc file from 'output_dir'/'output_filename'.bc

    the time series must has time as index, columns names as stations.
    the time series will be first converted into a equidistance timeseries
    with frequency specified in 'freq'. support [D, H,M,S]
    each columns-wise array will be converted into one bc timeseries.

    The time series has the quantity and unit as specified
    in 'quantity' and 'unit'.
    """
    time_unit = {"D": "days", "H": "hours", "M": "minutes", "S": "seconds"}

    df = df.resample(freq).ffill()
    time = df.index
    stations = df.columns

    with open(output_dir.joinpath(f"{output_filename}.bc"), "w") as f:
        f.write("[General]\n")
        f.write("\tfileVersion = 1.01\n")
        f.write("\tfileType = boundConds\n")
        for s in stations:
            d = df[s]
            f.write("\n")
            f.write("[forcing]\n")
            f.write(f"\tName = {d.name}\n")
            f.write("\tfunction = timeSeries\n")
            f.write("\ttimeInterpolation = linear\n")
            f.write(f"\tquantity = {quantity}\n")
            f.write(f"\tunit = {unit}\n")
            f.write("\tquantity = time\n")
            f.write(f"\tunit = {time_unit[freq]} since {time[0].date()}\n")
            f.write("\t0 0\n")
            for i, di in enumerate(d.values):
                f.write(f"\t{i} {di}\n")


def compute_meteo_forcings(
    df_meteo: pd.DataFrame = None,
    fill_value: float = 0.0,
    is_rate: bool = True,
    meteo_location: tuple = None,
    logger=logger,
) -> xr.DataArray:
    """
    Compute meteo forcings.

    Parameters
    ----------
    df_meteo : pd.DataFrame, optional
        pd.DataFrame containing the meteo timeseries values.
        If None, uses ``fill_value``.

        * Required variables: ["precip"]
    meteo_value : float, optional
        Constant value to use for global meteo if ``df_meteo`` is None and to
        fill in missing data in ``df_meteo``. By default 0.0 mm/day.
    is_rate : bool, optional
        Specify if the type of meteo data is direct "rainfall" (False)
        or "rainfall_rate" (True). By default True for "rainfall_rate".
        Note that Delft3DFM 1D2D Suite 2022.04 supports only "rainfall_rate".
        If rate, unit is expected to be in mm/day and else mm.
    meteo_location : tuple
        Global location for meteo timeseries
    logger
        Logger to log messages.

    Returns
    -------
    da_meteo : xr.DataArray
        xr.DataArray containing the meteo timeseries values. If None, uses ``df_meteo``.

        * Required variables if netcdf: [``precip``]
    """
    # Set units and type
    if is_rate:
        meteo_type = "rainfall_rate"
        meteo_unit = "mm/day"
    else:
        meteo_type = "rainfall"
        meteo_unit = "mm"

    # Timeseries boundary values

    logger.info("Preparing global (spatially uniform) timeseries.")
    # get data freq in seconds
    _TIMESTR = {"D": "days", "H": "hours", "T": "minutes", "S": "seconds"}
    dt = df_meteo.time[1] - df_meteo.time[0]
    freq = dt.resolution_string
    multiplier = 1
    if freq == "D":
        logger.warning(
            "time unit days is not supported by the current GUI version: 2022.04"
        )  # converting to hours as temporary solution
        # FIXME: day is converted to hours temporarily
        multiplier = 24
    if len(
        pd.date_range(df_meteo.iloc[0, :].time, df_meteo.iloc[-1, :].time, freq=dt)
    ) != len(df_meteo.time):
        logger.error("does not support non-equidistant time-series.")
    freq_name = _TIMESTR[freq]
    freq_step = getattr(dt.components, freq_name)
    meteo_times = np.array([(i * freq_step) for i in range(len(df_meteo.time))])
    if multiplier == 24:
        meteo_times = np.array(
            [(i * freq_step * multiplier) for i in range(len(df_meteo.time))]
        )
        freq_name = "hours"
    # instantiate xr.DataArray for global time series
    da_out = xr.DataArray(
        data=np.full((1, len(df_meteo)), df_meteo["precip"].values, dtype=np.float32),
        dims=["index", "time"],
        coords=dict(
            index=["global"],
            time=meteo_times,
            x=("index", meteo_location[0].values),
            y=("index", meteo_location[1].values),
        ),
        attrs=dict(
            function="TimeSeries",
            timeInterpolation="Linear",
            quantity=f"{meteo_type}",
            units=f"{meteo_unit}",
            time_unit=f"{freq_name} since {pd.to_datetime(df_meteo.time[0])}",
            # support only yyyy-mm-dd HH:MM:SS
        ),
    )
    # fill in na using default
    da_out = da_out.fillna(fill_value)
    da_out.name = f"{meteo_type}"
    da_out.dropna(dim="time")

    return da_out


def _standardize_forcing_timeindexes(da):
    """Standardize timeindexes frequency based on forcing DataArray"""
    _TIMESTR = {"D": "days", "H": "hours", "T": "minutes", "S": "seconds"}
    dt = pd.to_timedelta((da.time[1].values - da.time[0].values))
    freq = dt.resolution_string
    multiplier = 1
    if freq == "D":
        logger.warning(
            "time unit days is not supported by the current GUI version: 2022.04"
        )  # converting to hours as temporary solution # FIXME: day is converted to hours temporarily
        multiplier = 24
    if len(pd.date_range(da.time[0].values, da.time[-1].values, freq=dt)) != len(
        da.time
    ):
        logger.error("does not support non-equidistant time-series.")
    freq_name = _TIMESTR[freq]
    freq_step = getattr(dt.components, freq_name)
    bd_times = np.array([float(i * freq_step) for i in range(len(da.time))])
    if multiplier == 24:
        bd_times = np.array([(i * freq_step * multiplier) for i in range(len(da.time))])
        freq_name = "hours"
    return bd_times, freq_name


def compute_forcing_values_points(
    gdf: gpd.GeoDataFrame,
    da: xr.DataArray = None,
    forcing_value: float = 0.0,
    forcing_type: str = "lateral_discharge",
    forcing_unit: str = "m3/s",
    logger=logger,
):
    """
    Compute 1d forcing values.

    Used for 1D lateral point locations.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame of points to add 1D forcing data.

        * Required variables: ['geometry']
    da : xr.DataArray, optional
        xr.DataArray containing the forcing timeseries values.
        If None, uses a constant ``forcing_value`` for all forcings.

        * Required variables: ['forcing_type']

    forcing_value : float, optional
        Constant value to use for all forcings if ``da`` is None and to
        fill in missing data.
        By default 0.0 ``forcing_unit``
    forcing_type : {'lateral_discharge'}
        Type of forcing to use.
        For now only support 'lateral_discharge'.
        By default 'lateral_discharge'
    forcing_unit : {'m3/s'}
        Unit corresponding to ``forcing_type``.
        By default 'm3/s'
    logger
        Logger to log messages.
    """
    # TODO: harmonize for other point forcing #21
    # first process data based on either timeseries or constant
    # then update data based on either nodes or branches
    # Timeseries forcing values
    if da is not None:
        logger.info(f"Preparing 1D forcing type {forcing_type} from timeseries.")

        # get forcing data freq in seconds
        bd_times, freq_name = _standardize_forcing_timeindexes(da)

        # instantiate xr.DataArray for forcing data
        # NOTE only support points on branches
        da_out = xr.DataArray(
            data=da.data,
            dims=["index", "time"],
            coords=dict(
                index=gdf.index,
                time=bd_times,
                x=("index", gdf.geometry.x.values),
                y=("index", gdf.geometry.y.values),
                branchid=("index", gdf.branchid.values),
                chainage=("index", gdf.chainage.values),
            ),
            attrs=dict(
                function="TimeSeries",
                timeInterpolation="Linear",
                quantity=f"{forcing_type}",
                units=f"{forcing_unit}",
                time_unit=f"{freq_name} since {pd.to_datetime(da.time[0].values)}",  # support only yyyy-mm-dd HH:MM:SS
            ),
        )
        # fill in na using default
        da_out = da_out.fillna(forcing_value)

        # drop na in time
        da_out.dropna(dim="time")

        # add name
        da_out.name = f"{forcing_type}"
    else:
        logger.info(
            f"Using constant value {forcing_value} {forcing_unit} for all {forcing_type} forcings."
        )
        # instantiate xr.DataArray for bnd data with forcing_type directly
        da_out = xr.DataArray(
            data=np.full((len(gdf.index)), forcing_value, dtype=np.float32),
            dims=["index"],
            coords=dict(
                index=gdf.index,
                x=("index", gdf.geometry.x.values),
                y=("index", gdf.geometry.y.values),
                branchid=("index", gdf.branchid.values),
                chainage=("index", gdf.chainage.values),
            ),
            attrs=dict(
                function="constant",
                offset=0.0,
                factor=1.0,
                quantity=f"{forcing_type}",
                units=f"{forcing_unit}",
            ),
        )
        da_out.name = f"{forcing_type}"
    return da_out


def get_geometry_coords_for_polygons(gdf):
    """Gets xarray DataArray coordinates that describes polygon geometries.
    Inlcudes numcoordinates, xcoordinates and ycoordinates"""
    if gdf.geometry.type.iloc[0] == "Polygon":
        # Get the maximum number of coordinates for any polygon
        max_coords = gdf["geometry"].apply(lambda x: len(x.exterior.coords[:])).max()

        def get_xcoords(geom):
            coords = [xy[0] for xy in geom.exterior.coords[:]]
            return np.pad(
                coords,
                (0, max_coords - len(coords)),
                "constant",
                constant_values=np.nan,
            )

        def get_ycoords(geom):
            coords = [xy[1] for xy in geom.exterior.coords[:]]
            return np.pad(
                coords,
                (0, max_coords - len(coords)),
                "constant",
                constant_values=np.nan,
            )

        # Create the 2D arrays
        x_2d = np.vstack(gdf["geometry"].apply(get_xcoords))
        y_2d = np.vstack(gdf["geometry"].apply(get_ycoords))

        return dict(
            index=gdf.index,
            numcoordinates=np.arange(max_coords),
            xcoordinates=(("index", "numcoordinates"), x_2d),
            ycoordinates=(("index", "numcoordinates"), y_2d),
        )


def compute_forcing_values_polygon(
    gdf: gpd.GeoDataFrame,
    da: xr.DataArray = None,
    forcing_value: float = 0.0,
    forcing_type: str = "waterlevelbnd",
    forcing_unit: str = "m",
    logger=logger,
):
    """
    Compute 1d forcing values.

    Used for 1D lateral polygon locations.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame of polygons to add 1D forcing data.

        * Required variables: ['geometry']
    da : xr.DataArray, optional
        xr.DataArray containing the forcing timeseries values.
        If None, uses a constant ``forcing_value`` for all forcings.

        * Required variables: ['forcing_type']

    forcing_value : float, optional
        Constant value to use for all forcings if ``da`` is None and to
        fill in missing data.
        By default 0.0 ``forcing_unit``
    forcing_type : {'lateral_discharge'}
        Type of forcing to use.
        For now only support 'lateral_discharge'.
        By default 'lateral_discharge'
    forcing_unit : {'m3/s'}
        Unit corresponding to ``forcing_type``.
        By default 'm3/s'
    logger
        Logger to log messages.
    """
    # Timeseries forcing values
    if da is not None:
        logger.info(f"Preparing 1D forcing type {forcing_type} from timeseries.")

        # get forcing data time indes
        bd_times, freq_name = _standardize_forcing_timeindexes(da)

        # instantiate xr.DataArray for forcing data
        coords_dict = get_geometry_coords_for_polygons(gdf)
        # Prepare the data
        data_3d = np.tile(
            np.expand_dims(da.data, axis=-1), (1, 1, len(coords_dict["numcoordinates"]))
        )
        # Create the DataArray
        da_out = xr.DataArray(
            data=data_3d,
            dims=("index", "time", "numcoordinates"),
            coords={
                "index": coords_dict["index"],
                "numcoordinates": coords_dict["numcoordinates"],
                "xcoordinates": coords_dict["xcoordinates"],
                "ycoordinates": coords_dict["ycoordinates"],
                "time": bd_times,
            },
            attrs=dict(
                function="TimeSeries",
                timeInterpolation="Linear",
                quantity=f"{forcing_type}",
                units=f"{forcing_unit}",
                time_unit=f"{freq_name} since {pd.to_datetime(da.time[0].values)}",  # support only yyyy-mm-dd HH:MM:SS
            ),
        )
        # fill in na using default
        da_out = da_out.fillna(forcing_value)

        # drop na in time
        da_out.dropna(dim="time")

        # add name
        da_out.name = f"{forcing_type}"
    else:
        logger.info(
            f"Using constant value {forcing_value} {forcing_unit} for all {forcing_type} forcings."
        )
        # instantiate xr.DataArray for forcing data with forcing_type directly
        coords_dict = get_geometry_coords_for_polygons(gdf)
        data_3d = np.full(
            (len(coords_dict["index"]), len(coords_dict["numcoordinates"])),
            forcing_value,
            dtype=np.float32,
        )
        da_out = xr.DataArray(
            data=data_3d,
            coords={
                "index": coords_dict["index"],
                "numcoordinates": coords_dict["numcoordinates"],
                "xcoordinates": coords_dict["xcoordinates"],
                "ycoordinates": coords_dict["ycoordinates"],
            },
            attrs=dict(
                function="constant",
                offset=0.0,
                factor=1.0,
                quantity=f"{forcing_type}",
                units=f"{forcing_unit}",
            ),
        )
        da_out.name = f"{forcing_type}"

    return da_out.drop_duplicates(dim=...)
