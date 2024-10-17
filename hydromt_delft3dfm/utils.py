"""Utilities read/write functions for Delft3D-FM model."""

from os.path import join
from pathlib import Path
from typing import Dict, List, Tuple, Union

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from hydrolib.core.dflowfm import (
    Boundary,
    BranchModel,
    CrossDefModel,
    CrossLocModel,
    ExtModel,
    FMModel,
    ForcingModel,
    FrictionModel,
    Lateral,
    Meteo,
    PolyFile,
    StorageNodeModel,
    StructureModel,
)
from shapely.geometry import Point, Polygon

from . import gis_utils
from .workflows import boundaries

__all__ = [
    "read_branches_gui",
    "write_branches_gui",
    "read_crosssections",
    "write_crosssections",
    "read_friction",
    "write_friction",
    "read_structures",
    "write_structures",
    "read_manholes",
    "write_manholes",
    "read_1dboundary",
    "write_1dboundary",
    "read_2dboundary",
    "write_2dboundary",
    "read_meteo",
    "write_meteo",
]


def read_branches_gui(
    gdf: gpd.GeoDataFrame,
    fm_model: FMModel,
) -> gpd.GeoDataFrame:
    """
    Read branches.gui and add the properties to branches geodataframe.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the branches
    fm_model: hydrolib.core FMModel
        DflowFM model object from hydrolib

    Returns
    -------
    gdf_out: geopandas.GeoDataFrame
        gdf containing the branches updated with branches gui params
    """
    branchgui_model = BranchModel()
    filepath = fm_model.filepath.with_name(
        branchgui_model._filename() + branchgui_model._ext()
    )

    if not filepath.is_file():
        # Create df with all attributes from nothing;
        # all branches are considered as river
        df_gui = pd.DataFrame()
        df_gui["branchid"] = [b.branchid for b in gdf.itertuples()]
        df_gui["branchtype"] = "river"
        df_gui["manhole_up"] = ""
        df_gui["manhole_dn"] = ""

    else:
        # Create df with all attributes from gui
        branchgui_model = BranchModel(filepath)
        br_list = branchgui_model.branch
        df_gui = pd.DataFrame()
        df_gui["branchid"] = [b.name for b in br_list]
        df_gui["branchtype"] = [b.branchtype for b in br_list]
        df_gui["manhole_up"] = [b.sourcecompartmentname for b in br_list]
        df_gui["manhole_dn"] = [b.targetcompartmentname for b in br_list]

    # Adapt type and add close attribute
    # TODO: channel and tunnel types are not defined
    df_gui["branchtype"] = df_gui["branchtype"].replace(
        {
            0: "river",
            2: "pipe",
            1: "sewerconnection",
        }
    )

    # Merge the two df based on branchid
    df_gui = df_gui.drop_duplicates(subset="branchid")
    gdf_out = gdf.merge(df_gui, on="branchid", how="left", validate=None)

    return gdf_out


def write_branches_gui(
    gdf: gpd.GeoDataFrame,
    savedir: str,
) -> str:
    """
    write branches.gui file from branches geodataframe.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the branches
    savedir: str
        path to the directory where to save the file.
    manholes: geopandas.GeoDataFrame, optional
        If exists add manholes attributes

    Returns
    -------
    branchgui_fn: str
        relative filepath to branches_gui file.

    #TODO: branches.gui is written with a [general] section which is not recongnised by
    GUI. Improvement of the GUI is needed.
    #TODO: branches.gui has a column is custumised length written as bool, which is not
    recongnised by GUI. improvement of the hydrolib-core writer is needed.
    """
    if not all([col in gdf.columns for col in ["manhole_up", "manhole_dn"]]):
        gdf[["manhole_up", "manhole_dn"]] = ""

    branches = gdf[["branchid", "branchtype", "manhole_up", "manhole_dn"]]
    branches = branches.rename(
        columns={
            "branchid": "name",
            "manhole_up": "sourceCompartmentName",
            "manhole_dn": "targetCompartmentName",
        }
    )
    branches["branchtype"] = branches["branchtype"].replace(
        {"river": 0, "pipe": 2, "sewerconnection": 1}
    )
    branchgui_model = BranchModel(branch=branches.to_dict("records"))
    branchgui_fn = branchgui_model._filename() + branchgui_model._ext()
    branchgui_model.filepath = join(savedir, branchgui_fn)
    branchgui_model.save()

    return branchgui_fn


def read_crosssections(gdf: gpd.GeoDataFrame, fm_model: FMModel) -> gpd.GeoDataFrame:
    """
    Read crosssections from hydrolib-core crsloc and crsdef objects and add to branches.

    Also returns crosssections geodataframe.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the branches
    fm_model: hydrolib.core FMModel
        DflowFM model object from hydrolib

    Returns
    -------
    gdf_crs: geopandas.GeoDataFrame
        geodataframe copy of the cross-sections and data
    """

    def _list2Str(lst):
        if isinstance(lst, list):
            # apply conversion to list columns
            if isinstance(lst[0], float):
                return " ".join(["{}".format(i) for i in lst])
            elif isinstance(lst[0], str):
                return " ".join(["{}".format(i) for i in lst])
        else:
            return lst

    # Start with crsdef to create the crosssections attributes
    crsdef = fm_model.geometry.crossdeffile
    crs_dict = dict()
    for b in crsdef.definition:
        crs_dict[b.id] = b.__dict__
    df_crsdef = pd.DataFrame.from_dict(crs_dict, orient="index")
    df_crsdef = df_crsdef.drop("comments", axis=1)
    df_crsdef["crs_id"] = df_crsdef["id"]  # column to merge
    # convertion needed  for xyz/zw crossections
    # convert list to str ()
    df_crsdef = df_crsdef.applymap(lambda x: _list2Str(x))
    # except for frictionids
    if "frictionids" in df_crsdef.columns:
        df_crsdef["frictionids"] = df_crsdef["frictionids"].str.replace(
            " ", ";"
        )  # comma list sperated
    # convert float to int
    int_columns = list(
        set(df_crsdef.columns).intersection(("xyzcount", "sectioncount"))
    )
    df_crsdef[int_columns] = (
        df_crsdef[int_columns]
        .fillna(-999)
        .astype(int)
        .astype(object)
        .where(df_crsdef[int_columns].notnull())
    )
    # Rename to prepare for crossection geom
    _gdf_crsdef = df_crsdef.rename(
        columns={c: f"crsdef_{c}" for c in df_crsdef.columns if c != "crs_id"}
    )

    # Continue with locs to get the locations and branches id
    crsloc = fm_model.geometry.crosslocfile
    crsloc_dict = dict()
    for b in crsloc.crosssection:
        crsloc_dict[b.id] = b.__dict__
    df_crsloc = pd.DataFrame.from_dict(crsloc_dict, orient="index")
    df_crsloc = df_crsloc.drop("comments", axis=1)
    df_crsloc["crs_id"] = df_crsloc["definitionid"]  # column to merge
    # get crsloc geometry
    if df_crsloc.dropna(axis=1).columns.isin(["x", "y"]).all():
        df_crsloc["geometry"] = [
            Point(i.x, i.y) for i in df_crsloc[["x", "y"]].itertuples()
        ]
    else:
        _gdf = gdf.set_index("branchid")
        df_crsloc["geometry"] = [
            _gdf.loc[i.branchid, "geometry"].interpolate(i.chainage)
            for i in df_crsloc.itertuples()
        ]
    # Rename to prepare for crossection geom
    _gdf_crsloc = df_crsloc.rename(
        columns={
            c: f"crsloc_{c}"
            for c in df_crsloc.columns
            if c not in ["crs_id", "geometry"]
        }
    )

    # Combine def attributes with locs for crossection geom
    gdf_crs = _gdf_crsloc.merge(_gdf_crsdef, on="crs_id", how="outer")
    # use outer because some crsdefs are from structures,
    # therefore no crslocs associated
    gdf_crs = gpd.GeoDataFrame(gdf_crs, crs=gdf.crs)

    return gdf_crs


def write_crosssections(gdf: gpd.GeoDataFrame, savedir: str) -> Tuple[str, str]:
    """Write crosssections into hydrolib-core crsloc and crsdef objects.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the crosssections
    savedir: str
        path to the directory where to save the file.

    Returns
    -------
    crsdef_fn: str
        relative filepath to crsdef file.
    crsloc_fn: str
        relative filepath to crsloc file.
    """
    # crsdef
    # get crsdef from crosssections gpd
    gpd_crsdef = gdf[[c for c in gdf.columns if c.startswith("crsdef")]]
    gpd_crsdef = gpd_crsdef.rename(
        columns={c: c.removeprefix("crsdef_") for c in gpd_crsdef.columns}
    )
    gpd_crsdef = gpd_crsdef.drop_duplicates(subset="id")
    gpd_crsdef = gpd_crsdef.astype(object).replace(np.nan, None)
    crsdef = CrossDefModel(definition=gpd_crsdef.to_dict("records"))
    # fm_model.geometry.crossdeffile = crsdef

    crsdef_fn = crsdef._filename() + ".ini"
    crsdef.save(
        join(savedir, crsdef_fn),
        recurse=False,
    )

    # crsloc
    # get crsloc from crosssections gpd
    gpd_crsloc = gdf[[c for c in gdf.columns if c.startswith("crsloc")]]
    gpd_crsloc = gpd_crsloc.rename(
        columns={c: c.removeprefix("crsloc_") for c in gpd_crsloc.columns}
    )
    gpd_crsloc = gpd_crsloc.dropna(
        subset="id"
    )  # structures have crsdefs but no crslocs

    crsloc = CrossLocModel(crosssection=gpd_crsloc.to_dict("records"))

    crsloc_fn = crsloc._filename() + ".ini"
    crsloc.save(
        join(savedir, crsloc_fn),
        recurse=False,
    )

    return crsdef_fn, crsloc_fn


def read_friction(gdf: gpd.GeoDataFrame, fm_model: FMModel) -> gpd.GeoDataFrame:
    """
    Read friction files and add properties to branches geodataframe.

    Assumes cross-sections have been read before to contain per branch frictionid.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the crosssections
    fm_model: hydrolib.core FMModel
        DflowFM model object from hydrolib

    Returns
    -------
    gdf_out: geopandas.GeoDataFrame
        gdf containing the crosssections updated with the friction params
    """
    fric_list = fm_model.geometry.frictfile
    # Create dictionnaries with all attributes from fricfile
    # global only
    fricval = dict()
    frictype = dict()
    for i in range(len(fric_list)):
        for j in range(len(fric_list[i].global_)):
            fricval[fric_list[i].global_[j].frictionid] = (
                fric_list[i].global_[j].frictionvalue
            )
            frictype[fric_list[i].global_[j].frictionid] = (
                fric_list[i].global_[j].frictiontype
            )

    # Create friction value and type by replacing frictionid values with dict
    gdf_out = gdf.copy()
    if "crsdef_frictionid" in gdf_out:
        gdf_out["frictionvalue"] = gdf_out["crsdef_frictionid"]
        gdf_out["frictiontype"] = gdf_out["crsdef_frictionid"]
    if "crsdef_frictionids" in gdf_out:
        _do_not_support = (
            gdf_out["crsdef_frictionids"].str.split(";").apply(np.count_nonzero) > 1
        )
        gdf_out.loc[~_do_not_support, "crsdef_frictionid"] = gdf_out.loc[
            ~_do_not_support, "crsdef_frictionids"
        ].combine_first(gdf_out.loc[~_do_not_support, "crsdef_frictionids"])
        gdf_out.loc[~_do_not_support, "frictionvalue"] = gdf_out.loc[
            ~_do_not_support, "frictionvalue"
        ].combine_first(gdf_out.loc[~_do_not_support, "crsdef_frictionids"])
        gdf_out.loc[~_do_not_support, "frictiontype"] = gdf_out.loc[
            ~_do_not_support, "frictiontype"
        ].combine_first(gdf_out.loc[~_do_not_support, "crsdef_frictionids"])
    gdf_out["frictionvalue"] = gdf_out["frictionvalue"].replace(fricval)
    gdf_out["frictiontype"] = gdf_out["frictiontype"].replace(frictype)
    return gdf_out


def write_friction(gdf: gpd.GeoDataFrame, savedir: str) -> List[str]:
    """
    write friction files from crosssections geodataframe.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the crosssections
    savedir: str
        path to the directory where to save the file.

    Returns
    -------
    friction_fns: List of str
        list of relative filepaths to friction files.
    """
    # Do not support segmented frictions
    if "crsdef_frictionids" in gdf.columns:
        _do_not_support = (
            gdf["crsdef_frictionids"].str.split(";").apply(np.count_nonzero) > 1
        )
        gdf = gdf.loc[~_do_not_support]
        gdf["crsdef_frictionid"] = gdf["crsdef_frictionid"].fillna(
            gdf["crsdef_frictionids"]
        )

    frictions = gdf.copy()
    frictions = frictions.rename(columns={"crsdef_frictionid": "frictionid"})[
        ["frictionid", "frictionvalue", "frictiontype"]
    ]
    frictions = frictions.drop_duplicates().dropna(how="all")

    friction_fns = []
    # create a new friction
    for i, row in frictions.iterrows():
        if isinstance(row.frictionvalue, float) and not np.isnan(row.frictionvalue):
            fric_model = FrictionModel(global_=row.to_dict())
            fric_name = f"{row.frictionid}"
            fric_filename = f"{fric_model._filename()}_{fric_name}" + fric_model._ext()
            fric_model.filepath = join(savedir, fric_filename)
            fric_model.save(fric_model.filepath, recurse=False)
            # save relative path to mdu
            friction_fns.append(fric_filename)

    return friction_fns


def read_structures(branches: gpd.GeoDataFrame, fm_model: FMModel) -> gpd.GeoDataFrame:
    """
    Read structures into hydrolib-core structures objects.

    Returns structures geodataframe.

    Will drop compound structures.

    Parameters
    ----------
    branches: geopandas.GeoDataFrame
        gdf containing the branches
    fm_model: hydrolib.core FMModel
        DflowFM model object from hydrolib

    Returns
    -------
    gdf_structures: geopandas.GeoDataFrame
        geodataframe of the structures and data
    """
    structures_fns = fm_model.geometry.structurefile
    # Parse to dict and DataFrame
    structures_dict = dict()
    for structures_fn in structures_fns:
        structures = structures_fn.structure
        for st in structures:
            structures_dict[st.id] = st.__dict__
    df_structures = pd.DataFrame.from_dict(structures_dict, orient="index")

    if len(df_structures) > 0:
        # Drop comments
        df_structures = df_structures.drop(
            ["comments"],
            axis=1,
        )
        # Drop compound structures (only write but do not read it back)
        df_structures = df_structures[df_structures["type"] != "compound"]

        # Add geometry
        gdf_structures = gis_utils.get_gdf_from_branches(branches, df_structures)
    else:
        df_structures = None
    return gdf_structures


def write_structures(gdf: gpd.GeoDataFrame, savedir: str) -> str:
    """
    write structures into hydrolib-core structures objects.

    Will add compound structures.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the structures
    savedir: str
        path to the directory where to save the file.

    Returns
    -------
    structures_fn: str
        relative path to structures file.
    """
    # Add compound structures
    cmp_structures = gdf.groupby(["chainage", "branchid"])["id"].apply(list)
    for cmp_count, cmp_st in enumerate(cmp_structures, start=1):
        gdf = pd.concat(
            [
                gdf,
                pd.DataFrame(
                    index=[max(gdf.index) + 1],
                    data={
                        "id": [f"CompoundStructure_{cmp_count}"],
                        "name": [f"CompoundStructure_{cmp_count}"],
                        "type": ["compound"],
                        "numStructures": [len(cmp_st)],
                        "structureIds": [";".join(cmp_st)],
                    },
                ),
            ],
            axis=0,
        )

    # replace nan with None
    gdf = gdf.replace(np.nan, None)

    # Write structures
    structures = StructureModel(structure=gdf.to_dict("records"))

    structures_fn = structures._filename() + ".ini"
    structures.save(
        join(savedir, structures_fn),
        recurse=False,
    )

    return structures_fn


def read_manholes(gdf: gpd.GeoDataFrame, fm_model: FMModel) -> gpd.GeoDataFrame:
    """
    Read manholes from hydrolib-core storagenodes and network 1d nodes for locations.

    Returns manholes geodataframe.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the network1d nodes
    fm_model: hydrolib.core FMModel
        DflowFM model object from hydrolib

    Returns
    -------
    gdf_manholes: geopandas.GeoDataFrame
        geodataframe of the manholes and data
    """
    manholes = fm_model.geometry.storagenodefile
    # Parse to dict and DataFrame
    manholes_dict = dict()
    for b in manholes.storagenode:
        manholes_dict[b.id] = b.__dict__
    df_manholes = pd.DataFrame.from_dict(manholes_dict, orient="index")
    # replace 0e to 0d
    # # TODO: fix this when hydrolib-core fix issue
    # https://github.com/Deltares/HYDROLIB-core/issues/559
    df_manholes["id"] = df_manholes["id"].apply(
        lambda x: "0D" + x[2:] if isinstance(x, str) and x.startswith("0e") else x
    )
    df_manholes["name"] = df_manholes["name"].apply(
        lambda x: "0D" + x[2:] if isinstance(x, str) and x.startswith("0e") else x
    )
    df_manholes["manholeid"] = df_manholes["manholeid"].apply(
        lambda x: "0D" + x[2:] if isinstance(x, str) and x.startswith("0e") else x
    )
    df_manholes["nodeid"] = df_manholes["nodeid"].apply(
        lambda x: "0D" + x[2:] if isinstance(x, str) and x.startswith("0e") else x
    )

    # Drop variables
    df_manholes = df_manholes.drop(
        ["comments"],
        axis=1,
    )

    gdf_manholes = gpd.GeoDataFrame(
        df_manholes.merge(gdf, on="nodeid", how="left"), crs=gdf.crs
    )

    return gdf_manholes


def write_manholes(gdf: gpd.GeoDataFrame, savedir: str) -> str:
    """
    write manholes into hydrolib-core storage nodes objects.

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        gdf containing the manholes
    savedir: str
        path to the directory where to save the file.

    Returns
    -------
    storage_fn: str
        relative path to storage nodes file.
    """
    if gdf.columns.__contains__("numlevels"):
        gdf["numlevels"] = gdf["numlevels"].astype("Int64")

    storagenodes = StorageNodeModel(storagenode=gdf.to_dict("records"))

    storage_fn = storagenodes._filename() + ".ini"
    storagenodes.save(
        join(savedir, storage_fn),
        recurse=False,
    )

    return storage_fn


def _read_forcing_dataframe(
    df_forcing: pd.DataFrame,
    index_values: Union[np.ndarray, str],
    quantity: str,
) -> Tuple[np.ndarray, List, Dict, Dict]:
    """
    Read forcing dataframe and parse to xarray properties.

    Parameters
    ----------
    df_forcing: pd.DataFrame
        Dataframe with forcing data.
    index_values: np.ndarray or str
        Index values of the forcing data.
    quantity: str
        Name of quantity (eg 'waterlevel').

    Returns
    -------
    data: np.ndarray
        Data of the forcing.
    dims: list
        Dimensions of the forcing.
    coords: dict
        Coordinates of the forcing.
    bc: dict
        Attributes of the forcing.
    """
    # Initialise dataarray attributes bc
    bc = {"quantity": quantity}

    # Get data from forcing df
    # Check if all constant
    if np.all(df_forcing.function == "constant"):
        # Prepare data
        data = np.array([v[0][0] for v in df_forcing.datablock])
        data = data + df_forcing.offset.values * df_forcing.factor.values
        # Prepare dataarray properties
        dims = ["index"]
        coords = dict(index=index_values)
        bc["function"] = "constant"
        bc["units"] = df_forcing.quantityunitpair.iloc[0][0].unit
        # bc["factor"] = 1
        # bc["offset"] = 0
    # Check if all timeseries
    elif np.all(df_forcing.function == "timeseries"):
        # Prepare data
        data = list()
        for i in np.arange(len(df_forcing.datablock)):
            v = df_forcing.datablock.iloc[i]
            offset = df_forcing.offset.iloc[i]
            factor = df_forcing.factor.iloc[i]
            databl = [n[1] * factor + offset for n in v]
            data.append(databl)
        data = np.array(data)
        # Assume unique times
        times = np.array([n[0] for n in df_forcing.datablock.iloc[0]])
        # Prepare dataarray properties
        dims = ["index", "time"]
        coords = dict(index=index_values, time=times)
        bc["function"] = "timeseries"
        bc["timeinterpolation"] = df_forcing.timeinterpolation.iloc[0]
        bc["units"] = df_forcing.quantityunitpair.iloc[0][1].unit
        bc["time_unit"] = df_forcing.quantityunitpair.iloc[0][0].unit
        bc["factor"] = 1
        bc["offset"] = 0
    # Else not implemented yet
    else:
        raise NotImplementedError(
            "ForcingFile with several function for a single variable not implemented."
            f"Skipping reading forcing for variable {quantity}."
        )
    return data, dims, coords, bc


def read_1dboundary(
    df: pd.DataFrame, quantity: str, nodes: gpd.GeoDataFrame
) -> xr.DataArray:
    """
    Read for a specific quantity the external and forcing files and parse to xarray.

    # TODO: support external forcing for 2D.

    Parameters
    ----------
    df: pd.DataFrame
        External Model DataFrame filtered for quantity.
    quantity: str
        Name of quantity (eg 'waterlevel').
    nodes: gpd.GeoDataFrame
        Nodes locations of the boundary in df.

    Returns
    -------
    da_out: xr.DataArray
        External and focing values combined into a DataArray for variable quantity.
    """
    nodeids = df.nodeid.values
    nodeids = nodeids[nodeids != "nan"]
    # Assume one forcing file (hydromt writer) and read
    forcing = df.forcingfile.iloc[0]
    df_forcing = pd.DataFrame([f.__dict__ for f in forcing.forcing])
    # Filter for the current nodes, remove nans
    df_forcing = df_forcing[np.isin(df_forcing.name, nodeids)]

    data, dims, coords, bc = _read_forcing_dataframe(
        df_forcing,
        index_values=nodeids,
        quantity=quantity,
    )

    # Get nodeid coordinates
    node_geoms = nodes.set_index("nodeid").reindex(nodeids)
    # # get rid of missing geometries
    xs, ys = np.vectorize(
        lambda p: (np.nan, np.nan) if p is None else (p.xy[0][0], p.xy[1][0])
    )(node_geoms["geometry"])
    coords["x"] = ("index", xs)
    coords["y"] = ("index", ys)

    # Prep DataArray and add to forcing
    da_out = xr.DataArray(
        data=data,
        dims=dims,
        coords=coords,
        attrs=bc,
    )
    da_out.name = f"boundary1d_{quantity}bnd"

    return da_out


def write_1dboundary(forcing: Dict, savedir: str = None, ext_fn: str = None) -> Tuple:
    """
    Write 1dboundary ext and boundary files from forcing dict.

    Parameters
    ----------
    forcing: dict of xarray DataArray
        Dict of boundary DataArray for each variable
        Only forcing that starts with "boundary1d" is recognised.
    savedir: str, optional
        path to the directory where to save the file.
    ext_fn: str or Path, optional
        Path of the external forcing file (.ext) in which this function will append to.
    """

    # remove duplicated name and keep the last
    def _remove_old_forcing_based_on_name(forcing):
        data_names = [k for k, v in forcing.items()]
        data = [v for k, v in forcing.items()]
        seen_names = set()
        filtered_datanames = []
        filtered_data = []
        # Reverse the list
        data_names = data_names[::-1]
        data = data[::-1]
        for d_name, d in zip(data_names, data):
            if d["index"].values[0] not in seen_names:
                seen_names.add(d["index"].values[0])
                filtered_data.append(d)
                filtered_datanames.append(d_name)

        # Reverse the filtered list
        filtered_data = filtered_data[::-1]
        filtered_datanames = filtered_datanames[::-1]
        return {k: v for k, v in zip(filtered_datanames, filtered_data)}

    # filter for 1d boundary
    forcing = {
        key: forcing[key] for key in forcing.keys() if key.startswith("boundary1d")
    }
    forcing = _remove_old_forcing_based_on_name(forcing)

    if len(forcing) == 0:
        return

    extdict = list()
    bcdict = list()
    # Loop over forcing dict
    for name, da in forcing.items():
        for i in da.index.values:
            bc = da.attrs.copy()
            # Boundary
            ext = dict()
            ext["quantity"] = bc["quantity"]
            ext["nodeid"] = i
            extdict.append(ext)
            # Forcing
            bc["name"] = i
            if bc["function"] == "constant":
                # one quantityunitpair
                bc["quantityunitpair"] = [
                    {"quantity": bc["quantity"], "unit": bc["units"]}
                ]
                # only one value column (no times)
                bc["datablock"] = [[da.sel(index=i).values.item()]]
            else:
                # two quantityunitpair
                bc["quantityunitpair"] = [
                    {"quantity": "time", "unit": bc["time_unit"]},
                    {"quantity": bc["quantity"], "unit": bc["units"]},
                ]
                bc.pop("time_unit")
                # time/value datablock
                bc["datablock"] = [
                    [t, x] for t, x in zip(da.time.values, da.sel(index=i).values)
                ]
            bc.pop("quantity")
            bc.pop("units")
            bcdict.append(bc)

    # write forcing file
    forcing_model = ForcingModel(forcing=bcdict)
    forcing_fn = "boundarycondition1d.bc"
    forcing_model.save(join(savedir, forcing_fn), recurse=True)

    # add forcingfile to ext, note each node needs a forcingfile
    extdicts = []
    for ext in extdict:
        ext["forcingfile"] = forcing_fn
        extdicts.append(ext)

    # write external forcing file
    if ext_fn is not None:
        # write to external forcing file
        write_ext(
            extdicts, savedir, ext_fn=ext_fn, block_name="boundary", mode="append"
        )

    return forcing_fn, ext_fn


def read_1dlateral(
    df: pd.DataFrame,
    quantity: str = "lateral_discharge",
    branches: gpd.GeoDataFrame = None,
) -> xr.DataArray:
    """
    Read a 1D lateral from external and forcing files.

    Parsed to xarray.

    Parameters
    ----------
    df: pd.DataFrame
        External Model DataFrame filtered for quantity.
    quantity: str
        Name of quantity. Supports only "lateral_discharge".
    nodes: gpd.GeoDataFrame
        Nodes locations of the laterals in df.
        # not implemented
    branches: gpd.GeoDataFrame
        Branches on which the laterals in df are located.

    Returns
    -------
    da_out: xr.DataArray
        External and focing values combined into a DataArray for variable quantity.
    """
    # Assume one discharge (lateral specific) file (hydromt writer) and read
    forcing = df.discharge.iloc[0]
    df_forcing = pd.DataFrame([f.__dict__ for f in forcing.forcing])

    data, dims, coords, bc = _read_forcing_dataframe(
        df_forcing,
        index_values=df_forcing.name,
        quantity=quantity,
    )

    # Get lateral locations and update dimentions and coordinates
    if any(df.numcoordinates.values):
        # polygons
        _df = df[~df.numcoordinates.isna()]
        _data = data[~df.numcoordinates.isna()]
        # Update coords
        _df["geometry"] = _df.apply(
            lambda row: Polygon(zip(row["xcoordinates"], row["ycoordinates"])), axis=1
        )
        coords_dict = boundaries.get_geometry_coords_for_polygons(gpd.GeoDataFrame(_df))
        dims.append("numcoordinates")
        coords.update(coords_dict)
        # Updates the data
        _data = np.tile(
            np.expand_dims(_data, axis=-1),
            (1, 1, len(coords_dict["numcoordinates"])),
        )

        # Prep DataArray and add to forcing
        da_out = xr.DataArray(
            data=_data,
            dims=dims,
            coords=coords,
            attrs=bc,
        )
        da_out.name = "lateral1d_polygons"
    else:
        # points
        if any(df.nodeid.values):
            # TODO laterals on nodes #78
            pass
        elif any(df.branchid.values):
            _df = df[~df.branchid.isna()]
            _data = data[~df.numcoordinates.isnull()]
            # update coords
            _df["geometry"] = [
                branches.set_index("branchid")
                .loc[i.branchid, "geometry"]
                .interpolate(i.chainage)
                for i in _df.itertuples()
            ]
            # Updates the data
            coords["x"] = ("index", np.array([p.x for p in _df["geometry"].values]))
            coords["y"] = ("index", np.array([p.y for p in _df["geometry"].values]))
            coords["branchid"] = ("index", _df["branchid"].values)
            coords["chainage"] = ("index", _df["chainage"].values)
            # Prep DataArray and add to forcing
            da_out = xr.DataArray(
                data=_data,
                dims=dims,
                coords=coords,
                attrs=bc,
            )
            da_out.name = "lateral1d_points"

    return da_out


def write_1dlateral(
    forcing: Dict, savedir: str = None, ext_fn: str = None
) -> Union[None, Tuple]:
    """
    Write 1dlateral ext and bc files from forcing dict.

    Parameters
    ----------
    forcing: dict of xarray DataArray
        Dict of lateral DataArray for each variable
        Only forcing that starts with "lateral1d" is recognised.
    savedir: str, optional
        path to the directory where to save the file.
    ext_fn: str or Path, optional
        Path of the external forcing file (.ext) in which this function will append to.
    """
    # filter for 1d lateral
    forcing = {
        key: forcing[key] for key in forcing.keys() if key.startswith("lateral1d")
    }
    if len(forcing) == 0:
        return

    extdict = list()
    bcdict = list()
    # Loop over forcing dict
    for name, da in forcing.items():
        for i in da.index.values:
            if da.attrs["quantity"] == "lateral_discharge":
                # Lateral
                latid = f"{name}_{i}"
                # Ext
                ext = dict()
                ext["id"] = latid
                ext["name"] = latid
                ext["locationType"] = "1d"
                if "nodeid" in da.coords:
                    # TODO laterals on nodes #78
                    # ext["nodeid"] = da.sel(index=i).coords["nodeid"].item()
                    pass
                elif "branchid" in da.coords:  # for point laterals
                    ext["branchid"] = da.sel(index=i).coords["branchid"].item()
                    ext["chainage"] = da.sel(index=i).coords["chainage"].item()
                elif "numcoordinates" in da.coords:  # for polygon laterals
                    xs = da.sel(index=i)["xcoordinates"].values.tolist()
                    ys = da.sel(index=i)["ycoordinates"].values.tolist()
                    ext["xcoordinates"] = [x for x in xs if not np.isnan(x)]
                    ext["ycoordinates"] = [y for y in ys if not np.isnan(y)]
                    ext["numcoordinates"] = len(ext["ycoordinates"])
                else:
                    raise ValueError("Not supported.")
                extdict.append(ext)
                # Forcing
                bc = da.attrs.copy()
                bc["name"] = latid
                if bc["function"] == "constant":
                    # one quantityunitpair
                    bc["quantityunitpair"] = [
                        {"quantity": da.name, "unit": bc["units"]}
                    ]
                    # only one value column (no times)
                    bc["datablock"] = [[da.sel(index=i).values.item()]]
                else:
                    # two quantityunitpair
                    bc["quantityunitpair"] = [
                        {"quantity": "time", "unit": bc["time_unit"]},
                        {"quantity": da.name, "unit": bc["units"]},
                    ]
                    bc.pop("time_unit")
                    # time/value datablock
                    _d = da.sel(index=i).values
                    if len(_d.shape) == 1:
                        # point
                        bc["datablock"] = [
                            [t, x]
                            for t, x in zip(da.time.values, da.sel(index=i).values)
                        ]
                    else:
                        # polygon
                        bc["datablock"] = [
                            [t, x]
                            for t, x in zip(
                                da.time.values,
                                np.unique(
                                    da.sel(index=i).values, axis=1
                                ),  # get the unique value to reduce polygon dimention
                            )
                        ]
                bc.pop("quantity")
                bc.pop("units")
                bcdict.append(bc)

    # write forcing file
    forcing_model = ForcingModel(forcing=bcdict)
    forcing_fn = "lateral1d.bc"
    forcing_model.save(join(savedir, forcing_fn), recurse=True)

    # add forcingfile to ext, note forcing file is called discharge for lateral
    extdicts = []
    for ext in extdict:
        ext["discharge"] = forcing_fn
        extdicts.append(ext)

    # write external forcing file
    if ext_fn is not None:
        # write to external forcing file
        write_ext(extdicts, savedir, ext_fn=ext_fn, block_name="lateral", mode="append")

    return forcing_fn, ext_fn


def _forcingmodel_to_dataarray(forcing: ForcingModel):
    """Read forcingfile timeseries into dataarray."""
    df_forcing = pd.DataFrame([f.__dict__ for f in forcing.forcing])

    # Get data
    # Check if all constant, Assume only timeseries exist (hydromt writer) and read
    bc = {}
    if np.all(df_forcing.function == "constant"):
        # Prepare data
        data = np.array([v[0][0] for v in df_forcing.datablock])
        data = data + df_forcing.offset.values * df_forcing.factor.values
        # Prepare dataarray properties
        dims = ["index"]
        coords = dict(index=list(df_forcing.index))
        bc["function"] = "constant"
        bc["units"] = df_forcing.quantityunitpair.iloc[0][0].unit
        bc["quantity"] = df_forcing.quantityunitpair.iloc[0][1].quantity
        bc["factor"] = 1
        bc["offset"] = 0
    # Check if all timeseries
    elif np.all(df_forcing.function == "timeseries"):
        # Prepare data
        data = list()
        for i in np.arange(len(df_forcing.datablock)):
            v = df_forcing.datablock.iloc[i]
            offset = df_forcing.offset.iloc[i]
            factor = df_forcing.factor.iloc[i]
            databl = [n[1] * factor + offset for n in v]
            data.append(databl)
        data = np.array(data)
        # Assume unique times
        times = np.array([n[0] for n in df_forcing.datablock.iloc[0]])
        # Prepare dataarray properties
        dims = ["index", "time"]
        coords = dict(index=list(df_forcing.index), time=times)
        bc["function"] = "timeseries"
        bc["timeinterpolation"] = df_forcing.timeinterpolation.iloc[0]
        bc["units"] = df_forcing.quantityunitpair.iloc[0][1].unit
        bc["time_unit"] = df_forcing.quantityunitpair.iloc[0][0].unit
        bc["quantity"] = df_forcing.quantityunitpair.iloc[0][1].quantity
        bc["factor"] = 1
        bc["offset"] = 0
    # Else not implemented yet
    else:
        raise NotImplementedError(
            "ForcingFile with several function for a single variable "
            "not implemented yet. Skipping reading forcing."
        )

    # Prep DataArray and add to forcing
    da_out = xr.DataArray(
        data=data,
        dims=dims,
        coords=coords,
        attrs=bc,
    )
    return da_out


def read_2dboundary(df: pd.DataFrame, workdir: Path = Path.cwd()) -> xr.DataArray:
    """
    Read a 2d boundary forcing location and values, and parse to xarray.

    Parameters
    ----------
    df: pd.DataFrame
        External Model DataFrame filtered for 2d boundary.
    workdir: Path
        working directory, i.e. where the files are stored.
        Required for reading disk only file model. #FIXME: there might be better options

    Returns
    -------
    da_out: xr.DataArray
        External and forcing values combined into a DataArray with name starts with
        "boundary2d".
    """
    quantity = df.quantity.iloc[0]

    # Initialise dataarray attributes

    # location file
    # assume one location file has only one location (hydromt writer) and read
    locationfile = PolyFile(workdir.joinpath(df.locationfile.filepath))
    locationfile.objects[0].metadata.name
    pd.DataFrame([f.__dict__ for f in locationfile.objects[0].points])

    # Assume one forcing file (hydromt writer) and read
    forcing = df.forcingfile
    df_forcing = pd.DataFrame([f.__dict__ for f in forcing.forcing])

    data, dims, coords, bc = _read_forcing_dataframe(
        df_forcing,
        index_values=df_forcing.name.values,
        quantity=df.quantity,
    )

    bc["locationfile"] = df.locationfile.filepath.name

    # Assume one forcing file (hydromt writer) and read
    forcing = df.forcingfile.iloc[0]
    _da_forcing = _forcingmodel_to_dataarray(forcing)

    # Add location
    # one location file has only one location
    _df = df.copy()
    for i, df_i in df.iterrows():
        locationfile = PolyFile(workdir.joinpath(df_i.locationfile.filepath))
        x = [f.__dict__["x"] for f in locationfile.objects[0].points]
        y = [f.__dict__["y"] for f in locationfile.objects[0].points]
        from shapely import LineString

        line = LineString(list(zip(x, y)))
        _df.loc[i, "geometry"] = line

    # Update dataarray
    # prior
    data = _da_forcing.data
    dims = list(_da_forcing.dims)
    coords = {k: v.values for k, v in _da_forcing.coords.items()}
    bc = _da_forcing.attrs.copy()
    # post
    coords_dict = boundaries.get_geometry_coords_for_linestrings(gpd.GeoDataFrame(_df))
    dims.append("numcoordinates")
    coords["numcoordinates"] = coords_dict["numcoordinates"]
    coords["x"] = coords_dict["xcoordinates"]
    coords["y"] = coords_dict["ycoordinates"]
    data = np.tile(
        np.expand_dims(data, axis=-1),
        (1, 1, len(coords_dict["numcoordinates"])),
    )

    # Prep DataArray and add to forcing
    da_out = xr.DataArray(
        data=data,
        dims=dims,
        coords=coords,
        attrs=bc,
    )
    da_out.name = f"{quantity}"

    return da_out


def _write_ncdicts(ncdicts: Dict[str, xr.DataArray], savedir: str):
    """
    Save forcing dictionaries of xr.DataArray to netCDF conventions.

    Parameters
    ----------
    ncdicts: Dict[str, xr.DataArray]
        key: str, the identifier of the boundary, must match with pli file
        value: xr.DataArray, The DataArray that has certain forcing quantity as data.

        Required variables: ["x", "y", "time"]
    savedir : str
        Directory where the netCDF file will be saved.

    Returns
    -------
    None, but writes a netCDF file to 'savedir'.
    """
    for ncdict in ncdicts:
        # use bndid as output name
        bndid = list(ncdict.keys())[0]
        da = ncdict.get(bndid)

        # Extract values from the DataArray
        x_values = da["x"].values
        y_values = da["y"].values
        node_values = range(1, 1 + len(x_values))
        time_values = da.time.values
        da_values = da.data
        # Create location labels from x and y pairs
        location_labels = [f"{bndid}_{str(i).zfill(4)}" for i in node_values]
        location_strlen = max(
            map(len, location_labels)
        )  # Determine the max string length
        location_values = np.array(
            location_labels, dtype=("U1", location_strlen)
        )  # convert location labels to structured numpy array

        # Construct the Dataset
        ds = xr.Dataset(
            {
                "location": (("node", "strlen"), location_values),
                da.attrs["quantity"]: (["time", "node"], da_values),
                "x": ("node", x_values),
                "y": ("node", y_values),
            },
            coords={
                "time": time_values,
                "node": node_values,
            },
        )

        # Set attributes for variables and coordinates
        ds.location.attrs["cf_role"] = "timeseries_id"
        ds.time.attrs["units"] = da.attrs["time_unit"]
        ds.x.attrs.update(
            dict(
                axis="x",
                units="m",
                long_name="x coordinate of projection",
                standard_name="projection_x_coordinate",
            )
        )
        ds.y.attrs.update(
            dict(
                axis="y",
                units="m",
                long_name="y coordinate of projection",
                standard_name="projection_y_coordinate",
            )
        )
        ds[da.attrs["quantity"]].attrs.update(
            {
                "long_name": da.attrs["quantity"],
                "standard_name": da.attrs[
                    "quantity"
                ],  # Assuming the quantity attribute can be used as standard_name
                "units": da.attrs["units"],
                "_FillValue": -999.9,
            }
        )

        # Save to netCDF
        ds.to_netcdf(join(savedir, f"{bndid}.nc"))


def _create_pliobj_from_xy(xs: list, ys: list, name: str):
    """Create hydrolib-core pli objecti from list of x and y coordinates."""
    xs = [x for x in xs if not np.isnan(x)]
    ys = [y for y in ys if not np.isnan(y)]
    _points = [{"x": x, "y": y, "data": []} for x, y in zip(xs, ys)]
    pli_object = {
        "metadata": {"name": name, "n_rows": len(xs), "n_columns": 2},
        "points": _points,
    }
    return pli_object


def write_2dboundary(forcing: Dict, savedir: str, ext_fn: str = None) -> list[dict]:
    """
    Write 2 boundary forcings from forcing dict.

    Note! forcing file (.bc) and forcing locations (.pli) are written in this function.
    Use external forcing (.ext) file will be extended.

    Parameters
    ----------
    forcing: dict of xarray DataArray
        Dict of boundary DataArray.
        Only forcing that startswith "boundary" will be recognised.
    savedir: str, optional
        path to the directory where to save the file.
    ext_fn: str or Path, optional
        Path of the external forcing file (.ext) in which this function will append to.

    """
    _bcfilename = "boundary2d.bc"

    # filter for 2d boundary
    forcing = {
        key: forcing[key] for key in forcing.keys() if key.startswith("boundary2d")
    }
    if len(forcing) == 0:
        return

    extdicts = list()
    bcdicts = list()
    ncdicts = list()
    plidicts = list()
    # Loop over forcing dict
    for _, da in forcing.items():
        for i in da.index.values:
            # initialise for boundary
            bndid = da.sel(index=i).index.item()

            # drop invalid xy locations
            da_i = da.sel(index=i)
            da_i = da_i.sel(numcoordinates=~np.isnan(da_i["x"]) & ~np.isnan(da_i["y"]))

            # drop invalid timeseries
            da_i_tsvalid = da_i.dropna(dim="numcoordinates", how="all")

            # Ext file
            ext = dict(
                quantity=da_i.attrs["quantity"],
                locationfile=bndid + ".pli",
                _bcfilename=_bcfilename,
                _ncfilename=bndid + ".nc",
            )
            extdicts.append(ext)

            # Pli file
            pli = _create_pliobj_from_xy(da_i["x"].values, da_i["y"].values, name=bndid)
            plidicts.append(pli)

            # Bc file
            # check if only one support point
            if da_i_tsvalid.sizes["numcoordinates"] > 1:
                raise NotImplementedError(
                    "Timeseries at multiple support points are not yet implemented."
                )
            else:
                # flaten data for bc file
                datablock = da_i_tsvalid.sel(numcoordinates=0).drop_vars(
                    "numcoordinates"
                )
                datablock.attrs["name"] = bndid + "_0001"
            bc = datablock.attrs.copy()
            # get quantityunitpair
            if bc["function"] == "constant":
                # one quantityunitpair
                bc["quantityunitpair"] = [
                    {"quantity": datablock.name, "unit": bc["units"]}
                ]
                # only one value column (no times)
                bc["datablock"] = [[datablock.values.item()]]
            else:
                # two quantityunitpair
                bc["quantityunitpair"] = [
                    {"quantity": "time", "unit": bc["time_unit"]},
                    {"quantity": datablock.name, "unit": bc["units"]},
                ]
                bc.pop("time_unit")
                # time/value datablock
                bc["datablock"] = [
                    [t, x] for t, x in zip(datablock.time.values, datablock.values)
                ]

            bc.pop("quantity")
            bc.pop("units")
            bcdicts.append(bc)

            # Ncfiles
            ncdicts.append({bndid: da_i_tsvalid})

    # write polyfile (allow single pli per external forcing block)
    for _pliitem in plidicts:
        pli_model = PolyFile(objects=[_pliitem])
        pli_fn = _pliitem["metadata"]["name"] + ".pli"
        pli_model.save(join(savedir, pli_fn), recurse=True)

    # write forcing file (as bc)
    forcing_model = ForcingModel(forcing=bcdicts)
    forcing_fn = _bcfilename
    forcing_model.save(join(savedir, forcing_fn), recurse=True)

    # write forcing file (as netcdf)
    _write_ncdicts(ncdicts, savedir=savedir)

    # add forcingfile and locationfile
    extdicts_ = []
    for ext in extdicts:
        ext["forcingfile"] = ext["_bcfilename"]
        # ext["forcingfile"] = ext["_ncfilename"] # switch after hydrolib-core support
        extdicts_.append(ext)

    # write external forcing file
    if ext_fn is not None:
        # write to external forcing file
        write_ext(
            extdicts_, savedir, ext_fn=ext_fn, block_name="boundary", mode="append"
        )

    return forcing_fn, ext_fn


def read_meteo(df: pd.DataFrame, quantity: str) -> xr.DataArray:
    """
    Read for a specific quantity the external and forcing files and parse to xarray.

    Parameters
    ----------
    df: pd.DataFrame
        External Model DataFrame filtered for quantity.
    quantity: str
        Name of quantity (e.g. "rainfall_rate", "rainfall").

    Returns
    -------
    da_out: xr.DataArray
        External and focing values combined into a DataArray for variable quantity.
    """
    # Assume one forcing file (hydromt writer) and read
    forcing = df.forcingfile.iloc[0]
    df_forcing = pd.DataFrame([f.__dict__ for f in forcing.forcing])
    # Filter for the current nodes
    df_forcing = df_forcing[np.isin(df_forcing.name, "global")]

    data, dims, coords, bc = _read_forcing_dataframe(
        df_forcing,
        index_values="global",
        quantity=quantity,
    )
    # Do not apply to "global" meteo
    # coords["x"]
    # coords["y"]

    # Prep DataArray and add to forcing
    da_out = xr.DataArray(
        data=data,
        dims=dims,
        coords=coords,
        attrs=bc,
    )
    da_out.name = f"{quantity}"

    return da_out


def write_meteo(forcing: Dict, savedir: str, ext_fn: str = None) -> list[dict]:
    """
    Write 2d meteo forcing from forcing dict.

    Note! only forcing file (.bc) is written in this function.
    Use utils.write_ext() for writing external forcing (.ext) file.

    Parameters
    ----------
    forcing: dict of xarray DataArray
        Dict of boundary DataArray.
        Only forcing that startswith "meteo" will be recognised.
    savedir: str, optional
        path to the directory where to save the file.
    ext_fn: str or Path, optional
        Path of the external forcing file (.ext) in which this function will append to.

    """
    # filter for 2d meteo
    forcing = {key: forcing[key] for key in forcing.keys() if key.startswith("meteo")}
    if len(forcing) == 0:
        return

    extdicts = list()
    bcdict = list()
    # Loop over forcing dict
    for name, da in forcing.items():
        for i in da.index.values:
            bc = da.attrs.copy()
            # Meteo
            ext = dict()
            ext["quantity"] = bc["quantity"]
            ext["forcingFileType"] = "bcAscii"
            # FIXME: hardcoded, decide whether use bcAscii or netcdf in setup
            # Forcing
            bc["name"] = i
            if bc["function"] == "constant":
                # one quantityunitpair
                bc["quantityunitpair"] = [{"quantity": da.name, "unit": bc["units"]}]
                # only one value column (no times)
                bc["datablock"] = [[da.sel(index=i).values.item()]]
            else:
                # two quantityunitpair
                bc["quantityunitpair"] = [
                    {"quantity": "time", "unit": bc["time_unit"]},
                    {"quantity": da.name, "unit": bc["units"]},
                ]
                bc.pop("time_unit")
                # time/value datablock
                bc["datablock"] = [
                    [t, x] for t, x in zip(da.time.values, da.sel(index=i).values)
                ]
            bc.pop("quantity")
            bc.pop("units")
            bcdict.append(bc)

    forcing_model = ForcingModel(forcing=bcdict)
    forcing_fn = f"meteo_{forcing_model._filename()}.bc"
    forcing_model.save(join(savedir, forcing_fn), recurse=True)

    # add forcingfile to ext
    ext["forcingfile"] = forcing_fn
    extdicts.append(ext)

    # write external forcing file
    if ext_fn is not None:
        # write to external forcing file
        write_ext(extdicts, savedir, ext_fn=ext_fn, block_name="meteo", mode="append")

    return forcing_fn, ext_fn


def write_ext(
    extdicts: List,
    savedir: Path,
    ext_fn: str = None,
    block_name: str = "boundary",
    mode="append",
) -> str:
    """
    write external forcing file (.ext) from dictionary.

    Parameters
    ----------
    extdicts: list[dict]
        list of dictionary containing attributes of external forcing blocks
    savedir: str, Path
        path to the save location.
    ext_fn: str, optional
        filename to the external forcing file.
    block_name: str, optional
        name of the block in the external forcing file. Includes "boundary", "lateral"
        and "meteo".
    mode: str, optional
        "overwrite" or "append".
        By default, append.

    """
    # FIXME: requires change of working directory for the validator to work properly
    import os

    cwd = os.getcwd()
    os.chdir(savedir)

    # ensure correct operation
    if Path(savedir).joinpath(ext_fn).is_file():
        if mode == "append":
            ext_model = ExtModel(Path(savedir).joinpath(ext_fn))
    else:
        ext_model = ExtModel()

    for i in range(len(extdicts)):
        if block_name == "boundary":
            ext_model.boundary.append(Boundary(**{**extdicts[i]}))
        elif block_name == "lateral":
            ext_model.lateral.append(Lateral(**extdicts[i]))
        elif block_name == "meteo":
            ext_model.meteo.append(Meteo(**{**extdicts[i]}))
        else:
            pass

    # Save ext file
    ext_model.save(ext_fn, recurse=True)

    # Go back to working dir
    os.chdir(cwd)

    return ext_fn
