"""Implement Delft3D FM HydroMT plugin model class."""

import itertools
import logging
import os
from datetime import datetime, timedelta
from os.path import basename, dirname, isfile, join
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import geopandas as gpd
import hydromt
import numpy as np
import pandas as pd
import xarray as xr
import xugrid as xu
from hydrolib.core.dflowfm import FMModel, IniFieldModel
from hydrolib.core.dimr import DIMR, FMComponent, Start
from hydromt.models import MeshModel
from hydromt.workflows import create_mesh2d
from pyproj import CRS
from shapely.geometry import box

from . import DATADIR, gis_utils, mesh_utils, utils, workflows

__all__ = ["DFlowFMModel"]
logger = logging.getLogger(__name__)


class DFlowFMModel(MeshModel):
    """API for Delft3D FM models in HydroMT."""

    _NAME = "dflowfm"
    _CONF = "DFlowFM.mdu"
    _DATADIR = DATADIR
    _GEOMS = {}
    _API = {
        "crs": CRS,
        "config": Dict[str, Any],
        "region": gpd.GeoDataFrame,
        "geoms": Dict[str, gpd.GeoDataFrame],
        "maps": Dict[str, Union[xr.DataArray, xr.Dataset]],
        "mesh": Union[xu.UgridDataArray, xu.UgridDataset],
        "forcing": Dict[str, Union[xr.DataArray, xr.Dataset]],
        "results": Dict[str, Union[xr.DataArray, xr.Dataset]],
        "states": Dict[str, Union[xr.DataArray, xr.Dataset]],
    }
    _MAPS = {
        "elevtn": {
            "name": "bedlevel",
            "initype": "initial",
            "interpolation": "triangulation",
            "locationtype": "2d",
        },
        "waterlevel": {
            "name": "waterlevel",
            "initype": "initial",
            "interpolation": "mean",
            "locationtype": "2d",
            "averagingrelsize": 1.01,  # default
        },
        "waterdepth": {
            "name": "waterdepth",
            "initype": "initial",
            "interpolation": "mean",
            "locationtype": "2d",
            "averagingrelsize": 1.01,
        },
        "pet": {
            "name": "PotentialEvaporation",
            "initype": "initial",
            "interpolation": "triangulation",
            "locationtype": "2d",
        },
        "infiltcap": {
            "name": "InfiltrationCapacity",
            "initype": "initial",
            "interpolation": "triangulation",
            "locationtype": "2d",
        },
        "roughness_chezy": {
            "name": "frictioncoefficient",
            "initype": "parameter",
            "interpolation": "triangulation",
            "locationtype": "2d",
            "frictype": 0,
        },
        "roughness_manning": {
            "name": "frictioncoefficient",
            "initype": "parameter",
            "interpolation": "triangulation",
            "locationtype": "2d",
            "frictype": 1,
        },
        "roughness_walllawnikuradse": {
            "name": "frictioncoefficient",
            "initype": "parameter",
            "interpolation": "triangulation",
            "locationtype": "2d",
            "frictype": 2,
        },
        "roughness_whitecolebrook": {
            "name": "frictioncoefficient",
            "initype": "parameter",
            "interpolation": "triangulation",
            "locationtype": "2d",
            "frictype": 3,
        },
    }
    _FOLDERS = ["dflowfm", "geoms", "maps"]
    _CLI_ARGS = {"region": "setup_region"}
    _CATALOGS = join(_DATADIR, "parameters_data.yml")

    def __init__(
        self,
        root: Union[str, Path],
        mode: str = "w",
        config_fn: str = None,
        data_libs: List[str] = [],  # yml
        crs: Union[int, str] = None,
        dimr_fn: str = None,
        network_snap_offset=25,
        snap_newbranches_to_branches_at_snapnodes=True,
        openwater_computation_node_distance=40,
        logger=logger,
    ):
        """Initialize the DFlowFMModel.

        Parameters
        ----------
        root : str or Path
            The model root location.
        mode : {'w','r','r+'}
            Write/read/append mode.
            Default is "w".
        config_fn : str, optional
            The D-Flow FM model configuration file (.mdu).
            If None, default configuration file is used.
            Default is None.
        data_libs : list of str, optional
            List of data catalog yaml files.
            Default is None.
        crs : EPSG code, int
            EPSG code of the model.
        dimr_fn: str, optional
            Path to the dimr configuration file.
            If None, default dimr configuration file is used.
            Default is None.
        network_snap_offset: float, optional
            Global option for generation of the mesh1d network. Snapping tolerance to
            automatically connecting branches.
            By default 25 m.
        snap_newbranches_to_branches_at_snapnodes: bool, optional
            Global option for generation of the mesh1d network.
            By default True.
        openwater_computation_node_distance: float, optional
            Global option for generation of the mesh1d network. Distance to generate
            mesh1d nodes for open water system (rivers, channels). By default 40 m.
        logger
            The logger used to log messages.
        """
        if not isinstance(root, (str, Path)):
            raise ValueError("The 'root' parameter should be a of str or Path.")

        super().__init__(
            root=root,
            mode=mode,
            config_fn=config_fn,
            data_libs=data_libs,
            logger=logger,
        )

        # model specific
        self._branches = None
        self._dimr = None
        self._dimr_fn = "dimr_config.xml" if dimr_fn is None else dimr_fn
        self._dfmmodel = None
        self._config_fn = (
            join("dflowfm", self._CONF) if config_fn is None else config_fn
        )  # FIXME Xiaohan config needs to be derived from dimr_fn if dimr_fn exsit
        self.data_catalog.from_yml(self._CATALOGS)

        self.config

        # Global options for generation of the mesh1d network
        self._network_snap_offset = network_snap_offset
        self._snap_newbranches_to_branches_at_snapnodes = (
            snap_newbranches_to_branches_at_snapnodes
        )
        self._openwater_computation_node_distance = openwater_computation_node_distance
        self._res = None

        # crs
        self._crs = CRS.from_user_input(crs) if crs else None
        self._check_crs()

    def setup_region(self, region):
        """HYDROMT CORE METHOD NOT USED FOR DFlowFMModel."""
        raise ValueError(
            "setup_region() method not implemented for DFlowFMModel."
            "The region will be set in the methods preparing the mesh: "
            "[setup_mesh2d, setup_rivers, setup_rivers_from_dem, setup_channels,"
            "setup_pipes]"
            "Pass the region argument to these methods directly and not in the"
            "command line."
        )

    def setup_channels(
        self,
        region: dict,
        channels_fn: str,
        channels_defaults_fn: str = "channels_defaults",
        channel_filter: str = None,
        friction_type: str = "Manning",
        friction_value: float = 0.023,
        crosssections_fn: str = None,
        crosssections_type: str = None,
        spacing: float = np.inf,
        snap_offset: float = 0.0,
        maxdist: float = 1.0,
        allow_intersection_snapping: bool = True,
    ):
        """Prepare the 1D channels and adds to branches 1D network.

        Adds model layers:

        * **channels** geom: 1D channels vector
        * **branches** geom: 1D branches vector

        Parameters
        ----------
        region : dict, optional
            Dictionary describing region of interest for extracting 1D channels, e.g.:

            * {'bbox': [xmin, ymin, xmax, ymax]}

            * {'geom': 'path/to/polygon_geometry'}
            Note that only crs=4326 is supported for 'bbox'.
        channels_fn : str
            Name of data source for channelsparameters, see data/data_sources.yml.
            Note only the lines that intersect with the region polygon will be used.

            * Optional variables: [branchid, branchtype, branchorder, material,
            friction_type, friction_value]
        channels_defaults_fn : str, optional
            Path to a csv file containing all defaults values per 'branchtype'.
            Default is None.
        channel_filter: str, optional
            Keyword in branchtype column of channels_fn used to filter river lines.
            If None all lines in channels_fn are used (default).
        friction_type : str, optional
            Type of friction to use. One of ["Manning", "Chezy", "wallLawNikuradse",
            "WhiteColebrook", "StricklerNikuradse", "Strickler", "deBosBijkerk"].
            By default "Manning".
        friction_value : float, optional.
            Friction value.
            Units corresponding to [friction_type] are ["Chézy C [m 1/2 /s]",
            "Manning n [s/m 1/3 ]", "Nikuradse k_n [m]", "Nikuradse k_n [m]",
            "Nikuradse k_n [m]", "Strickler k_s [m 1/3 /s]",
            "De Bos-Bijkerk γ [1/s]"]
            By default 0.023.
        crosssections_fn : str or Path, optional
            Name of data source for crosssections, see data/data_sources.yml.
            If ``crosssections_type`` = "xyzpoints"

            * Required variables: [crsid, order, z]
            If ``crosssections_type`` = "points"

            * Required variables: [crsid, order, z]
            By default None, crosssections will be set from branches
        crosssections_type : str, optional
            Type of crosssections read from crosssections_fn.
            One of ["xyzpoints", "point"].
            By default None.
        snap_offset: float, optional
            Snapping tolerance to automatically connecting branches.
            By default 0.0, no snapping is applied.
        maxdist: float, optional
            Maximum distance allowed for crosssections to be applied on branches.
            Only used for `crosssections_type` = point.
            By default 1.0.
        allow_intersection_snapping: bool, optional
            Switch to choose whether snapping of multiple branch ends are allowed when
            ``snap_offset`` is used.
            By default True.

        See Also
        --------
        dflowfm._setup_branches
        """
        self.logger.info("Preparing 1D channels.")

        # filter for allowed columns
        br_type = "channel"
        _allowed_columns = [
            "geometry",
            "branchid",
            "branchtype",
            "branchorder",
            "material",
            "shape",
            "diameter",
            "width",
            "t_width",
            "height",
            "bedlev",
            "closed",
            "friction_type",
            "friction_value",
        ]

        # Read data and filter within region
        region = workflows.parse_region_geometry(region, self.crs)
        gdf_br = self.data_catalog.get_geodataframe(
            channels_fn, geom=region, buffer=0, predicate="intersects"
        )
        # Read defaults table
        defaults = self.data_catalog.get_dataframe(channels_defaults_fn)
        # Prepare branches
        channels, channel_nodes = workflows.prepare_branches(
            gdf_br,
            params=defaults,
            br_type=br_type,
            dst_crs=self.crs,
            id_start=len(self.branches) + 1,
            filter=channel_filter,
            spacing=spacing,
            snap_offset=snap_offset,
            allow_intersection_snapping=allow_intersection_snapping,
            allowed_columns=_allowed_columns,
            logger=self.logger,
        )
        # Prepare friction and crosssections
        channels = workflows.prepare_default_friction_and_crosssection(
            channels,
            br_type=br_type,
            friction_type=friction_type,
            friction_value=friction_value,
            logger=self.logger,
        )

        # setup crosssections
        if crosssections_type is None:
            crosssections_type = "branch"
            # TODO: maybe assign a specific one for river, like branch_river
        assert {crosssections_type}.issubset({"xyzpoints", "branch", "point"})
        crosssections = self._setup_crosssections(
            branches=channels,
            region=region,
            crosssections_fn=crosssections_fn,
            crosssections_type=crosssections_type,
            maxdist=maxdist,
        )

        # add crosssections to exisiting ones and update geoms
        self.logger.debug("Adding crosssections vector to geoms.")
        crosssections = workflows.add_crosssections(
            self.geoms.get("crosssections"), crosssections
        )
        self.set_geoms(crosssections, "crosssections")

        # setup geoms
        self.logger.debug("Adding branches and branch_nodes vector to geoms.")
        self.set_geoms(channels, "channels")
        self.set_geoms(channel_nodes, "channel_nodes")

        # add to branches geoms
        branches = workflows.add_branches(
            self.mesh_datasets.get("mesh1d"),
            self.branches,
            channels,
            self._snap_newbranches_to_branches_at_snapnodes,
            self._network_snap_offset,
        )
        workflows.validate_branches(branches)
        self.set_branches(branches)

        # update mesh
        # Convert self.branches to xugrid
        mesh1d, network1d = workflows.mesh1d_network1d_from_branches(
            self.opensystem,
            self.closedsystem,
            self._openwater_computation_node_distance,
        )
        self.set_mesh(network1d, grid_name="network1d", overwrite_grid=True)
        self.set_mesh(mesh1d, grid_name="mesh1d", overwrite_grid=True)

    def setup_rivers_from_dem(
        self,
        region: dict,
        hydrography_fn: str,
        river_geom_fn: str = None,
        rivers_defaults_fn: str = "rivers_defaults",
        rivdph_method="gvf",
        rivwth_method="geom",
        river_upa=25.0,
        river_len=1000,
        min_rivwth=50.0,
        min_rivdph=1.0,
        rivbank=True,
        rivbankq=25,
        segment_length=3e3,
        smooth_length=10e3,
        friction_type: str = "Manning",
        friction_value: float = 0.023,
        constrain_rivbed=True,
        constrain_estuary=True,
        **kwargs,  # for workflows.get_river_bathymetry method
    ) -> None:
        """
        Set the all river parameters from hydrograph and dem maps.

        River cells are based on the `river_mask_fn` raster file if
        `rivwth_method='mask'`, or if `rivwth_method='geom'` the rasterized segments
        buffered with half a river width ("rivwth" [m]) if that attribute is found in
        `river_geom_fn`.

        If a river segment geometry file `river_geom_fn` with bedlevel column
        ("zb" [m+REF]) or a river depth ("rivdph" [m]) in combination with
        `rivdph_method='geom'` is provided, this attribute is used directly.

        Otherwise, a river depth is estimated based on bankfull discharge
        ("qbankfull" [m3/s]) attribute taken from the nearest river segment in
        `river_geom_fn` or `qbankfull_fn` upstream river boundary points if provided.

        The river depth is relative to the bankfull elevation profile if `rivbank=True`
        (default), which is estimated as the `rivbankq` elevation percentile [0-100] of
        cells neighboring river cells. This option requires the flow direction
        ("flwdir") and upstream area ("uparea") maps to be set using the
        hydromt.flw.flwdir_from_da method. If `rivbank=False` the depth is simply
        subtracted from the elevation of river cells.

        Missing river width and river depth values are filled by propagating valid
        values downstream and using the constant minimum values `min_rivwth` and
        `min_rivdph` for the remaining missing values.

        Updates model layer:

        * **dep** map: combined elevation/bathymetry [m+ref]

        Adds model layers

        * **rivmsk** map: map of river cells (not used by DFlowFM)
        * **rivers** geom: geometry of rivers (not used by DFlowFM)

        Parameters
        ----------
        region : dict, optional
            Dictionary describing region of interest for extracting 1D rivers, e.g.:

            * {'bbox': [xmin, ymin, xmax, ymax]}

            * {'geom': 'path/to/polygon_geometry'}
        hydrography_fn : str
            Hydrography data to derive river shape and characteristics from.

            * Required variables: ['elevtn']
            * Optional variables: ['flwdir', 'uparea']
        river_geom_fn : str, optional
            Line geometry with river attribute data.

            * Required variable for direct bed level burning: ['zb']
            * Required variable for direct river depth burning: ['rivdph']
            (only in combination with rivdph_method='geom')
            * Variables used for river depth estimates: ['qbankfull', 'rivwth']
        rivers_defaults_fn : str Path
            Path to a csv file containing all defaults values per 'branchtype'.
        rivdph_method : {'gvf', 'manning', 'powlaw'}
            River depth estimate method, by default 'gvf'
        rivwth_method : {'geom', 'mask'}
            Derive the river with from either the `river_geom_fn` (geom) or
            `river_mask_fn` (mask; default) data.
        river_upa : float, optional
            Minimum upstream area threshold for rivers [km2], by default 25.0
        river_len: float, optional
            Mimimum river length within the model domain threshhold [m], by default
            1000 m.
        min_rivwth, min_rivdph: float, optional
            Minimum river width [m] (by default 50.0) and depth [m] (by default 1.0)
        rivbank: bool, optional
            If True (default), approximate the reference elevation for the river depth
            based on the river bankfull elevation at cells neighboring river cells.
            Otherwise use the elevation of the local river cell as reference level.
        rivbankq : float, optional
            quantile [1-100] for river bank estimation, by default 25
        segment_length : float, optional
            Approximate river segment length [m], by default 5e3
        smooth_length : float, optional
            Approximate smoothing length [m], by default 10e3
        friction_type : str, optional
            Type of friction tu use. One of ["Manning", "Chezy", "wallLawNikuradse",
            "WhiteColebrook", "StricklerNikuradse", "Strickler", "deBosBijkerk"].
            By default "Manning".
        friction_value : float, optional.
            Units corresponding to [friction_type] are ["Chézy C [m 1/2 /s]",
            "Manning n [s/m 1/3 ]", "Nikuradse k_n [m]", "Nikuradse k_n [m]",
            "Nikuradse k_n [m]", "Strickler k_s [m 1/3 /s]", "De Bos-Bijkerk γ [1/s]"]
            Friction value. By default 0.023.
        constrain_estuary : bool, optional
            If True (default) fix the river depth in estuaries based on the upstream
            river depth.
        constrain_rivbed : bool, optional
            If True (default) correct the river bed level to be hydrologically correct,
            i.e. sloping downward in downstream direction.

        See Also
        --------
        workflows.get_river_bathymetry

        Raises
        ------
        ValueError

        """
        self.logger.info("Preparing river shape from hydrography data.")
        # parse region argument
        region = workflows.parse_region_geometry(region, self.crs)

        # read data
        ds_hydro = self.data_catalog.get_rasterdataset(
            hydrography_fn,
            geom=region,
        )
        if isinstance(ds_hydro, xr.DataArray):
            ds_hydro = ds_hydro.to_dataset()

        # read river line geometry data
        gdf_riv = None
        if river_geom_fn is not None:
            gdf_riv = self.data_catalog.get_geodataframe(
                river_geom_fn, geom=region, buffer=1
            ).to_crs(ds_hydro.raster.crs)

        # check if flwdir and uparea in ds_hydro
        if "flwdir" not in ds_hydro.data_vars:
            da_flw = hydromt.flw.d8_from_dem(ds_hydro["elevtn"])
        else:
            da_flw = ds_hydro["flwdir"]
        flwdir = hydromt.flw.flwdir_from_da(da_flw, ftype="d8")
        if "uparea" not in ds_hydro.data_vars:
            da_upa = xr.DataArray(
                dims=ds_hydro["elevtn"].raster.dims,
                coords=ds_hydro["elevtn"].raster.coords,
                data=flwdir.upstream_area(unit="km2"),
                name="uparea",
            )
            da_upa.raster.set_nodata(-9999)
            ds_hydro["uparea"] = da_upa

        # get river shape and bathymetry
        if friction_type == "Manning":
            kwargs.update(manning=friction_value)
        elif rivdph_method == "gvf":
            raise ValueError(
                "rivdph_method 'gvf' requires friction_type='Manning'."
                "Use 'geom' or 'powlaw' instead."
            )
        gdf_riv, _ = workflows.get_river_bathymetry(
            ds_hydro,
            flwdir=flwdir,
            gdf_riv=gdf_riv,
            gdf_qbf=None,
            rivdph_method=rivdph_method,
            rivwth_method=rivwth_method,
            river_upa=river_upa,
            river_len=river_len,
            min_rivdph=min_rivdph,
            min_rivwth=min_rivwth,
            rivbank=rivbank,
            rivbankq=rivbankq,
            segment_length=segment_length,
            smooth_length=smooth_length,
            constrain_estuary=constrain_estuary,
            constrain_rivbed=constrain_rivbed,
            logger=self.logger,
            **kwargs,
        )
        # Rename river properties column and reproject
        rm_dict = {"rivwth": "width", "rivdph": "height", "zb": "bedlev"}
        gdf_riv = gdf_riv.rename(columns=rm_dict).to_crs(self.crs)

        # filter for allowed columns
        br_type = "river"
        _allowed_columns = [
            "geometry",
            "branchid",
            "branchtype",
            "branchorder",
            "material",
            "shape",
            "diameter",
            "width",
            "t_width",
            "height",
            "bedlev",
            "closed",
            "friction_type",
            "friction_value",
        ]

        # assign default attributes
        # Read defaults table
        defaults = self.data_catalog.get_dataframe(rivers_defaults_fn)
        # Make sure default shape is rectangle
        defaults["shape"] = np.array(["rectangle"])

        # Prepare branches
        rivers, river_nodes = workflows.prepare_branches(
            gdf_riv,
            params=defaults,
            br_type=br_type,
            dst_crs=self.crs,
            id_start=len(self.branches) + 1,
            allowed_columns=_allowed_columns,
            logger=self.logger,
        )
        # Prepare friction
        branches = workflows.prepare_default_friction_and_crosssection(
            rivers,
            br_type=br_type,
            friction_type=friction_type,
            friction_value=friction_value,
            logger=self.logger,
        )

        # setup crosssections
        crosssections = self._setup_crosssections(
            branches=rivers,
            crosssections_fn=None,
            crosssections_type="branch",
        )

        # add crosssections to exisiting ones and update geoms
        self.logger.debug("Adding crosssections vector to geoms.")
        crosssections = workflows.add_crosssections(
            self.geoms.get("crosssections"), crosssections
        )
        self.set_geoms(crosssections, "crosssections")

        # setup geoms #TODO do we still need channels?
        self.logger.debug("Adding rivers and river_nodes vector to geoms.")
        self.set_geoms(rivers, "rivers")
        self.set_geoms(river_nodes, "rivers_nodes")

        # add to branches geoms
        branches = workflows.add_branches(
            self.mesh_datasets.get("mesh1d"),
            self.branches,
            rivers,
            self._snap_newbranches_to_branches_at_snapnodes,
            self._network_snap_offset,
        )
        workflows.validate_branches(branches)
        self.set_branches(branches)

        # update mesh
        # Convert self.branches to xugrid
        mesh1d, network1d = workflows.mesh1d_network1d_from_branches(
            self.opensystem,
            self.closedsystem,
            self._openwater_computation_node_distance,
        )
        self.set_mesh(network1d, grid_name="network1d", overwrite_grid=True)
        self.set_mesh(mesh1d, grid_name="mesh1d", overwrite_grid=True)

    def setup_rivers(
        self,
        region: dict,
        rivers_fn: str,
        rivers_defaults_fn: str = "rivers_defaults",
        river_filter: str = None,
        friction_type: str = "Manning",
        friction_value: float = 0.023,
        crosssections_fn: Union[int, list] = None,
        crosssections_type: Union[int, list] = None,
        snap_offset: float = 0.0,
        maxdist: float = 1.0,
        allow_intersection_snapping: bool = True,
    ):
        """Prepare the 1D rivers and adds to 1D branches.

        1D rivers must contain valid geometry, friction and crosssections.

        The river geometry is read from ``rivers_fn``. If defaults attributes
        [branchorder, spacing, material, shape, width, t_width, height, bedlev, closed]
        are not present in ``rivers_fn``, they are added from defaults values in
        ``rivers_defaults_fn``. For branchid and branchtype, they are created on the fly
        if not available in rivers_fn ("river" for Type and "river_{i}" for Id).

        Friction attributes are either taken from ``rivers_fn`` or filled in using
        ``friction_type`` and ``friction_value`` arguments.
        Note for now only branch friction or global friction is supported.

        Crosssections are read from ``crosssections_fn`` based on the
        ``crosssections_type``. If there is no ``crosssections_fn`` values are derived
        at the centroid of each river line based on defaults. If there are multiple
        types of crossections, specify them as lists.

        Adds/Updates model layers:

        * **rivers** geom: 1D rivers vector
        * **branches** geom: 1D branches vector
        * **crosssections** geom: 1D crosssection vector

        Parameters
        ----------
        region : dict, optional
            Dictionary describing region of interest for extracting 1D rivers, e.g.:

            * {'bbox': [xmin, ymin, xmax, ymax]}

            * {'geom': 'path/to/polygon_geometry'}
        rivers_fn : str
            Name of data source for rivers parameters, see data/data_sources.yml.
            Note only the lines that intersect with the region polygon will be used.

            * Optional variables: [branchid, branchtype, branchorder, material,
            friction_type, friction_value]
        rivers_defaults_fn : str Path
            Path to a csv file containing all defaults values per 'branchtype'.
            By default None.
        river_filter: str, optional
            Keyword in branchtype column of rivers_fn used to filter river lines.
            If None all lines in rivers_fn are used (default).
        friction_type : str, optional
            Type of friction to use. One of ["Manning", "Chezy", "wallLawNikuradse",
            "WhiteColebrook", "StricklerNikuradse", "Strickler", "deBosBijkerk"].
            By default "Manning".
        friction_value : float, optional.
            Units corresponding to [friction_type] are ["Chézy C [m 1/2 /s]",
            "Manning n [s/m 1/3 ]", "Nikuradse k_n [m]", "Nikuradse k_n [m]",
            "Nikuradse k_n [m]", "Strickler k_s [m 1/3 /s]", "De Bos-Bijkerk γ [1/s]"]
            Friction value. By default 0.023.
        crosssections_fn : str, Path, or a list of str or Path, optional
            Name of data source for crosssections, see data/data_sources.yml. One or a
            list corresponding to ``crosssections_type`` .
            If ``crosssections_type`` = "branch"
            crosssections_fn should be None
            If ``crosssections_type`` = "xyz"

            * Required variables: [crsid, order, z]
            If ``crosssections_type`` = "point"

            * Required variables: [crsid, shape, shift]
            By default None, crosssections will be set from branches
        crosssections_type : str, or a list of str, optional
            Type of crosssections read from crosssections_fn. One or a list of
            ["branch", "xyz", "point"].
            By default None.
        snap_offset: float, optional
            Snapping tolerance to automatically connecting branches.
            By default 0.0, no snapping is applied.
        maxdist: float, optional
            Maximum distance allowed for crosssections to be applied on branches.
            Only used for `crosssections_type` = point.
            By default 1.0.
        allow_intersection_snapping: bool, optional
            Switch to choose whether snapping of multiple branch ends are allowed when
            ``snap_offset`` is used.
            By default True.

        See Also
        --------
        dflowfm._setup_branches
        dflowfm._setup_crosssections
        """
        self.logger.info("Preparing 1D rivers.")
        # filter for allowed columns
        br_type = "river"
        _allowed_columns = [
            "geometry",
            "branchid",
            "branchtype",
            "branchorder",
            "material",
            "shape",
            "width",
            "t_width",
            "height",
            "bedlev",
            "closed",
            "friction_type",
            "friction_value",
        ]

        # Read data and filter within region
        region = workflows.parse_region_geometry(region, self.crs)
        gdf_br = self.data_catalog.get_geodataframe(
            rivers_fn, geom=region, buffer=0, predicate="intersects"
        )
        # Read defaults table
        defaults = self.data_catalog.get_dataframe(rivers_defaults_fn)
        # Prepare branches
        rivers, river_nodes = workflows.prepare_branches(
            gdf_br,
            params=defaults,
            br_type=br_type,
            dst_crs=self.crs,
            id_start=len(self.branches) + 1,
            filter=river_filter,
            snap_offset=snap_offset,
            allow_intersection_snapping=allow_intersection_snapping,
            allowed_columns=_allowed_columns,
            logger=self.logger,
        )
        # Prepare friction and crosssections
        rivers = workflows.prepare_default_friction_and_crosssection(
            rivers,
            br_type=br_type,
            friction_type=friction_type,
            friction_value=friction_value,
            logger=self.logger,
        )

        # setup crosssections
        crosssections_type, crosssections_fn = workflows.init_crosssections_options(
            crosssections_type, crosssections_fn
        )
        for crs_fn, crs_type in zip(crosssections_fn, crosssections_type):
            crosssections = self._setup_crosssections(
                branches=rivers,
                region=region,
                crosssections_fn=crs_fn,
                crosssections_type=crs_type,
                maxdist=maxdist,
            )
            crosssections = workflows.add_crosssections(
                self.geoms.get("crosssections"), crosssections
            )
            # setup geoms for crosssections
            self.set_geoms(crosssections, "crosssections")

        # setup branch orders
        # for crossection type yz or xyz, always use branchorder = -1,
        # because no interpolation can be applied.
        # TODO: change to lower case is needed
        _overwrite_branchorder = self.geoms["crosssections"][
            self.geoms["crosssections"]["crsdef_type"].str.contains("yz")
        ]["crsdef_branchid"].tolist()
        if len(_overwrite_branchorder) > 0:
            rivers.loc[
                rivers["branchid"].isin(_overwrite_branchorder), "branchorder"
            ] = -1

        # setup geoms for rivers and river_nodes
        self.logger.debug("Adding rivers and river_nodes vector to geoms.")
        self.set_geoms(rivers, "rivers")
        self.set_geoms(river_nodes, "rivers_nodes")

        # setup branches
        branches = workflows.add_branches(
            self.mesh_datasets.get(
                "mesh1d"
            ),  # FIXME Xiaohan self.get_mesh("mesh1d") gives an error, is that desired?
            self.branches,
            rivers,
            self._snap_newbranches_to_branches_at_snapnodes,
            self._network_snap_offset,
        )
        workflows.validate_branches(branches)
        self.set_branches(branches)

        # update mesh
        # Convert self.branches to xugrid
        mesh1d, network1d = workflows.mesh1d_network1d_from_branches(
            self.opensystem,
            self.closedsystem,
            self._openwater_computation_node_distance,
        )
        self.set_mesh(network1d, grid_name="network1d", overwrite_grid=True)
        self.set_mesh(mesh1d, grid_name="mesh1d", overwrite_grid=True)

    def setup_pipes(
        self,
        region: dict,
        pipes_fn: str,
        pipes_defaults_fn: str = "pipes_defaults",
        pipe_filter: Union[str, None] = None,
        spacing: float = np.inf,
        friction_type: str = "WhiteColebrook",
        friction_value: float = 0.003,
        crosssections_shape: str = "circle",
        crosssections_value: Union[int, list] = 0.5,
        dem_fn: Union[str, None] = None,
        pipes_depth: float = 2.0,
        pipes_invlev: float = -2.5,
        snap_offset: float = 0.0,
        allow_intersection_snapping: bool = True,
    ):
        """Prepare the 1D pipes and adds to 1D branches.

        Note that 1D manholes must also be set-up  when setting up 1D pipes.

        1D pipes must contain valid geometry, friction and crosssections.

        The pipe geometry is read from ``pipes_fn``.

        If branchtype is present in ``pipes_fn``, it is possible to filter pipe geometry
        using an additional filter specificed in``pipe_filter``. If defaults attributes
        ["branchorder"] are not present in ``pipes_fn``, they are added from defaults
        values in ``pipes_defaults_fn``. For branchid and branchtype, if they are not
        present in ``pipes_fn``, they are created on the fly ("pipe" for branchtype and
        "pipe_{i}" for branchid).

        The pipe geometry can be processed using splitting based on ``spacing``.

        Friction attributes ["branchid", "frictionvalue"] are either taken from
        ``pipes_fn`` or filled in using ``friction_type`` and ``friction_value``
        arguments.
        Note for now only branch friction or global friction is supported.

        Crosssections definition attributes ["shape", "diameter", "width", "height",
        "closed"] are either taken from ``pipes_fn`` or filled in using
        ``crosssections_shape`` and ``crosssections_value``.

        Crosssections location attributes ["invlev_up", "invlev_dn"] are either taken
        from ``pipes_fn``, or derived from ``dem_fn`` minus a fixed depth ``pipe_depth``
        [m], or from a constant ``pipe_invlev`` [m asl] (not recommended! should be
        edited before a model run).

        Adds/Updates model layers:

        * **pipes** geom: 1D pipes vector
        * **branches** geom: 1D branches vector
        * **crosssections** geom: 1D crosssection vector

        Parameters
        ----------
        region : dict, optional
            Dictionary describing region of interest for extracting 1D pipes, e.g.:

            * {'bbox': [xmin, ymin, xmax, ymax]}

            * {'geom': 'path/to/polygon_geometry'}
            Note that only manholes within the region are kept.
        pipes_fn : str
            Name of data source for pipes parameters, see data/data_sources.yml.
            Note only the lines that are within the region polygon will be used.

            * Optional variables: [branchid, branchtype, branchorder, spacing, branchid,
            frictionvalue, shape, diameter, width, height, closed, invlev_up, invlev_dn]
            #TODO: material table is used for friction which is not implemented
        pipes_defaults_fn : str Path
            Path to a csv file containing all defaults values per "branchtype"'".
        pipe_filter: str, optional
            Keyword in branchtype column of pipes_fn used to filter pipe lines.
            If None all lines in pipes_fn are used (default).
        spacing: float, optional
            Spacing value in meters to split the long pipelines lines into shorter
            pipes. By default inf - no splitting is applied.
        friction_type : str, optional
            Type of friction to use. One of ["Manning", "Chezy", "wallLawNikuradse",
            "WhiteColebrook", "StricklerNikuradse", "Strickler", "deBosBijkerk"].
            By default "WhiteColeBrook".
        friction_value : float, optional.
            Units corresponding to ''friction_type'' are ["Chézy C [m 1/2 /s]",
            "Manning n [s/m 1/3 ]", "Nikuradse k_n [m]", "Nikuradse k_n [m]",
            "Nikuradse k_n [m]", "Strickler k_s [m 1/3 /s]", "De Bos-Bijkerk γ [1/s]"]
            Friction value. By default 0.003.
        crosssections_shape : str, optional
            Shape of pipe crosssections. Either "circle" (default) or "rectangle".
        crosssections_value : int or list of int, optional
            Crosssections parameter value.
            If ``crosssections_shape`` = "circle", expects a diameter
            (default with 0.5 m) [m]
            If ``crosssections_shape`` = "rectangle", expects a list with
            [width, height] (e.g. [1.0, 1.0]) [m]. closed rectangle by default.
        dem_fn: str, optional
            Name of data source for dem data. Used to derive default invert levels
            values (DEM - pipes_depth - pipes diameter/height).

            * Required variables: [elevtn]
        pipes_depth: float, optional
            Depth of the pipes underground [m] (default 2.0 m). Used to derive defaults
            invert levels values (DEM - pipes_depth - pipes diameter/height).
        pipes_invlev: float, optional
            Constant default invert levels of the pipes [m asl] (default -2.5 m asl).
            This method is recommended to be used together with the dem method to fill
            remaining nan values. It slone is not a recommended method.
        snap_offset: float, optional
            Snapping tolenrance to automatically connecting branches. Tolenrance must be
            smaller than the shortest pipe length.
            By default 0.0, no snapping is applied.
        allow_intersection_snapping: bool, optional
            Switch to choose whether snapping of multiple branch ends are allowed when
            ``snap_offset`` is used.
            By default True.

        See Also
        --------
        dflowfm._setup_branches
        dflowfm._setup_crosssections
        """
        self.logger.info("Preparing 1D pipes.")

        # filter for allowed columns
        br_type = "pipe"
        _allowed_columns = [
            "geometry",
            "branchid",
            "branchtype",
            "branchorder",
            "material",
            "spacing",
            "branchid",
            "frictionvalue",
            "shape",  # circle or rectangle
            "diameter",  # circle
            "width",  # rectangle
            "height",  # rectangle
            "invlev_up",
            "invlev_dn",
        ]

        # Read data and filter within region
        region = workflows.parse_region_geometry(region, self.crs)
        gdf_br = self.data_catalog.get_geodataframe(
            pipes_fn, geom=region, buffer=0, predicate="intersects"
        )
        # Read defaults table
        defaults = self.data_catalog.get_dataframe(pipes_defaults_fn)
        # Prepare branches
        pipes, pipe_nodes = workflows.prepare_branches(
            gdf_br,
            params=defaults,
            br_type=br_type,
            dst_crs=self.crs,
            id_start=len(self.branches) + 1,
            filter=pipe_filter,
            spacing=spacing,
            snap_offset=snap_offset,
            allow_intersection_snapping=allow_intersection_snapping,
            allowed_columns=_allowed_columns,
            logger=self.logger,
        )
        # Prepare friction and crosssections
        pipes = workflows.prepare_default_friction_and_crosssection(
            pipes,
            br_type=br_type,
            friction_type=friction_type,
            friction_value=friction_value,
            crosssections_shape=crosssections_shape,
            crosssections_value=crosssections_value,
            logger=self.logger,
        )
        # filter extra time for geting clipped pipes within the region (better match)
        # remove the index name to avoid "ValueError: cannot insert branchid,
        # already exists" in geopandas>=1
        pipes.index.name = None
        pipes = gpd.sjoin(pipes, region, predicate="within")

        # setup crosssections
        # setup invert levels
        # 1. check if invlev up and dn are fully filled in (nothing needs to be done)
        if "invlev_up" and "invlev_dn" in pipes.columns:
            inv = pipes[["invlev_up", "invlev_dn"]]
            if inv.isnull().sum().sum() > 0:  # nodata values in pipes for invert levels
                fill_invlev = True
                self.logger.info(
                    f"{pipes_fn} data has {inv.isnull().sum().sum()} no data values"
                    "for invert levels. Will be filled using dem_fn or"
                    f"default value {pipes_invlev}"
                )
            else:
                fill_invlev = False
        else:
            fill_invlev = True
            self.logger.info(
                f"{pipes_fn} does not have columns [invlev_up, invlev_dn]."
                "Invert levels will be generated from dem_fn or"
                f"default value {pipes_invlev}"
            )
        # 2. filling use dem_fn + pipe_depth
        if fill_invlev and dem_fn is not None:
            dem = self.data_catalog.get_rasterdataset(
                dem_fn,
                geom=region,
                buffer=self._network_snap_offset,
                variables=["elevtn"],
            )
            pipes = workflows.invert_levels_from_dem(
                gdf=pipes, dem=dem, depth=pipes_depth
            )
            if pipes[["invlev_up", "invlev_dn"]].isnull().sum().sum() > 0:
                fill_invlev = True
            else:
                fill_invlev = False
        # 3. filling use pipes_invlev
        if fill_invlev and pipes_invlev is not None:
            self.logger.warning(
                "!Using a constant up and down invert levels for all pipes."
                "May cause issues when running the delft3dfm model.!"
            )
            df_inv = pd.DataFrame(
                data={
                    "branchtype": ["pipe"],
                    "invlev_up": [pipes_invlev],
                    "invlev_dn": [pipes_invlev],
                }
            )
            pipes = workflows.update_data_columns_attributes(
                pipes, df_inv, brtype="pipe"
            )

        # TODO: check that geometry lines are properly oriented from up to dn
        # when deriving invert levels from dem

        # Update crosssections object
        crosssections = self._setup_crosssections(
            pipes,
            crosssections_type="branch",
            midpoint=False,
        )
        # add crosssections to exisiting ones and update geoms
        self.logger.debug("Adding crosssections vector to geoms.")
        crosssections = workflows.add_crosssections(
            self.geoms.get("crosssections"), crosssections
        )
        self.set_geoms(crosssections, "crosssections")

        # setup geoms
        self.logger.debug("Adding pipes and pipe_nodes vector to geoms.")
        self.set_geoms(pipes, "pipes")
        self.set_geoms(pipe_nodes, "pipe_nodes")  # TODO: for manholes

        # add to branches
        branches = workflows.add_branches(
            self.mesh_datasets.get("mesh1d"),
            self.branches,
            pipes,
            self._snap_newbranches_to_branches_at_snapnodes,
            self._network_snap_offset,
        )
        workflows.validate_branches(branches)
        self.set_branches(branches)

        # update mesh
        # Convert self.branches to xugrid
        mesh1d, network1d = workflows.mesh1d_network1d_from_branches(
            self.opensystem,
            self.closedsystem,
            self._openwater_computation_node_distance,
        )
        self.set_mesh(network1d, grid_name="network1d", overwrite_grid=True)
        self.set_mesh(mesh1d, grid_name="mesh1d", overwrite_grid=True)

    def _setup_crosssections(
        self,
        branches,
        region: gpd.GeoDataFrame = None,
        crosssections_fn: str = None,
        crosssections_type: str = "branch",
        midpoint=True,
        maxdist=1.0,
    ) -> gpd.GeoDataFrame:
        """Prepare 1D crosssections from branches, points and xyz.

        # TODO to be extended also from dem data for rivers/channels?

        Crosssection must only be used after friction has been setup.

        Crosssections are read from ``crosssections_fn``.
        Crosssection types of this file is read from ``crosssections_type``

        If ``crosssections_fn`` is not defined, default method is
        ``crosssections_type`` = 'branch', meaning that branch attributes will be used
        to derive regular crosssections.
        Crosssections are derived at branches mid points if ``midpoints`` is True,
        else at both upstream and downstream extremities of branches if False.

        Adds/Updates model layers:

        * **crosssections** geom: 1D crosssection vector

        Parameters
        ----------
        branches : gpd.GeoDataFrame
            geodataframe of the branches to apply crosssections.

            * Required variables: [branchid, branchtype, branchorder]
            If ``crosssections_type`` = "branch"
                if shape = 'circle': 'diameter'
                if shape = 'rectangle': 'width', 'height', 'closed'
                if shape = 'trapezoid': 'width', 't_width', 'height', 'closed'
            * Optional variables: [material, friction_type, friction_value]
        region : gpd.GeoDataFrame, optional
            geodataframe of the region of interest for extracting crosssections_fn,
            by default None
        crosssections_fn : str Path, optional
            Name of data source for crosssections, see data/data_sources.yml.
            Note that for point crossections, only ones within the snap_network_offset
            will be used.
            If ``crosssections_type`` = "xyz"
            Note that only points within the region + 1000m buffer will be read.

            * Required variables: crsid, order, z
            * Optional variables:
            If ``crosssections_type`` = "point"

            * Required variables: crsid, shape, shift, closed
            * Optional variables:
                if shape = 'rectangle': 'width', 'height'
                if shape = 'trapezoid': 'width', 't_width', 'height'
                if shape = 'yz': 'yzcount','ycoordinates','zcoordinates'
                if shape = 'zw': 'numlevels', 'levels', 'flowwidths','totalwidths',
                    'fricitonid', 'frictiontype', 'frictionvalue'
                if shape = 'zwRiver': Not Supported
                Note that list input must be strings seperated by a whitespace ''.
            By default None, crosssections will be set from branches
        crosssections_type : {'branch', 'xyz', 'point'}
            Type of crosssections read from crosssections_fn. One of
            ['branch', 'xyz', 'point'].
            By default `branch`.
        maxdist: float, optional
            Maximum distance allowed for crosssections to be applied on branches.
            Only used for `crosssections_type` = point.
            By default 1.0.

        Returns
        -------
        gdf_cs : gpd.GeoDataFrame
            geodataframe of the new cross-sections

        Raise:
        ------
        NotImplementedError: if ``crosssection_type`` is not recongnised.
        """
        # setup crosssections
        if crosssections_fn is None and crosssections_type == "branch":
            # TODO: set a seperate type for rivers because other branch types
            # might require upstream/downstream
            # TODO: check for required columns
            # read crosssection from branches
            self.logger.info("Preparing crossections from branch.")
            gdf_cs = workflows.set_branch_crosssections(branches, midpoint=midpoint)

        elif crosssections_type == "xyz":
            # Read the crosssection data
            gdf_cs = self.data_catalog.get_geodataframe(
                crosssections_fn,
                geom=region,
                buffer=1000,
                predicate="contains",
            )

            # check if feature valid
            if len(gdf_cs) == 0:
                self.logger.warning(
                    f"No {crosssections_fn} 1D xyz crosssections found within domain"
                )
                return None
            valid_attributes = gis_utils.check_gpd_attributes(
                gdf_cs, required_columns=["crsid", "order", "z"]
            )
            if not valid_attributes:
                self.logger.error(
                    "Required attributes [crsid, order, z] in xyz crosssections"
                    "do not exist"
                )
                return None

            # assign id
            id_col = "crsid"
            gdf_cs.index = gdf_cs[id_col]
            gdf_cs.index.name = id_col

            # reproject to model crs
            gdf_cs.to_crs(self.crs)

            # set crsloc and crsdef attributes to crosssections
            self.logger.info(f"Preparing 1D xyz crossections from {crosssections_fn}")
            gdf_cs = workflows.set_xyz_crosssections(branches, gdf_cs)

        elif crosssections_type == "point":
            # Read the crosssection data
            gdf_cs = self.data_catalog.get_geodataframe(
                crosssections_fn,
                geom=region,
                buffer=100,
                predicate="contains",
            )

            # check if feature valid
            if len(gdf_cs) == 0:
                self.logger.warning(
                    f"No {crosssections_fn} 1D point crosssections found within domain"
                )
                return None
            valid_attributes = gis_utils.check_gpd_attributes(
                gdf_cs, required_columns=["crsid", "shape", "shift"]
            )
            if not valid_attributes:
                self.logger.error(
                    "Required attributes [crsid, shape, shift] in point crosssections"
                    "do not exist"
                )
                return None

            # assign id
            id_col = "crsid"
            gdf_cs.index = gdf_cs[id_col]
            gdf_cs.index.name = id_col

            # reproject to model crs
            gdf_cs.to_crs(self.crs)

            # set crsloc and crsdef attributes to crosssections
            self.logger.info(f"Preparing 1D point crossections from {crosssections_fn}")
            gdf_cs = workflows.set_point_crosssections(
                branches, gdf_cs, maxdist=maxdist
            )
        else:
            raise NotImplementedError(
                f"Method {crosssections_type} is not implemented."
            )

        return gdf_cs

    def setup_manholes(
        self,
        manholes_fn: str = None,
        manholes_defaults_fn: str = "manholes_defaults",
        bedlevel_shift: float = -0.5,
        dem_fn: str = None,
        snap_offset: float = 1e-3,
    ):
        """
        Prepare the 1D manholes to pipes or tunnels.

        Can only be used after all branches are setup.

        The manholes are generated based on a set of standards specified in
        ``manholes_defaults_fn`` (default)  and can be overwritten with manholes
        read from ``manholes_fn``.

        Use ``manholes_fn`` to set the manholes from a dataset of point locations.
        Only locations within the model region are selected. They are snapped to the
        model network nodes locations within a max distance defined in ``snap_offset``.

        Manhole attributes ["area", "streetstoragearea", "storagetype", "streetlevel"]
        are either taken from ``manholes_fn`` or filled in using defaults in
        ``manholes_defaults_fn``.
        Manhole attribute ["bedlevel"] is always generated from invert levels of the
        pipe/tunnel network plus a shift defined in ``bedlevel_shift``. This is needed
        for numerical stability.
        Manhole attribute ["streetlevel"]  can also be overwriten with values dervied
        from "dem_fn".

        #TODO probably needs another parameter to apply different sampling method for
        the manholes, e.g. min within 2 m radius.

        Adds/Updates model layers:

        * **manholes** geom: 1D manholes vector

        Parameters
        ----------
        manholes_fn: str Path, optional
            Path or data source name for manholes see data/data_sources.yml.
            Note only the points that are within the region polygon will be used.

            * Optional variables: ["area", "streetstoragearea", "storagetype",
            "streetlevel"]
        manholes_defaults_fn : str Path, optional
            Path to a csv file containing all defaults values per "branchtype".
            Use multiple rows to apply defaults per ["shape", "diameter"/"width"] pairs.
            By default `hydrolib.hydromt_delft3dfm.data.manholes.manholes_defaults.csv`
            is used.

            * Allowed variables: ["area", "streetlevel", "streeStorageArea",
            "storagetype"]
        dem_fn: str, optional
            Name of data source for dem data. Used to derive default invert levels
            values (DEM - pipes_depth - pipes diameter/height).

            * Required variables: [elevtn]
        bedlevel_shift: float, optional
            Shift applied to lowest pipe invert levels to derive manhole bedlevels [m]
            (default -0.5 m, meaning bedlevel = pipe invert - 0.5m).
        snap_offset: float, optional
            Snapping tolenrance to automatically connecting manholes to network nodes.
            By default 0.001. Use a higher value if large number of user manholes are
            missing.
        """
        # geom columns for manholes
        # id = storage node id, considered identical to manhole id
        # when using single compartment manholes
        _allowed_columns = [
            "geometry",
            "id",
            "name",
            "manholeid",
            "nodeid",
            "area",
            "bedlevel",
            "streetlevel",
            "streetstoragearea",
            "storagetype",
            "usetable",
        ]

        # generate manhole locations and bedlevels
        self.logger.info("generating manholes locations and bedlevels. ")
        manholes, branches = workflows.generate_manholes_on_branches(
            self.branches,
            bedlevel_shift=bedlevel_shift,
            use_branch_variables=["diameter", "width"],
            id_prefix="manhole_",
            id_suffix="_generated",
            logger=self.logger,
        )
        # FIXME Xiaohan: why do we need set_branches here? Because of branches.gui
        # --> add a high level write_gui files same level as write_mesh
        self.set_branches(branches)

        # add manhole attributes from defaults
        defaults = self.data_catalog.get_dataframe(manholes_defaults_fn)

        # add defaults
        manholes = workflows.update_data_columns_attributes(manholes, defaults)

        # read user manhole
        if manholes_fn:
            self.logger.info(f"reading manholes street level from file {manholes_fn}. ")
            # read
            gdf_manhole = self.data_catalog.get_geodataframe(
                manholes_fn,
                geom=self.region,
                buffer=self._network_snap_offset,
                predicate="contains",
            )
            # reproject
            if gdf_manhole.crs != self.crs:
                gdf_manhole = gdf_manhole.to_crs(self.crs)
            # filter for allowed columns
            allowed_columns = set(_allowed_columns).intersection(gdf_manhole.columns)
            self.logger.debug(
                f'filtering for allowed columns:{",".join(allowed_columns)}'
            )
            gdf_manhole = gpd.GeoDataFrame(
                gdf_manhole[list(allowed_columns)], crs=gdf_manhole.crs
            )
            # replace generated manhole using user manholes
            self.logger.debug("overwriting generated manholes using user manholes.")
            manholes = hydromt.gis_utils.nearest_merge(
                manholes, gdf_manhole, max_dist=snap_offset, overwrite=True
            )

        # generate manhole streetlevels from dem
        if dem_fn is not None:
            self.logger.info("overwriting manholes street level from dem. ")
            dem = self.data_catalog.get_rasterdataset(
                dem_fn,
                geom=self.region,
                variables=["elevtn"],
                buffer=self._network_snap_offset,
            )
            # reproject of manholes is done in sample method
            manholes["_streetlevel_dem"] = dem.raster.sample(manholes).values
            manholes["_streetlevel_dem"].fillna(manholes["streetlevel"], inplace=True)
            manholes["streetlevel"] = manholes["_streetlevel_dem"]
            self.logger.debug(
                f'street level mean is {np.mean(manholes["streetlevel"])}'
            )

        # internal administration
        # drop duplicated manholeid
        self.logger.debug("dropping duplicated manholeid")
        manholes.drop_duplicates(subset="manholeid")
        # add nodeid to manholes
        network1d_nodes = mesh_utils.network1d_nodes_geodataframe(
            self.mesh_datasets["network1d"]
        )
        manholes = hydromt.gis_utils.nearest_merge(
            manholes, network1d_nodes, max_dist=0.1, overwrite=False
        )
        # add additional required columns
        manholes["id"] = manholes["nodeid"]
        # id of the storage nodes id, identical to manholeid
        # when single compartment manholes are used
        manholes["name"] = manholes["manholeid"]
        manholes["usetable"] = False

        # validate
        if manholes[_allowed_columns].isna().any().any():
            self.logger.error(
                "manholes contain no data."
                "Use manholes_defaults_fn to apply no data filling."
            )

        # setup geoms
        self.logger.debug("Adding manholes vector to geoms.")
        self.set_geoms(manholes, "manholes")

    def setup_retentions(
        self,
        retentions_fn: str = None,
        retention_defaults_fn: str = "retentions_defaults",
        snap_offset: float = 1.0,
    ):
        """
        Prepare the 1D retentions branches.

        Can only be used after all branches are setup.
        The retentions are read from ``retentions_fn``.

        Use ``retentions_fn`` to set the retentions from a dataset of point locations.
        Only locations within the model region are selected. They are snapped to the
        model network nodes within a max distance defined in ``snap_offset``.
        #TODO: allow branch and chainage
        retention attributes ["numLevels", "levels", "storageArea", "interpolate"]
        are either taken from ``retentions_fn`` or filled in using defaults in
        ``retention_defaults_fn``.

        Adds/Updates model layers:
        * **retentions** geom: 1D retentions vector

        Parameters
        ----------
        retentions_fn: str Path, optional
            Path or data source name for retentions_fn, see data/data_sources.yml.
            Note: only the points that are within the region polygon will be used.
            * Optional variables: ["numLevels", "levels", "storageArea", "interpolate"]
        retention_defaults_fn : str Path, optional
            Path to a csv file containing all defaults values per "branchtype".
            By default
            `hydrolib.hydromt_delft3dfm.data.storages.retentions_defaults.csv` is used.
            * Allowed variables: ["numLevels", "levels", "storageArea", "interpolate"]
        snap_offset: float, optional
            Snapping tolenrance to automatically connecting retentions to network nodes.
            By default 0.001. Use a higher value if needed.
        """
        _allowed_columns = [
            "geometry",
            "id",
            "name",
            "nodeid",
            "numlevels",
            "levels",
            "storagearea",
            "interpolate",
            "usetable",
        ]

        # read user manhole
        if retentions_fn:
            self.logger.info(f"reading retentions from {retentions_fn}. ")
            # read
            gdf_retentions = self.data_catalog.get_geodataframe(
                retentions_fn,
                geom=self.region,
                buffer=self._network_snap_offset,
                predicate="contains",
            )
            # reproject
            if gdf_retentions.crs != self.crs:
                gdf_retentions = gdf_retentions.to_crs(self.crs)
            # set index
            if "id" not in gdf_retentions:
                gdf_retentions["id"] = [
                    f"retention_{i}" for i in range(len(gdf_retentions))
                ]
            else:
                gdf_retentions["id"] = gdf_retentions["id"].astype(str)
            # filter for allowed columns
            allowed_columns = set(_allowed_columns).intersection(gdf_retentions.columns)
            self.logger.debug(
                f'filtering for allowed columns:{",".join(allowed_columns)}'
            )
            gdf_retentions = gpd.GeoDataFrame(
                gdf_retentions[list(allowed_columns)], crs=gdf_retentions.crs
            )
            # add defaults
            gdf_retentions["branchtype"] = "river"
            defaults = self.data_catalog.get_dataframe(retention_defaults_fn)
            gdf_retentions = workflows.update_data_columns_attributes(
                gdf_retentions, defaults
            )
            if len(gdf_retentions) == 0:
                self.logger.error("No retention ponds were setup.")
                return None
        else:
            raise ValueError("Path to 'retentions_fn' missing.")

        if len(gdf_retentions) > 0:
            self.logger.info(f"Process {len(gdf_retentions)} retention ponds.")
            # add nodeid to retentions
            network1d_nodes = mesh_utils.network1d_nodes_geodataframe(
                self.mesh_datasets["network1d"]
            )
            retentions = hydromt.gis_utils.nearest_merge(
                gdf_retentions, network1d_nodes, max_dist=snap_offset, overwrite=False
            )
            # drop not snapped
            dropped_retentions = retentions[retentions["nodeid"].isna()]
            retentions = retentions[~retentions["nodeid"].isna()]
            self.logger.info(f"Drop unsnapped {list(dropped_retentions.id)}")

            # drop duplicated nodeid
            self.logger.debug("Dropping duplicated retentions")
            retentions.drop_duplicates(subset=["nodeid"])

            if len(retentions) == 0:
                self.logger.error("No retention ponds left after processing.")
                return None

            # add additional required columns
            retentions["name"] = retentions["id"]
            retentions["usetable"] = True
            retentions["numlevels"] = retentions["levels"].apply(
                lambda x: len(x.split())
            )

            retentions = gpd.GeoDataFrame(
                retentions[list(_allowed_columns)], crs=gdf_retentions.crs
            )

            # setup geoms
            self.logger.debug("Adding retentions vector to geoms.")
            self.set_geoms(retentions, "retentions")

    def setup_1dboundary(
        self,
        boundaries_geodataset_fn: str = None,
        boundaries_timeseries_fn: str = None,
        boundary_value: float = -2.5,
        branch_type: str = "river",
        boundary_type: str = "waterlevel",
        boundary_unit: str = "m",
        boundary_locs: str = "downstream",
        snap_offset: float = 1.0,
    ):
        """
        Prepare the 1D ``boundary_type`` boundaries using timeseries or a constant.

        Boundaries are prepared for a specific ``branch_type`` at the ``boundary_locs``
        locations. E.g. 'waterlevel' boundaries for 'downstream''river' branches.

        The values can either be a constant using ``boundary_value`` (default) or
        timeseries read from ``boundaries_geodataset_fn``.

        Use ``boundaries_geodataset_fn`` to set the boundary values from a dataset
        of point location timeseries. Only locations within the possible model boundary
        locations + snap_offset are used. They are snapped to the model boundary
        locations within a max distance defined in ``snap_offset``. If
        ``boundaries_geodataset_fn`` has missing values, the constant ``boundary_value``
        will be used.

        The dataset/timeseries are clipped to the model time based on the model config
        tstart and tstop entries.

        Adds/Updates model layers:

        * **boundary1d_{boundary_type}bnd_{branch_type}** forcing: 1D boundaries
            DataArray

        Parameters
        ----------
        boundaries_geodataset_fn : str, Path
            Path or data source name for geospatial point timeseries file.
            This can either be a netcdf file with geospatial coordinates
            or a combined point location file with a timeseries data csv file
            which can be setup through the data_catalog yml file.

            * Required variables if netcdf: ['discharge', 'waterlevel'] depending on
                ``boundary_type``
            * Required coordinates if netcdf: ['time', 'index', 'y', 'x']

            * Required variables if a combined point location file: ['index'] with type
                int
            * Required index types if a time series data csv file: int
            NOTE: Require equidistant time series
        boundaries_timeseries_fn: str, Path
            Path to tabulated timeseries csv file with time index in first column
            and location IDs in the first row,
            see :py:meth:`hydromt.open_timeseries_from_table`, for details.
            NOTE: tabulated timeseries files can only in combination with point location
            coordinates be set as a geodataset in the data_catalog yml file.
            NOTE: Require equidistant time series
        boundary_value : float, optional
            Constant value to use for all boundaries if ``boundaries_geodataset_fn`` is
            None and to fill in missing data. By default -2.5 m.
        branch_type: str
            Type of branch to apply boundaries on. One of ["river", "pipe"].
        boundary_type : str, optional
            Type of boundary tu use. One of ["waterlevel", "discharge"].
            By default "waterlevel".
        boundary_unit : str, optional.
            Unit corresponding to [boundary_type].
            If ``boundary_type`` = "waterlevel"
                Allowed unit is [m]
            if ''boundary_type`` = "discharge":
               Allowed unit is [m3/s]
            By default m.
        boundary_locs:
            Boundary locations to consider. One of ["upstream", "downstream", "both"].
            Only used for river waterlevel which can be upstream, downstream or both.
            By default "downstream".
            For the others, it is automatically derived from branch_type and
            boundary_type.
        snap_offset : float, optional
            Snapping tolerance to automatically applying boundaries at the correct
            network nodes. By default 0.1, a small snapping is applied to avoid
            precision errors.
        """
        self.logger.info(f"Preparing 1D {boundary_type} boundaries for {branch_type}.")

        # 1. get potential boundary locations based on branch_type and boundary_type
        boundaries_branch_type = workflows.select_boundary_type(
            self.boundaries, branch_type, boundary_type, boundary_locs
        )

        # 2. read boundary from user data
        _, da_bnd = self._read_forcing_geodataset(
            boundaries_geodataset_fn, boundary_type, snap_offset
        )

        # 3. Derive DataArray with boundary values at boundary locations
        # in boundaries_branch_type
        da_out = workflows.compute_boundary_values(
            boundaries=boundaries_branch_type,
            da_bnd=da_bnd,
            boundary_value=boundary_value,
            boundary_type=boundary_type,
            boundary_unit=boundary_unit,
            snap_offset=snap_offset,
            logger=self.logger,
        )

        # 4. set boundaries
        self.set_forcing(da_out, name=f"boundary1d_{da_out.name}_{branch_type}")
        # FIXME: this format cannot be read back due to lack of branch type info
        # from model files

    def _read_forcing_geodataset(
        self,
        forcing_geodataset_fn: Union[str, Path],
        forcing_name: str = "discharge",
        region_buffer=0.0,
    ):
        """Read forcing geodataset."""
        refdate, tstart, tstop = self.get_model_time()  # time slice

        if (
            forcing_geodataset_fn is not None
            and self.data_catalog[forcing_geodataset_fn].data_type == "GeoDataset"
        ):
            da = self.data_catalog.get_geodataset(
                forcing_geodataset_fn,
                geom=self.region.buffer(region_buffer),  # buffer region
                variables=[forcing_name],
                time_tuple=(tstart, tstop),
            )
            # error if time mismatch
            if np.logical_and(
                pd.to_datetime(da.time.values[0]) == pd.to_datetime(tstart),
                pd.to_datetime(da.time.values[-1]) == pd.to_datetime(tstop),
            ):
                pass
            else:
                self.logger.error(
                    "Forcing has different start and end time."
                    + " Please check the forcing file. Support yyyy-mm-dd HH:MM:SS. "
                )
            # reproject if needed and convert to location
            if da.vector.crs != self.crs:
                da = da.vector.to_crs(self.crs)
                # TODO update after HydroMT release >0.9.0
            # get geom
            gdf = da.vector.to_gdf(reducer=np.mean)
        elif (
            forcing_geodataset_fn is not None
            and self.data_catalog[forcing_geodataset_fn].data_type == "GeoDataFrame"
        ):
            gdf = self.data_catalog.get_geodataframe(
                forcing_geodataset_fn,
                geom=self.region.buffer(self._network_snap_offset),
            )
            # reproject
            if gdf.crs != self.crs:
                gdf = gdf.to_crs(self.crs)
            da = None
        else:
            gdf = None
            da = None
        return gdf, da

    def setup_1dlateral_from_points(
        self,
        laterals_geodataset_fn: str = None,
        lateral_value: float = 0.0,
        snap_offset: float = 1.0,
        branch_type: str = "river",
    ):
        """
        Prepare the 1D lateral discharge from geodataset of point geometries.

        E.g. '1' m3/s for all lateral locations.

        Use ``laterals_geodataset_fn`` to set the lateral values from a geodataset
        of point locations.
        Support also geodataframe of point locations in combination of `lateral_value`.

        Only locations that are snapped to the network of `branch_type`
        within a max distance defined in ``snap_offset`` are used.

        The discharge can either be a constant using ``lateral_value`` (default) or
        a timeseries read from ``laterals_geodataset_fn``.
        If the timeseries has missing values, constant ``lateral_value`` will be used.

        The timeseries are clipped to the model time based on the model config
        tstart and tstop entries.

        Adds/Updates model layers:
            ** lateral1d_points** forcing: DataArray with points coordinates.

        Parameters
        ----------
        laterals_geodataset_fn : str, Path
            Path or data source name for geospatial point location file.
            * Required variables if geodataset is provided ['lateral_discharge']
            NOTE: Require equidistant time series
        lateral_value : float, optional
            Constant value, used if ``laterals_geodataset_fn`` is a geodataframe,
            or for filling in missing data.
            By default 0 [m3/s].
        snap_offset : float, optional
            Snapping tolerance to snap boundaries to the correct network nodes.
            By default 0.1, a small snapping is applied to avoid precision errors.
        branch_type: str, optional
            Type of branch to apply laterals on. One of ["river", "pipe"].
            If None, all branches are used.
            By defalt None.
        """
        self.logger.info(f"Preparing 1D laterals for {branch_type}.")
        network_by_branchtype = self.staticgeoms[f"{branch_type}s"]

        # 1. read lateral geodataset and snap to network
        gdf_laterals, da_lat = self._read_forcing_geodataset(
            laterals_geodataset_fn, "lateral_discharge", snap_offset
        )

        # snap laterlas to selected branches
        gdf_laterals = workflows.snap_geom_to_branches_and_drop_nonsnapped(
            branches=network_by_branchtype.set_index("branchid"),
            geoms=gdf_laterals,
            snap_offset=snap_offset,
        )

        if len(gdf_laterals) == 0:
            return None

        # 2. Compute lateral dataarray
        da_out = workflows.compute_forcing_values_points(
            gdf=gdf_laterals,
            da=da_lat,
            forcing_value=lateral_value,
            forcing_type="lateral_discharge",
            forcing_unit="m3/s",
            logger=self.logger,
        )

        # 3. set laterals
        self.set_forcing(da_out, name="lateral1d_points")

    def setup_1dlateral_from_polygons(
        self,
        laterals_geodataset_fn: str = None,
        lateral_value: float = -2.5,
    ):
        """
        Prepare the 1D lateral discharge from geodataset of polygons.

        E.g. '1' m3/s for all lateral locations.

        Use ``laterals_geodataset_fn`` to set the lateral values from a geodatasets
        of polygons.
        Support also geodataframe of polygons in combination of `lateral_value`.

        The discharge can either be a constant using ``lateral_value`` (default) or
        a timeseries read from ``laterals_geodataset_fn``.
        If the timeseries has missing values, constant ``lateral_value`` will be used.

        The timeseries are clipped to the model time based on the model config
        tstart and tstop entries.

        Adds/Updates model layers:
            * ** lateral1d_polygons** forcing: DataArray with polygon coordinates.

        Parameters
        ----------
        laterals_geodataset_fn : str, Path
            Path or data source name for geospatial point location file.
            * Required variables if geodataset is provided ['lateral_discharge']
            NOTE: Require equidistant time series
        lateral_value : float, optional
            Constant value, used if ``laterals_geodataset_fn`` is a geodataframe,
            or for filling in missing data.
            By default 0 [m3/s].
        """
        self.logger.info("Preparing 1D laterals for polygons.")

        # 1. read lateral geodataset
        gdf_laterals, da_lat = self._read_forcing_geodataset(
            laterals_geodataset_fn, "lateral_discharge"
        )

        if len(gdf_laterals) == 0:
            return None

        # 2. Compute lateral dataarray
        da_out = workflows.compute_forcing_values_polygon(
            gdf=gdf_laterals,
            da=da_lat,
            forcing_value=lateral_value,
            forcing_type="lateral_discharge",
            forcing_unit="m3/s",
            logger=self.logger,
        )

        # 3. set laterals
        self.set_forcing(da_out, name="lateral1d_polygons")

    def setup_bridges(
        self,
        bridges_fn: Optional[str] = None,
        bridges_defaults_fn: Optional[str] = "1D_bridges_defaults",
        bridge_filter: Optional[str] = None,
        snap_offset: Optional[float] = None,
    ):
        """Prepare bridges, including bridge locations and bridge crossections.

        The bridges are read from ``bridges_fn`` and if any missing, filled with
        information provided in ``bridges_defaults_fn``.

        When reading ``bridges_fn``, only locations within the region will be read.
        Read locations are then filtered for value specified in ``bridge_filter`` on
        the column "structure_type".
        Remaining locations are snapped to the existing network within a max distance
        defined in ``snap_offset`` and will be dropped if not snapped.

        A default rectangle bridge profile can be found in ``bridges_defaults_fn``
        as an example.

        Structure attributes ['structure_id', 'structure_type'] are either taken from
        data or generated in the script.
        Structure attributes ['shape', 'diameter', 'width', 't_width', 'height',
        'closed', 'shift', 'length', 'pillarwidth', 'formfactor', 'friction_type',
        'friction_value', 'allowedflowdir',  'inletlosscoeff', 'outletlosscoeff']
        are either taken from data,  or in case of missing read from defaults.

        Adds/Updates model layers:

        * **bridges** geom: 1D bridges vector

        Parameters
        ----------
        bridges_fn: str Path, optional
            Path or data source name for bridges, see data/data_sources.yml.
            Note only the points that are within the region polygon will be used.

            * Optional variables: ['structure_id', 'structure_type', 'shape',
            'diameter', 'width', 't_width', 'height', 'closed', 'shift', 'length',
            'pillarwidth', 'formfactor', 'friction_type', 'friction_value',
            'allowedflowdir',  'inletlosscoeff', 'outletlosscoeff']

        bridges_defaults_fn : str Path, optional
            Path to a csv file containing all defaults values per "structure_type".
            By default `hydrolib.hydromt_delft3dfm.data.bridges.bridges_defaults.csv`
            is used. This file describes a minimum rectangle bridge profile.

            * Allowed variables: ['structure_type', 'shape', 'diameter', 'width',
                't_width', 'height', 'closed', 'shift', 'length', 'pillarwidth',
                'formfactor', 'friction_type', 'friction_value', 'allowedflowdir',
                'inletlosscoeff', 'outletlosscoeff']

        river_filter: str, optional
            Keyword in "structure_type" column of ``bridges_fn`` used to filter bridge
            features. If None all features are used (default).

        snap_offset: float, optional
            Snapping tolenrance to automatically snap bridges to network and add
            ['branchid', 'chainage'] attributes.
            By default None. In this case, global variable "network_snap_offset"
            will be used.

        See Also
        --------
        workflows.prepare_1dstructures
        """
        snap_offset = self._network_snap_offset if snap_offset is None else snap_offset
        _st_type = "bridge"
        _allowed_columns = [
            "id",
            "name",
            "type",
            "branchid",
            "chainage",
            "allowedflowdir",
            "pillarwidth",
            "formfactor",
            "csdefid",
            "shift",
            "inletlosscoeff",
            "outletlosscoeff",
            "frictiontype",
            "friction",
            "length",
        ]

        # read data
        gdf_bridges = self.data_catalog.get_geodataframe(
            bridges_fn, geom=self.region, buffer=0, predicate="contains"
        )
        defaults = self.data_catalog.get_dataframe(bridges_defaults_fn)

        # setup general 1d structures
        bridges = workflows.prepare_1dstructures(
            self.branches,
            gdf_bridges,
            params=defaults,
            filter=bridge_filter,
            snap_offset=snap_offset,
            st_type=_st_type,
        )

        # setup crossection definitions
        crosssections = bridges[[c for c in bridges.columns if c.startswith("crsdef")]]
        crosssections = workflows.add_crosssections(
            self.geoms.get("crosssections"), crosssections
        )
        self.set_geoms(crosssections, "crosssections")

        # setup bridges (note the allowed columns are different per structure)
        bridges.columns = bridges.columns.str.lower()
        bridges.rename(
            columns={
                "structure_id": "id",
                "structure_name": "name",
                "structure_type": "type",
                "structure_branchid": "branchid",
                "structure_chainage": "chainage",
                "crsdef_id": "csdefid",
                "friction_type": "frictiontype",
                "friction_value": "friction",
            },
            inplace=True,
        )
        # filter for allowed columns
        allowed_columns = set(_allowed_columns).intersection(bridges.columns)
        allowed_columns.update({"geometry"})
        bridges = gpd.GeoDataFrame(bridges[list(allowed_columns)], crs=bridges.crs)
        self.set_geoms(bridges, "bridges")

    def setup_culverts(
        self,
        culverts_fn: Optional[str] = None,
        culverts_defaults_fn: Optional[str] = "1D_culverts_defaults",
        culvert_filter: Optional[str] = None,
        snap_offset: Optional[float] = None,
    ):
        """Prepare culverts, including locations and crossections.

        Note that only subtype culvert is supported, i.e. inverted siphon is not
        supported.

        The culverts are read from ``culverts_fn`` and if any missing, filled with
        information provided in ``culverts_defaults_fn``.

        When reading ``culverts_fn``, only locations within the region will be read.
        Read locations are then filtered for value specified in ``culvert_filter`` on
        the column "structure_type" .
        Remaining locations are snapped to the existing network within a max distance
        defined in ``snap_offset`` and will be dropped if not snapped.

        A default ``culverts_defaults_fn`` that defines a circle culvert profile can be
        found in dflowfm.data.culverts as an example.

        Structure attributes ['structure_id', 'structure_type'] are either taken from
        data or generated in the script.
        Structure attributes ['shape', 'diameter', 'width', 't_width', 'height',
        'closed', 'leftlevel', 'rightlevel', 'length','valveonoff',
        'valveopeningheight', 'numlosscoeff', 'relopening', 'losscoeff',
        'friction_type', 'friction_value', 'allowedflowdir',  'inletlosscoeff',
        'outletlosscoeff'] are either taken from data,  or in case of missing read
        from defaults.

        Adds/Updates model layers:

        * **culverts** geom: 1D culverts vector

        Parameters
        ----------
        culverts_fn: str Path, optional
            Path or data source name for culverts, see data/data_sources.yml.
            Note only the points that are within the region polygon will be used.

            * Optional variables: ['structure_id', 'structure_type', 'shape',
                'diameter', 'width', 't_width', 'height', 'closed', 'leftlevel',
                'rightlevel', 'length', 'valveonoff', 'valveopeningheight',
                'numlosscoeff', 'relopening', 'losscoeff', 'friction_type',
                'friction_value', 'allowedflowdir',  'inletlosscoeff',
                'outletlosscoeff']

        culverts_defaults_fn : str Path, optional
            Path to a csv file containing all defaults values per "structure_type".
            By default `hydrolib.hydromt_delft3dfm.data.culverts.culverts_defaults.csv`
            is used.
            This file describes a default circle culvert profile.

            * Allowed variables: ['structure_type', 'shape', 'diameter', 'width',
                't_width', 'height', 'closed', 'leftlevel', 'rightlevel', 'length',
                'valveonoff', 'valveopeningheight', 'numlosscoeff', 'relopening',
                'losscoeff', 'friction_type', 'friction_value', 'allowedflowdir',
                'inletlosscoeff', 'outletlosscoeff']

        culvert_filter: str, optional
            Keyword in "structure_type" column of ``culverts_fn`` used to filter
            culvert features. If None all features are used (default).

        snap_offset: float, optional
            Snapping tolenrance to automatically snap culverts to network and add
            ['branchid', 'chainage'] attributes.
            By default None. In this case, global variable "network_snap_offset" will
            be used.

        See Also
        --------
        workflows.prepare_1dstructures
        """
        snap_offset = self._network_snap_offset if snap_offset is None else snap_offset
        _st_type = "culvert"
        _allowed_columns = [
            "id",
            "name",
            "type",
            "branchid",
            "chainage",
            "allowedflowdir",
            "leftlevel",
            "rightlevel",
            "csdefid",
            "length",
            "inletlosscoeff",
            "outletlosscoeff",
            "valveonoff",
            "valveopeningheight",
            "numlosscoeff",
            "relopening",
            "losscoeff",
            "bedfrictiontype",
            "bedfriction",
        ]

        # read data
        gdf_culverts = self.data_catalog.get_geodataframe(
            culverts_fn, geom=self.region, buffer=0, predicate="contains"
        )
        defaults = self.data_catalog.get_dataframe(culverts_defaults_fn)

        # setup general 1d structures
        culverts = workflows.prepare_1dstructures(
            self.branches,
            gdf_culverts,
            params=defaults,
            filter=culvert_filter,
            snap_offset=snap_offset,
            st_type=_st_type,
        )

        # setup crossection definitions
        crosssections = culverts[
            [c for c in culverts.columns if c.startswith("crsdef")]
        ]
        crosssections = workflows.add_crosssections(
            self.geoms.get("crosssections"), crosssections
        )
        self.set_geoms(crosssections, "crosssections")

        # setup culverts (note the allowed columns are different per structure)
        culverts.columns = culverts.columns.str.lower()
        culverts.rename(
            columns={
                "structure_id": "id",
                "structure_name": "name",
                "structure_type": "type",
                "structure_branchid": "branchid",
                "structure_chainage": "chainage",
                "crsdef_id": "csdefid",
                "friction_type": "bedfrictiontype",
                "friction_value": "bedfriction",
            },
            inplace=True,
        )
        # filter for allowed columns
        allowed_columns = set(_allowed_columns).intersection(culverts.columns)
        allowed_columns.update({"geometry"})
        culverts = gpd.GeoDataFrame(culverts[list(allowed_columns)], crs=culverts.crs)
        self.set_geoms(culverts, "culverts")

    def setup_mesh2d(
        self,
        region: dict,
        res: Optional[float] = None,
    ) -> xu.UgridDataset:
        """Create a 2D unstructured mesh according UGRID conventions.

        Grids are read according to UGRID conventions. A 2D unstructured mesh
        will be created as 2D rectangular grid from a geometry (geom_fn) or bbox.
        If an existing 2D mesh is given, then no new mesh will be generated but an
        extent can be extracted using the `bounds` argument of region.

        # Note that:
        # (1) Refinement of the mesh is a seperate setup function, however an existing
        # grid with refinement (mesh_fn) can already be read.
        # (2) If mesh also has 1D, 1D2Dlinks are created in a separate setup function.
        # (3) At mesh border, cells that intersect with geometry border will be kept.
        # (4) Only existing mesh with only 2D grid can be read. So 1D2D network files
        # are not supported as mesh2d_fn.

        Adds/Updates model layers:

        * **grid_name** mesh topology: add grid_name 2D topology to mesh object

        Parameters
        ----------
        region : dict
            Dictionary describing region of interest, bounds can be provided for type
            'mesh'.
        CRS for 'bbox' and 'bounds' should be 4326; e.g.:

            * {'bbox': [xmin, ymin, xmax, ymax]}

            * {'geom': 'path/to/polygon_geometry'}

            * {'mesh': 'path/to/2dmesh_file'}

            * {'mesh': 'path/to/2dmesh_file', 'bounds': [xmin, ymin, xmax, ymax]}
        res: float
            Resolution used to generate 2D mesh [unit of the CRS], required if region
            is not based on 'mesh'.

        Returns
        -------
        mesh2d : xu.UgridDataset
            Generated mesh2d.

        """  # noqa: E501
        # Create the 2dmesh
        mesh2d = create_mesh2d(
            region=region,
            res=res,
            crs=self.crs,
        )
        # Add to self.mesh
        self.set_mesh(mesh2d, grid_name="mesh2d", overwrite_grid=False)
        # update res
        self._res = res

    def setup_mesh2d_refine(
        self,
        polygon_fn: Optional[str] = None,
        sample_fn: Optional[str] = None,
        steps: Optional[int] = 1,
    ):
        """
        Refine the 2d mesh.

        Refinement are done within the geometry based on polygon `polygon_fn` or
        raster samples `sample_fn`.

        The number of refinement is defined by `steps` if `polygon_fn` is used.

        Note that this function can only be applied on an existing regular rectangle
        2d mesh.

        # FIXME: how to identify if the mesh is uniform rectangle 2d mesh?
        regular irregular?

        Adds/Updates model layers:

        * **2D mesh** geom: By any changes in 2D grid

        Parameters
        ----------
        polygon_fn : str Path, optional
            Path to a polygon or MultiPolygon used to refine the 2D mesh
        sample_fn  : str Path, optional
            Path to a raster sample file used to refine the 2D mesh.
            The value of each sample point is the number of steps to refine the mesh.
            Allow only single values. The resolution of the raster should be the same
            as the desired end result resolution.

            * Required variable: ['steps']
        steps : int, optional
            Number of steps in the refinement when `polygon_fn' is used.
            By default 1, i.e. no refinement is applied.

        """
        if "mesh2d" not in self.mesh_names:
            logger.error(
                "2d mesh is not available, use setup_mesh2d before refinement."
            )
            return

        if polygon_fn is not None:
            self.logger.info(f"reading geometry from file {polygon_fn}. ")
            # read
            gdf = self.data_catalog.get_geodataframe(
                polygon_fn, geom=self.region, buffer=0, predicate="contains"
            )
            # reproject
            if gdf.crs != self.crs:
                gdf = gdf.to_crs(self.crs)

        elif sample_fn is not None:
            self.logger.info(f"reading samples from file {sample_fn}. ")
            # read
            da = self.data_catalog.get_rasterdataset(
                sample_fn,
                geom=self.region,
                buffer=0,
                predicate="contains",
                variables=["steps"],
                single_var_as_array=True,
            ).astype(
                np.float64
            )  # float64 is needed by mesh kernel to convert into c double
            # reproject
            if da.raster.crs != self.crs:
                self.logger.warning(
                    "Sample grid has a different resolution than model."
                    "Reprojecting with nearest but some information might be lost."
                )
                da = da.raster.reproject(self.crs, method="nearest")

        # refine
        mesh2d, res = workflows.mesh2d_refine(
            mesh2d=self.get_mesh("mesh2d"),
            res=self._res,
            gdf_polygon=gdf if polygon_fn is not None else None,
            da_sample=da if sample_fn is not None else None,
            steps=steps,
            logger=self.logger,
        )

        # set mesh2d
        self.set_mesh(mesh2d, grid_name="mesh2d", overwrite_grid=True)
        # update res
        self._res = res

    def setup_link1d2d(
        self,
        link_direction: Optional[str] = "1d_to_2d",
        link_type: Optional[str] = "embedded",
        polygon_fn: Optional[str] = None,
        branch_type: Optional[str] = None,
        max_length: Union[float, None] = np.inf,
        dist_factor: Union[float, None] = 2.0,
        **kwargs,
    ):
        """
        Generate 1d2d links that link mesh1d and mesh2d according UGRID conventions.

        1d2d links are added to allow water exchange between 1d and 2d for a 1d2d model.
        They can only be added if both mesh1d and mesh2d are present. By default,
        1d_to_2d links are generated for the entire mesh1d except boundary locations.

        When ''polygon_fn'' is specified, only links within the polygon will be added.
        When ''branch_type'' is specified, only 1d branches matching the specified type
        will be used for generating 1d2d link.

        # TODO: This option should also allows more customised setup for pipes and
        tunnels: 1d2d links will also be generated at boundary locations.

        Parameters
        ----------
        link_direction : str, optional
            Direction of the links: ["1d_to_2d", "2d_to_1d"].
            Default to 1d_to_2d.
        link_type : str, optional
            Type of the links to be generated: ["embedded", "lateral"].
            Only used when ''link_direction'' = '2d_to_1d'.
            Default to None.
        polygon_fn: str Path, optional
            Source name of raster data in data_catalog.
            Default to None.
        branch_type: str, Optional
            Type of branch to be used for 1d: ["river","pipe","channel", "tunnel"].
            When ''branch_type'' = "pipe" or "tunnel" are specified, 1d2d links will
            also be generated at boundary locations.
            Default to None.
            Add 1d2d links for the all branches at non-boundary locations.
        max_length : Union[float, None], optional
            Max allowed edge length for generated links.
            Only used when ''link_direction'' = '2d_to_1d'  and
            ''link_type'' = 'lateral'.
            Defaults to infinity.
        dist_factor : Union[float, None], optional:
            Factor to determine which links are kept.
            Only used when ''link_direction'' = '2d_to_1d'  and
            ''link_type'' = 'lateral'.
            Defaults to 2.0.
            Links with an intersection distance larger than 2 times the center to edge
            distance of the cell, are removed.

        See Also
        --------
        workflows.links1d2d_add_links_1d_to_2d
        workflows.links1d2d_add_links_2d_to_1d_embedded
        workflows.links1d2d_add_links_2d_to_1d_lateral
        """
        # check existing network
        if "mesh1d" not in self.mesh_names or "mesh2d" not in self.mesh_names:
            self.logger.error(
                "cannot setup link1d2d: either mesh1d or mesh2d or both do not exist"
            )
            return None

        # FIXME: question - how to seperate if the user wants to update the
        # entire 1d2d links object or simply wants to add another set of links? #1
        # TODO: would be nice in hydrolib to allow clear of subset of 1d2d links
        # for specific branches

        # check input
        if polygon_fn is not None:
            within = self.data_catalog.get_geodataframe(polygon_fn).geometry
            self.logger.info(f"adding 1d2d links only within polygon {polygon_fn}")
        else:
            within = None

        if branch_type is not None:
            branchids = self.branches[
                self.branches.branchtype == branch_type
            ].branchid.to_list()  # use selective branches
            self.logger.info(f"adding 1d2d links for {branch_type} branches.")
        else:
            branchids = None  # use all branches
            self.logger.warning(
                "adding 1d2d links for all branches at non boundary locations."
            )

        # setup 1d2d links
        if link_direction == "1d_to_2d":
            self.logger.info("setting up 1d_to_2d links.")
            # recompute max_length based on the diagonal distance of the max mesh area
            max_length = np.sqrt(self.mesh_grids["mesh2d"].area.max()) * np.sqrt(2)
            link1d2d = workflows.links1d2d_add_links_1d_to_2d(
                self.mesh, branchids=branchids, within=within, max_length=max_length
            )

        elif link_direction == "2d_to_1d":
            if link_type == "embedded":
                self.logger.info("setting up 2d_to_1d embedded links.")

                link1d2d = workflows.links1d2d_add_links_2d_to_1d_embedded(
                    self.mesh, branchids=branchids, within=within
                )
            elif link_type == "lateral":
                self.logger.info("setting up 2d_to_1d lateral links.")
                link1d2d = workflows.links1d2d_add_links_2d_to_1d_lateral(
                    self.mesh,
                    branchids=branchids,
                    within=within,
                    max_length=max_length,
                    dist_factor=dist_factor,
                )
            else:
                self.logger.error(f"link_type {link_type} is not recognised.")

        else:
            self.logger.error(f"link_direction {link_direction} is not recognised.")

        # Add link1d2d to xu Ugrid mesh
        if len(link1d2d["link1d2d"]) == 0:
            self.logger.warning("No 1d2d links were generated.")
        else:
            self.set_link1d2d(link1d2d)

    def setup_maps_from_rasterdataset(
        self,
        raster_fn: str,
        variables: Optional[list] = None,
        fill_method: Optional[str] = None,
        reproject_method: Optional[str] = "nearest",
        interpolation_method: Optional[str] = "triangulation",
        locationtype: Optional[str] = "2d",
        name: Optional[str] = None,
        split_dataset: Optional[bool] = True,
    ) -> None:
        """
        Add data variable(s) from ``raster_fn`` to maps object.

        If raster is a dataset, all variables will be added unless ``variables`` list
        is specified.

        Adds model layers:

        * **raster.name** maps: data from raster_fn

        Parameters
        ----------
        raster_fn: str
            Source name of raster data in data_catalog.
        variables: list, optional
            List of variables to add to maps from raster_fn. By default all.
            Available variables: ['elevtn', 'waterlevel', 'waterdepth', 'pet',
            'infiltcap', 'roughness_chezy', 'roughness_manning',
            'roughness_walllawnikuradse', 'roughness_whitecolebrook']
        fill_method : str, optional
            If specified, fills no data values using fill_nodata method. Available
            methods are ['linear', 'nearest', 'cubic', 'rio_idw'].
        reproject_method : str, optional
            CRS reprojection method from rasterio.enums.Resampling. By default nearest.
            Available methods: [ 'nearest', 'bilinear', 'cubic', 'cubic_spline',
            'lanczos', 'average', 'mode', 'gauss', 'max', 'min', 'med', 'q1', 'q3',
            'sum', 'rms']
        interpolation_method : str, optional
            Interpolation method for DFlow-FM. By default mean for waterlevel and
            waterdepth, and triangulation for all other variables. When methods other
            than 'triangulation' are used, the relative search cell size will be
            estimated based on resolution of the raster.
            Available methods: ['triangulation', 'mean', 'nearestNb', 'max', 'min',
            'invDist', 'minAbs', 'median']
        locationtype : str, optional
            LocationType in initial fields. Either 2d (default), 1d or all.
        name: str, optional
            Variable name, only in case data is of type DataArray or if a Dataset is
            added as is (split_dataset=False).
        split_dataset: bool, optional
            If data is a xarray.Dataset, either add it as a Dataset to maps or split it
            into a xarray.DataArrays per variable.
            Default to True.
        """
        # check for name when split_dataset is False
        if split_dataset is False and name is None:
            self.logger.error("name must be specified when split_dataset = False")

        # Call super method from HydroMT Core
        variables = super().setup_maps_from_rasterdataset(
            raster_fn=raster_fn,
            variables=variables,
            fill_method=fill_method,
            reproject_method=reproject_method,
            name=name,
            split_dataset=split_dataset,
        )

        for var in variables:
            da = self.maps[var]
            da.where(da == da.raster.nodata, -999.0)
            self.set_maps(da, var)

        allowed_methods = [
            "triangulation",
            "mean",
            "nearestNb",
            "max",
            "min",
            "invDist",
            "minAbs",
            "median",
        ]
        if not np.isin(interpolation_method, allowed_methods):
            raise ValueError(
                f"Interpolation method {interpolation_method} not allowed."
                f"Select from {allowed_methods}"
            )
        if not np.isin(locationtype, ["2d", "1d", "all"]):
            raise ValueError(
                f"Locationtype {locationtype} not allowed."
                "Select from ['2d', '1d', 'all']"
            )

        for var in variables:
            self.__set_map_parameters_based_on_variable(
                var, locationtype, interpolation_method
            )

    def setup_maps_from_raster_reclass(
        self,
        raster_fn: str,
        reclass_table_fn: str,
        reclass_variables: list,
        fill_method: Optional[str] = None,
        reproject_method: Optional[str] = "nearest",
        interpolation_method: Optional[str] = "triangulation",
        locationtype: Optional[str] = "2d",
        name: Optional[str] = None,
        split_dataset: Optional[bool] = True,
        **kwargs,
    ) -> None:
        """
        Add data variable(s) to maps by reclassifying values from ``raster_fn``.

        Reclassification is done bycombining values in ``raster_mapping_fn`` to
        spatial layer ``raster_fn``.

        The ``mapping_variables`` rasters are first created by mapping variables values
        from ``raster_mapping_fn`` to value in the ``raster_fn`` grid.

        Adds model layers:

        * **mapping_variables** maps: data from raster_mapping_fn spatially
            distributed with raster_fn

        Parameters
        ----------
        raster_fn: str
            Source name of raster data in data_catalog. Should be a DataArray. Else use
            **kwargs to select variables/time_tuple in
            hydromt.data_catalog.get_rasterdataset method
        reclass_table_fn: str
            Source name of mapping table of raster_fn in data_catalog. Make sure the
            data type is consistant for a ``reclass_variables`` including nodata.
            For example, for roughness, it is common that the data type is float,
            then use no data value as -999.0.
        reclass_variables: list
            List of mapping_variables from raster_mapping_fn table to add to mesh.
            Index column should match values in raster_fn.
            Available variables: ['elevtn', 'waterlevel', 'waterdepth', 'pet',
            'infiltcap', 'roughness_chezy', 'roughness_manning',
            'roughness_walllawnikuradse', 'roughness_whitecolebrook']
        fill_method : str, optional
            If specified, fills no data values using fill_nodata method. Available
            methods are {'linear', 'nearest', 'cubic', 'rio_idw'}.
        reproject_method : str, optional
            CRS reprojection method from rasterio.enums.Resampling. By default nearest.
            Available methods: ['nearest', 'bilinear', 'cubic', 'cubic_spline',
            'lanczos', 'average', 'mode', 'gauss', 'max', 'min', 'med', 'q1', 'q3',
            'sum', 'rms']
        interpolation_method : str, optional
            Interpolation method for DFlow-FM. By default triangulation. Except for
            waterlevel and waterdepth then the default is mean.
            When methods other than 'triangulation', the relative search cell size will
            be estimated based on resolution of the raster.
            Available methods: ['triangulation', 'mean', 'nearestNb', 'max', 'min',
            'invDist', 'minAbs', 'median']
        locationtype : str, optional
            LocationType in initial fields. Either 2d (default), 1d or all.
        name: str, optional
            Variable name, only in case data is of type DataArray or if a Dataset is
            added as is (split_dataset=False).
        split_dataset: bool, optional
            If data is a xarray.Dataset, either add it as is to maps or split it into
            several xarray.DataArrays.
            Default to True.
        """
        # check for name when split_dataset is False
        if split_dataset is False and name is None:
            self.logger.error("name must be specified when split_dataset = False")

        # Call super method
        reclass_variables = super().setup_maps_from_raster_reclass(
            raster_fn=raster_fn,
            reclass_table_fn=reclass_table_fn,
            reclass_variables=reclass_variables,
            fill_method=fill_method,
            reproject_method=reproject_method,
            name=name,
            split_dataset=split_dataset,
            **kwargs,
        )

        allowed_methods = [
            "triangulation",
            "mean",
            "nearestNb",
            "max",
            "min",
            "invDist",
            "minAbs",
            "median",
        ]
        if not np.isin(interpolation_method, allowed_methods):
            raise ValueError(
                f"Interpolation method {interpolation_method} not allowed."
                f"Select from {allowed_methods}"
            )
        if not np.isin(locationtype, ["2d", "1d", "all"]):
            raise ValueError(
                f"Locationtype {locationtype} not allowed."
                "Select from ['2d', '1d', 'all']"
            )
        for var in reclass_variables:
            self.__set_map_parameters_based_on_variable(
                var, locationtype, interpolation_method
            )

    def __set_map_parameters_based_on_variable(
        self, var: str, locationtype: str, interpolation_method: str
    ) -> None:
        """Set map parameters by updating user inputs to default self._MAP."""
        if var in self._MAPS:
            self._MAPS[var]["locationtype"] = locationtype
            self._MAPS[var]["interpolation"] = interpolation_method
            if interpolation_method != "triangulation":
                # adjust relative search cell size for averaging methods
                if self.maps[var].raster.res[0] > self.res:
                    relsize = np.round(
                        np.abs(self.maps[var].raster.res[0]) / self.res * np.sqrt(2)
                        + 0.05,
                        2,
                    )
                else:
                    relsize = 1.01
                self._MAPS[var]["averagingrelsize"] = relsize

    def setup_2dboundary(
        self,
        boundaries_fn: str = None,
        boundaries_timeseries_fn: str = None,
        boundary_value: float = 0.0,
        boundary_type: str = "waterlevel",
        tolerance: float = 3.0,
    ):
        """
        Prepare the 2D boundaries from line geometries.

        The values can either be a spatially-uniform constant using ``boundaries_fn``
        and ``boundary_value`` (default), or spatially-varying timeseries using
        ``boundaries_fn`` and  ``boundaries_timeseries_fn``
        The ``boundary_type`` can either be "waterlevel" or "discharge".

        If ``boundaries_timeseries_fn`` has missing values, the constant
        ``boundary_value`` will be used.

        The dataset/timeseries are clipped to the model region (see below note), and
        model time based on the model config tstart and tstop entries.

        Note that:
        (1) Only line geometry that are contained within the distance of
            ``tolenrance`` to grid cells are allowed.
        (2) Because of the above, this function must be called before the mesh
            refinement. #FIXME: check this after deciding on mesh refinement being
            a workflow or function
        (3) when using constant boundary, the output forcing will be written to time
            series with constant values.

        Adds/Updates model layers:

        * **boundary2d_{boundary_name}** forcing: 2D boundaries DataArray

        Parameters
        ----------
        boundaries_fn: str Path
            Path or data source name for line geometry file.

            * Required variables if a combined time series data csv file:
                ["boundary_id"] with type int
        boundaries_timeseries_fn: str, Path
            Path to tabulated timeseries csv file with time index in first column
            and location index with type int in the first row, matching "boundary_id"
            in ``boundaries_fn`.
            see :py:meth:`hydromt.get_dataframe`, for details.
            NOTE: Require equidistant time series
        boundary_value : float, optional
            Constant value to use for all boundaries, and to
            fill in missing data. By default 0.0 m.
        boundary_type : str, optional
            Type of boundary tu use. One of ["waterlevel", "discharge"].
            By default "waterlevel".
        tolerance: float, optional
            Search tolerance factor between boundary polyline and grid cells.
            Unit: in cell size units (i.e., not meters)
            By default, 3.0

        Raises
        ------
        AssertionError
            if "boundary_id" in "boundaries_fn" does not match the columns of
            ``boundaries_timeseries_fn``.

        """
        self.logger.info("Preparing 2D boundaries.")

        if boundary_type == "waterlevel":
            boundary_unit = "m"
        if boundary_type == "discharge":
            boundary_unit = "m3/s"

        _mesh = self.mesh_grids["mesh2d"]
        _mesh_region = gpd.GeoDataFrame(
            geometry=_mesh.to_shapely(dim=_mesh.face_dimension)
        ).unary_union
        _boundary_region = _mesh_region.buffer(tolerance * self.res).difference(
            _mesh_region
        )  # region where 2d boundary is allowed
        _boundary_region = gpd.GeoDataFrame(
            {"geometry": [_boundary_region]},
            crs=self.crs,
        )

        refdate, tstart, tstop = self.get_model_time()  # time slice

        # 1. read boundary geometries
        if boundaries_fn is not None:
            gdf_bnd = self.data_catalog.get_geodataframe(
                boundaries_fn,
                geom=_boundary_region,
                crs=self.crs,
                predicate="contains",
            )
            if len(gdf_bnd) == 0:
                self.logger.error(
                    "Boundaries are not found. Check if the boundary are outside of"
                    "recognisable boundary region (cell size * tolerance to the mesh)."
                )
            # preprocess
            gdf_bnd = gdf_bnd.explode(index_parts=True)
            # set index
            if "boundary_id" not in gdf_bnd:
                gdf_bnd["boundary_id"] = [
                    f"2dboundary_{i}" for i in range(len(gdf_bnd))
                ]
            else:
                gdf_bnd["boundary_id"] = gdf_bnd["boundary_id"].astype(str)
        else:
            gdf_bnd = None
        # 2. read timeseries boundaries
        if boundaries_timeseries_fn is not None:
            self.logger.info("reading timeseries boundaries")
            df_bnd = self.data_catalog.get_dataframe(
                boundaries_timeseries_fn, time_tuple=(tstart, tstop)
            )  # could not use open_geodataset due to line geometry
            # error if time mismatch or wrong parsing of dates
            if np.dtype(df_bnd.index).type != np.datetime64:
                raise ValueError(
                    "Dates in boundaries_timeseries_fn were not parsed correctly. "
                    "Update the source kwargs in the DataCatalog based on the driver"
                    "function arguments (eg pandas.read_csv for csv driver)."
                )
            if (df_bnd.index[-1] - df_bnd.index[0]) < (tstop - tstart):
                raise ValueError(
                    "Time in boundaries_timeseries_fn were shorter than model sim time."
                    "Update the source kwargs in the DataCatalog based on the driver"
                    "function arguments (eg pandas.read_csv for csv driver)."
                )
            if gdf_bnd is not None:
                # check if all boundary_id are in df_bnd
                assert all(
                    [bnd in df_bnd.columns for bnd in gdf_bnd["boundary_id"].unique()]
                ), "Not all boundary_id are in df_bnd"
        else:
            # default timeseries
            d_bnd = {bnd_id: np.nan for bnd_id in gdf_bnd["boundary_id"].unique()}
            d_bnd.update(
                {
                    "time": pd.date_range(
                        start=pd.to_datetime(tstart),
                        end=pd.to_datetime(tstop),
                        freq="D",
                    )
                }
            )
            df_bnd = pd.DataFrame(d_bnd).set_index("time")

        # 4. Derive DataArray with boundary values at boundary locations
        # in boundaries_branch_type
        da_out_dict = workflows.compute_2dboundary_values(
            boundaries=gdf_bnd,
            df_bnd=df_bnd.reset_index(),
            boundary_value=boundary_value,
            boundary_type=boundary_type,
            boundary_unit=boundary_unit,
            logger=self.logger,
        )

        # 5. set boundaries
        for da_out_name, da_out in da_out_dict.items():
            self.set_forcing(da_out, name=f"boundary2d_{da_out_name}")

        # adjust parameters
        self.set_config("geometry.openboundarytolerance", tolerance)

    def setup_rainfall_from_constant(
        self,
        constant_value: float,
    ):
        """
        Prepare constant 2D daily rainfall_rate timeseries based on ``constant_value``.

        Adds/Updates model layers:

        * **meteo_{meteo_type}** forcing: DataArray

        Parameters
        ----------
        constant_value: float
            Constant value for the rainfall_rate timeseries in mm/day.
        """
        self.logger.info("Preparing rainfall meteo forcing from uniform timeseries.")

        refdate, tstart, tstop = self.get_model_time()  # time slice
        meteo_location = (
            self.region.centroid.x,
            self.region.centroid.y,
        )  # global station location

        df_meteo = pd.DataFrame(
            {
                "time": pd.date_range(
                    start=pd.to_datetime(tstart), end=pd.to_datetime(tstop), freq="D"
                ),
                "precip": constant_value,
            }
        )

        # 3. Derive DataArray with meteo values
        da_out = workflows.compute_meteo_forcings(
            df_meteo=df_meteo,
            fill_value=constant_value,
            is_rate=True,
            meteo_location=meteo_location,
            logger=self.logger,
        )

        # 4. set meteo forcing
        self.set_forcing(da_out, name=f"meteo_{da_out.name}")

        # 5. set meteo in mdu
        self.set_config("external_forcing.rainfall", 1)

    def setup_rainfall_from_uniform_timeseries(
        self,
        meteo_timeseries_fn: Union[str, Path],
        fill_value: float = 0.0,
        is_rate: bool = True,
    ):
        """
        Prepare spatially uniform 2D rainfall forcings from ``meteo_timeseries_fn``.

        For now only support global  (spatially uniform) timeseries.

        If ``meteo_timeseries_fn`` has missing values or shorter than model simulation
        time, the constant ``fill_value`` will be used, e.g. 0.

        The dataset/timeseries are clipped to the model time based on the model config
        tstart and tstop entries.

        Adds/Updates model layers:
            * **meteo_{meteo_type}** forcing: DataArray

        Parameters
        ----------
        meteo_timeseries_fn: str, Path
            Path or data source name to tabulated timeseries csv file with time index
            in first column.

            * Required variables : ['precip']

            see :py:meth:`hydromt.get_dataframe`, for details.
            NOTE:
            Require equidistant time series
            If ``is_rate`` = True, unit is expected to be in mm/day and else mm.
        fill_value : float, optional
            Constant value to use to fill in missing data. By default 0.
        is_rate : bool, optional
            Specify if the type of meteo data is direct "rainfall" (False) or
            "rainfall_rate" (True).
            By default True for "rainfall_rate".
            Note that Delft3DFM 1D2D Suite 2022.04 supports only "rainfall_rate".

        """
        self.logger.info("Preparing rainfall meteo forcing from uniform timeseries.")

        refdate, tstart, tstop = self.get_model_time()  # time slice
        meteo_location = (
            self.region.centroid.x,
            self.region.centroid.y,
        )  # global station location

        # get meteo timeseries
        df_meteo = self.data_catalog.get_dataframe(
            meteo_timeseries_fn, variables=["precip"], time_tuple=(tstart, tstop)
        )
        # error if time mismatch or wrong parsing of dates
        if np.dtype(df_meteo.index).type != np.datetime64:
            raise ValueError(
                "Dates in meteo_timeseries_fn were not parsed correctly. "
                "Update the source kwargs in the DataCatalog based on the driver"
                "function arguments (eg pandas.read_csv for csv driver)."
            )
        if (df_meteo.index[-1] - df_meteo.index[0]) < (tstop - tstart):
            self.logger.warning(
                "Time in meteo_timeseries_fn were shorter than model simulation time. "
                "Will fill in using fill_value."
            )
            dt = df_meteo.index[1] - df_meteo.index[0]
            t_index = pd.DatetimeIndex(pd.date_range(start=tstart, end=tstop, freq=dt))
            df_meteo = df_meteo.reindex(t_index).fillna(fill_value)
        df_meteo["time"] = df_meteo.index

        # 3. Derive DataArray with meteo values
        da_out = workflows.compute_meteo_forcings(
            df_meteo=df_meteo,
            fill_value=fill_value,
            is_rate=is_rate,
            meteo_location=meteo_location,
            logger=self.logger,
        )

        # 4. set meteo forcing
        self.set_forcing(da_out, name=f"meteo_{da_out.name}")

        # 5. set meteo in mdu
        self.set_config("external_forcing.rainfall", 1)

    # ## I/O
    def read(self):
        """
        Read the complete model schematization and configuration from file.

        # FIXME: where to read crs?.
        """
        self.logger.info(f"Reading model data from {self.root}")
        self.read_dimr()
        self.read_config()
        self.read_mesh()
        self.read_maps()
        self.read_geoms()  # needs mesh so should be done after
        self.read_forcing()
        self._check_crs()

    def write(self):  # complete model
        """Write the complete model schematization and configuration to file."""
        self.logger.info(f"Writing model data to {self.root}")
        # if in r, r+ mode, only write updated components
        if not self._write:
            self.logger.warning("Cannot write in read-only mode")
            return

        if self._maps:
            self.write_maps()
        if self._geoms:
            self.write_geoms()
        if self._mesh is not None or not self.branches.empty:
            self.write_mesh()
        if self._forcing:
            self.write_forcing()
        if self.config:  # dflowfm config, should always be last!
            self.write_config()
        if self.dimr:  # dimr config, should always be last after dflowfm config!
            self.write_dimr()
        self.write_data_catalog()

    def read_config(self) -> None:
        """Use Hydrolib-core reader and return to dictionnary."""
        # Read via init_dfmmodel
        if self._dfmmodel is None:
            self.init_dfmmodel()
        # Convert to full dictionnary without hydrolib-core objects
        cf_dict = dict()
        for k, v in self._dfmmodel.__dict__.items():
            if v is None or k == "filepath":
                cf_dict[k] = v
            else:
                ci_dict = dict()
                for ki, vi in v.__dict__.items():
                    if ki == "frictfile" and isinstance(vi, list):  # list of filepath
                        ci_dict[ki] = ";".join([str(vj.filepath) for vj in vi])
                    elif ki != "comments":
                        if hasattr(vi, "filepath"):
                            # need to change the filepath object to path
                            ci_dict[ki] = vi.filepath
                        else:
                            ci_dict[ki] = vi
                cf_dict[k] = ci_dict
        self._config = cf_dict

    def write_config(self) -> None:
        """From config dict to Hydrolib MDU."""
        # Not sure if this is worth it compared to just calling write_config from super
        # advantage is the validator but the whole model is then read
        # when initialising FMModel
        self._assert_write_mode()

        cf_dict = self._config.copy()
        # Need to switch to dflowfm folder for files to be found and properly added
        mdu_fn = cf_dict.pop("filepath", None)
        mdu_fn = Path(join(self.root, self._config_fn))
        cwd = os.getcwd()
        os.chdir(dirname(mdu_fn))
        mdu = FMModel(**cf_dict)
        # add filepath
        mdu.filepath = mdu_fn
        # temporarily remove sediment section to avoid error in Delft3D FM 1D2D 2024.03
        # https://issuetracker.deltares.nl/browse/FM1D2D-3047
        del mdu.sediment
        # write
        mdu.save(recurse=False)
        # Go back to working dir
        os.chdir(cwd)

    def read_maps(self) -> Dict[str, Union[xr.Dataset, xr.DataArray]]:
        """Read maps from initialfield and parse to dict of xr.DataArray."""
        self._assert_read_mode()
        # Read initial fields
        inifield_model = self.dfmmodel.geometry.inifieldfile
        if inifield_model:
            # seperate 1d and 2d
            # inifield_model_1d = [
            #     i for i in inifield_model.initial if "1d" in i.locationtype
            # ] # not supported yet
            inifield_model_2dinitial = [
                i for i in inifield_model.initial if "2d" in i.locationtype
            ]
            inifield_model_2dparameter = [
                i for i in inifield_model.parameter if "2d" in i.locationtype
            ]
            inifield_model_2d = inifield_model_2dinitial + inifield_model_2dparameter
        else:
            inifield_model_2d = []

        if any(inifield_model_2d):
            # Loop over initial / parameter to read the geotif
            inilist = inifield_model_2d

            if len(inilist) > 0:
                # DFM map names
                rm_dict = dict()
                for v in self._MAPS:
                    rm_dict[self._MAPS[v]["name"]] = v
                for inidict in inilist:
                    _fn = inidict.datafile.filepath
                    # Bug: when initialising IniFieldModel hydrolib-core
                    # does not parse correclty the relative path
                    # For now re-update manually....
                    if not isfile(_fn):
                        _fn = join(self.root, "maps", _fn.name)
                    inimap = hydromt.io.open_raster(_fn)
                    name = inidict.quantity
                    # Need to get branchid from config
                    if name == "frictioncoefficient":
                        frictype = self.get_config("physics.uniffricttype", fallback=1)
                        fricname = [
                            n
                            for n in self._MAPS
                            if self._MAPS[n].get("frictype", None) == frictype
                        ]
                        rm_dict[name] = fricname[0]
                    # Check if name in self._MAPS to update properties
                    if name in rm_dict:
                        # update all keywords
                        if "comments" in inidict.__dict__:
                            inidict.__dict__.pop("comments")
                        self._MAPS[rm_dict[name]].update(inidict)
                        # Update default interpolation method
                        if inidict.interpolationmethod == "averaging":
                            interpmethod = inidict.averagingtype
                        else:
                            interpmethod = inidict.interpolationmethod
                        self._MAPS[rm_dict[name]]["interpolation"] = interpmethod
                        # Rename to HydroMT name
                        name = rm_dict[name]
                    # Add to maps
                    inimap.name = name
                    self.set_maps(inimap, name)

    def write_maps(self) -> None:
        """Write maps as tif files in maps folder and update initial fields."""
        if len(self._maps) == 0:
            self.logger.debug("No maps data found, skip writing.")
            return
        self._assert_write_mode()
        # Global parameters
        mapsroot = join(self.root, "maps")
        inilist = []
        paramlist = []
        self.logger.info(f"Writing maps files to {mapsroot}")

        def _prepare_inifields(da_dict, da):
            # Write tif files
            name = da_dict["name"]
            type = da_dict["initype"]
            interp_method = da_dict["interpolation"]
            locationtype = da_dict["locationtype"]
            _fn = join(mapsroot, f"{name}.tif")
            if da.raster.nodata is None or np.isnan(da.raster.nodata):
                da.raster.set_nodata(-999)
            da.raster.to_raster(_fn)
            self.logger.info(f"Writing file {mapsroot}/{name}.tif")
            # Prepare dict
            if interp_method == "triangulation":
                inidict = {
                    "quantity": name,
                    "dataFile": f"../maps/{name}.tif",
                    "dataFileType": "GeoTIFF",
                    "interpolationMethod": interp_method,
                    "operand": da_dict.get("oprand", "O"),
                    "locationType": locationtype,
                }
            else:
                inidict = {
                    "quantity": name,
                    "dataFile": f"../maps/{name}.tif",
                    "dataFileType": "GeoTIFF",
                    "interpolationMethod": "averaging",
                    "operand": da_dict.get("oprand", "O"),
                    "averagingType": interp_method,
                    "averagingRelSize": da_dict.get("averagingrelsize"),
                    "locationType": locationtype,
                }
            if type == "initial":
                inilist.append(inidict)
            elif type == "parameter":
                paramlist.append(inidict)

        # Only write maps that are listed in self._MAPS, rename tif on the fly
        # TODO raise value error if both waterdepth and waterlevel are given in maps
        for name, ds in self._maps.items():
            if isinstance(ds, xr.DataArray):
                if name in self._MAPS:
                    _prepare_inifields(self._MAPS[name], ds)
                    # update config if friction
                    if "frictype" in self._MAPS[name]:
                        self.set_config(
                            "physics.uniffricttype", self._MAPS[name]["frictype"]
                        )
                    # update config if infiltration
                    if name == "infiltcap":
                        self.set_config("grw.infiltrationmodel", 2)
                else:
                    self.logger.error(
                        f"Could not write map to model: {name} not recognized"
                    )
            elif isinstance(ds, xr.Dataset):
                for v in ds.data_vars:
                    if v in self._MAPS:
                        _prepare_inifields(self._MAPS[v], ds[v])
                        # update config if frcition
                        if self._MAPS[v] == "frictype":
                            self.set_config(
                                "physics.uniffricttype", self._MAPS[name]["frictype"]
                            )
                        # update config if infiltration
                        if v == "infiltcap":
                            self.set_config("grw.infiltrationmodel", 2)
                    else:
                        self.logger.error(
                            f"Could not write map to model: {v} not found in map {name}"
                        )
        # Assign initial fields to model and write
        inifield_model = IniFieldModel(initial=inilist, parameter=paramlist)
        # Bug: when initialising IniFieldModel hydrolib-core does not parse correclty
        # the relative path
        # For now re-update manually....
        for i in range(len(inifield_model.initial)):
            path = Path(f"../maps/{inifield_model.initial[i].datafile.filepath.name}")
            inifield_model.initial[i].datafile.filepath = path
        for i in range(len(inifield_model.parameter)):
            path = Path(f"../maps/{inifield_model.parameter[i].datafile.filepath.name}")
            inifield_model.parameter[i].datafile.filepath = path
        # Write inifield file
        inifield_model_filename = inifield_model._filename() + ".ini"
        fm_dir = dirname(join(self.root, self._config_fn))
        inifield_model.save(
            join(fm_dir, inifield_model_filename),
            recurse=False,
        )
        # save filepath in the config
        self.set_config("geometry.inifieldfile", inifield_model_filename)

    def read_geoms(self) -> None:  # FIXME: gives an error when only 2D model.
        """
        Read model geometries files at <root>/<geoms> and add to geoms property.

        For branches / boundaries etc... the reading of hydrolib-core objects happens
        in read_mesh. There the geoms geojson copies are re-set based on dflowfm files
        content.
        """
        self._assert_read_mode()
        super().read_geoms(fn="geoms/region.geojson")

        if self.dfmmodel.geometry.crosslocfile is not None:
            # Read cross-sections and friction
            # Add crosssections properties, should be done before friction
            # Branches are needed do derive locations,
            # self.branches should start the read if not done yet
            self.logger.info("Reading cross-sections files")
            crosssections = utils.read_crosssections(self.branches, self.dfmmodel)

            # Add friction properties from roughness files
            # self.logger.info("Reading friction files")
            crosssections = utils.read_friction(crosssections, self.dfmmodel)
            self.set_geoms(crosssections, "crosssections")

        # Read manholes
        if self.dfmmodel.geometry.storagenodefile is not None:
            self.logger.info("Reading manholes file")
            network1d_nodes = mesh_utils.network1d_nodes_geodataframe(
                self.mesh_datasets["network1d"]
            )
            manholes = utils.read_manholes(network1d_nodes, self.dfmmodel)
            self.set_geoms(manholes, "manholes")

        # Read structures
        if self.dfmmodel.geometry.structurefile is not None:
            self.logger.info("Reading structures file")
            structures = utils.read_structures(self.branches, self.dfmmodel)
            for st_type in structures["type"].unique():
                self.set_geoms(structures[structures["type"] == st_type], f"{st_type}s")

    def write_geoms(self, write_mesh_gdf=True) -> None:
        """Write model geometries to a GeoJSON file at <root>/<geoms>."""
        self._assert_write_mode()

        # Optional: also write mesh_gdf object
        if write_mesh_gdf:
            for name, gdf in self.mesh_gdf.items():
                self.set_geoms(gdf, name)

        # Write geojson equivalent of all objects.
        # Note that these files are not directly used when updating the model
        super().write_geoms(fn="geoms/{name}.geojson")

        # Write dfm files
        savedir = dirname(join(self.root, self._config_fn))

        # Write cross-sections (inc. friction)
        if "crosssections" in self._geoms:
            # Crosssections
            gdf_crs = self.geoms["crosssections"]
            self.logger.info("Writting cross-sections files crsdef and crsloc")
            crsdef_fn, crsloc_fn = utils.write_crosssections(gdf_crs, savedir)
            self.set_config("geometry.crossdeffile", crsdef_fn)
            self.set_config("geometry.crosslocfile", crsloc_fn)

            # Friction
            self.logger.info("Writting friction file(s)")
            friction_fns = utils.write_friction(gdf_crs, savedir)
            self.set_config("geometry.frictfile", ";".join(friction_fns))

        # Write structures
        # Manholes
        if "manholes" in self._geoms:
            self.logger.info("Writting manholes file.")
            storage_fn = utils.write_manholes(
                self.geoms["manholes"],
                savedir,
            )
            self.set_config("geometry.storagenodefile", storage_fn)

        # Write structures
        existing_structures = [st for st in ["bridges", "culverts"] if st in self.geoms]
        if len(existing_structures) > 0:
            # combine all structures
            structures = []
            for st in existing_structures:
                structures.append(self.geoms.get(st).to_dict("records"))
            structures = list(itertools.chain.from_iterable(structures))
            structures = pd.DataFrame(structures).replace(np.nan, None)
            # write
            self.logger.info("Writting structures file.")
            structures_fn = utils.write_structures(
                structures,
                savedir,
            )
            self.set_config("geometry.structurefile", structures_fn)

        # write hydromt
        # Optional: also write mesh_gdf object
        if write_mesh_gdf:
            for name, gdf in self.mesh_gdf.items():
                self.set_geoms(gdf, name)

        # Write geojson equivalent of all objects.
        # NOTE these files are not used for model update.
        # convert any list in geoms to strings
        def convert_lists_to_strings(df):
            for column_name in df.columns:
                if df[column_name].apply(isinstance, args=(list,)).any():
                    df[column_name] = df[column_name].apply(
                        lambda x: " ".join(f"{x}") if isinstance(x, list) else x
                    )
            return df

        for name in self.geoms:
            self.set_geoms(convert_lists_to_strings(self.geoms[name]), name)

        super().write_geoms(fn="geoms/{name}.geojson")

    def read_forcing(
        self,
    ) -> None:  # FIXME reading of forcing should include boundary, lateral and meteo
        """Read forcing at <root/?/> and parse to dict of xr.DataArray."""
        self._assert_read_mode()
        # Read external forcing
        ext_model = self.dfmmodel.external_forcing.extforcefilenew
        if ext_model is not None:
            # boundary
            if len(ext_model.boundary) > 0:
                df_ext = pd.DataFrame([f.__dict__ for f in ext_model.boundary])
                # 1d boundary
                df_ext_1d = df_ext.loc[~df_ext.nodeid.isna(), :]
                if len(df_ext_1d) > 0:
                    # Forcing data arrays to prepare for each quantity
                    forcing_names = np.unique(df_ext_1d.quantity).tolist()
                    # Loop over forcing names to build data arrays
                    for name in forcing_names:
                        # Get the dataframe corresponding to the current variable
                        df = df_ext_1d[df_ext_1d.quantity == name]
                        # Get the corresponding nodes gdf
                        network1d_nodes = mesh_utils.network1d_nodes_geodataframe(
                            self.mesh_datasets["network1d"]
                        )
                        node_geoms = network1d_nodes[
                            np.isin(network1d_nodes["nodeid"], df.nodeid.values)
                        ]
                        da_out = utils.read_1dboundary(
                            df, quantity=name, nodes=node_geoms
                        )
                        # Add to forcing
                        self.set_forcing(da_out)
                # 2d boundary
                df_ext_2d = df_ext.loc[df_ext.nodeid.isna(), :]
                if len(df_ext_2d) > 0:
                    for _, df in df_ext_2d.iterrows():
                        da_out = utils.read_2dboundary(
                            df, workdir=self.dfmmodel.filepath.parent
                        )
                        # Add to forcing
                        self.set_forcing(da_out)
            # lateral
            if len(ext_model.lateral) > 0:
                df_ext = pd.DataFrame([f.__dict__ for f in ext_model.lateral])
                da_out = utils.read_1dlateral(
                    df_ext, branches=self.branches
                )  # TODO extend support to get laterals on nodes #78
                # Add to forcing
                self.set_forcing(da_out)
            # meteo
            if len(ext_model.meteo) > 0:
                df_ext = pd.DataFrame([f.__dict__ for f in ext_model.meteo])
                # Forcing dataarrays to prepare for each quantity
                forcing_names = np.unique(df_ext.quantity).tolist()
                # Loop over forcing names to build data arrays
                for name in forcing_names:
                    # Get the dataframe corresponding to the current variable
                    df = df_ext[df_ext.quantity == name]
                    da_out = utils.read_meteo(df, quantity=name)
                    # Add to forcing
                    self.set_forcing(da_out)
            # TODO lateral

    def write_forcing(self) -> None:
        """Write forcing into hydrolib-core ext and forcing models."""
        if len(self._forcing) == 0:
            self.logger.debug("No forcing data found, skip writing.")
        else:
            self._assert_write_mode()
            self.logger.info("Writting forcing files.")
            savedir = dirname(join(self.root, self._config_fn))
            # create new external forcing file
            ext_fn = "bnd.ext"
            Path(join(savedir, ext_fn)).unlink(missing_ok=True)
            # populate external forcing file
            utils.write_1dboundary(self.forcing, savedir, ext_fn=ext_fn)
            utils.write_2dboundary(self.forcing, savedir, ext_fn=ext_fn)
            utils.write_1dlateral(self.forcing, savedir, ext_fn=ext_fn)
            utils.write_meteo(self.forcing, savedir, ext_fn=ext_fn)
            self.set_config("external_forcing.extforcefilenew", ext_fn)

    def read_mesh(self):
        """Read network file with Hydrolib-core and extract mesh/branches info."""
        self._assert_read_mode()

        # Read mesh
        # hydrolib-core convention
        network = self.dfmmodel.geometry.netfile.network
        # FIXME: crs info is not available in dfmmodel, so get it from region.geojson
        # Cannot use read_geoms yet because for some some geoms
        # (crosssections, manholes) mesh needs to be read first...
        region_fn = join(self.root, "geoms", "region.geojson")
        if (not self._crs) and isfile(region_fn):
            crs = gpd.read_file(region_fn).crs
            self._crs = crs

        crs = self.crs

        # convert to xugrid
        mesh = mesh_utils.mesh_from_hydrolib_network(network, crs=crs)
        # set mesh
        self._mesh = mesh

        # update resolution
        if "mesh2d" in self.mesh_names:
            if self._res is None:
                self._res = np.max(np.diff(self.mesh_grids["mesh2d"].node_x))

        # creates branches geometry from network1d
        if "network1d" in self.mesh_names:
            network1d_dataset = self.mesh_datasets["network1d"]
            # Create the branches GeoDataFrame (from geom)
            # network1d_geometry = self.mesh_gdf["network1d"] this returns the network
            branches = mesh_utils.network1d_geoms_geodataframe(network1d_dataset)
            # branches["branchtype"] = network1d_dataset["network1d_branch_type"]
            # might support in the future
            # https://github.com/Deltares/HYDROLIB-core/issues/561

            # Add branchtype, properties from branches.gui file
            self.logger.info("Reading branches GUI file")
            branches = utils.read_branches_gui(branches, self.dfmmodel)

            # Set branches
            self.set_branches(branches)

    def write_mesh(self, write_gui=True):
        """Write 1D branches and 2D mesh at <root/dflowfm/fm_net.nc>."""
        self._assert_write_mode()
        savedir = join(self.root, "dflowfm")
        mesh_filename = "fm_net.nc"

        # write mesh
        # HydroMT convention - FIXME hydrolib does not seem to read the 1D and links
        # part of the mesh
        # super().write_mesh(fn=join(savedir, mesh_filename))

        # write with hydrolib-core
        # Note: hydrolib-core writes more information including attributes and
        # converts some variables using start_index
        # FIXME: does not write crs that is recongnised by Delft3D FM GUI.
        # check dfm_tools/meshkernel_helpers.py#L82

        network = mesh_utils.hydrolib_network_from_mesh(self.mesh)
        network.to_file(Path(join(savedir, mesh_filename)))

        # save relative path to mdu
        self.set_config("geometry.netfile", mesh_filename)

        # other mesh1d related geometry TODO update
        if "mesh1d" in self.mesh_names and write_gui:
            self.logger.info("Writting branches.gui file")
            if "manholes" in self.geoms:
                utils.write_branches_gui(self.branches, savedir)

    def read_states(self):
        """Read states at <root/?/> and parse to dict of xr.DataArray."""
        return self._states
        # raise NotImplementedError()

    def write_states(self):
        """Write states at <root/?/> in model ready format."""
        pass
        # raise NotImplementedError()

    def read_results(self):
        """Read results at <root/?/> and parse to dict of xr.DataArray."""
        return self._results
        # raise NotImplementedError()

    def write_results(self):
        """Write results at <root/?/> in model ready format."""
        pass
        # raise NotImplementedError()

    @property
    def crs(self):
        """Return model crs."""
        # return pyproj.CRS.from_epsg(self.get_config("global.epsg", fallback=4326))
        if self._crs is None:
            # try to read it from mesh usingMeshModel method
            self._crs = super().crs
        return self._crs

    @property
    def bounds(self) -> Tuple:
        """Return model mesh bounds."""
        return self.region.total_bounds

    @property
    def region(self) -> gpd.GeoDataFrame:
        """Return geometry of region of the model area of interest."""
        # First tries in geoms
        if "region" in self.geoms:
            region = self.geoms["region"]
        # Else derives from mesh or branches
        else:
            if self.mesh is not None:
                bounds = self.mesh.ugrid.total_bounds
                crs = self.crs
            elif not self.branches.empty:
                bounds = self.branches.total_bounds
                crs = self.branches.crs
            else:
                # Finally raise error assuming model is empty
                raise ValueError(
                    "Could not derive region from geoms, or mesh. Model may be empty."
                )
            region = gpd.GeoDataFrame(geometry=[box(*bounds)], crs=crs)
            self.set_geoms(region, "region")

        return region

    @property
    def dfmmodel(self):
        """Hydrolib-core FMModel object."""
        if self._dfmmodel is None:
            self.init_dfmmodel()
        return self._dfmmodel

    def init_dfmmodel(self):
        """Initialise the hydrolib-core FMModel object."""
        # create a new MDU-Model
        mdu_fn = Path(join(self.root, self._config_fn))
        if isfile(mdu_fn) and self._read:
            self.logger.info(f"Reading mdu file at {mdu_fn}")
            self._dfmmodel = FMModel(filepath=mdu_fn)
        else:  # use hydrolib template
            self._assert_write_mode()
            self.logger.info("Initialising empty mdu file")
            self._dfmmodel = FMModel()
            self._dfmmodel.filepath = mdu_fn

    @property
    def dimr(self):
        """DIMR file object."""
        if not self._dimr:
            self.read_dimr()
        return self._dimr

    def read_dimr(self, dimr_fn: Optional[str] = None) -> None:
        """Read DIMR from file and else create from hydrolib-core."""
        if dimr_fn is None:
            dimr_fn = join(self.root, self._dimr_fn)
        # if file exist, read
        if isfile(dimr_fn) and self._read:
            self.logger.info(f"Reading dimr file at {dimr_fn}")
            dimr = DIMR(filepath=Path(dimr_fn))
        # else initialise
        else:
            self._assert_write_mode()
            self.logger.info("Initialising empty dimr file")
            dimr = DIMR()
        self._dimr = dimr

    def write_dimr(self, dimr_fn: Optional[str] = None):
        """Write the dmir file.

        In write mode, updates first the FMModel component.
        """
        # force read
        self.dimr
        if dimr_fn is not None:
            self._dimr.filepath = join(self.root, dimr_fn)
        else:
            self._dimr.filepath = join(self.root, self._dimr_fn)

        if not self._read:
            # Updates the dimr file first before writing
            self.logger.info("Adding dflowfm component to dimr config")

            # update component
            components = self._dimr.component
            if len(components) != 0:
                components = []
            fmcomponent = FMComponent(
                name="dflowfm",
                workingdir="dflowfm",
                inputfile=basename(self._config_fn),
                model=self.dfmmodel,
            )
            components.append(fmcomponent)
            self._dimr.component = components
            # update control
            controls = self._dimr.control
            if len(controls) != 0:
                controls = []
            control = Start(name="dflowfm")
            controls.append(control)
            self._dimr.control = control

        # write
        self.logger.info(f"Writing model dimr file to {self._dimr.filepath}")
        self.dimr.save(recurse=False)

    @property
    def branches(self):
        """
        Return the branches (gpd.GeoDataFrame object) representing the 1D network.

        Contains several "branchtype" for : channel, river, pipe, tunnel.
        """
        if self._branches is None and self._read:
            self.read_mesh()
        if self._branches is None:
            self._branches = gpd.GeoDataFrame()
        return self._branches

    def set_branches(self, branches: gpd.GeoDataFrame):
        """Update the branches object as well as the linked geoms."""
        # Check if "branchtype" col in new branches
        if "branchtype" in branches.columns:
            self._branches = branches
        else:
            self.logger.error(
                "'branchtype' column absent from the new branches, could not update."
            )

        # Update channels/pipes in geoms
        _ = self.set_branches_component(name="river")
        _ = self.set_branches_component(name="channel")
        _ = self.set_branches_component(name="pipe")

        # update geom
        self.logger.debug("Adding branches vector to geoms.")
        self.set_geoms(branches, "branches")

        self.logger.debug("Updating branches in network.")

    def set_branches_component(self, name: str):
        """Extract component name from branches and add it to geoms."""
        gdf_comp = self.branches[self.branches["branchtype"] == name]
        if gdf_comp.index.size > 0:
            self.set_geoms(gdf_comp, name=f"{name}s")
        return gdf_comp

    @property
    def rivers(self):
        """Extract rivers from branches."""
        if "rivers" in self.geoms:
            gdf = self.geoms["rivers"]
        else:
            gdf = self.set_branches_component("river")
        return gdf

    @property
    def channels(self):
        """Extract channels from branches."""
        if "channels" in self.geoms:
            gdf = self.geoms["channels"]
        else:
            gdf = self.set_branches_component("channel")
        return gdf

    @property
    def pipes(self):
        """Extract pipes from branches."""
        if "pipes" in self.geoms:
            gdf = self.geoms["pipes"]
        else:
            gdf = self.set_branches_component("pipe")
        return gdf

    @property
    def opensystem(self):
        """Open system branches (river, channel)."""
        if len(self.branches) > 0:
            gdf = self.branches[self.branches["branchtype"].isin(["river", "channel"])]
        else:
            gdf = gpd.GeoDataFrame()
        return gdf

    @property
    def closedsystem(self):
        """Closed system branches (pipe, tunnel)."""
        if len(self.branches) > 0:
            gdf = self.branches[self.branches["branchtype"].isin(["pipe", "tunnel"])]
        else:
            gdf = gpd.GeoDataFrame()
        return gdf

    @property
    def boundaries(self):
        """1D boundary locations."""
        if "boundaries" not in self.geoms:
            self.set_geoms(
                workflows.get_boundaries_with_nodeid(
                    self.branches,
                    mesh_utils.network1d_nodes_geodataframe(
                        self.mesh_datasets["network1d"]
                    ),
                ),
                "boundaries",
            )
        return self.geoms["boundaries"]

    def get_model_time(self):
        """
        Return (refdate, tstart, tstop) tuple.

        It is parsed from model reference datem start and end time.
        """
        refdate = datetime.strptime(str(self.get_config("time.refdate")), "%Y%m%d")
        tstart = refdate + timedelta(seconds=float(self.get_config("time.tstart")))
        tstop = refdate + timedelta(seconds=float(self.get_config("time.tstop")))
        return refdate, tstart, tstop

    @property
    def res(self):
        """Resolution of the mesh2d."""
        if self._res is not None:
            return self._res

    def set_mesh(
        self,
        data: Union[xu.UgridDataArray, xu.UgridDataset],
        name: Optional[str] = None,
        grid_name: Optional[str] = None,
        overwrite_grid: Optional[bool] = False,
    ) -> None:
        """Add data to mesh.

        All layers of mesh have identical spatial coordinates in Ugrid conventions.

        Parameters
        ----------
        data: xugrid.UgridDataArray or xugrid.UgridDataset
            new layer to add to mesh
        name: str, optional
            Name of new object layer, this is used to overwrite the name of
            a UgridDataArray.
        grid_name: str, optional
            Name of the mesh grid to add data to. If None, inferred from data.
            Can be used for renaming the grid.
        overwrite_grid: bool, optional
            If True, overwrite the grid with the same name as the grid in self.mesh.
        """
        # Check if new grid_name
        if grid_name not in self.mesh_names:
            new_grid = True
        else:
            new_grid = False

        # First call super method to add mesh data
        super().set_mesh(
            data=data,
            name=name,
            grid_name=grid_name,
            overwrite_grid=overwrite_grid,
        )

        # check if 1D and 2D and and 1D2D links and overwrite
        # then send warning that setup_link1d2d should be run again
        if overwrite_grid and "link1d2d" in self.mesh.data_vars:
            if grid_name == "mesh1d" or grid_name == "mesh2d":
                # TODO check if warning is enough or if we should remove to be sure?
                self.logger.warning(
                    f"{grid_name} grid was updated in self.mesh. "
                    "Re-run setup_link1d2d method to update the model 1D2D links."
                )

        # update related geoms if necessary: region - boundaries
        # the region is done in HydroMT Core
        if overwrite_grid or new_grid:
            # 1D boundaries
            if grid_name == "mesh1d":
                self.set_geoms(
                    workflows.get_boundaries_with_nodeid(
                        self.branches,
                        mesh_utils.network1d_nodes_geodataframe(
                            self.mesh_datasets["network1d"]
                        ),
                    ),
                    "boundaries",
                )

    def set_link1d2d(
        self,
        link1d2d: xr.Dataset,
    ):
        """
        Add or replace the link1d2d in the model mesh.

        Parameters
        ----------
        link1d2d: xr.Dataset
            link1d2d dataset with variables: [link1d2d, link1d2d_ids,
            link1d2d_long_names, link1d2d_contact_type]
        """
        # Check if link1d2d already in self.mesh
        # FIXME current implementation of below does not support updating partial
        # 1d2d links. Either document or adapt. #1
        if "link1d2d" in self.mesh.data_vars:
            self.logger.info("Overwriting existing link1d2d in self.mesh.")
            self._mesh = self._mesh.drop_vars(
                [
                    "link1d2d",
                    "link1d2d_id",
                    "link1d2d_long_name",
                    "link1d2d_contact_type",
                ]
            )

        # Add link1d2d to mesh
        self._mesh = self._mesh.merge(link1d2d)

    def _model_has_2d(self):
        """Check if model has 2D mesh part."""
        if "mesh2d" in self.mesh_names:
            return True
        else:
            return False

    def _model_has_1d(self):
        """Check if model has 1D mesh part."""
        if "mesh1d" in self.mesh_names:
            return True
        else:
            return False

    def _check_crs(self):
        """Check if model crs is defined."""
        if self.crs is None:
            if self._read:
                self.logger.warning(
                    "Could not derive CRS from reading the mesh file."
                    "Please define the CRS in the [global] init attributes before"
                    "setting up the model."
                )
            else:
                raise ValueError(
                    "CRS is not defined. Please define the CRS in the [global] init"
                    "attributes before setting up the model."
                )
        else:
            self.logger.info(f"project crs: {self.crs.to_epsg()}")
