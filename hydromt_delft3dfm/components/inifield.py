"""Delft3D FM inifield component."""

import logging
from os.path import dirname, isfile, join
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from hydrolib.core.dflowfm import IniFieldModel
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import SpatialDatasetsComponent
from hydromt.readers import open_raster

__all__ = ["IniFieldComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class IniFieldComponent(SpatialDatasetsComponent):
    """
    Manage the Delft3D-FM ini field files for model initial fields.

    It extends the base SpatialDatasetsComponent from hydromt and consists of a
    dictionary of xarray Datasets.
    """

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

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "maps/{name}.tif",
        region_component: str | None = None,
    ):
        """Initialize IniFieldComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            by default "maps/{name}.tif", i.e. one file per dataset in
            the data dictionary.
        region_component : str, optional
            The name of the region component to use as reference for this component's
            region. If None, the region will be set to the union of all datasets in
            the data dictionary.
        """
        super().__init__(
            model=model,
            filename=filename,
            region_component=region_component,
        )

    ### I/O methods ###
    @hydromt_step
    def read(self) -> None:
        """Read maps from initialfield and parse to dict of xr.DataArray."""
        self.root.is_reading_mode()

        # Read initial fields
        inifield_model = self.model.dfmmodel.geometry.inifieldfile
        if inifield_model:
            # separate 1d and 2d
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
                        _fn = join(self.root.path, "maps", _fn.name)
                    inimap = open_raster(_fn)
                    name = inidict.quantity
                    # Need to get branchid from mdu
                    if name == "frictioncoefficient":
                        frictype = self.model.mdu.get_value(
                            "physics.uniffricttype", fallback=1
                        )
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
                    self.set(inimap, name)

    @hydromt_step
    def write(self) -> None:
        """Write maps as tif files in maps folder and update initial fields."""
        self.root._assert_write_mode()
        if self._data is None:
            logger.debug("No maps data found, skip writing.")
            return

        # Global parameters
        mapsroot = join(self.root.path, "maps")
        Path(mapsroot).mkdir(parents=True, exist_ok=True)
        inilist = []
        paramlist = []
        logger.info(f"Writing maps files to {mapsroot}")

        def _prepare_inifields(da_dict, da):
            # Write tif files
            name = da_dict["name"]
            initype = da_dict["initype"]
            interp_method = da_dict["interpolation"]
            locationtype = da_dict["locationtype"]
            _fn = join(mapsroot, f"{name}.tif")
            if da.raster.nodata is None or np.isnan(da.raster.nodata):
                da.raster.set_nodata(-999)
            da.raster.to_raster(_fn)
            logger.info(f"Writing file {mapsroot}/{name}.tif")
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
            if initype == "initial":
                inilist.append(inidict)
            elif initype == "parameter":
                paramlist.append(inidict)

        # Only write maps that are listed in self._MAPS, rename tif on the fly
        # TODO raise value error if both waterdepth and waterlevel are given in maps
        for name, ds in self._data.items():
            if isinstance(ds, xr.DataArray):
                if name in self._MAPS:
                    _prepare_inifields(self._MAPS[name], ds)
                    # update mdu if friction
                    if "frictype" in self._MAPS[name]:
                        self.model.mdu.set(
                            "physics.uniffricttype", self._MAPS[name]["frictype"]
                        )
                    # update mdu if infiltration
                    if name == "infiltcap":
                        self.model.mdu.set("grw.infiltrationmodel", 2)
                else:
                    logger.error(f"Could not write map to model: {name} not recognized")
            elif isinstance(ds, xr.Dataset):
                for v in ds.data_vars:
                    if v in self._MAPS:
                        _prepare_inifields(self._MAPS[v], ds[v])
                        # update mdu if friction
                        if self._MAPS[v] == "frictype":
                            self.model.mdu.set(
                                "physics.uniffricttype", self._MAPS[v]["frictype"]
                            )
                        # update mdu if infiltration
                        if v == "infiltcap":
                            self.model.mdu.set("grw.infiltrationmodel", 2)
                    else:
                        logger.error(
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
        fm_dir = dirname(join(self.root.path, self.model.mdu._filename))
        inifield_model.save(
            join(fm_dir, inifield_model_filename),
            recurse=False,
        )
        # save filepath in the config
        self.model.mdu.set("geometry.inifieldfile", inifield_model_filename)

    ### Add data methods ###
    @hydromt_step
    def add_raster_data_from_rasterdataset(
        self,
        raster_filename: str | Path | xr.Dataset,
        variables: list = None,
        fill_method: str | None = None,
        reproject_method: str | None = "nearest",
        interpolation_method: str | None = "triangulation",
        locationtype: str = "2d",
        name: str | None = None,
        split_dataset: bool = True,
    ) -> None:
        """
        Add data variable(s) from ``raster_fn`` to maps object.

        If raster is a dataset, all variables will be added unless ``variables`` list
        is specified.

        Adds model layers:

        * **raster.name** maps: data from raster_fn

        Parameters
        ----------
        raster_filename: str
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
            raise ValueError("name must be specified when split_dataset = False")

        # Call super method from HydroMT Core
        variables = super().add_raster_data_from_rasterdataset(
            raster_filename=raster_filename,
            variables=variables,
            fill_method=fill_method,
            reproject_method=reproject_method,
            name=name,
            split_dataset=split_dataset,
        )

        for var in variables:
            da = self.data[var]
            da = da.where(da != da.raster.nodata, -999.0)
            da.raster.set_nodata(-999.0)
            self.set(da, var)

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

    @hydromt_step
    def add_raster_data_from_raster_reclass(
        self,
        raster_filename: str | Path | xr.DataArray,
        reclass_table_filename: str | Path | pd.DataFrame,
        reclass_variables: list,
        fill_method: str | None = None,
        reproject_method: str | None = "nearest",
        interpolation_method: str | None = "triangulation",
        locationtype: str = "2d",
        name: str | None = None,
        split_dataset: bool = True,
    ) -> None:
        """
        Add data variable(s) to inifield by reclassifying values from raster.

        Reclassification is done bycombining values in ``raster_reclass_table_filename``
        to spatial layer ``raster_filename``.

        The ``mapping_variables`` rasters are first created by mapping variables values
        from ``raster_reclass_table_filename`` to value in the ``raster_filename`` grid.
        Adds model layers:

        * **mapping_variables** maps: data from raster_mapping_fn spatially
            distributed with raster_fn

        Parameters
        ----------
        raster_filename: str
            Source name of raster data in data_catalog. Should be a DataArray. Else use
            **kwargs to select variables/time_range in
            hydromt.data_catalog.get_rasterdataset method
        reclass_table_filename: str
            Source name of mapping table of raster_filename in data_catalog. Make sure
            the data type is consistant for a ``reclass_variables`` including nodata.
            For example, for roughness, it is common that the data type is float,
            then use no data value as -999.0.
        reclass_variables: list
            List of mapping_variables from raster_reclass_table_filename table to add to
            mesh. Index column should match values in raster_filename.
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
            raise ValueError("name must be specified when split_dataset = False")

        # Call super method
        reclass_variables = super().add_raster_data_from_raster_reclass(
            raster_filename=raster_filename,
            reclass_table_filename=reclass_table_filename,
            reclass_variables=reclass_variables,
            fill_method=fill_method,
            reproject_method=reproject_method,
            name=name,
            split_dataset=split_dataset,
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
                if self.data[var].raster.res[0] > self.model.mesh.res:
                    relsize = np.round(
                        np.abs(self.data[var].raster.res[0])
                        / self.model.mesh.res
                        * np.sqrt(2)
                        + 0.05,
                        2,
                    )
                else:
                    relsize = 1.01
                self._MAPS[var]["averagingrelsize"] = relsize
