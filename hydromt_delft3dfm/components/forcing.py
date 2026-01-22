"""DFlowFM Forcing component."""

import logging
from os.path import dirname, join
from pathlib import Path

import numpy as np
import pandas as pd
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import SpatialDatasetsComponent

from hydromt_delft3dfm import mesh_utils, utils

__all__ = ["DFlowFMForcingComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DFlowFMForcingComponent(SpatialDatasetsComponent):
    """
    Manage the Delft3D-FM forcing files for model forcings.

    This class is used to manage spatial datasets related to model forcings.
    The forcing component data stored in the ``data`` property is a dictionary of
    xarray DataArrays or Datasets.
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "bnd.ext",
        region_component: str | None = None,
    ):
        """Initialize DFlowFMForcingComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            by default "bnd.ext", i.e. one file per dataset in the data dictionary.
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
        """Read forcing at <root/?/> and parse to dict of xr.DataArray."""
        self.root.is_reading_mode()
        # Read external forcing
        ext_model = self.model.dfmmodel.external_forcing.extforcefilenew

        if ext_model is None:
            logger.warning("No external forcing file found.")
            return

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
                        self.model.mesh.mesh_datasets["network1d"]
                    )
                    node_geoms = network1d_nodes[
                        np.isin(network1d_nodes["nodeid"], df.nodeid.values)
                    ]
                    da_out = utils.read_1dboundary(df, quantity=name, nodes=node_geoms)
                    # Add to forcing
                    self.set(da_out)
            # 2d boundary
            df_ext_2d = df_ext.loc[df_ext.nodeid.isna(), :]
            if len(df_ext_2d) > 0:
                for _, df in df_ext_2d.iterrows():
                    da_out = utils.read_2dboundary(
                        df, workdir=self.model.dfmmodel.filepath.parent
                    )
                    # Add to forcing
                    self.set(da_out)
        # lateral
        if len(ext_model.lateral) > 0:
            df_ext = pd.DataFrame([f.__dict__ for f in ext_model.lateral])
            da_out = utils.read_1dlateral(
                df_ext, branches=self.model.branches
            )  # TODO extend support to get laterals on nodes #78
            # Add to forcing
            self.set(da_out)
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
                self.set(da_out)
        # TODO lateral

    @hydromt_step
    def write(self) -> None:
        """Write forcing into hydrolib-core ext and forcing models."""
        if len(self._data) == 0:
            logger.debug("No forcing data found, skip writing.")
        else:
            self.root.is_writing_mode()
            logger.info("Writing forcing files.")
            savedir = dirname(join(self.root.path, self.model.mdu._filename))
            # create new external forcing file
            ext_fn = "bnd.ext"
            Path(join(savedir, ext_fn)).unlink(missing_ok=True)
            # populate external forcing file
            utils.write_1dboundary(self.data, savedir, ext_fn=ext_fn)
            utils.write_2dboundary(self.data, savedir, ext_fn=ext_fn)
            utils.write_1dlateral(self.data, savedir, ext_fn=ext_fn)
            utils.write_meteo(self.data, savedir, ext_fn=ext_fn)
            self.model.mdu.set("external_forcing.extforcefilenew", ext_fn)
