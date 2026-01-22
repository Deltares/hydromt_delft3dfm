"""Implement Delft3D FM 2D3D HydroMT plugin model class."""

import logging
from pathlib import Path

from hydromt.model import Model
from pyproj import CRS

from hydromt_delft3dfm.components import (
    MDUComponent,
)

__all__ = ["DFlowFMM2D3DModel"]
__hydromt_eps__ = ["DFlowFMM2D3DModel"]  # core entrypoints
logger = logging.getLogger(f"hydromt.{__name__}")


class DFlowFMM2D3DModel(Model):
    """API for Delft3D-FM 2D3D models in HydroMT."""

    name: str = "dflowfm_2d3d"

    def __init__(
        self,
        root: str | Path,
        mode: str = "w",
        mdu_filename: str = None,
        data_libs: list[str] = [],  # yml
        crs: int | str = None,
    ):
        """Initialize the DFlowFM1D2DModel.

        Parameters
        ----------
        root : str or Path
            The model root location.
        mode : {'w','r','r+'}
            Write/read/append mode.
            Default is "w".
        mdu_filename : str, optional
            The D-Flow FM model configuration file (.mdu).
            If None, default mdu file is used.
            Default is None.
        data_libs : list of str, optional
            List of data catalog yaml files.
            Default is None.
        crs : EPSG code, int, optional
            EPSG code of the model.
        """
        if not isinstance(root, (str, Path)):
            raise ValueError("The 'root' parameter should be a of str or Path.")

        # FIXME Xiaohan mdu needs to be derived from dimr_fn if dimr_fn exists
        components = {
            "mdu": MDUComponent(self, filename=str(mdu_filename)),
        }
        super().__init__(
            root=root,
            components=components,
            mode=mode,
            data_libs=data_libs,
        )

        # crs
        self._crs = CRS.from_user_input(crs) if crs else None
