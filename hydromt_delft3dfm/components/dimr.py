"""Custom Delft3D-FM config component."""

import logging
from os.path import basename, isfile, join
from pathlib import Path

from hydrolib.core.dimr import DIMR, FMComponent, Start
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import ModelComponent

__all__ = ["DIMRComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DIMRComponent(ModelComponent):
    """Manage the dflowfm DIMR file for model simulations/settings.

    ``DIMRComponent`` data is stored in a hydrolib core DIMR object. The component
    is used to prepare and update model simulations/settings of the delft3d model.
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "dimr_config.xml",
    ):
        """Manage DIMR configuration file for dflowfm simulations/settings.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            A path relative to the root where the configuration files will
            be read and written if user does not provide a path themselves.
            By default 'dimr_config.xml'.
        """
        self._data = None
        self._filename: str = filename

        super().__init__(model=model)

    @property
    def data(self):
        """Model dimr values."""
        if self._data is None:
            self._initialize()
        if self._data is None:
            raise RuntimeError("Could not load data for model dimr component")
        else:
            return self._data

    def _initialize(self, skip_read=False) -> None:
        """Initialize the model dimr."""
        if self._data is None:
            if not skip_read:
                # no check for read mode here
                # model dimr is read if in read-mode and it exists
                # default dimr if in write-mode
                self.read()
            else:
                self._data = DIMR()

    @hydromt_step
    def read(self) -> None:
        """
        Read the dflowfm DIMR file from <root/filename>.

        If the file does not exist, the hydrolib-core template will be used.
        """
        self._initialize(skip_read=True)

        dimr_fn = join(self.root.path, self._filename)
        # if file exist, read
        if isfile(dimr_fn) and self.root.is_reading_mode():
            logger.info(f"Reading dimr file at {dimr_fn}")
            dimr = DIMR(filepath=Path(dimr_fn))
        # else initialise
        else:
            self.root.is_writing_mode()
            logger.info("Initialising empty dimr file")
            dimr = DIMR()
        self._data = dimr

    @hydromt_step
    def write(self):
        """Write the dflowfm DIMR to a file."""
        # Check if dict is empty
        if not self.data:
            logger.info("No dflowfm dimr data found, skip writing.")
            return

        self._data.filepath = join(self.root.path, self._filename)

        if not self.root.is_reading_mode():
            # Updates the dimr file first before writing
            logger.info("Adding dflowfm component to dimr config")

            # update component
            components = self._data.component
            if len(components) != 0:
                components = []
            fmcomponent = FMComponent(
                name="dflowfm",
                workingdir="dflowfm",
                inputfile=basename(self.model.mdu._filename),
                model=self.model.dfmmodel,
            )
            components.append(fmcomponent)
            self._data.component = components
            # update control
            controls = self._data.control
            if len(controls) != 0:
                controls = []
            control = Start(name="dflowfm")
            controls.append(control)
            self._data.control = control

        # write
        logger.info(f"Writing model dimr file to {self._data.filepath}")
        self.data.save(recurse=False)
