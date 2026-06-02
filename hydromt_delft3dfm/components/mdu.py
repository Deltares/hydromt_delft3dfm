"""Custom Delft3D FM config component."""

import datetime as dt
import logging
import os
from os.path import dirname, join
from pathlib import Path
from typing import cast

from hydrolib.core.dflowfm import FMModel
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import ConfigComponent

__all__ = ["MDUComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class MDUComponent(ConfigComponent):
    """Manage the dflowfm MDU file for model simulations/settings.

    ``MDUComponent`` data is stored in a dictionary. The component
    is used to prepare and update model simulations/settings of the dflowfm model.
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "dflowfm/DFlowFM.mdu",
    ):
        """Manage MDU configuration file for dflowfm simulations/settings.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            A path relative to the root where the configuration files will
            be read and written if user does not provide a path themselves.
            By default 'dflowfm/DFlowFM.mdu'.
        """
        super().__init__(
            model,
            filename=filename,
        )

    def _initialize(self, skip_read=False) -> None:
        """Initialize the model config."""
        if self._data is None:
            self._data = {}
            if not skip_read:
                # no check for read mode here
                # model config is read if in read-mode and it exists
                # default config if in write-mode
                self.read()

    def get_model_time(self):
        """
        Return (refdate, tstart, tstop) tuple.

        It is parsed from model startdatetime/stopdatetime, or from the refdate/tunit/
        tstart/tstop if not available.
        """

        def parse_time_to_datetime(key: str) -> dt.datetime:
            date_str = str(self.get_value(key))
            # expected format is yyyymmddhhmmss, but hhmmss may be omitted
            # (default:000000).
            if len(date_str) == 8:
                date_str += "000000"
            date_dt = dt.datetime.strptime(date_str, "%Y%m%d%H%M%S")
            return date_dt

        startdatetime_str = self.get_value("time.startdatetime", "")
        stopdatetime_str = self.get_value("time.stopdatetime", "")
        if startdatetime_str == "" or stopdatetime_str == "":
            logger.debug("get_model_time(): fallback to refdate/tstart/tstop")
            refdate = parse_time_to_datetime(key="time.refdate")
            tunit = self.get_value("time.tunit")
            tstart = float(self.get_value("time.tstart"))
            tstop = float(self.get_value("time.tstop"))
            if tunit.lower() == "s":
                tstart = refdate + dt.timedelta(seconds=tstart)
                tstop = refdate + dt.timedelta(seconds=tstop)
            elif tunit.lower() == "m":
                tstart = refdate + dt.timedelta(minutes=tstart)
                tstop = refdate + dt.timedelta(minutes=tstop)
            elif tunit.lower() == "h":
                tstart = refdate + dt.timedelta(hours=tstart)
                tstop = refdate + dt.timedelta(hours=tstop)
            elif tunit.lower() == "d":
                tstart = refdate + dt.timedelta(days=tstart)
                tstop = refdate + dt.timedelta(days=tstop)
            else:
                raise ValueError(f"tunit='{tunit}' not supported by get_model_time()")
        else:
            logger.debug("get_model_time(): from startdatetime/stopdatetime")
            tstart = parse_time_to_datetime(key="time.startdatetime")
            tstop = parse_time_to_datetime(key="time.stopdatetime")

        return tstart, tstop

    @hydromt_step
    def read(self) -> None:
        """
        Read the dflowfm MDU file from <root/filename>.

        If the file does not exist, the hydrolib template will be used.
        """
        self._initialize(skip_read=True)

        # Read via init_dfmmodel
        if self.model._dfmmodel is None:
            self.model.init_dfmmodel()
        # Convert to full dictionary without hydrolib-core objects
        cf_dict = dict()
        for k, v in self.model._dfmmodel.__dict__.items():
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
        self._data = cf_dict

    @hydromt_step
    def write(self):
        """Write the dflowfm MDU to a file."""
        # Check if dict is empty
        if not self.data:
            logger.info("No dflowfm config data found, skip writing.")
            return

        # Not sure if this is worth it compared to just calling write_config from super
        # advantage is the validator but the whole model is then read
        # when initialising FMModel
        self.root.is_writing_mode()

        cf_dict = self.data.copy()
        # Need to switch to dflowfm folder for files to be found and properly added
        _ = cf_dict.pop("filepath", None)
        mdu_fn = Path(join(self.root.path, self._filename))
        mdu_fn.parent.mkdir(parents=True, exist_ok=True)
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

    def test_equal(self, other: ConfigComponent) -> tuple[bool, dict[str, str]]:
        """Test if two components are equal.

        Parameters
        ----------
        other : ModelComponent
            The component to compare against.

        Returns
        -------
        tuple[bool, Dict[str, str]]
            True if the components are equal, and a dict with the associated errors per
            property checked.
        """
        errors = {}
        if not isinstance(other, self.__class__):
            errors["__class__"] = f"other does not inherit from {self.__class__}."
        eq = len(errors) == 0
        if not eq:
            return eq, errors
        other_config = cast(ConfigComponent, other)

        # not enough details in python recursion
        errors.update(**_check_equal(self.data, other_config.data, "mdu"))
        return len(errors) == 0, errors


def _check_equal(a, b, name="") -> dict[str, str]:
    """Recursive test of model components.

    Returns dict with component name and associated error message.
    """
    errors = {}
    try:
        assert isinstance(b, type(a)), "property types do not match"
        if isinstance(a, dict):
            for key in a:
                assert key in b, f"{key} missing"
                errors.update(**_check_equal(a[key], b[key], f"{name}.{key}"))
        else:
            assert a == b, "values not equal"
    except AssertionError as e:
        errors.update({name: e})
    return errors
