"""hydroMT plugin for Delft3D FM models."""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import pandas as pd
from packaging.version import Version

try:
    __version__ = version("hydromt_delft3dfm")
except PackageNotFoundError:
    # Package is not installed
    __version__ = "unknown"

DATADIR = Path(__file__).parent / "data"

from hydromt_delft3dfm.dflowfm import DFlowFMModel

__all__ = ["DFlowFMModel"]

# Avoid pandas FutureWarning from pandas 2.2 onwards: "Downcasting behavior
#  in `replace` is deprecated and will be removed in a future version"
#  by avoiding silent downcasting and converting dtypes manually.
# Setting this option generates a warning since pandas 3.0.0
#  "Pandas4Warning: future.no_silent_downcasting' is deprecated"
#  so only do this in pandas 3.
#  Also be sure to remove this when allowing for pandas 4:
#  https://github.com/Deltares/hydromt_delft3dfm/issues/354
pd_version = version("pandas")
if Version(pd_version) < Version("3.0.0"):
    pd.set_option("future.no_silent_downcasting", True)
