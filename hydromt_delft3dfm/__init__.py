"""hydroMT plugin for Delft3D FM models."""

from pathlib import Path

import pandas as pd

__version__ = "0.3.1.dev"
DATADIR = Path(__file__).parent / "data"

from hydromt_delft3dfm.dflowfm import DFlowFMModel

__all__ = ["DFlowFMModel"]

# Avoid pandas FutureWarning from pandas 2.2 onwards: "Downcasting behavior
#  in `replace` is deprecated and will be removed in a future version"
#  by avoiding silent downcasting and converting dtypes manually.
# TODO: setting this option generates a warning since pandas 3.0.0
#  "Pandas4Warning: future.no_silent_downcasting' is deprecated"
#  so be sure to remove it or silence it when supporting pandas 3 in
#  https://github.com/Deltares/hydromt_delft3dfm/issues/229
pd.set_option("future.no_silent_downcasting", True)
