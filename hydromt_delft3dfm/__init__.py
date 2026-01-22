"""hydroMT plugin for Delft3D FM models."""

from pathlib import Path

__version__ = "0.3.1.dev"
DATADIR = Path(__file__).parent / "data"

from hydromt_delft3dfm.models.dflowfm_1d2d import *
from hydromt_delft3dfm.models.dflowfm_2d3d import *
