"""hydroMT plugin for Delft3D-FM models."""

from os.path import abspath, dirname, join

__version__ = "0.1.3.dev"

DATADIR = join(dirname(abspath(__file__)), "data")

from .dflowfm import *
