"""

Miscellaneous tools for the analysis of single cell genomics data in python. A not comprehensive port of the R package `scmisc <https://github.com/ddiez/scmisc>`_


"""

from ._version import __version__
from . import plotting as pl
from ._data import *
from ._doublet import *
from ._utils import *
