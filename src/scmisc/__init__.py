"""

Miscellaneous tools for the analysis of single cell genomics data in python. A not comprehensive port of the R package `scmisc <https://github.com/ddiez/scmisc>`_


"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("scmisc")
except PackageNotFoundError:
  # package is not installed
  pass

from . import atac
from . import tools
from . import plotting as pl
from ._data import *
from ._doublet import *
from ._tools import *
from ._utils import *
