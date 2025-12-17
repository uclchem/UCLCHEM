"""The UCLCHEM python module is divided into several parts.
`model` contains the functions for running chemical models under different physics.
`analysis` contains functions for reading and plotting output files as well as investigating the chemistry.
`advanced` provides access to heating and cooling mechanism controls and advanced solver parameters.
`tests` contains functions for testing the code.
"""

from . import advanced as advanced
from . import analysis as analysis

# This following contains the virtual submodule `functional`, which allows for calling the new API in the legacy format.
from . import model as model
from . import plot as plot
from . import tests as tests
from . import utils as utils
from . import version as version
