"""The UCLCHEM python module is divided into three parts.
`model` contains the functions for running chemical models under different physics.
`analysis` contains functions for reading and plotting output files as well as investigating the chemistry.
`tests` contains functions for testing the code.
"""

from . import analysis as analysis
from . import model as model
from . import plot as plot
from . import tests as tests
from . import utils as utils
from . import version as version
from . import trialmodel as trialmodel
