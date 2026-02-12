"""
UCLCHEM Constants Module

This module provides access to canonical constants defined in the Fortran backend.
Values are read directly from the compiled f2py_constants module to ensure
Python and Fortran stay in sync without code generation.

The values in this module reflect the current compiled state of UCLCHEM.
If you run MakeRates to change the network size or coolants, you must
reinstall UCLCHEM for these constants to update:
    1. Run makerates: python -m uclchem.makerates.makerates user_settings.yaml
    2. Reinstall: pip install . --force-reinstall --no-deps
    3. Verify: python -c "from uclchem import constants; print(constants.n_species)"
"""

# Import canonical values from compiled Fortran module
from uclchemwrap import f2py_constants

# Import PHYSICAL_PARAMETERS from its canonical source
# This is defined in makerates to avoid circular dependency
from uclchem.makerates.io_functions import PHYSICAL_PARAMETERS

# Read canonical values from Fortran
n_species = int(f2py_constants.nspec)
n_reactions = int(f2py_constants.nreac)
N_PHYSICAL_PARAMETERS = int(f2py_constants.n_physics_params)
NCOOLANTS = int(f2py_constants.ncoolants)
N_DVODE_STATS = int(f2py_constants.n_dvode_stats)
N_TOTAL_LEVELS = int(f2py_constants.n_total_levels)
N_SE_STATS_PER_COOLANT = int(f2py_constants.n_se_stats_per_coolant)

# DVODE solver statistics names
DVODE_STAT_NAMES = [
    "ISTATE",
    "HU",
    "HCUR",
    "TCUR",
    "TOLSF",
    "NST",
    "NFE",
    "NJE",
    "NQU",
    "NQCUR",
    "IMXER",
    "LENRW",
    "LENIW",
    "NLU",
    "NNI",
    "NCFN",
    "NETF",
    "CPU_TIME",
]

# SE solver statistics names (per coolant: convergence flag, iterations, max relative change)
SE_STAT_NAMES = []
for i in range(NCOOLANTS):
    SE_STAT_NAMES.extend(
        [
            f"COOLANT_{i:02d}_CONVERGED",
            f"COOLANT_{i:02d}_ITERATIONS",
            f"COOLANT_{i:02d}_MAX_REL_CHANGE",
        ]
    )

# Validate consistency
if len(PHYSICAL_PARAMETERS) != N_PHYSICAL_PARAMETERS:
    raise RuntimeError(
        f"PHYSICAL_PARAMETERS length ({len(PHYSICAL_PARAMETERS)}) does not match "
        f"N_PHYSICAL_PARAMETERS from Fortran ({N_PHYSICAL_PARAMETERS}). "
        "This indicates a build inconsistency. Please run MakeRates and reinstall."
    )

# User-configurable constants (not from Fortran)
TIMEPOINTS = 2000  # Number of timepoints for Fortran interface

# Default parameter dictionary
# These are default values for model parameters, not network structure constants
default_param_dictionary = {
    "initialtemp": 10.0,
    "initialdens": 100.0,
    "finaldens": 100000.0,
    "currenttime": 0.0,
    "finaltime": 5000000.0,
    "radfield": 1.0,
    "zeta": 1.0,
    "rout": 0.05,
    "rin": 0.0,
    "baseav": 2.0,
    "points": 1,
    "bm0": 1.0,
    "freezefactor": 1.0,
    "parcelstoppingmode": 0,  # Default: never stop (0=never, 1=stop all when outermost reaches max, 2=stop each individually)
    "freefall": True,
    "freefallfactor": 1.0,
    "desorb": True,
    "h2desorb": True,
    "crdesorb": True,
    "uvdesorb": True,
    "thermdesorb": True,
    "instantsublimation": True,
    "cosmicrayattenuation": True,
    "ionmodel": "L",
    "improvedh2crpdissociation": True,
    "outputfile": None,
    "columnfile": None,
    "ratefile": None,
    "fluxfile": None,
    "writestep": 1,
    "abundsavefile": None,
    "abundloadfile": None,
    "metallicity": 1.0,
    "ion": 2,
    "fh": 0.5,
    "fhe": 0.1,
    "fc": 0.000177,
    "fo": 0.000334,
    "fn": 6.18e-05,
    "fs": 3.51e-06,
    "fmg": 2.256e-06,
    "fsi": 1.78e-06,
    "fcl": 3.39e-08,
    "fp": 7.78e-08,
    "ffe": 2.01e-07,
    "ff": 3.6e-08,
    "fd": 0.0,
    "fli": 0.0,
    "fna": 0.0,
    "fpah": 0.0,
    "f15n": 0.0,
    "f13c": 0.0,
    "f18o": 0.0,
    "reltol": 1e-08,
    "abstol_factor": 1e-14,
    "abstol_min": 1e-25,
    "mxstep": 10000,
    "ebmaxh2": 1210.0,
    "ebmaxcr": 1210.0,
    "ebmaxuvcr": 10000.0,
    "epsilon": 0.01,
    "uv_yield": 0.03,
    "phi": 100000.0,
    "uvcreff": 0.001,
    "omega": 0.5,
}
