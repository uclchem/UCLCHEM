"""UCLCHEM Constants Module.

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

# Canonical definition of physical parameters
PHYSICAL_PARAMETERS = [
    "Time",
    "Density",
    "gasTemp",
    "dustTemp",
    "Av",
    "radfield",
    "zeta",
    "dstep",
    "parcel_radius",
]

# Read canonical values from Fortran
n_species = int(f2py_constants.nspec)
n_reactions = int(f2py_constants.nreac)
N_PHYSICAL_PARAMETERS = int(f2py_constants.n_physics_params)
NCOOLANTS = int(f2py_constants.ncoolants)
N_DVODE_STATS = int(f2py_constants.n_dvode_stats)
N_TOTAL_LEVELS = int(f2py_constants.n_total_levels)
N_SE_STATS_PER_COOLANT = int(f2py_constants.n_se_stats_per_coolant)

# DVODE solver statistics names
# Note: Stats are now written for EVERY solver attempt (including retries)
# TRAJECTORY_INDEX links solver stats to trajectory timesteps
DVODE_STAT_NAMES = [
    "TRAJECTORY_INDEX",  # Links to trajectory timestep (dtime)
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
    msg = (
        f"PHYSICAL_PARAMETERS length ({len(PHYSICAL_PARAMETERS)}) does not match "
        f"N_PHYSICAL_PARAMETERS from Fortran ({N_PHYSICAL_PARAMETERS}). "
        "This indicates a build inconsistency. Please run MakeRates and reinstall."
    )
    raise RuntimeError(msg)

# User-configurable constants (not from Fortran)
TIMEPOINTS = 2000  # Number of timepoints for Fortran interface

CENTIMETERS_PER_PARSEC = 3.086e18  # parsec in cgs
SECONDS_PER_YEAR = 3.15569e7

SPEED_OF_LIGHT_CGS = 2.99792458e10  # speed of light cm/s
PLANCK_CONSTANT_CGS = 6.62606896e-27  # Planck constant erg*s

# Physical constants matching collapse.f90
HYDROGEN_MASS_CGS = 1.6736e-24  # hydrogen mass in g
BOLTZMANN_CONSTANT_CGS = 1.38e-16  # Boltzmann constant in erg/K
GRAVITATIONAL_CONSTANT_CGS = 6.67e-8  # gravitational constant in cgs

# Default parameter dictionary
# These are default values for model parameters, not network structure constants
default_param_dictionary = {
    "initialtemp": 10.0,
    "initialdens": 100.0,
    "finaldens": 100000.0,
    "currenttime": 0.0,
    "finaltime": 5000000.0,
    "endatfinaldensity": False,
    "writetimestepinfo": False,
    "radfield": 1.0,
    "zeta": 1.0,
    "r_out": 0.05,
    "r_in": 0.0,
    "baseav": 2.0,
    "points": 1,
    "bm0": 1.0,
    "freezefactor": 1.0,
    "parcelstoppingmode": 0,  # Default: never stop (0=never, 1=stop all when outermost reaches max, 2=stop each individually)
    "freefall": False,
    "freefallfactor": 1.0,
    "desorb": True,
    "h2desorb": False,
    "crdesorb": True,
    "uvdesorb": True,
    "thermdesorb": True,
    "instantsublimation": False,
    "cosmicrayattenuation": False,
    "ionmodel": "L",
    "improvedh2crpdissociation": False,
    "outputfile": None,
    "columnfile": None,
    "rateConstantFile": None,
    "ratesFile": None,
    "heatingFile": None,
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
    "abstol_ice_factor": 1e-10,
    "abstol_ice_min": 1e-20,
    "reltol_phys": 1e-4,
    "abstol_phys_factor": 1e-4,
    "abstol_t_min": 0.01,
    "abstol_nh_min": 1.0,
    "mxstep": 10000,
    "ebmaxh2": 1210.0,
    "ebmaxcr": 1210.0,
    "ebmaxuvcr": 10000.0,
    "epsilon": 0.01,
    "uv_yield": 0.03,
    "phi": 100000.0,
    "uvcreff": 0.001,
    "omega": 0.5,
    # 1D radiative transfer defaults
    "enable_radiative_transfer": False,
    "density_scale_radius": 0.05,
    "density_power_index": 2.0,
    "lum_star": 1000000.0,
    "temp_star": 45000.0,
    # Advanced surface chemistry parameters
    "h2encounterdesorption": True,
    "hencounterdesorption": False,
    "edendothermicityfactor": 0.0,
    "h2stickingcoeffbyh2coverage": False,
    "hstickingcoeffbyh2coverage": False,
    "hdiffusionbarrier": -1.0,
    "usecustomdiffusionbarriers": True,
    "separatediffanddesorbprefactor": True,
    "usetstprefactors": False,  # Set this one to True and add the intertias from dijkhuis25.
    "usecustomprefactors": False,
    "useminissaleicechemdesefficiency": False,
    "maxgraintemp": 150.0,
    "parameterizeh2form": 2,
    # Heating/cooling ODE cache tolerances
    # tempDot is recomputed only when |ΔT| > min(abstol, reltol * T)
    "heating_temp_abstol": 1.0,
    "heating_temp_reltol": 1.0e-2,
    # Coolant/validation tolerances
    # freq_rel_tol default is auto-computed at makerates time from LAMDA file deviations
    "freq_rel_tol": float(getattr(f2py_constants, "suggested_freq_rel_tol", 0.1)),
    "pop_rel_tol": 0.1,
}

default_elements_to_check: list[str] = ["H", "N", "C", "O"]
