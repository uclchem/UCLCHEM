# This file was machine generated with Makerates on 2024-07-02 10:35:59.575116
# This file contains the default magic numbers that ensure that fortran and
# python are in sync. If you adjust anything here, you must:
# 1. rerun makerates (this puts the magic numbers in fortran)
# 2. reinstall the python (this compiles fortran and the python wrapper)
# 3. rerun the tests (this ensures that everything is working)

PHYSICAL_PARAMETERS = [
    "age",
    "density",
    "gasTemp",
    "dustTemp",
    "Av",
    "radfield",
    "zeta",
    "dstep",
]
N_PHYSICAL_PARAMETERS = len(PHYSICAL_PARAMETERS)
TIMEPOINTS = 500
MAX_SPECIES = 335
