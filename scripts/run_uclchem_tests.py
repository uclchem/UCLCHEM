"""Test some basic UCLCHEM models.

This should be run from the UCLCHEM root directory.

"""

import logging
from pathlib import Path
from time import perf_counter

import uclchem

logging.basicConfig(level=logging.DEBUG)


if __name__ == "__main__":
    settings = uclchem.advanced.GeneralSettings()

    output_dir = Path("examples") / "test-output"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    print("Running test models...")
    # set a parameter dictionary for static model
    out_species = ["OH", "OCS", "CO", "CS", "CH3OH"]
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
        "outputFile": output_dir / "static-full.dat",
        "abundSaveFile": output_dir / "startstatic.dat",
        "reltol": 1e-6,
        "abstol_factor": 1e-12,
        "abstol_min": 1e-20,
        "writeTimestepInfo": True,
    }

    start = perf_counter()
    uclchem.functional.cloud(param_dict=params, out_species=out_species)
    stop = perf_counter()
    print(f"Static model in {stop - start:.1f} seconds")

    # change to collapsing phase1 params
    params["finalTime"] = 5e6
    params["freefall"] = True
    params["endAtFinalDensity"] = True
    params["initialDens"] = 1e2
    params["abundSaveFile"] = output_dir / "startcollapse.dat"
    params["outputFile"] = output_dir / "phase1-full.dat"
    params["columnFile"] = output_dir / "phase1-column.dat"
    start = perf_counter()
    uclchem.functional.cloud(param_dict=params, out_species=out_species)
    stop = perf_counter()
    print(f"Phase 1 in {stop - start:.1f} seconds")

    # finally, run phase 2 from the phase 1 model.
    params["initialDens"] = 1e5
    params["freezeFactor"] = 0.0
    params["thermdesorb"] = True
    params["endAtFinalDensity"] = False
    params["freefall"] = False
    params["finalTime"] = 1e6
    # Use default parameters to avoid convergence issues
    params["abundLoadFile"] = output_dir / "startcollapse.dat"
    params["outputFile"] = output_dir / "phase2-full.dat"
    params.pop("columnFile")
    start = perf_counter()
    uclchem.functional.prestellar_core(
        3, 300.0, param_dict=params, out_species=out_species
    )
    stop = perf_counter()
    print(f"Phase 2 in {stop - start:.1f} seconds")
