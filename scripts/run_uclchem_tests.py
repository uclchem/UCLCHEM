"""Test some basic UCLCHEM models.

This should be run from the UCLCHEM root directory.

"""

from pathlib import Path
from time import perf_counter

import uclchem

if __name__ == "__main__":
    uclchem.utils.configure_logging(level="INFO")

    out_dir = Path("examples/test-output")
    out_dir.mkdir(parents=True, exist_ok=True)
    save_file = out_dir / "models.h5"
    if save_file.exists():
        save_file.unlink()

    print("Running test models...")

    # Static cloud model (no collapse)
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
        "outputFile": str(out_dir / "static-full.dat"),
    }

    start = perf_counter()
    static_model = uclchem.model.Cloud(param_dict=params)
    stop = perf_counter()
    print(f"Static model in {stop - start:.1f} seconds")
    static_model.save_model(file=str(save_file), name="static", overwrite=True)

    # Phase 1: collapsing cloud to build up initial abundances
    params = {
        "endAtFinalDensity": True,
        "freefall": True,
        "initialDens": 1e2,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5e6,
        "outputFile": str(out_dir / "phase1-full.dat"),
    }

    start = perf_counter()
    phase1_model = uclchem.model.Cloud(param_dict=params)
    stop = perf_counter()
    print(f"Phase 1 in {stop - start:.1f} seconds")
    phase1_model.save_model(file=str(save_file), name="phase1", overwrite=True)

    # Phase 2: prestellar core warm-up using phase 1 abundances in memory
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
        "outputFile": str(out_dir / "phase2-full.dat"),
        "reltol": 1e-6,
        "abstol_factor": 1e-8,
        # MUST match gas — freeze-out moves species between phases
        "abstol_ice_factor": 1e-8,
        "abstol_min": 1e-20,
        "abstol_ice_min": 1e-25,
        "mxstep": 100_000,
        "writeTimeStepInfo": True,
    }

    start = perf_counter()
    phase2_model = uclchem.model.PrestellarCore(
        temp_indx=3,
        max_temperature=300.0,
        param_dict=params,
        previous_model=phase1_model,
    )
    stop = perf_counter()
    print(f"Phase 2 in {stop - start:.1f} seconds")
    phase2_model.save_model(file=str(save_file), name="phase2", overwrite=True)
