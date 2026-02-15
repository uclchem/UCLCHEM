"""
Test multi-stage model runs with DISK-BASED I/O only.

This test uses outputFile, abundSaveFile, abundLoadFile, and columnFile
to ensure all model stages work with Fortran disk I/O.
"""

import shutil
import tempfile
from pathlib import Path

import pandas as pd
import pytest

try:
    import uclchem

    uclchem_imported = True
except ImportError:
    uclchem_imported = False


def test_import_uclchem():
    if not uclchem_imported:
        pytest.fail(
            "uclchem module could not be imported, "
            "make sure your environment is loaded and UCLCHEM is installed."
        )


@pytest.fixture(scope="function")
def test_output_directory(request):
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


def test_static_model_disk(test_output_directory):
    """Test static cloud model with disk output"""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
        "outputFile": str(test_output_directory / "static-full.dat"),
        "abundSaveFile": str(test_output_directory / "startstatic.dat"),
    }
    return_code = uclchem.functional.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Static model returned with nonzero exit code {return_code[0]}"

    # Verify output files exist
    assert Path(
        params["outputFile"]
    ).exists(), f"Output file not created: {params['outputFile']}"
    assert Path(
        params["abundSaveFile"]
    ).exists(), f"Abundance save file not created: {params['abundSaveFile']}"

    # Verify finalTime is respected (within 10% tolerance)
    output_df = pd.read_csv(params["outputFile"], sep=",", skipinitialspace=True)
    max_time = output_df["Time"].max()
    assert (
        max_time <= 1.1 * params["finalTime"]
    ), f"Model exceeded finalTime tolerance: {max_time:.2e} > {1.1 * params['finalTime']:.2e}"
    assert (
        max_time >= 0.9 * params["finalTime"]
    ), f"Model stopped too early: {max_time:.2e} < {0.9 * params['finalTime']:.2e}"


def test_collapse_hotcore_disk(test_output_directory):
    """Test collapse -> hot core chained models with disk I/O"""
    # Stage 1: Collapse with endAtFinalDensity=True
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "finalDens": 1e6,
        "finalTime": 1e5,
        "abundSaveFile": str(test_output_directory / "startstage1.dat"),
        "outputFile": str(test_output_directory / "stage1-full.dat"),
        "columnFile": str(test_output_directory / "stage1-column.dat"),
    }
    return_code = uclchem.functional.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Stage 1 returned with nonzero exit code {return_code[0]}"

    # Verify output files exist
    assert Path(
        params["outputFile"]
    ).exists(), f"Output file not created: {params['outputFile']}"
    assert Path(
        params["abundSaveFile"]
    ).exists(), f"Abundance save file not created: {params['abundSaveFile']}"
    assert Path(
        params["columnFile"]
    ).exists(), f"Column file not created: {params['columnFile']}"

    # Verify endAtFinalDensity=True behavior: stops at finalTime OR finalDens
    output_df = pd.read_csv(params["outputFile"], skipinitialspace=True)
    max_time = output_df["Time"].max()
    max_density = output_df["Density"].max()
    assert (
        max_time <= 1.1 * params["finalTime"] or max_density >= 0.9 * params["finalDens"]
    ), f"Collapse should stop at finalTime OR finalDens: time={max_time:.2e}, density={max_density:.2e}"

    # Stage 2: Hot core
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
        "outputFile": str(test_output_directory / "stage2-full.dat"),
        "abundLoadFile": str(test_output_directory / "startstage1.dat"),
    }
    return_code = uclchem.functional.prestellar_core(
        3, 300.0, param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Stage 2 returned with nonzero exit code {return_code[0]}"

    # Verify output files exist
    assert Path(
        params["outputFile"]
    ).exists(), f"Output file not created: {params['outputFile']}"
    assert Path(
        params["abundLoadFile"]
    ).exists(), f"Abundance load file doesn't exist: {params['abundLoadFile']}"

    # Verify finalTime is respected (within 10% tolerance)
    output_df = pd.read_csv(params["outputFile"], skipinitialspace=True)
    max_time = output_df["Time"].max()
    assert (
        max_time <= 1.1 * params["finalTime"]
    ), f"Hot core exceeded finalTime tolerance: {max_time:.2e} > {1.1 * params['finalTime']:.2e}"
    assert (
        max_time >= 0.9 * params["finalTime"]
    ), f"Hot core stopped too early: {max_time:.2e} < {0.9 * params['finalTime']:.2e}"


def test_cshock_disk(test_output_directory):
    """Test C-shock model with disk output"""
    # Pre-shock cloud
    param_dict = {
        "endAtFinalDensity": False,
        "freefall": True,
        "initialDens": 1e2,
        "finalDens": 1e4,
        "initialTemp": 10.0,
        "finalTime": 6.0e6,
        "rout": 0.1,
        "baseAv": 1.0,
        "abundSaveFile": str(test_output_directory / "start_cshock.dat"),
        "outputFile": str(test_output_directory / "pre_cshock.dat"),
    }
    return_code = uclchem.functional.cloud(param_dict=param_dict)
    assert (
        return_code[0] == 0
    ), f"Pre-cshock cloud returned with nonzero exit code {return_code[0]}"

    # Verify output files exist
    assert Path(
        param_dict["outputFile"]
    ).exists(), f"Output file not created: {param_dict['outputFile']}"
    assert Path(
        param_dict["abundSaveFile"]
    ).exists(), f"Abundance save file not created: {param_dict['abundSaveFile']}"

    # Verify finalTime is respected (within 10% tolerance)
    output_df = pd.read_csv(param_dict["outputFile"], skipinitialspace=True)
    max_time = output_df["Time"].max()
    assert (
        max_time <= 1.1 * param_dict["finalTime"]
    ), f"Pre-shock cloud exceeded finalTime: {max_time:.2e} > {1.1 * param_dict['finalTime']:.2e}"

    # C-shock
    param_dict["initialDens"] = 1e4
    param_dict["finalTime"] = 1e6
    param_dict["abundLoadFile"] = str(test_output_directory / "start_cshock.dat")
    param_dict["outputFile"] = str(test_output_directory / "cshock.dat")
    param_dict.pop("abundSaveFile", None)  # Remove abundSaveFile for second stage
    return_code = uclchem.functional.cshock(shock_vel=40, param_dict=param_dict)
    assert (
        return_code[0] == 0
    ), f"C-shock returned with nonzero exit code {return_code[0]}"

    # Verify output files exist
    assert Path(
        param_dict["outputFile"]
    ).exists(), f"Output file not created: {param_dict['outputFile']}"
    assert Path(
        param_dict["abundLoadFile"]
    ).exists(), f"Abundance load file doesn't exist: {param_dict['abundLoadFile']}"

    # Verify finalTime is respected (shock models may stop early)
    output_df = pd.read_csv(param_dict["outputFile"], skipinitialspace=True)
    max_time = output_df["Time"].max()
    assert (
        max_time <= 1.1 * param_dict["finalTime"]
    ), f"C-shock exceeded finalTime: {max_time:.2e} > {1.1 * param_dict['finalTime']:.2e}"


def test_endAtFinalDensity_validation_disk(test_output_directory):
    """Test that endAtFinalDensity=True raises error without freefall for Cloud with disk I/O"""
    params = {
        "endAtFinalDensity": True,  # Invalid without freefall
        "freefall": False,
        "initialDens": 1e4,
        "finalTime": 1e5,
        "outputFile": str(test_output_directory / "invalid-test.dat"),
    }

    # Should raise ValueError for Cloud without freefall
    with pytest.raises(ValueError, match="parcelStoppingMode != 0 can only be used"):
        uclchem.functional.cloud(
            param_dict=params,
            out_species=["OH", "CO"],
        )


def main():
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
