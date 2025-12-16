"""
Test multi-stage model runs with DISK-BASED I/O only.

This test uses outputFile, abundSaveFile, abundLoadFile, and columnFile
to ensure all model stages work with Fortran disk I/O.
"""

import shutil
import tempfile
from pathlib import Path

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


@pytest.fixture(scope="module")
def common_output_directory(request):
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture(scope="module", autouse=True)
def reset_output_mode():
    """Reset OUTPUT_MODE at start of module to allow disk tests"""
    import uclchem.model as model

    original_mode = model.OUTPUT_MODE
    model.OUTPUT_MODE = ""
    yield
    model.OUTPUT_MODE = original_mode


def test_static_model_disk(common_output_directory):
    """Test static cloud model with disk output"""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
        "outputFile": common_output_directory / "static-full.dat",
        "abundSaveFile": common_output_directory / "startstatic.dat",
    }
    return_code = uclchem.model.functional.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Static model returned with nonzero exit code {return_code[0]}"


def test_collapse_hotcore_disk(common_output_directory):
    """Test collapse -> hot core chained models with disk I/O"""
    # Stage 1: Collapse
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "abundSaveFile": common_output_directory / "startstage1.dat",
        "outputFile": common_output_directory / "stage1-full.dat",
        "columnFile": common_output_directory / "stage1-column.dat",
    }
    return_code = uclchem.model.functional.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Stage 1 returned with nonzero exit code {return_code[0]}"

    # Stage 2: Hot core
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
        "outputFile": common_output_directory / "stage2-full.dat",
        "abundLoadFile": common_output_directory / "startstage1.dat",
    }
    return_code = uclchem.model.functional.prestellar_core(
        3, 300.0, param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Stage 2 returned with nonzero exit code {return_code[0]}"


def test_cshock_disk(common_output_directory):
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
        "abundSaveFile": common_output_directory / "start_cshock.dat",
        "outputFile": common_output_directory / "pre_cshock.dat",
    }
    return_code = uclchem.model.functional.cloud(param_dict=param_dict)
    assert (
        return_code[0] == 0
    ), f"Pre-cshock cloud returned with nonzero exit code {return_code[0]}"

    # C-shock
    param_dict["initialDens"] = 1e4
    param_dict["finalTime"] = 1e6
    param_dict["abundLoadFile"] = common_output_directory / "start_cshock.dat"
    param_dict["outputFile"] = common_output_directory / "cshock.dat"
    param_dict.pop("abundSaveFile")
    return_code = uclchem.model.functional.cshock(shock_vel=40, param_dict=param_dict)
    assert (
        return_code[0] == 0
    ), f"C-shock returned with nonzero exit code {return_code[0]}"


def main():
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
