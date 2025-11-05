"""
Test multi-stage model runs with IN-MEMORY return modes only.

This test uses return_array and return_dataframe with starting_chemistry
to ensure all model stages work with Python in-memory arrays.
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
    """Reset OUTPUT_MODE at start of module to allow memory tests"""
    import uclchem.model as model

    original_mode = model.OUTPUT_MODE
    model.OUTPUT_MODE = ""
    yield
    model.OUTPUT_MODE = original_mode


def test_static_model_return_array(common_output_directory):
    """Test static cloud model with return_array"""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
    }
    physics, chemistry, rates, abundances_start, return_code = (
        uclchem.model.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_array=True,
        )
    )
    assert return_code == 0, (
        f"Static model returned with nonzero exit code {return_code}"
    )


def test_static_model_return_dataframe(common_output_directory):
    """Test static cloud model with return_dataframe"""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
    }
    physics, chemistry, rates, abundances_start, return_code = (
        uclchem.model.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_dataframe=True,
        )
    )
    assert return_code == 0, (
        f"Static model returned with nonzero exit code {return_code}"
    )


def test_collapse_hotcore_return_array(common_output_directory):
    """Test collapse -> hot core chained models with return_array"""
    # Stage 1: Collapse
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }
    physics, chemistry, rates, abundances_start, return_code = (
        uclchem.model.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_array=True,
        )
    )
    assert return_code == 0, (
        f"Stage 1 returned with nonzero exit code {return_code}"
    )

    # Stage 2: Hot core using starting_chemistry
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    physics, chemistry, rates, abundances_start, return_code = (
        uclchem.model.hot_core(
            3,
            300.0,
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_array=True,
            starting_chemistry=abundances_start,
        )
    )
    assert return_code == 0, (
        f"Stage 2 returned with nonzero exit code {return_code}"
    )


def test_collapse_hotcore_return_dataframe(common_output_directory):
    """Test collapse -> hot core chained models with return_dataframe"""
    # Stage 1: Collapse
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }
    physics, chemistry, rates, abundances_start, return_code = (
        uclchem.model.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_dataframe=True,
        )
    )
    assert return_code == 0, (
        f"Stage 1 returned with nonzero exit code {return_code}"
    )

    # Stage 2: Hot core using starting_chemistry
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    physics, chemistry, rates, abundances_start, return_code = (
        uclchem.model.hot_core(
            3,
            300.0,
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_dataframe=True,
            starting_chemistry=abundances_start,
        )
    )
    assert return_code == 0, (
        f"Stage 2 returned with nonzero exit code {return_code}"
    )


def test_cshock_return_dataframe(common_output_directory):
    """Test C-shock model with return_dataframe and starting_chemistry"""
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
    }
    (
        df_stage1_physics,
        df_stage1_chemistry,
        rates,
        final_abundances,
        return_code,
    ) = uclchem.model.cloud(param_dict=param_dict, return_dataframe=True)
    assert return_code == 0, (
        f"Pre-cshock cloud returned with nonzero exit code {return_code}"
    )

    # C-shock with starting_chemistry
    param_dict["initialDens"] = 1e4
    param_dict["finalTime"] = 1e6
    (
        df_stage2_physics,
        df_stage2_chemistry,
        rates,
        dissipation_time,
        final_abundances,
        return_code,
    ) = uclchem.model.cshock(
        shock_vel=40,
        param_dict=param_dict,
        return_dataframe=True,
        starting_chemistry=final_abundances,
    )
    assert return_code == 0, (
        f"C-shock returned with nonzero exit code {return_code}"
    )


def main():
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
