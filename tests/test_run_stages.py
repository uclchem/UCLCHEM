import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

# gc.set_debug(gc.DEBUG_LEAK)

try:
    import uclchem

    uclchem_imported = True
except ImportError:
    uclchem_imported = False


# Define a test to check if uclchem was imported successfully
def test_import_uclchem():
    if not uclchem_imported:
        pytest.fail(
            "uclchem module could not be imported, make sure your environment is loaded and UCLCHEM is installed correctly."
        )


@pytest.fixture(scope="module")
def common_output_directory(request):
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    # Clean up the temporary directory after all tests in the module
    shutil.rmtree(temp_dir, ignore_errors=True)


# Test for running the static model (Stage 1)
def test_static_model_on_disk(common_output_directory):
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
    return_code = uclchem.model.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Static model returned with nonzero exit code {return_code[0]}"


def test_static_model_return_array(common_output_directory):
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
    }
    (
        physics,
        chemistry,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_array=True,
    )
    assert (
        return_code == 0
    ), f"Static model returned with nonzero exit code {return_code}"
    np.save(common_output_directory / "physics_static_array.npy", physics)
    np.save(common_output_directory / "chemistry_static_array.npy", chemistry)


def test_static_model_return_dataframe(common_output_directory):
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
    }
    (
        physics,
        chemistry,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_dataframe=True,
    )
    assert (
        return_code == 0
    ), f"Static model returned with nonzero exit code {return_code}"
    physics.to_csv(common_output_directory / "physics_static_array.csv")
    chemistry.to_csv(common_output_directory / "chemistry_static_array.csv")


# Test for running on disk
def test_collapse_hotcore_on_disk(common_output_directory):
    # Stage 1
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "abundSaveFile": common_output_directory / "startstage1.dat",
        "outputFile": common_output_directory / "stage1-full.dat",
        "columnFile": common_output_directory / "stage1-column.dat",
    }
    return_code = uclchem.model.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Stage 1 returned with nonzero exit code {return_code[0]}"

    # Stage 2
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
    return_code = uclchem.model.hot_core(
        3, 300.0, param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        return_code[0] == 0
    ), f"Stage 2 returned with nonzero exit code {return_code[0]}"


def test_collapse_hotcore_return_array(common_output_directory):
    # STAGE 1
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }
    (
        physics,
        chemistry,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_array=True,
    )
    assert return_code == 0, f"Stage 1 returned with nonzero exit code {return_code}"
    # STAGE 2
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    (
        physics,
        chemistry,
        abundances_start,
        return_code,
    ) = uclchem.model.hot_core(
        3,
        300.0,
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_array=True,
        starting_chemistry=abundances_start,
    )
    assert return_code == 0, f"Stage 2 returned with nonzero exit code {return_code}"


def test_collapse_hotcore_return_dataframe(common_output_directory):
    # STAGE 1
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }
    (
        physics,
        chemistry,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_dataframe=True,
    )
    assert return_code == 0, f"Stage 1 returned with nonzero exit code {return_code}"
    # STAGE 2
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    (
        physics,
        chemistry,
        abundances_start,
        return_code,
    ) = uclchem.model.hot_core(
        3,
        300.0,
        param_dict=params,
        out_species=[
            "OH",
            "OCS",
            "CO",
            "CS",
            "CH3OH",
        ],
        return_dataframe=True,
        starting_chemistry=abundances_start,
    )
    assert return_code == 0, f"Stage 2 returned with nonzero exit code {return_code}"


def test_cshock_return_dataframe(common_output_directory):
    # STAGE 1 - cshock
    param_dict = {
        "endAtFinalDensity": False,  # stop at finalTime
        "freefall": True,  # increase density in freefall
        "initialDens": 1e2,  # starting density
        "finalDens": 1e4,  # final density
        "initialTemp": 10.0,  # temperature of gas
        "finalTime": 6.0e6,  # final time
        "rout": 0.1,  # radius of cloud in pc
        "baseAv": 1.0,  # visual extinction at cloud edge.
    }
    df_stage1_physics, df_stage1_chemistry, final_abundances, return_code = (
        uclchem.model.cloud(
            param_dict=param_dict,
            return_dataframe=True,
        )
    )
    assert (
        return_code == 0
    ), f"Stage 1 pre-cshock returned with nonzero exit code {return_code}"
    # STAGE 2 - cshock
    param_dict["initialDens"] = 1e4
    param_dict["finalTime"] = 1e6
    (
        df_stage2_physics,
        df_stage2_chemistry,
        dissipation_time,
        final_abundances,
        return_code,
    ) = uclchem.model.cshock(
        shock_vel=40,
        param_dict=param_dict,
        return_dataframe=True,
        starting_chemistry=final_abundances,
    )
    assert (
        return_code == 0
    ), f"Stage 2 cshock returned with nonzero exit code {return_code}"


# jshock is super slow, so disable it for now:
# def test_jshock_return_dataframe(common_output_directory):
#     # STAGE 1 - jshock
#     param_dict = {
#     "endAtFinalDensity": False,
#     "freefall": False,
#     "initialDens": 1e2,
#     "finalDens": 1e4,
#     "initialTemp": 10.0,
#     "finalTime": 6.0e6,
#     "rout": 0.1,
#     "baseAv": 1.0,
#     "reltol": 1e-12
#     }
#     df_stage1_physics, df_stage1_chemistry, final_abundances, return_code = uclchem.model.cloud(param_dict=param_dict, return_dataframe=True, )
#     assert return_code == 0, f"Stage 1 pre-jshock returned with nonzero exit code {return_code}"
#     param_dict["initialDens"] = 1e4
#     param_dict["freefall"] = False
#     param_dict["reltol"] = 1e-12
#     shock_vel = 10.0
#     df_jshock_physics, df_jshock_chemistry, final_abundances, return_code = uclchem.model.jshock(shock_vel=shock_vel,param_dict=param_dict, return_dataframe=True, starting_chemistry=final_abundances, timepoints=2000)
#     assert return_code == 0, f"Stage 2 jshock returned with nonzero exit code {return_code}"


def main():
    # Run the tests using pytest
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
