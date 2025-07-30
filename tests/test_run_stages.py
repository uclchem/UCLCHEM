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
    cloud = uclchem.model.Cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        cloud.success_flag == 0
    ), f"Static model returned with nonzero exit code {cloud.success_flag}"


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
    #(
    #    physics,
    #    chemistry,
    #    rates,
    #    abundances_start,
    #    return_code,
    #) =
    cloud = uclchem.model.Cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        cloud.success_flag == 0
    ), f"Static model returned with nonzero exit code {cloud.success_flag}"
    np.save(common_output_directory / "physics_static_array.npy", cloud.physics_array)
    np.save(common_output_directory / "chemistry_static_array.npy", cloud.chemical_abun_array)


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
    cloud = uclchem.model.Cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_dataframe=True,
    )
    physics, chemistry = cloud.get_dataframes(joined=False)
    assert (
        cloud.success_flag == 0
    ), f"Static model returned with nonzero exit code {cloud.success_flag}"
    physics.to_csv(common_output_directory / "physics_static_array.csv")
    chemistry.to_csv(common_output_directory / "chemistry_static_array.csv")


# Test for running on disk
def test_collapse_prestellarcore_on_disk(common_output_directory):
    # Stage 1
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "abundSaveFile": common_output_directory / "startstage1.dat",
        "outputFile": common_output_directory / "stage1-full.dat",
        "columnFile": common_output_directory / "stage1-column.dat",
    }
    cloud = uclchem.model.Cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        cloud.success_flag == 0
    ), f"Stage 1 returned with nonzero exit code {cloud.success_flag}"

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
    p_core = uclchem.model.PrestellarCore(
        3, 300.0, param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    assert (
        p_core.success_flag == 0
    ), f"Stage 2 returned with nonzero exit code {p_core.success_flag}"


def test_collapse_prestellarcore_return_array(common_output_directory):
    # STAGE 1
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }
    cloud = uclchem.model.Cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
    )
    assert cloud.success_flag == 0, f"Stage 1 returned with nonzero exit code {cloud.success_flag}"
    # STAGE 2
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    p_core = uclchem.model.PrestellarCore(
        3,
        300.0,
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        previous_model=cloud,
    )
    assert p_core.success_flag == 0, f"Stage 2 returned with nonzero exit code {p_core.success_flag}"
    np.save(common_output_directory / "physics_prestellarcore_array.npy", p_core.physics_array)
    np.save(common_output_directory / "chemistry_prestellarcore_array.npy", p_core.chemical_abun_array)


def test_collapse_prestellarcore_return_dataframe(common_output_directory):
    # STAGE 1
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }
    cloud = uclchem.model.Cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
    )
    assert cloud.success_flag == 0, f"Stage 1 returned with nonzero exit code {cloud.success_flag}"
    # STAGE 2
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    p_core = uclchem.model.PrestellarCore(
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
        previous_model=cloud,
    )
    physics, chemistry = p_core.get_dataframes(joined=False)
    assert p_core.success_flag == 0, f"Stage 2 returned with nonzero exit code {p_core.success_flag}"
    physics.to_csv(common_output_directory / "physics_static_array.csv")
    chemistry.to_csv(common_output_directory / "chemistry_static_array.csv")


def test_cshock_return_array(common_output_directory):
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
    shock_start = uclchem.model.Cloud(
        param_dict=param_dict,
    )
    assert shock_start.success_flag==0, f"Stage 1 pre-cshock returned with nonzero exit code {shock_start.success_flag}"
    # STAGE 2 - cshock
    param_dict["initialDens"] = 1e4
    param_dict["finalTime"] = 1e6
    cshock = uclchem.model.CShock(
        shock_vel=40,
        param_dict=param_dict,
        previous_model=shock_start,
    )
    assert cshock.success_flag == 0, f"Stage 2 cshock returned with nonzero exit code {cshock.success_flag}"
    np.save(common_output_directory / "physics_prestellarcore_array.npy", cshock.physics_array)
    np.save(common_output_directory / "chemistry_prestellarcore_array.npy", cshock.chemical_abun_array)


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
    shock_start = uclchem.model.Cloud(
            param_dict=param_dict,
        )
    assert shock_start.success_flag==0, f"Stage 1 pre-cshock returned with nonzero exit code {shock_start.success_flag}"
    # STAGE 2 - cshock
    param_dict["initialDens"] = 1e4
    param_dict["finalTime"] = 1e6
    cshock = uclchem.model.CShock(
        shock_vel=40,
        param_dict=param_dict,
        previous_model=shock_start,
    )
    physics, chemistry = cshock.get_dataframes(joined=False)
    assert (cshock.success_flag == 0), f"Stage 2 cshock returned with nonzero exit code {cshock.success_flag}"
    physics.to_csv(common_output_directory / "physics_static_array.csv")
    chemistry.to_csv(common_output_directory / "chemistry_static_array.csv")



''' jshock is super slow, so disable it for now:
def test_jshock_return_dataframe(common_output_directory):
    # STAGE 1 - jshock
    param_dict = {
    "endAtFinalDensity": False,
    "freefall": False,
    "initialDens": 1e2,
    "finalDens": 1e4,
    "initialTemp": 10.0,
    "finalTime": 6.0e6,
    "rout": 0.1,
    "baseAv": 1.0,
    "reltol": 1e-12
    }
    df_stage1_physics, df_stage1_chemistry, final_abundances, return_code = uclchem.model.cloud(param_dict=param_dict, return_dataframe=True, )
    assert return_code == 0, f"Stage 1 pre-jshock returned with nonzero exit code {return_code}"
    param_dict["initialDens"] = 1e4
    param_dict["freefall"] = False
    param_dict["reltol"] = 1e-12
    shock_vel = 10.0
    df_jshock_physics, df_jshock_chemistry, final_abundances, return_code = uclchem.model.jshock(shock_vel=shock_vel,param_dict=param_dict, return_dataframe=True, starting_chemistry=final_abundances, timepoints=2000)
    assert return_code == 0, f"Stage 2 jshock returned with nonzero exit code {return_code}"
'''

def main():
    # Run the tests using pytest
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
