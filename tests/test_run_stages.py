import gc
import os
import shutil
import tempfile
from pathlib import Path
from time import perf_counter

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


# # Test for ODE element conservation
# def test_element_conservation():
#     result = uclchem.tests.test_ode_conservation()
#     for key, value in result.items():
#         assert abs(value) < 1e-12, f"{key} not conserved with total rate of change {value:.2e}"


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
    start = perf_counter()
    print("Running the model")
    return_code = uclchem.model.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    print("Finished running!")
    print(return_code)
    stop = perf_counter()
    elapsed_time = stop - start
    assert (
        return_code[0] == 0
    ), f"Static model returned with nonzero exit code {return_code[0]}"

    # Test for running the static model (Stage 1)


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
    start = perf_counter()
    print("Running the model")
    (
        physics,
        chemistry,
        physics_start,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_array=True,
    )
    print("Finished running!")
    stop = perf_counter()
    elapsed_time = stop - start
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
    start = perf_counter()
    print("Running the model")
    (
        physics,
        chemistry,
        physics_start,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_dataframe=True,
    )
    print("Finished running!")
    stop = perf_counter()
    elapsed_time = stop - start
    assert (
        return_code == 0
    ), f"Static model returned with nonzero exit code {return_code}"
    physics.to_csv(common_output_directory / "physics_static_array.csv")
    chemistry.to_csv(common_output_directory / "chemistry_static_array.csv")


# Test for running Stage 1


def test_stage1_on_disk(common_output_directory):
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "abundSaveFile": common_output_directory / "startstage1.dat",
        "outputFile": common_output_directory / "stage1-full.dat",
        "columnFile": common_output_directory / "stage1-column.dat",
    }

    start = perf_counter()
    return_code = uclchem.model.cloud(
        param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    stop = perf_counter()
    elapsed_time = stop - start
    print(return_code)
    assert (
        return_code[0] == 0
    ), f"Stage 1 returned with nonzero exit code {return_code[0]}"


def test_stage1_return_array(common_output_directory):
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }

    start = perf_counter()
    (
        physics,
        chemistry,
        physics_start,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_array=True,
    )
    stop = perf_counter()
    elapsed_time = stop - start
    print(return_code)
    assert return_code == 0, f"Stage 1 returned with nonzero exit code {return_code}"
    np.save("stage1_phys_array.npy", physics_start)
    np.save("stage1_abund_array.npy", abundances_start)
    np.save("physics_1_array.npy", physics)
    np.save("abund_1_array.npy", chemistry)


def test_stage1_return_dataframe(common_output_directory):
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
    }

    start = perf_counter()
    (
        physics,
        chemistry,
        physics_start,
        abundances_start,
        return_code,
    ) = uclchem.model.cloud(
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_array=True,
    )
    stop = perf_counter()
    elapsed_time = stop - start
    print(return_code)
    assert return_code == 0, f"Stage 1 returned with nonzero exit code {return_code}"
    np.save("stage1_phys_df.npy", physics_start)
    np.save("stage1_abund_df.npy", abundances_start)


# Test for running Stage 2
def test_stage2_on_disk(common_output_directory):
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
    start = perf_counter()
    return_code = uclchem.model.hot_core(
        3, 300.0, param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"]
    )
    stop = perf_counter()
    elapsed_time = stop - start
    assert (
        return_code[0] == 0
    ), f"Stage 2 returned with nonzero exit code {return_code[0]}"


# Test for running Stage 2
def test_stage2_return_array(common_output_directory):
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    start_physics = np.asfortranarray(np.load("stage1_phys_array.npy"))
    start_abundances = np.asfortranarray(np.load("stage1_abund_array.npy"))
    start = perf_counter()
    (
        physics,
        chemistry,
        physics_start,
        abundances_start,
        return_code,
    ) = uclchem.model.hot_core(
        3,
        300.0,
        param_dict=params,
        out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
        return_array=True,
        starting_physics=start_physics,
        starting_chemistry=start_abundances,
    )
    stop = perf_counter()
    elapsed_time = stop - start
    assert return_code == 0, f"Stage 2 returned with nonzero exit code {return_code}"


# Test for running Stage 2
def test_stage2_return_dataframe(common_output_directory):
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
    start_physics = np.asfortranarray(np.load("stage1_phys_df.npy"))
    start_abundances = np.asfortranarray(np.load("stage1_abund_df.npy"))
    start = perf_counter()
    (
        physics,
        chemistry,
        physics_start,
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
        starting_physics=start_physics,
        starting_chemistry=start_abundances,
    )
    stop = perf_counter()
    elapsed_time = stop - start
    assert return_code == 0, f"Stage 2 returned with nonzero exit code {return_code}"


def main():
    import uclchem

    # Run the tests using pytest
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()


# if __name__ == "__main__":
#     datapath = Path("data")
#     datapath.mkdir(exist_ok=True, parents=True)
#     print("static 1")
#     test_static_model_on_disk(datapath)
#     print("static 2")
#     test_static_model_return_array(datapath)
#     print("static 3")
#     test_static_model_return_dataframe(datapath)
#     print("Stage 1 1")
#     test_stage1_on_disk(datapath)
#     print("Stage 1 2")
#     test_stage1_return_array(datapath)
#     print("Stage 1 3")
#     test_stage1_return_dataframe(datapath)
#     print("Stage 2")
#     test_stage2_on_disk(datapath)
#     test_stage2_return_array(datapath)
#     test_stage2_return_dataframe(datapath)
