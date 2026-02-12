"""
Unit tests to verify that users cannot mix in-memory and write-to-disk modes.

Core principle: Never cross-pollinate between Fortran-based disk I/O
and Python in-memory arrays/dataframes.

These tests check the pre_flight_checklist function in model.py to ensure:
1. Users cannot specify file parameters with return_array/return_dataframe
2. Users cannot switch between disk and memory modes in the same session
3. starting_chemistry requires return_array or return_dataframe
"""

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

try:
    import uclchem

    uclchem_imported = True
except ImportError:
    uclchem_imported = False


def test_import_uclchem():
    """Verify uclchem was imported successfully"""
    if not uclchem_imported:
        pytest.fail(
            "uclchem module could not be imported, make sure your environment is loaded and UCLCHEM is installed correctly."
        )


@pytest.fixture(scope="function")
def temp_output_directory():
    """Create a temporary directory for each test"""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def basic_params():
    """Basic parameter dictionary for testing"""
    return {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 1.0e3,
    }


# Test 1: Files cannot be specified with return_array
def test_return_array_with_file_raises_error(basic_params):
    """Test that specifying any *File parameter with return_array raises RuntimeError"""
    params = basic_params.copy()
    params["outputFile"] = "test_output.dat"

    with pytest.raises(
        RuntimeError,
        match="return_array or return_dataframe cannot be used if any output or input file is specified",
    ):
        uclchem.functional.cloud(param_dict=params, return_array=True)


# Test 2: Files cannot be specified with return_dataframe
def test_return_dataframe_with_file_raises_error(basic_params):
    """Test that specifying any *File parameter with return_dataframe raises RuntimeError"""
    params = basic_params.copy()
    params["abundSaveFile"] = "test_abund.dat"

    with pytest.raises(
        RuntimeError,
        match="return_array or return_dataframe cannot be used if any output or input file is specified",
    ):
        uclchem.functional.cloud(param_dict=params, return_dataframe=True)


# Test 3: starting_chemistry can now be used with any mode (no longer restricted)
def test_starting_chemistry_with_memory_mode(basic_params):
    """Test that starting_chemistry works with return_array/return_dataframe"""
    params = basic_params.copy()
    dummy_abundances = np.ones(uclchem.constants.n_species) * 1e-20
    # Ensure that common species have sensible defaults (use model lookup for indices):
    species = uclchem.model.get_species_names()
    species_idx = {name: i for i, name in enumerate(species)}
    # Set common species where available
    dummy_abundances[species_idx["H2"]] = 0.45
    dummy_abundances[species_idx["H"]] = 0.1
    dummy_abundances[species_idx["HE"]] = 0.1
    dummy_abundances[species_idx["C"]] = 1e-7
    dummy_abundances[species_idx["O"]] = 3e-7

    # This should work fine now
    result = uclchem.functional.cloud(
        param_dict=params, starting_chemistry=dummy_abundances, return_array=True
    )
    assert result[-1] == 0  # success_flag should be 0


# Test 4: return_rates requires memory mode
def test_return_rates_with_file_raises_error(basic_params, temp_output_directory):
    """Test that return_rates with file output raises error"""
    params = basic_params.copy()
    params["outputFile"] = temp_output_directory / "test_output.dat"

    with pytest.raises(
        RuntimeError,
        match="return_array or return_dataframe cannot be used if any output or input file is specified",
    ):
        uclchem.functional.cloud(param_dict=params, return_rates=True)


# Test 5: Cannot run memory mode after disk mode
# TODO: DEPRECATED - Mode mixing checks no longer needed since IO is handled by model objects directly, not in Fortran
# def test_disk_then_memory_mode_raises_error(
#     basic_params, temp_output_directory, reset_output_mode
# ):
#     """Test that running disk-based model then in-memory model raises AssertionError"""
#     # First run a disk-based model
#     params_disk = basic_params.copy()
#     params_disk["outputFile"] = temp_output_directory / "test1.dat"
#     result = uclchem.functional.cloud(param_dict=params_disk)
#     assert result[0] == 0
#
#     # Now try to run an in-memory model - should fail
#     params_memory = basic_params.copy()
#     with pytest.raises(
#         AssertionError,
#         match="Cannot run an in memory based model after running a disk based one",
#     ):
#         uclchem.functional.cloud(param_dict=params_memory, return_array=True)


# Test 6: Cannot run disk mode after memory mode
# TODO: DEPRECATED - Mode mixing checks no longer needed since IO is handled by model objects directly, not in Fortran
# def test_memory_then_disk_mode_raises_error(
#     basic_params, temp_output_directory, reset_output_mode
# ):
#     """Test that running in-memory model then disk-based model raises AssertionError"""
#     # First run an in-memory model
#     params_memory = basic_params.copy()
#     physics, chemistry, rates, heating, abundances, return_code = (
#         uclchem.functional.cloud(
#             param_dict=params_memory, return_array=True, return_rates=True
#         )
#     )
#     assert return_code == 0
#
#     # Now try to run a disk-based model - should fail
#     params_disk = basic_params.copy()
#     params_disk["outputFile"] = temp_output_directory / "test2.dat"
#     with pytest.raises(
#         AssertionError,
#         match="Cannot run a disk based model after running an in memory one",
#     ):
#         uclchem.functional.cloud(param_dict=params_disk)


# Test 7: Multiple memory models succeed
def test_multiple_memory_models_succeed(basic_params):
    """Test that running multiple in-memory models in sequence works"""
    params = basic_params.copy()

    # Run first in-memory model
    physics1, chemistry1, rates1, heating1, abundances1, return_code1 = (
        uclchem.functional.cloud(param_dict=params, return_array=True, return_rates=True)
    )
    assert return_code1 == 0

    # Run second in-memory model - should succeed
    physics2, chemistry2, rates2, heating2, abundances2, return_code2 = (
        uclchem.functional.cloud(param_dict=params, return_dataframe=True)
    )
    assert return_code2 == 0


# Test 8: Multiple disk models succeed
def test_multiple_disk_models_succeed(basic_params, temp_output_directory):
    """Test that running multiple disk-based models in sequence works"""
    # Run first disk-based model
    params1 = basic_params.copy()
    params1["outputFile"] = temp_output_directory / "test1.dat"
    result1 = uclchem.functional.cloud(param_dict=params1)
    assert result1[0] == 0

    # Run second disk-based model - should succeed
    params2 = basic_params.copy()
    params2["outputFile"] = temp_output_directory / "test2.dat"
    result2 = uclchem.functional.cloud(param_dict=params2)
    assert result2[0] == 0


# Test 9: Chained models work with starting_chemistry in memory
def test_chained_models_in_memory(basic_params):
    """Test Stage 1 -> Stage 2 workflow using in-memory arrays with starting_chemistry"""
    # Stage 1: Cloud collapse
    params_stage1 = {
        "endAtFinalDensity": False,
        "freefall": True,
        "initialDens": 1e2,
        "finalDens": 1e6,
        "initialTemp": 10.0,
        "finalTime": 6.0e5,
        "rout": 0.1,
        "baseAv": 1.0,
    }
    _, _, _, _, final_abundances, result1 = uclchem.functional.cloud(
        param_dict=params_stage1, return_dataframe=True
    )
    assert result1 == 0

    # Stage 2: Hot core using starting_chemistry
    params_stage2 = {
        "initialDens": 1e6,
        "finalTime": 1e5,
        "freefall": False,
        "freezeFactor": 0.0,
    }
    _, _, _, _, final_abundances2, result2 = uclchem.functional.prestellar_core(
        temp_indx=3,
        max_temperature=300.0,
        param_dict=params_stage2,
        return_dataframe=True,
        starting_chemistry=final_abundances,
    )
    assert result2 == 0


# Test 10: Cannot mix disk and memory in chained models
# TODO: DEPRECATED - Mode mixing checks no longer needed since IO is handled by model objects directly, not in Fortran
# def test_cannot_mix_disk_and_memory_in_chain(
#     basic_params, temp_output_directory, reset_output_mode
# ):
#     """Test that you cannot start with disk mode then switch to memory mode in a chain"""
#     # Stage 1: Cloud collapse (disk mode)
#     params_stage1 = {
#         "endAtFinalDensity": False,
#         "freefall": True,
#         "initialDens": 1e2,
#         "finalDens": 1e6,
#         "initialTemp": 10.0,
#         "finalTime": 6.0e5,
#         "rout": 0.1,
#         "baseAv": 1.0,
#         "outputFile": temp_output_directory / "stage1-full.dat",
#     }
#     result1 = uclchem.functional.cloud(param_dict=params_stage1)
#     assert result1[0] == 0
#
#     # Stage 2: Try to use memory mode - should fail
#     params_stage2 = {
#         "initialDens": 1e6,
#         "finalTime": 1e5,
#         "freefall": False,
#     }
#     with pytest.raises(
#         AssertionError,
#         match="Cannot run an in memory based model after running a disk based one",
#     ):
#         uclchem.functional.prestellar_core(
#             temp_indx=3,
#             max_temperature=300.0,
#             param_dict=params_stage2,
#             return_dataframe=True,
#         )


def main():
    """Run the tests using pytest"""
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
