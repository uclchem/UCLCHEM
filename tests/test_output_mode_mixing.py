"""
Unit tests to verify that users cannot mix in-memory and write-to-disk modes.

These tests check the pre_flight_checklist function in model.py to ensure:
1. Users cannot specify output/input files when using return_array or return_dataframe
2. Users cannot run disk-based models after in-memory models in the same session
3. Users cannot run in-memory models after disk-based models in the same session
4. starting_chemistry can only be used with return_array or return_dataframe
5. return_rates requires return_array or return_dataframe
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


@pytest.fixture(scope="function")
def reset_output_mode():
    """Reset the OUTPUT_MODE global variable between tests"""
    import uclchem.model as model

    original_mode = model.OUTPUT_MODE
    model.OUTPUT_MODE = ""
    yield
    model.OUTPUT_MODE = original_mode


@pytest.fixture
def basic_params():
    """Basic parameter dictionary for testing"""
    return {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 1.0e3,  # Shorter time for faster tests
    }


class TestFileSpecificationWithInMemoryMode:
    """Test that output/input files cannot be specified with return_array/return_dataframe"""

    def test_return_array_with_outputFile_raises_error(
        self, basic_params, reset_output_mode
    ):
        """Test that specifying outputFile with return_array raises RuntimeError"""
        params = basic_params.copy()
        params["outputFile"] = "test_output.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cloud(
                param_dict=params,
                return_array=True,
            )

    def test_return_dataframe_with_outputFile_raises_error(
        self, basic_params, reset_output_mode
    ):
        """Test that specifying outputFile with return_dataframe raises RuntimeError"""
        params = basic_params.copy()
        params["outputFile"] = "test_output.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cloud(
                param_dict=params,
                return_dataframe=True,
            )

    def test_return_array_with_abundSaveFile_raises_error(
        self, basic_params, reset_output_mode
    ):
        """Test that specifying abundSaveFile with return_array raises RuntimeError"""
        params = basic_params.copy()
        params["abundSaveFile"] = "test_abund.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cloud(
                param_dict=params,
                return_array=True,
            )

    def test_return_dataframe_with_abundLoadFile_raises_error(
        self, basic_params, reset_output_mode
    ):
        """Test that specifying abundLoadFile with return_dataframe raises RuntimeError"""
        params = basic_params.copy()
        params["abundLoadFile"] = "test_abund.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cloud(
                param_dict=params,
                return_dataframe=True,
            )

    def test_return_array_with_columnFile_raises_error(
        self, basic_params, reset_output_mode
    ):
        """Test that specifying columnFile with return_array raises RuntimeError"""
        params = basic_params.copy()
        params["columnFile"] = "test_column.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cloud(
                param_dict=params,
                return_array=True,
            )

    def test_return_dataframe_with_multiple_files_raises_error(
        self, basic_params, reset_output_mode
    ):
        """Test that specifying multiple *File parameters with return_dataframe raises RuntimeError"""
        params = basic_params.copy()
        params["outputFile"] = "test_output.dat"
        params["abundSaveFile"] = "test_abund.dat"
        params["columnFile"] = "test_column.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cloud(
                param_dict=params,
                return_dataframe=True,
            )


class TestStartingChemistryRestrictions:
    """Test that starting_chemistry can only be used with return_array or return_dataframe"""

    def test_starting_chemistry_without_return_mode_raises_error(
        self, basic_params, reset_output_mode
    ):
        """Test that using starting_chemistry without return_array/return_dataframe raises AssertionError"""
        params = basic_params.copy()
        dummy_abundances = np.zeros(uclchem.constants.n_species)

        with pytest.raises(
            AssertionError,
            match="starting_chemistry can only be used with return_array or return_dataframe set to True",
        ):
            uclchem.model.cloud(
                param_dict=params,
                starting_chemistry=dummy_abundances,
            )

    def test_starting_chemistry_with_return_array_succeeds(
        self, basic_params, reset_output_mode
    ):
        """Test that using starting_chemistry with return_array works"""
        params = basic_params.copy()
        dummy_abundances = np.zeros(uclchem.constants.n_species)

        # This should not raise an error
        physics, chemistry, rates, abundances_start, return_code = uclchem.model.cloud(
            param_dict=params,
            return_array=True,
            starting_chemistry=dummy_abundances,
        )
        assert return_code == 0

    def test_starting_chemistry_with_return_dataframe_succeeds(
        self, basic_params, reset_output_mode
    ):
        """Test that using starting_chemistry with return_dataframe works"""
        params = basic_params.copy()
        dummy_abundances = np.zeros(uclchem.constants.n_species)

        # This should not raise an error
        physics, chemistry, rates, abundances_start, return_code = uclchem.model.cloud(
            param_dict=params,
            return_dataframe=True,
            starting_chemistry=dummy_abundances,
        )
        assert return_code == 0


class TestReturnRatesRestrictions:
    """Test that return_rates can only be used with return_array or return_dataframe"""

    def test_return_rates_without_return_mode_raises_error(
        self, basic_params, reset_output_mode, temp_output_directory
    ):
        """Test that return_rates with outputFile but no return_array/dataframe raises error"""
        params = basic_params.copy()
        params["outputFile"] = temp_output_directory / "test_output.dat"

        # return_rates=True triggers the file check even without return_array/dataframe
        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cloud(
                param_dict=params,
                return_rates=True,
            )

    def test_return_rates_with_return_array_succeeds(
        self, basic_params, reset_output_mode
    ):
        """Test that return_rates with return_array works and returns rates"""
        params = basic_params.copy()

        physics, chemistry, rates, abundances_start, return_code = uclchem.model.cloud(
            param_dict=params,
            return_array=True,
            return_rates=True,
        )
        assert return_code == 0
        assert rates is not None
        assert rates.shape[2] == uclchem.constants.n_reactions

    def test_return_rates_with_return_dataframe_succeeds(
        self, basic_params, reset_output_mode
    ):
        """Test that return_rates with return_dataframe works and returns rates DataFrame"""
        params = basic_params.copy()

        physics, chemistry, rates, abundances_start, return_code = uclchem.model.cloud(
            param_dict=params,
            return_dataframe=True,
            return_rates=True,
        )
        assert return_code == 0
        assert rates is not None
        assert len(rates.columns) == uclchem.constants.n_reactions


class TestSessionModeConsistency:
    """
    Test that users cannot switch between in-memory and disk modes within a session.

    NOTE: These tests are isolated with reset_output_mode fixture to avoid affecting other tests.
    In a real session without restarting the kernel, switching modes would fail.
    """

    def test_disk_then_memory_mode_raises_error(
        self, basic_params, temp_output_directory, reset_output_mode
    ):
        """Test that running disk-based model then in-memory model raises AssertionError"""
        # First run a disk-based model
        params_disk = basic_params.copy()
        params_disk["outputFile"] = temp_output_directory / "test1.dat"
        result = uclchem.model.cloud(param_dict=params_disk)
        assert result[0] == 0

        # Now try to run an in-memory model - should fail
        params_memory = basic_params.copy()
        with pytest.raises(
            AssertionError,
            match="Cannot run an in memory based model after running a disk based one",
        ):
            uclchem.model.cloud(
                param_dict=params_memory,
                return_array=True,
            )

    def test_memory_then_disk_mode_raises_error(
        self, basic_params, temp_output_directory, reset_output_mode
    ):
        """Test that running in-memory model then disk-based model raises AssertionError"""
        # First run an in-memory model
        params_memory = basic_params.copy()
        physics, chemistry, rates, abundances, return_code = uclchem.model.cloud(
            param_dict=params_memory,
            return_array=True,
        )
        assert return_code == 0

        # Now try to run a disk-based model - should fail
        params_disk = basic_params.copy()
        params_disk["outputFile"] = temp_output_directory / "test2.dat"
        with pytest.raises(
            AssertionError,
            match="Cannot run a disk based model after running an in memory one",
        ):
            uclchem.model.cloud(param_dict=params_disk)

    def test_multiple_memory_models_succeed(self, basic_params, reset_output_mode):
        """Test that running multiple in-memory models in sequence works"""
        params = basic_params.copy()

        # Run first in-memory model
        physics1, chemistry1, rates1, abundances1, return_code1 = uclchem.model.cloud(
            param_dict=params,
            return_array=True,
        )
        assert return_code1 == 0

        # Run second in-memory model - should succeed
        physics2, chemistry2, rates2, abundances2, return_code2 = uclchem.model.cloud(
            param_dict=params,
            return_dataframe=True,
        )
        assert return_code2 == 0

    def test_multiple_disk_models_succeed(
        self, basic_params, temp_output_directory, reset_output_mode
    ):
        """Test that running multiple disk-based models in sequence works"""
        # Run first disk-based model
        params1 = basic_params.copy()
        params1["outputFile"] = temp_output_directory / "test1.dat"
        result1 = uclchem.model.cloud(param_dict=params1)
        assert result1[0] == 0

        # Run second disk-based model - should succeed
        params2 = basic_params.copy()
        params2["outputFile"] = temp_output_directory / "test2.dat"
        result2 = uclchem.model.cloud(param_dict=params2)
        assert result2[0] == 0


class TestChainedModelsWorkflow:
    """Test the typical workflow from notebooks 2a and 2b for chained models"""

    def test_chained_models_on_disk(
        self, basic_params, temp_output_directory, reset_output_mode
    ):
        """Test Stage 1 -> Stage 2 workflow using disk files (like notebook 2a)"""
        # Stage 1: Cloud collapse
        params_stage1 = {
            "endAtFinalDensity": False,
            "freefall": True,
            "initialDens": 1e2,
            "finalDens": 1e6,
            "initialTemp": 10.0,
            "finalTime": 6.0e5,  # Shorter for testing
            "rout": 0.1,
            "baseAv": 1.0,
            "abundSaveFile": temp_output_directory / "startcollapse.dat",
            "outputFile": temp_output_directory / "stage1-full.dat",
        }
        result1 = uclchem.model.cloud(param_dict=params_stage1)
        assert result1[0] == 0

        # Stage 2: Hot core using saved abundances
        params_stage2 = {
            "initialDens": 1e6,
            "finalTime": 1e5,  # Shorter for testing
            "freefall": False,
            "freezeFactor": 0.0,
            "abundLoadFile": temp_output_directory / "startcollapse.dat",
            "outputFile": temp_output_directory / "stage2-full.dat",
        }
        result2 = uclchem.model.hot_core(
            temp_indx=3,
            max_temperature=300.0,
            param_dict=params_stage2,
        )
        assert result2[0] == 0

    def test_chained_models_in_memory(self, basic_params, reset_output_mode):
        """Test Stage 1 -> Stage 2 workflow using in-memory arrays (like notebook 2b)"""
        # Stage 1: Cloud collapse
        params_stage1 = {
            "endAtFinalDensity": False,
            "freefall": True,
            "initialDens": 1e2,
            "finalDens": 1e6,
            "initialTemp": 10.0,
            "finalTime": 6.0e5,  # Shorter for testing
            "rout": 0.1,
            "baseAv": 1.0,
        }
        (
            df_stage1_physics,
            df_stage1_chemistry,
            df_stage1_rates,
            final_abundances,
            result1,
        ) = uclchem.model.cloud(
            param_dict=params_stage1,
            return_dataframe=True,
        )
        assert result1 == 0

        # Stage 2: Hot core using in-memory abundances
        params_stage2 = {
            "initialDens": 1e6,
            "finalTime": 1e5,  # Shorter for testing
            "freefall": False,
            "freezeFactor": 0.0,
        }
        (
            df_stage2_physics,
            df_stage2_chemistry,
            df_stage2_rates,
            final_abundances2,
            result2,
        ) = uclchem.model.hot_core(
            temp_indx=3,
            max_temperature=300.0,
            param_dict=params_stage2,
            return_dataframe=True,
            starting_chemistry=final_abundances,
        )
        assert result2 == 0

    def test_cannot_mix_disk_and_memory_in_chain(
        self, basic_params, temp_output_directory, reset_output_mode
    ):
        """Test that you cannot start with disk mode then switch to memory mode in a chain"""
        # Stage 1: Cloud collapse (disk mode)
        params_stage1 = {
            "endAtFinalDensity": False,
            "freefall": True,
            "initialDens": 1e2,
            "finalDens": 1e6,
            "initialTemp": 10.0,
            "finalTime": 6.0e5,
            "rout": 0.1,
            "baseAv": 1.0,
            "outputFile": temp_output_directory / "stage1-full.dat",
        }
        result1 = uclchem.model.cloud(param_dict=params_stage1)
        assert result1[0] == 0

        # Stage 2: Try to use memory mode - should fail
        params_stage2 = {
            "initialDens": 1e6,
            "finalTime": 1e5,
            "freefall": False,
        }
        with pytest.raises(
            AssertionError,
            match="Cannot run an in memory based model after running a disk based one",
        ):
            uclchem.model.hot_core(
                temp_indx=3,
                max_temperature=300.0,
                param_dict=params_stage2,
                return_dataframe=True,
            )


class TestDifferentModelTypes:
    """Test mode mixing restrictions apply to all model types (cloud, hot_core, cshock, jshock, collapse)"""

    def test_hot_core_respects_mode_restrictions(self, basic_params, reset_output_mode):
        """Test that hot_core also enforces mode restrictions"""
        params = basic_params.copy()
        params["outputFile"] = "test_output.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.hot_core(
                temp_indx=3,
                max_temperature=300.0,
                param_dict=params,
                return_array=True,
            )

    def test_cshock_respects_mode_restrictions(self, basic_params, reset_output_mode):
        """Test that cshock also enforces mode restrictions"""
        params = basic_params.copy()
        params["outputFile"] = "test_output.dat"
        params["initialDens"] = 1e4

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.cshock(
                shock_vel=40,
                param_dict=params,
                return_dataframe=True,
            )

    def test_collapse_respects_mode_restrictions(self, basic_params, reset_output_mode):
        """Test that collapse also enforces mode restrictions"""
        params = basic_params.copy()
        params["outputFile"] = "test_output.dat"

        with pytest.raises(
            RuntimeError,
            match="return_array or return_dataframe cannot be used if any output of input file is specified",
        ):
            uclchem.model.collapse(
                collapse="BE1.1",
                physics_output=None,
                param_dict=params,
                return_array=True,
            )


def main():
    """Run the tests using pytest"""
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
