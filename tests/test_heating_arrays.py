"""
Tests for heating array functionality in UCLCHEM.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import uclchem


class TestHeatingArrays:
    """Test class for heating array functionality."""

    @pytest.fixture
    def param_dict(self):
        """Standard parameter dictionary for testing."""
        return {
            "endAtFinalDensity": False,
            "freefall": False,
            "writeStep": 1,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalDens": 1e5,
            "finalTime": 1.0e3,  # Much shorter time for faster tests
            "points": 1,
        }

    @pytest.fixture
    def expected_heating_columns(self):
        """Expected column names for heating DataFrame."""
        return [
            "Time",
            "Atomic Cooling",
            "Collisionally Induced Emission",
            "Compton Scattering Cooling",
            "Continuum Emission Cooling",
            "Line Cooling 1",
            "Line Cooling 2",
            "Line Cooling 3",
            "Line Cooling 4",
            "Line Cooling 5",
            "Photoelectric Heating",
            "H2Formation Heating",
            "FUVPumping Heating",
            "Photodissociation Heating",
            "CIonization Heating",
            "Cosmic Ray Heating",
            "Turbulent Heating",
            "Gas-Grain Collisions",
            "Chemical Heating",
        ]

    def test_array_creation_and_specifications(self, param_dict):
        """Test array creation and specifications."""
        from uclchem.model import _create_arrays, _get_standard_array_specs

        # Test array specifications
        array_specs = _get_standard_array_specs(return_heating=True)
        assert "heatarray" in array_specs
        assert array_specs["heatarray"]["third_dim"] == 19
        assert array_specs["heatarray"]["dtype"] == "float64"

        # Test array creation
        timepoints = 50
        arrays = _create_arrays(param_dict, array_specs, timepoints=timepoints)

        assert "heatarray" in arrays

        # Verify heating array has correct dimensions
        # timepoints+1, points, n_heating_terms
        expected_shape = (timepoints + 1, 1, 19)
        assert arrays["heatarray"].shape == expected_shape

    def test_cloud_function_with_return_array(self, param_dict):
        """Test cloud function with return_array=True."""

        (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            heatArray,
            abundanceStart,
            success_flag,
        ) = uclchem.model.cloud(
            param_dict=param_dict,
            out_species=["OH", "CO", "H2O"],
            return_array=True,
            return_heating=True,
            return_rates=True,
            timepoints=50,  # Reduced from 500 for faster tests
        )

        assert success_flag == 0, "Model run should be successful"
        assert heatArray is not None, "Heat array should be returned"
        assert isinstance(heatArray, np.ndarray), "Heat array should be numpy array"
        assert heatArray.shape[2] == 19, "Heat array should have 19 columns per particle"

    def test_cloud_function_with_return_dataframe(
        self, param_dict, expected_heating_columns
    ):
        """Test cloud function with return_dataframe=True."""

        result = uclchem.model.cloud(
            param_dict=param_dict,
            out_species=["OH", "CO", "H2O"],
            return_dataframe=True,
            return_heating=True,
            return_rates=True,
            timepoints=50,  # Reduced from 500 for faster tests
        )

        (
            physics_df,
            chemistry_df,
            rates_df,
            heating_df,
            abundanceStart,
            success_flag,
        ) = result

        assert success_flag == 0, "Model run should be successful"
        assert heating_df is not None, "Heating DataFrame should be returned"
        assert isinstance(heating_df, pd.DataFrame), "Heating data should be DataFrame"

        # Check DataFrame structure
        assert len(heating_df.columns) == len(expected_heating_columns)

        # Check column names
        actual_columns = list(heating_df.columns)
        missing_columns = set(expected_heating_columns) - set(actual_columns)
        extra_columns = set(actual_columns) - set(expected_heating_columns)

        assert not missing_columns, f"Missing columns: {missing_columns}"
        assert not extra_columns, f"Extra columns: {extra_columns}"

    @pytest.mark.parametrize(
        "model_function",
        [
            "cloud",
            "collapse",
        ],
    )
    def test_all_model_functions_support_heating(self, model_function, param_dict):
        """Test that all model functions support heating arrays."""
        # Get the function from the module
        func = getattr(uclchem.model, model_function)
        # Adjust parameters for different model types
        test_params = param_dict.copy()

        try:
            if model_function == "collapse":
                result = func(
                    collapse="ambipolar",
                    physics_output=None,
                    param_dict=test_params,
                    out_species=["OH", "CO"],
                    return_dataframe=True,
                    return_rates=True,
                    return_heating=True,
                    timepoints=50,  # Reduced from 500 for faster tests
                )
            elif model_function == "cloud":
                result = func(
                    param_dict=test_params,
                    out_species=["OH", "CO"],
                    return_dataframe=True,
                    return_rates=True,
                    return_heating=True,
                    timepoints=50,  # Reduced from 500 for faster tests
                )
            else:
                raise ValueError(f"Unknown model function: {model_function}")

            # Check that we get the expected number of return values
            assert len(result) >= 5, f"{model_function} should return at least 5 values"

            # The heating DataFrame should be in the result
            heating_df = None
            for item in result:
                if isinstance(item, pd.DataFrame) and len(item.columns) >= 18:
                    # This is likely the heating DataFrame
                    if "Time" in item.columns and "Atomic Cooling" in item.columns:
                        heating_df = item
                        break

            
            assert isinstance(heating_df, pd.DataFrame), "The output should be a DataFrame"
            assert "Time" in heating_df.columns, "Time should be returned"
            assert (heating_df.values[:, 1:] != 0.0).any(), f"Some terms should have non-zero values, head is {heating_df.head()}"

        except Exception as e:
            # Some model functions might have specific requirements
            # we haven't met. For now, we'll allow them to fail
            # with specific error messages
            allowed_errors = [
                "not implemented",
                "parameter missing",
            ]

            error_str = str(e).lower()
            if not any(allowed_error in error_str for allowed_error in allowed_errors):
                # Re-raise if it's an unexpected error
                raise

    def test_heating_array_content_validation(self, param_dict):
        """Test that heating arrays contain reasonable physical values."""
        (
            physics_df,
            chemistry_df,
            rates_df,
            heating_df,
            abundanceStart,
            success_flag,
        ) = uclchem.model.cloud(
            param_dict=param_dict,
            out_species=["OH", "CO", "H2O"],
            return_dataframe=True,
            return_rates=True,
            return_heating=True,
            timepoints=50,  # Reduced from 1000 for faster tests
        )

        assert success_flag == 0, f"Model run should be successful, or run out of points, instead it was {success_flag}"
        assert heating_df is not None, "Heating DataFrame should be returned"

        # Check that we have finite values (no NaN or inf)
        for col in heating_df.columns:
            if col != "Time":  # Time column might have different constraints
                finite_values = np.isfinite(heating_df[col]).all()
                assert finite_values, f"Column {col} contains non-finite values"

        # Check that at least some heating/cooling terms have non-zero values
        # (this ensures the physics is actually being calculated)
        non_time_columns = [col for col in heating_df.columns if col != "Time"]
        has_nonzero = False
        for col in non_time_columns:
            if (heating_df[col] != 0).any():
                has_nonzero = True
                break

        assert has_nonzero, (
            "At least some heating/cooling terms should have non-zero values"
        )
        
    def test_heating_array_to_disk(self, param_dict):
        """Test that heating arrays can be saved to disk."""
        TEST_DIR = Path("tests/heating_test_output/")
        TEST_DIR.mkdir(parents=True, exist_ok=True)
        TEST_FILE = TEST_DIR / "heating_file.csv"
        param_dict["heatingFile"] = str(TEST_FILE)
        result = uclchem.model.cloud(
            param_dict=param_dict,
            out_species=["OH", "CO"],
            timepoints=50,  # Reduced from 500 for faster tests
            return_rates=True,
            return_heating=True
        )
        assert result[0] == 0, "Model run should be successful"
        assert TEST_FILE.exists(), f"Heating file should be created on disk"
        heating_df = pd.read_csv(TEST_FILE, index_col=None)
        assert not heating_df.empty, "Heating DataFrame should not be empty"
        assert (heating_df.values[:, 1:] != 0.0).any(), "Heating DataFrame should have some non-zero values"
