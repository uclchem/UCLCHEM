"""Integration test for new surface chemistry features.

Runs a cloud model and verifies the new physics produces different/better results.
"""

import shutil
import tempfile

import numpy as np
import pytest

import uclchem


def finite_difference(x, y):
    slope = np.diff(y) / np.diff(x)
    midpoints = np.diff(x) / 2 + x[:-1]
    return midpoints, slope


@pytest.fixture
def temp_output_dir():
    """Create temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp(prefix="uclchem_analysis_test_")
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


def test_analysis_matches(temp_output_dir):
    """Test that the dy and reaction rates match, and that they match the actual
    slope of abundances calculated using finite differences.
    """
    param_dict = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1.0e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e6,
        "points": 1,  # Explicitly set to 0D mode to avoid state pollution from 1D tests
        "enable_radiative_transfer": False,  # Explicitly disable to avoid state pollution from 1D tests
    }

    result = uclchem.model.Cloud(param_dict=param_dict)

    # Basic checks
    assert result is not None, "Model failed to run"
    result.check_error()

    # Get dataframe
    physics_df, abundances_df, rate_constants_df = result.get_dataframes(
        joined=False, with_rate_constants=True
    )
    network = uclchem.makerates.network.Network.from_csv()
    dy, reaction_rates_df = uclchem.analysis.rate_constants_to_dy_and_rates(
        physics_df,
        abundances_df,
        rate_constants_df,
        network=network,
    )

    species_to_check = ["H2", "OH", "CH3", "#CH3", "#CH4", "#H2O", "#H", "@CO2"]
    for species in species_to_check:
        fd_midpoints, fd_slope = finite_difference(
            physics_df["Time"], abundances_df[species]
        )
        # Convert from change in abundances per year to change in abundances per second,
        # same unit of time as what is given in rate_constant_df
        fd_slope /= uclchem.constants.SECONDS_PER_YEAR

        production, destruction = uclchem.analysis.get_production_and_destruction(
            dataframe=reaction_rates_df, species=species
        )
        assert np.allclose(
            dy[species], production.sum(axis=1) - destruction.sum(axis=1), atol=1e-30
        ), (
            f"dy and sum of production and destruction rates do not match for species {species}"
        )

        # This atol might be too tight to get right, needs to be tuned.
        assert np.allclose(dy[species].iloc[-10:-1], fd_slope[-9:], atol=1e-20), (
            f"dy and finite differences slope do not match for species {species}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
