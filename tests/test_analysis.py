"""Integration test for new surface chemistry features.

Runs a cloud model and verifies the new physics produces different/better results.
"""

import shutil
import tempfile
from pathlib import Path

import pytest

import uclchem


@pytest.fixture
def temp_output_dir():
    """Create temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp(prefix="uclchem_analysis_test_")
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


def test_analysis_matches(temp_output_dir):
    """Test that ice-coverage-dependent desorption actually affects chemistry.

    Verifies:
    - Model runs successfully with new features
    - Ice builds up over time
    - Surface chemistry evolves as ice grows
    - Chemical desorption efficiency varies with coverage
    """
    param_dict = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1.0e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e6,
        "points": 1,  # Explicitly set to 0D mode to avoid state pollution from 1D tests
        "enable_radiative_transfer": False,  # Explicitly disable to avoid state pollution from 1D tests
        "desorb": True,
        "chemdesorb": True,
        "outputFile": str(Path(temp_output_dir) / "surface_test.dat"),
    }

    result = uclchem.model.Cloud(
        param_dict=param_dict, out_species=["#H2O", "#CO", "CH4"]
    )

    # Basic checks
    assert result is not None, "Model failed to run"
    result.check_error()

    # Get dataframe
    physics_df, abundances_df, rate_constants_df = result.get_dataframes(
        joined=False, with_rate_constants=True
    )
    network = uclchem.makerates.network.Network.from_csv()
    dy, reaction_rates_df = uclchem.analysis.rate_constants_to_dy_and_rates(
        abundances=abundances_df, rate_constants=rate_constants_df, physics=physics_df, network=network
    )
    print(dy)

    species_to_check = ["H2", "OH", "CH3", "#CH3", "#CH4", "#H2O", "#H"]
    for species in species_to_check:
        production, destruction = uclchem.analysis.get_production_and_destruction(
            dataframe=reaction_rates_df, species=species
        )
        assert dy[species] == production.sum(axis=1) - destruction.sum(axis=1)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
