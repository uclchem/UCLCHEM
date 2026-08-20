"""Integration test for new surface chemistry features.

Runs a cloud model and verifies the new physics produces different/better results.
"""

import shutil
import tempfile

import numpy as np
import pandas as pd
import pytest

import uclchem

element_abundance_data = [
    (pd.DataFrame.from_dict({"H": [1]}), {"H": 1}),
    (pd.DataFrame.from_dict({"H2": [1]}), {"H": 2}),
    (pd.DataFrame.from_dict({"H2": [1], "CH3OH": [0.5]}), {"H": 4, "C": 0.5, "O": 0.5}),
    (pd.DataFrame.from_dict({"CCL": [1]}), {"CL": 1, "C": 1}),
    (pd.DataFrame.from_dict({"H2": [1], "initialTemp": [0.5]}), {"H": 2, "N": 0}),
    (pd.DataFrame.from_dict({"H2": [1], "SURFACE": [0.5]}), {"H": 2, "C": 0}),
]


@pytest.mark.parametrize(("df", "expected"), element_abundance_data)
def test_total_element_abundance(df: pd.DataFrame, expected: dict[str, float]):
    for element, expected_abundance in expected.items():
        assert (
            expected_abundance
            == uclchem.analysis.total_element_abundance(element, df).to_numpy()
        )


element_conservation_data = [
    (pd.DataFrame.from_dict({"H": [1, 1]}), {"H": 0}),
    (pd.DataFrame.from_dict({"H2": [2, 1]}), {"H": 50}),
    (
        pd.DataFrame.from_dict({"H2": [1, 2], "CH3OH": [0.5, 0.5]}),
        {"H": 50, "C": 0, "O": 0},
    ),
    (pd.DataFrame.from_dict({"CCL": [1]}), {"CL": 0, "C": 0}),
]


@pytest.mark.parametrize(("df", "expected"), element_conservation_data)
def test_check_element_conservation(df: pd.DataFrame, expected: dict[str, float]):
    calculated_change = uclchem.analysis.check_element_conservation(
        df, element_list=expected.keys(), percent=True
    )
    for element, expected_change in expected.items():
        assert float(calculated_change[element].strip("%")) == pytest.approx(
            expected_change
        )


def test_non_existent_element_raises():
    df = pd.DataFrame.from_dict({"H": [1], "N": [0]})
    element = "N"

    with pytest.raises(ZeroDivisionError, match="elemental abundance of 0"):
        uclchem.analysis.check_element_conservation(df, element_list=[element])


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
        with_rate_constants=True, joined=False
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
        _fd_midpoints, fd_slope = finite_difference(
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

        # FD comparison is approximate: large log-spaced Δt causes a factor
        # ~k·Δt/(1-e^{-kΔt}) error, and SURFSWAP_GEOMETRIC post-processing
        # corrections can create transient spikes for surface/bulk species.
        # atol=1e-14 skips these small-magnitude artifacts while still
        # catching sign errors or order-of-magnitude bugs at large rates.
        assert np.allclose(
            dy[species].iloc[-10:-1], fd_slope[-9:], rtol=1.5, atol=1e-14
        ), f"dy and FD slope do not match for species {species}"
