"""
Test DVODE solver statistics collection and exposure through the Python API.
"""

import numpy as np
import pytest

try:
    import uclchem
    from uclchem.constants import DVODE_STAT_NAMES, N_DVODE_STATS

    uclchem_imported = True
except ImportError:
    uclchem_imported = False


def test_import_uclchem():
    if not uclchem_imported:
        pytest.fail(
            "uclchem module could not be imported, "
            "make sure your environment is loaded and UCLCHEM is installed."
        )


def test_stats_array_populated():
    """Test that stats_array is populated after a model run."""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e5,  # Short run for speed
    }
    model = uclchem.model.Cloud(param_dict=params)

    assert model.stats_array is not None, "stats_array should not be None"
    assert (
        model.stats_array.shape[2] == N_DVODE_STATS
    ), f"stats_array should have {N_DVODE_STATS} columns, got {model.stats_array.shape[2]}"

    # Check that some stats are non-zero (the solver must have done some work)
    # NST (number of steps) is column 5 (0-indexed)
    nst_values = model.stats_array[:, 0, 5]
    assert np.any(
        nst_values > 0
    ), "NST (number of steps) should be non-zero for at least some timesteps"

    # NFE (number of f evaluations) is column 6 (0-indexed)
    nfe_values = model.stats_array[:, 0, 6]
    assert np.any(
        nfe_values > 0
    ), "NFE (number of f evaluations) should be non-zero for at least some timesteps"

    # CPU_TIME is column 17 (0-indexed)
    cpu_values = model.stats_array[:, 0, 17]
    assert np.any(
        cpu_values > 0
    ), "CPU_TIME should be non-zero for at least some timesteps"


def test_stats_dataframe_columns():
    """Test that get_dataframes with_stats returns correct columns."""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e5,
    }
    model = uclchem.model.Cloud(param_dict=params)

    # Test joined DataFrame
    df = model.get_dataframes(with_stats=True)
    for stat_name in DVODE_STAT_NAMES:
        assert (
            stat_name in df.columns
        ), f"Column {stat_name} should be in joined DataFrame"

    # Test separate DataFrames
    result = model.get_dataframes(joined=False, with_stats=True)
    stats_df = result[-1]  # stats_df is the last element when with_stats=True
    assert (
        list(stats_df.columns) == DVODE_STAT_NAMES
    ), f"stats_df columns should match DVODE_STAT_NAMES, got {list(stats_df.columns)}"


def test_stats_reasonable_values():
    """Test that stat values are within reasonable ranges."""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e5,
    }
    model = uclchem.model.Cloud(param_dict=params)
    df = model.get_dataframes(with_stats=True)

    # Filter to rows where the solver actually ran (non-zero timesteps)
    active_rows = df[df["Time"] > 0]

    if len(active_rows) > 0:
        # ISTATE should be 1 or 2 for successful integrations
        istate_vals = active_rows["ISTATE"]
        assert istate_vals.min() >= -5, "ISTATE should not be extremely negative"

        # NST (steps) should be positive when solver ran
        nst_vals = active_rows["NST"]
        assert nst_vals.max() > 0, "NST should be positive for some timesteps"

        # NFE should be at least as large as NST (at least one f evaluation per step)
        nfe_vals = active_rows["NFE"]
        assert nfe_vals.max() >= nst_vals.max(), "NFE should be >= NST"


def test_stats_with_sequential_model():
    """Test that stats work with SequentialModel (chained stages)."""
    config = {
        "Cloud": {
            "param_dict": {
                "endAtFinalDensity": True,
                "freefall": True,
                "initialDens": 1e2,
                "finalDens": 1e4,
                "initialTemp": 10.0,
                "finalTime": 1.0e7,
            }
        },
    }

    seq_model = uclchem.model.SequentialModel(config)

    # Get the first (and only) stage's model
    cloud_model = seq_model.models[0]["Model"]

    assert (
        cloud_model.stats_array is not None
    ), "stats_array should exist on chained model"
    assert cloud_model.stats_array.shape[2] == N_DVODE_STATS


def test_functional_api_stats():
    """Test that functional API returns stats when requested."""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e5,
    }

    # Test return_array mode with stats
    result = uclchem.functional.cloud(
        param_dict=params,
        out_species=["CO"],
        return_array=True,
        return_stats=True,
    )
    # Result is: (phys, chem, rates, heat, stats, next_abun, flag)
    stats_array = result[4]
    assert stats_array is not None, "stats_array should not be None"
    assert stats_array.shape[2] == N_DVODE_STATS

    # Test return_dataframe mode with stats
    result_df = uclchem.functional.cloud(
        param_dict=params,
        out_species=["CO"],
        return_dataframe=True,
        return_stats=True,
    )
    # Result is: (phys_df, chem_df, rates_df, heat_df, stats_df, next_abun, flag)
    stats_df = result_df[4]
    assert stats_df is not None, "stats_df should not be None"
    assert list(stats_df.columns) == DVODE_STAT_NAMES
