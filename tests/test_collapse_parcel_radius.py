"""
Tests for collapse model parcel_radius output.

Verifies that parcel_radius is:
  - present in the physics DataFrame for all models
  - non-zero and monotonically ordered across parcels for collapse models
  - zero for non-collapse (cloud) models
"""

import numpy as np
import pytest

import uclchem

COLLAPSE_PARAMS = {
    # "initialDens": 1e4,
    # "finalDens": 1e6,
    "initialTemp": 10.0,
    "endAtFinalDensity": True,
}

CLOUD_PARAMS = {
    "freefall": False,
    "initialDens": 1e4,
    "finalDens": 1e5,
    "finalTime": 5e5,
}


@pytest.mark.parametrize("collapse_mode", ["BE1.1", "BE4", "filament", "ambipolar"])
def test_parcel_radius_in_collapse_physics(collapse_mode):
    """parcel_radius column is present and non-zero for collapse models."""
    physics, chemistry, rates, heating, abundances, return_code = (
        uclchem.functional.collapse(
            collapse=collapse_mode,
            param_dict=COLLAPSE_PARAMS,
            out_species=["CO", "OH"],
            return_dataframe=True,
        )
    )

    assert return_code == 0, f"Collapse ({collapse_mode}) failed with code {return_code}"
    assert "parcel_radius" in physics.columns, (
        f"'parcel_radius' column missing from physics DataFrame for collapse={collapse_mode}"
    )

    radii = physics["parcel_radius"].values
    assert np.any(radii > 0), (
        f"All parcel_radius values are zero for collapse={collapse_mode}; expected non-zero radii"
    )


def test_parcel_radius_decreases_during_collapse():
    """For a collapse model with points=1, the single parcel moves inward over time."""
    physics, chemistry, rates, heating, abundances, return_code = (
        uclchem.functional.collapse(
            collapse="BE4",
            param_dict=COLLAPSE_PARAMS,
            out_species=["CO"],
            return_dataframe=True,
        )
    )

    assert return_code == 0, f"Collapse failed with code {return_code}"
    radii = physics["parcel_radius"].values
    # Radius should decrease (or at least not increase) over time
    assert radii[0] >= radii[-1], (
        f"Parcel radius did not decrease during collapse: r_initial={radii[0]:.4f}, r_final={radii[-1]:.4f}"
    )


def test_parcel_radius_zero_for_cloud():
    """parcel_radius is zero for a non-collapse cloud model."""
    physics, chemistry, rates, heating, abundances, return_code = (
        uclchem.functional.cloud(
            param_dict=CLOUD_PARAMS,
            out_species=["CO", "OH"],
            return_dataframe=True,
            timepoints=5000,
        )
    )

    assert return_code == 0, f"Cloud model failed with code {return_code}"
    assert "parcel_radius" in physics.columns, (
        "'parcel_radius' column missing from physics DataFrame for cloud model"
    )

    radii = physics["parcel_radius"].values
    assert np.all(radii == 0.0), (
        f"Expected all parcel_radius to be 0 for cloud model, got max={radii.max():.4e}"
    )


def test_parcel_radius_initial_value_matches_rout():
    """Initial parcel radius should be approximately rout (outermost shell = rout)."""
    rout = 0.05  # default rout in parsec
    params = {**COLLAPSE_PARAMS, "rout": rout}

    physics, chemistry, rates, heating, abundances, return_code = (
        uclchem.functional.collapse(
            collapse="BE4",
            param_dict=params,
            out_species=["CO"],
            return_dataframe=True,
        )
    )

    assert return_code == 0, f"Collapse failed with code {return_code}"
    initial_radius = physics["parcel_radius"].iloc[0]
    assert abs(initial_radius - rout) < 0.1 * rout, (
        f"Initial parcel_radius {initial_radius:.4f} deviates more than 10% from rout={rout}"
    )


if __name__ == "__main__":
    pytest.main(["-v", __file__])
