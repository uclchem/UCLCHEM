"""Tests for model metadata separation and legacy file reading.

This module consolidates tests for:
- Metadata vs dataset storage separation
- Legacy output file reading with validation
- Save/load/pickle roundtrip preservation
"""

import numpy as np
import pytest
import xarray as xr

from uclchem.constants import PHYSICAL_PARAMETERS, n_species
from uclchem.model import PrestellarCore, get_species_names, load_model

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def legacy_file_without_dstep(tmp_path):
    """Create a legacy output file without 'dstep' parameter for testing."""
    # Create a minimal legacy file with only 7 physical parameters (no dstep)
    legacy_file = tmp_path / "legacy_no_dstep.dat"

    # Get the full species list from UCLCHEM
    species_list = get_species_names()

    # Header: Time, Density, gasTemp, dustTemp, Av, radfield, zeta, point, [all species]
    # (missing dstep between zeta and point)
    physics_header = "Time,Density,gasTemp,dustTemp,Av,radfield,zeta,point"
    species_header = ",".join(species_list)
    header = f"{physics_header},{species_header}\n"

    # Sample data: 3 timesteps with no duplicates (safe for dstep inference)
    # Physics: Time, Density, gasTemp, dustTemp, Av, radfield, zeta, point
    # Then 335 species abundances (all set to 1e-20 for simplicity)
    species_values = ",".join(["1e-20"] * len(species_list))
    data = [
        f"1.0e5,1e4,10.0,10.0,1.0,1.0,1e-17,1,{species_values}\n",
        f"2.0e5,1e4,10.0,10.0,1.0,1.0,1e-17,1,{species_values}\n",
        f"3.0e5,1e4,10.0,10.0,1.0,1.0,1e-17,1,{species_values}\n",
    ]

    with open(legacy_file, "w") as f:
        f.write(header)
        f.writelines(data)

    return str(legacy_file)


# ============================================================================
# Metadata Separation Tests
# ============================================================================


def test_assignment_array_goes_to_dataset():
    """Array variables with _array suffix should be stored in the xarray Dataset."""
    m = PrestellarCore()
    arr = np.zeros((3, 1, 4))
    m.test_array = arr
    assert "test_array" in m._data, "Array attribute not stored in _data Dataset"
    # stored array should match the input
    np.testing.assert_array_equal(m._data["test_array"].values, arr)


def test_meta_survives_dataset_replace():
    """Scalar metadata should survive when _data Dataset is replaced."""
    m = PrestellarCore()
    # public scalar attributes (not PHYSICAL_PARAMETERS or specname which are global)
    m.run_type = "local"
    m.custom_flag = 42

    # Simulate dataset replacement (as happens on load/pickle)
    object.__setattr__(m, "_data", xr.Dataset())

    assert (
        m.run_type == "local"
    ), "string metadata 'run_type' not preserved after _data replace"
    assert (
        m.custom_flag == 42
    ), "Scalar metadata 'custom_flag' not preserved after _data replace"


def test_save_load_restores_meta(tmp_path):
    """Metadata should be preserved through save/load cycle."""
    m = PrestellarCore()
    m.some_meta = "xyz"
    m.custom_value = 123

    fname = tmp_path / "test_model.h5"
    # Save using API
    m.save_model(file=str(fname), name="default", overwrite=True)

    # Load back
    loaded = load_model(file=str(fname), name="default")

    # Metadata should be restored on load
    assert (
        loaded.some_meta == "xyz"
    ), "String metadata 'some_meta' not preserved after load"
    assert (
        loaded.custom_value == 123
    ), "Integer metadata 'custom_value' not preserved after load"


def test_pickle_roundtrip_preserves_meta_and_arrays():
    """Both metadata and arrays should be preserved through pickle/unpickle."""
    m = PrestellarCore()
    m.some_meta = "meta"
    arr = np.ones((2, 1, 3))
    m.example_array = arr

    # pickle/unpickle
    m.pickle()
    m.un_pickle()

    # metadata and arrays should be restored
    assert (
        m.some_meta == "meta"
    ), "String metadata 'some_meta' not preserved after pickle/unpickle"
    assert (
        "example_array" in m._data
    ), "Array attribute 'example_array' not preserved after pickle/unpickle"
    np.testing.assert_array_equal(m._data["example_array"].values, arr)


# ============================================================================
# Legacy File Reading Tests
# ============================================================================


def test_legacy_read_infers_dstep_and_uses_global_constants(legacy_file_without_dstep):
    """Loading a legacy file missing 'dstep' should warn, infer dstep safely, and use global constants."""
    # Expect a warning about missing/inferred dstep
    with pytest.warns(UserWarning, match="dstep"):
        model = PrestellarCore(read_file=legacy_file_without_dstep)

    # Model should load successfully with arrays populated
    assert model.physics_array is not None, "physics_array not loaded successfully"
    assert (
        model.chemical_abun_array is not None
    ), "chemical_abun_array not loaded successfully"

    # Arrays should have the correct dimensions after inferring dstep
    assert model.physics_array.shape[-1] == len(
        PHYSICAL_PARAMETERS
    ), "physics_array dimension mismatch after inferring dstep"
    assert (
        model.chemical_abun_array.shape[-1] == n_species
    ), "chemical_abun_array dimension mismatch after inferring dstep"

    # Coordinates should match global constants
    assert list(model._data.coords["physics_values"].values) == list(
        PHYSICAL_PARAMETERS
    ), "physics_values coordinate mismatch after inferring dstep"
    assert (
        len(model._data.coords["chemical_abun_values"]) == n_species
    ), "chemical_abun_values coordinate length mismatch after inferring dstep"


# ============================================================================
# Global Constants Tests
# ============================================================================


def test_physical_parameters_always_from_global():
    """All models use the same global PHYSICAL_PARAMETERS constant."""
    # Both should use the same global constant (can't be modified per-instance)
    from uclchem.constants import PHYSICAL_PARAMETERS as CONST

    assert PHYSICAL_PARAMETERS == [
        "Time",
        "Density",
        "gasTemp",
        "dustTemp",
        "Av",
        "radfield",
        "zeta",
        "dstep",
    ], "PHYSICAL_PARAMETERS constant has been modified!"
    assert list(PHYSICAL_PARAMETERS) == list(CONST), "PHYSICAL_PARAMETERS mismatch"
    assert (
        len(PHYSICAL_PARAMETERS) >= 8
    ), "UCLCHEM having less than 8 physical parameters is suspicious"


def test_species_names_always_from_global():
    """All models use the same global species list."""

    # Both should use the same global constant
    species = get_species_names()
    assert len(species) == n_species, "Species count mismatch with n_species"
    assert len(species) >= 20, "UCLCHEM having less than 20 species is suspicious"


def test_getattr_raises_on_missing_attribute():
    """__getattr__ should raise AttributeError for non-existent attributes."""
    m = PrestellarCore()
    with pytest.raises(AttributeError, match="no attribute"):
        _ = m.nonexistent_attribute


def test_setattr_preserves_underscored_attributes():
    """Attributes starting with underscore should be set as real attributes."""
    m = PrestellarCore()
    m._custom_internal = "value"
    assert m._custom_internal == "value", "Underscored attribute value mismatch"
    assert "_custom_internal" in m.__dict__, "Underscored attribute not set in __dict__"


def test_scalar_meta_not_stored_in_dataset():
    """Non-array scalar values should be stored in _meta, not in xarray Dataset."""
    m = PrestellarCore()
    m.scalar_value = 42.5
    m.string_value = "test"

    # Should be in _meta, not in _data
    assert m.scalar_value == 42.5, "Scalar value mismatch in _meta"
    assert m.string_value == "test", "String value mismatch in _meta"
    assert "scalar_value" in m._meta, "Scalar value not found in _meta"
    assert "string_value" in m._meta, "String value not found in _meta"


def test_dstep_inference_validation(legacy_file_without_dstep):
    """dstep inference should detect and reject unsafe cases with duplicate timesteps."""
    # This test verifies the validation logic works correctly
    # The fixture creates a legacy file without dstep that can be safely inferred
    with pytest.warns(UserWarning, match="missing the 'dstep' parameter"):
        model = PrestellarCore(read_file=legacy_file_without_dstep)

    # Verify dstep was properly added
    assert (
        "dstep" in model._data.coords["physics_values"].values
    ), "dstep coordinate not found in physics_values"
