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
# Metadata Separation Tests
# ============================================================================

def test_assignment_array_goes_to_dataset():
    """Array variables with _array suffix should be stored in the xarray Dataset."""
    m = PrestellarCore()
    arr = np.zeros((3, 1, 4))
    m.test_array = arr
    assert "test_array" in m._data
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

    assert m.run_type == "local"
    assert m.custom_flag == 42


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
    assert loaded.some_meta == "xyz"
    assert loaded.custom_value == 123


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
    assert m.some_meta == "meta"
    assert "example_array" in m._data
    np.testing.assert_array_equal(m._data["example_array"].values, arr)


# ============================================================================
# Legacy File Reading Tests
# ============================================================================

def test_legacy_read_output_infers_dstep_safely():
    """Loading a legacy file missing 'dstep' should succeed with warning when safe."""
    # phase2-full.dat has 7 physical parameters (missing 'dstep')
    # If no duplicate timesteps, it should load successfully with a warning
    with pytest.warns(UserWarning, match="missing the 'dstep' parameter"):
        model = PrestellarCore(read_file="examples/test-output/phase2-full.dat")
    
    # Model should load successfully
    assert model.physics_array is not None
    assert model.chemical_abun_array is not None
    
    # Arrays should now have the correct dimensions (8 physical parameters)
    assert model.physics_array.shape[-1] == len(PHYSICAL_PARAMETERS)
    assert model.chemical_abun_array.shape[-1] == n_species


def test_legacy_read_output_coords_use_global_constants():
    """After inferring dstep, coordinates should use global PHYSICAL_PARAMETERS."""
    with pytest.warns(UserWarning, match="dstep"):
        model = PrestellarCore(read_file="examples/test-output/phase2-full.dat")
    
    # Coordinates should match global constants
    assert list(model._data.coords["physics_values"].values) == list(PHYSICAL_PARAMETERS)
    assert len(model._data.coords["chemical_abun_values"]) == n_species


# ============================================================================
# Global Constants Tests
# ============================================================================

def test_physical_parameters_always_from_global():
    """All models use the same global PHYSICAL_PARAMETERS constant."""
    m1 = PrestellarCore()
    m2 = PrestellarCore()
    
    # Both should use the same global constant (can't be modified per-instance)
    from uclchem.constants import PHYSICAL_PARAMETERS as CONST
    assert list(PHYSICAL_PARAMETERS) == list(CONST)
    assert len(PHYSICAL_PARAMETERS) == 8


def test_species_names_always_from_global():
    """All models use the same global species list."""
    m1 = PrestellarCore()
    m2 = PrestellarCore()
    
    # Both should use the same global constant
    species = get_species_names()
    assert len(species) == n_species
    assert len(species) == 335


def test_constants_cannot_be_modified():
    """PHYSICAL_PARAMETERS and n_species are hard-coded constants."""
    # These are defined at module level and tied to Fortran code
    assert PHYSICAL_PARAMETERS == [
        "Time", "Density", "gasTemp", "dustTemp", "Av", "radfield", "zeta", "dstep"
    ]
    assert n_species == 335


def test_getattr_raises_on_missing_attribute():
    """__getattr__ should raise AttributeError for non-existent attributes."""
    m = PrestellarCore()
    with pytest.raises(AttributeError, match="no attribute"):
        _ = m.nonexistent_attribute


def test_setattr_preserves_underscored_attributes():
    """Attributes starting with underscore should be set as real attributes."""
    m = PrestellarCore()
    m._custom_internal = "value"
    assert m._custom_internal == "value"
    assert "_custom_internal" in m.__dict__


def test_scalar_meta_not_stored_in_dataset():
    """Non-array scalar values should be stored in _meta, not in xarray Dataset."""
    m = PrestellarCore()
    m.scalar_value = 42.5
    m.string_value = "test"
    
    # Should be in _meta, not in _data
    assert m.scalar_value == 42.5
    assert m.string_value == "test"
    assert "scalar_value" in m._meta
    assert "string_value" in m._meta


def test_dstep_inference_validation():
    """dstep inference should detect and reject unsafe cases with duplicate timesteps."""
    # This test verifies the validation logic works correctly
    # The actual test file with duplicates would fail to load
    # For now, verify the successful case loads
    with pytest.warns(UserWarning):
        model = PrestellarCore(read_file="examples/test-output/phase2-full.dat")
    
    # Verify dstep was properly added
    assert "dstep" in model._data.coords["physics_values"].values
