"""Tests for the HDF5 storage backend (h5py-based save/load).

This branch replaces the xarray-based netCDF storage with direct h5py read/write.
Tests cover:
- Basic save/load roundtrip for model data
- Metadata and param_dict preservation
- Multiple models in a single file
- Overwrite behavior (overwrite=True vs overwrite=False)
- Error handling (missing file, missing model name)
- HDF5 file structure verification
- String/numeric dtype roundtrip
- SequentialRunner save/load with new List[Dict] format
- _write_array / _read_array low-level helpers
- from_file classmethod
"""

import json
import warnings

import h5py
import numpy as np
import pytest
import xarray as xr

try:
    from uclchem.model import (
        AbstractModel,
        Cloud,
        SequentialRunner,
        _read_array,
        load_model,
    )

    uclchem_imported = True
except ImportError:
    uclchem_imported = False


pytestmark = pytest.mark.skipif(not uclchem_imported, reason="uclchem not installed")


# ============================================================================
# Fixtures
# ============================================================================


# Default params used across fixtures
_DEFAULT_PARAMS = {
    "parcelStoppingMode": 0,
    "freefall": False,
    "writeStep": 1,
    "initialDens": 1e4,
    "initialTemp": 10.0,
    "finalDens": 1e5,
    "finalTime": 1.0e5,
}


@pytest.fixture(scope="module")
def cloud_model():
    """Create a single Cloud model to reuse across read-only tests.

    Uses a short finalTime to keep the test fast while still producing
    meaningful physics and chemistry arrays.

    WARNING: Do not mutate this model (e.g. setting attributes or calling
    save_model which is destructive). Use `fresh_cloud_model` instead.
    """
    model = Cloud(param_dict=dict(_DEFAULT_PARAMS), out_species=["H", "N", "C", "O"])
    return model


@pytest.fixture(scope="module")
def saved_model_file(cloud_model, tmp_path_factory):
    """Save the module-scoped cloud model to a temp file ONCE and
    return (path, original arrays).

    Because save_model is destructive (drops vars from _data), we snapshot the
    arrays before saving so tests can compare against the originals.
    """
    # Snapshot arrays before save_model mutates _data
    orig_physics = cloud_model.physics_array.copy()
    orig_chem = cloud_model.chemical_abun_array.copy()
    orig_param_dict = dict(cloud_model._param_dict)
    orig_success_flag = cloud_model.success_flag

    tmpdir = tmp_path_factory.mktemp("saved")
    fpath = str(tmpdir / "saved.h5")
    cloud_model.save_model(file=fpath, name="default", overwrite=True)
    return fpath, orig_physics, orig_chem, orig_param_dict, orig_success_flag


@pytest.fixture
def fresh_cloud_model():
    """Create a fresh Cloud model (costly — use only when mutation is needed)."""
    return Cloud(param_dict=dict(_DEFAULT_PARAMS), out_species=["H", "N", "C", "O"])


# ============================================================================
# Basic save / load roundtrip
# ============================================================================


class TestSaveLoadRoundtrip:
    """Verify that model data survives a save→load cycle."""

    def test_basic_roundtrip(self, saved_model_file):
        """Physics and chemistry arrays are preserved through save/load."""
        fpath, orig_physics, orig_chem, _, _ = saved_model_file
        loaded = load_model(file=fpath, name="default")

        np.testing.assert_array_equal(
            orig_physics,
            loaded.physics_array,
            err_msg="physics_array not preserved through save/load",
        )
        np.testing.assert_array_equal(
            orig_chem,
            loaded.chemical_abun_array,
            err_msg="chemical_abun_array not preserved through save/load",
        )

    def test_roundtrip_preserves_param_dict(self, saved_model_file):
        """The _param_dict should be identical after load."""
        fpath, _, _, orig_param_dict, _ = saved_model_file
        loaded = load_model(file=fpath, name="default")

        assert loaded._param_dict == orig_param_dict, (
            "_param_dict mismatch after save/load"
        )

    def test_roundtrip_preserves_metadata(self, fresh_cloud_model, tmp_path):
        """Custom scalar metadata should survive save/load."""
        model = fresh_cloud_model
        model.custom_meta_string = "test_value"
        model.custom_meta_int = 42

        fpath = str(tmp_path / "meta_test.h5")
        model.save_model(file=fpath, name="test_meta", overwrite=True)
        loaded = load_model(file=fpath, name="test_meta")

        assert loaded.custom_meta_string == "test_value"
        assert loaded.custom_meta_int == 42

    def test_from_file_classmethod(self, saved_model_file):
        """Cloud.from_file() should produce identical results to load_model()."""
        fpath, _, _, _, _ = saved_model_file
        loaded_via_classmethod = Cloud.from_file(file=fpath, name="default")
        loaded_via_function = load_model(file=fpath, name="default")

        np.testing.assert_array_equal(
            loaded_via_classmethod.physics_array,
            loaded_via_function.physics_array,
        )
        assert type(loaded_via_classmethod) is type(loaded_via_function)


# ============================================================================
# Multiple models in one file
# ============================================================================


class TestMultipleModelsInFile:
    """Test storing and loading multiple models from a single HDF5 file."""

    def test_two_models_different_names(self, tmp_path):
        """Two models saved under different names should both be loadable."""
        model_a = Cloud(
            param_dict=dict(_DEFAULT_PARAMS), out_species=["H", "N", "C", "O"]
        )
        fpath = str(tmp_path / "multi.h5")
        model_a.save_model(file=fpath, name="model_a", overwrite=True)

        model_b = Cloud(
            param_dict=dict(_DEFAULT_PARAMS, initialTemp=50.0),
            out_species=["H", "N", "C", "O"],
        )
        model_b.save_model(file=fpath, name="model_b", overwrite=True)

        loaded_a = load_model(file=fpath, name="model_a")
        loaded_b = load_model(file=fpath, name="model_b")

        # Both should be loadable and have correct initial temp
        assert loaded_a._param_dict["initialtemp"] == 10.0
        assert loaded_b._param_dict["initialtemp"] == 50.0


# ============================================================================
# Overwrite behavior
# ============================================================================


class TestOverwrite:
    """Test overwrite=True and overwrite=False semantics."""

    def test_overwrite_false_warns_and_does_not_replace(self, tmp_path):
        """overwrite=False should warn and keep existing data unchanged."""
        model = Cloud(param_dict=dict(_DEFAULT_PARAMS), out_species=["H", "N", "C", "O"])
        orig_physics = model.physics_array.copy()
        fpath = str(tmp_path / "ow.h5")
        model.save_model(file=fpath, name="dup", overwrite=True)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model2 = Cloud(
                param_dict=dict(_DEFAULT_PARAMS, initialTemp=50.0),
                out_species=["H", "N", "C", "O"],
            )
            model2.save_model(file=fpath, name="dup", overwrite=False)

        # Should have issued a warning
        assert any("already exists" in str(w.message) for w in caught), (
            "Expected warning about existing model when overwrite=False"
        )

        # Model should still be loadable (old data intact)
        loaded = load_model(file=fpath, name="dup")
        np.testing.assert_array_equal(loaded.physics_array, orig_physics)

    def test_overwrite_true_replaces_existing(self, tmp_path):
        """overwrite=True should replace the model data with new data."""
        params_v1 = {
            "parcelStoppingMode": 0,
            "freefall": False,
            "writeStep": 1,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalDens": 1e5,
            "finalTime": 1.0e5,
        }
        model_v1 = Cloud(param_dict=params_v1, out_species=["H", "N", "C", "O"])

        fpath = str(tmp_path / "overwrite.h5")
        model_v1.save_model(file=fpath, name="model", overwrite=True)

        # Create a second model with different params
        params_v2 = dict(params_v1, initialTemp=50.0)
        model_v2 = Cloud(param_dict=params_v2, out_species=["H", "N", "C", "O"])
        model_v2.save_model(file=fpath, name="model", overwrite=True)

        loaded = load_model(file=fpath, name="model")
        # The loaded model should have the v2 initial temp
        assert loaded._param_dict["initialtemp"] == 50.0, (
            "Overwritten model should have new initialTemp"
        )


# ============================================================================
# Error handling
# ============================================================================


class TestErrorHandling:
    """Test error paths in load_model."""

    def test_load_nonexistent_file_raises(self, tmp_path):
        """Loading from a file that doesn't exist should raise FileNotFoundError."""
        fake_path = str(tmp_path / "does_not_exist.h5")
        with pytest.raises(FileNotFoundError):
            load_model(file=fake_path, name="default")

    def test_load_nonexistent_model_name_raises(self, saved_model_file):
        """Loading a model name that doesn't exist should raise an Exception."""
        fpath, _, _, _, _ = saved_model_file
        with pytest.raises(Exception, match="was not found"):
            load_model(file=fpath, name="nonexistent_name")


# ============================================================================
# HDF5 file structure verification
# ============================================================================


class TestHDF5Structure:
    """Verify the internal HDF5 layout produced by save_model."""

    def test_coords_subgroup_exists(self, saved_model_file):
        """A _coords subgroup should exist inside the model group."""
        fpath, _, _, _, _ = saved_model_file
        with h5py.File(fpath, "r") as f:
            assert "_coords" in f["default"], "_coords subgroup missing"

    def test_datasets_have_dims_attr(self, saved_model_file):
        """Every dataset should have a _dims attribute."""
        fpath, _, _, _, _ = saved_model_file
        with h5py.File(fpath, "r") as f:
            group = f["default"]
            for ds_name in group:
                if ds_name == "_coords":
                    # Check coord datasets
                    for coord_name in group["_coords"]:
                        ds = group["_coords"][coord_name]
                        assert "_dims" in ds.attrs, (
                            f"Coord dataset '{coord_name}' missing _dims attr"
                        )
                else:
                    ds = group[ds_name]
                    assert "_dims" in ds.attrs, f"Dataset '{ds_name}' missing _dims attr"


# ============================================================================
# _write_array / _read_array low-level tests
# ============================================================================


class TestWriteReadArray:
    """Test the static _write_array and module-level _read_array helpers."""

    def test_numeric_array_roundtrip(self, tmp_path):
        """Numeric arrays should roundtrip exactly."""
        fpath = str(tmp_path / "arrays.h5")
        original = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        xr_var = xr.Variable(["row", "col"], original)

        with h5py.File(fpath, "w") as f:
            grp = f.create_group("test")
            AbstractModel._write_array(grp, "numeric", xr_var)

        with h5py.File(fpath, "r") as f:
            result = _read_array(f["test"], "numeric")
            np.testing.assert_array_equal(result.values, original)
            assert list(result.dims) == ["row", "col"]

    def test_unicode_string_roundtrip(self, tmp_path):
        """Unicode string arrays should be converted to bytes on write and back on read."""
        fpath = str(tmp_path / "strings.h5")
        original = np.array(["hello", "world", "test"])
        xr_var = xr.Variable(["items"], original)

        with h5py.File(fpath, "w") as f:
            grp = f.create_group("test")
            AbstractModel._write_array(grp, "strings", xr_var)

        # Verify stored as bytes
        with h5py.File(fpath, "r") as f:
            raw = f["test"]["strings"][()]
            assert raw.dtype.kind == "S", f"Expected byte string dtype, got {raw.dtype}"

        # Verify read back as str
        with h5py.File(fpath, "r") as f:
            result = _read_array(f["test"], "strings")
            np.testing.assert_array_equal(result.values, original)
            assert result.values.dtype.kind == "U", "Read-back strings should be Unicode"

    def test_scalar_in_1d_array_roundtrip(self, tmp_path):
        """A 1-element array (like JSON blobs) should roundtrip."""
        fpath = str(tmp_path / "scalar.h5")
        original = np.array([json.dumps({"key": "value"})])
        xr_var = xr.Variable(["dim_0"], original)

        with h5py.File(fpath, "w") as f:
            grp = f.create_group("test")
            AbstractModel._write_array(grp, "json_blob", xr_var)

        with h5py.File(fpath, "r") as f:
            result = _read_array(f["test"], "json_blob")
            loaded = json.loads(result.values[0])
            assert loaded == {"key": "value"}


# ============================================================================
# SequentialRunner save / load
# ============================================================================


class TestSequentialRunner:
    """Test SequentialRunner with the new List[Dict] format."""

    @pytest.fixture
    def sequential_model(self):
        """Build a two-stage sequential model (Cloud -> Cloud)."""
        config = [
            {
                "Cloud": {
                    "param_dict": {
                        "parcelStoppingMode": 1,
                        "freefall": True,
                        "initialDens": 1e2,
                        "finalDens": 1e4,
                        "initialTemp": 10.0,
                        "finalTime": 1.0e5,
                    }
                },
            },
            {
                "Cloud": {
                    "param_dict": {
                        "parcelStoppingMode": 0,
                        "freefall": False,
                        "initialDens": 1e4,
                        "finalDens": 1e4,
                        "initialTemp": 10.0,
                        "finalTime": 1.0e5,
                    }
                },
            },
        ]
        return SequentialRunner(config)

    def test_sequential_model_creates_multiple_stages(self, sequential_model):
        """SequentialRunner with two stages should produce two model entries."""
        assert len(sequential_model.models) == 2

    def test_sequential_model_save_load(self, sequential_model, tmp_path):
        """SequentialRunner save/load roundtrip should preserve all stages."""
        fpath = str(tmp_path / "sequential.h5")
        sequential_model.save_model(file=fpath, name="seq")

        # Each stage should be loadable
        for i, model_entry in enumerate(sequential_model.models):
            model_type = model_entry["Model_Type"]
            model_order = model_entry["Model_Order"]
            group_name = f"seq_{model_order}_{model_type}"
            loaded = load_model(file=fpath, name=group_name)

            np.testing.assert_array_equal(
                loaded.physics_array,
                model_entry["Model"].physics_array,
                err_msg=f"Stage {i} physics_array mismatch after save/load",
            )

    def test_sequential_model_naming_convention(self, sequential_model, tmp_path):
        """Saved group names should follow {prefix}_{order}_{type} convention."""
        fpath = str(tmp_path / "naming.h5")
        sequential_model.save_model(file=fpath, name="test")

        with h5py.File(fpath, "r") as f:
            groups = list(f.keys())

        # With two Cloud stages (orders 0 and 1):
        assert "test_0_Cloud" in groups, f"Expected 'test_0_Cloud' in {groups}"
        assert "test_1_Cloud" in groups, f"Expected 'test_1_Cloud' in {groups}"

    def test_sequential_allows_duplicate_model_types(self):
        """The new List[Dict] format should allow two Cloud stages."""
        config = [
            {
                "Cloud": {
                    "param_dict": {
                        "parcelStoppingMode": 1,
                        "freefall": True,
                        "initialDens": 1e2,
                        "finalDens": 1e4,
                        "initialTemp": 10.0,
                        "finalTime": 1.0e5,
                    }
                },
            },
            {
                "Cloud": {
                    "param_dict": {
                        "parcelStoppingMode": 0,
                        "freefall": False,
                        "initialDens": 1e4,
                        "finalDens": 1e4,
                        "initialTemp": 10.0,
                        "finalTime": 1.0e5,
                    }
                },
            },
        ]
        seq = SequentialRunner(config)
        # Both stages should be Cloud
        assert seq.models[0]["Model_Type"] == "Cloud"
        assert seq.models[1]["Model_Type"] == "Cloud"


# ============================================================================
# Save to new file (creates file)
# ============================================================================


class TestSaveCreatesFile:
    """Verify save_model creates a new file when it doesn't exist."""

    def test_save_is_non_destructive(self, tmp_path):
        """save_model should not mutate the model's internal _data dataset."""
        model = Cloud(param_dict=dict(_DEFAULT_PARAMS), out_species=["H", "N", "C", "O"])
        physics_before = model.physics_array.copy()
        chem_before = model.chemical_abun_array.copy()
        vars_before = set(model._data.variables)

        fpath = str(tmp_path / "non_destructive.h5")
        model.save_model(file=fpath, name="default")

        # Arrays should be unchanged
        np.testing.assert_array_equal(model.physics_array, physics_before)
        np.testing.assert_array_equal(model.chemical_abun_array, chem_before)

        # Dataset variables should be the same (no synthetic vars added, no scalars removed)
        vars_after = set(model._data.variables)
        assert vars_before == vars_after, (
            f"save_model mutated _data variables: removed={vars_before - vars_after}, "
            f"added={vars_after - vars_before}"
        )

    def test_save_twice_produces_same_result(self, tmp_path):
        """Saving the same model twice to different files should produce identical loads."""
        model = Cloud(param_dict=dict(_DEFAULT_PARAMS), out_species=["H", "N", "C", "O"])

        fpath1 = str(tmp_path / "first.h5")
        fpath2 = str(tmp_path / "second.h5")
        model.save_model(file=fpath1, name="default")
        model.save_model(file=fpath2, name="default")

        loaded1 = load_model(file=fpath1, name="default")
        loaded2 = load_model(file=fpath2, name="default")

        np.testing.assert_array_equal(loaded1.physics_array, loaded2.physics_array)
        np.testing.assert_array_equal(
            loaded1.chemical_abun_array, loaded2.chemical_abun_array
        )
        assert loaded1._param_dict == loaded2._param_dict


# ============================================================================
# Coordinate preservation
# ============================================================================


class TestCoordinatePreservation:
    """Verify that xarray coordinates survive the save/load cycle."""

    def test_coords_roundtrip(self, tmp_path):
        """Coordinates in the xarray Dataset should be preserved."""
        model = Cloud(param_dict=dict(_DEFAULT_PARAMS), out_species=["H", "N", "C", "O"])
        # Snapshot coords before destructive save
        orig_coords = {k: v.values.copy() for k, v in model._data.coords.items()}

        fpath = str(tmp_path / "coords.h5")
        model.save_model(file=fpath, name="default", overwrite=True)
        loaded = load_model(file=fpath, name="default")

        for coord_name, orig_values in orig_coords.items():
            assert coord_name in loaded._data.coords, (
                f"Coordinate '{coord_name}' missing after load"
            )
            np.testing.assert_array_equal(
                orig_values,
                loaded._data.coords[coord_name].values,
                err_msg=f"Coordinate '{coord_name}' values mismatch",
            )


# ============================================================================
# Data integrity checks
# ============================================================================


class TestDataIntegrity:
    """Deeper checks on the data arrays after roundtrip."""

    def test_success_flag_preserved(self, saved_model_file):
        """success_flag should be preserved through save/load."""
        fpath, _, _, _, orig_flag = saved_model_file
        loaded = load_model(file=fpath, name="default")
        assert loaded.success_flag == orig_flag

    def test_model_type_stored_in_attributes(self, saved_model_file):
        """The model_type should be stored in the attributes_dict JSON."""
        fpath, _, _, _, _ = saved_model_file
        with h5py.File(fpath, "r") as f:
            raw = f["default"]["attributes_dict"][()]
            # Convert bytes to string if needed
            if raw.dtype.kind == "S":
                raw = raw.astype(str)
            attrs = json.loads(raw[0] if raw.ndim > 0 else raw)
            assert "model_type" in attrs
            assert attrs["model_type"] == "Cloud"

    def test_param_dict_stored_as_json(self, saved_model_file):
        """_param_dict should be stored as a JSON string in the HDF5 file."""
        fpath, _, _, _, _ = saved_model_file
        with h5py.File(fpath, "r") as f:
            raw = f["default"]["_param_dict"][()]
            if raw.dtype.kind == "S":
                raw = raw.astype(str)
            params = json.loads(raw[0] if raw.ndim > 0 else raw)
            assert isinstance(params, dict)
            assert "initialdens" in params or "initialDens" in params

    def test_array_precision(self, tmp_path):
        model = Cloud(param_dict=_DEFAULT_PARAMS)

        fp64_path = tmp_path / "fp64_precision.h5"
        fp32_path = tmp_path / "fp32_precision.h5"
        fp16_path = tmp_path / "fp16_precision.h5"

        model.save_model(fp64_path, array_dtype=np.float64)
        model.save_model(fp32_path, array_dtype=np.float32)
        model.save_model(fp16_path, array_dtype=np.float16)

        fp64_size = fp64_path.stat().st_size
        fp32_size = fp32_path.stat().st_size
        fp16_size = fp16_path.stat().st_size
        assert fp32_size < fp64_size and fp16_size < fp32_size, (
            f"Expected smaller file size when saving arrays at half precision, but got sizes {fp64_size} (full), {fp32_size} (half), and {fp16_size} (quarter) in Bytes"
        )

        fp64_model = load_model(fp64_path)
        assert fp64_model.chemical_abun_array.dtype == np.float64

        fp32_model = load_model(fp32_path)
        assert fp32_model.chemical_abun_array.dtype == np.float32

        assert np.allclose(
            fp64_model.chemical_abun_array,
            fp32_model.chemical_abun_array,
            atol=0,
        )


# ============================================================================
# Engine parameter removal
# ============================================================================


class TestEngineParameterRemoved:
    """Verify that the old 'engine' parameter is no longer accepted."""

    def test_save_model_rejects_engine_kwarg(self, fresh_cloud_model, tmp_path):
        """save_model should not accept an 'engine' parameter."""
        fpath = str(tmp_path / "engine.h5")
        with pytest.raises(TypeError):
            fresh_cloud_model.save_model(file=fpath, name="default", engine="h5netcdf")

    def test_load_model_rejects_engine_kwarg(self, saved_model_file):
        """load_model should not accept an 'engine' parameter."""
        fpath, _, _, _, _ = saved_model_file
        with pytest.raises(TypeError):
            load_model(file=fpath, name="default", engine="h5netcdf")


# ============================================================================
# Entry point
# ============================================================================


def main():
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
