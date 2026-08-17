from pathlib import Path

import numpy as np
import pytest

from uclchem.model import Cloud, load_model


def test_out_species_on_oo_model_raises():
    with pytest.raises(TypeError):
        Cloud(out_species=["CO"], run_type="external")


def test_get_final_abundances_for_species():
    model = Cloud()
    species = ["CO", "H2O", "#CH3"]
    final_abundances = model.get_final_abundances_of_species(species)

    assert len(final_abundances) == len(species)

    _phys_df, chem_df = model.get_dataframes(joined=False)
    final_abundances_from_df = [chem_df[spec].iloc[-1] for spec in species]
    assert np.all(final_abundances_from_df == final_abundances)


def test_model_precision_saving_storage_size(tmp_path: Path):
    model = Cloud()
    model_df = model.get_joined_dataframes()

    dtypes = ["fp32", "fp64", "fp128"]
    paths = [tmp_path / f"model_{dtype}.h5" for dtype in dtypes]
    for dtype, path in zip(dtypes, paths, strict=True):
        model.save_model(file=path, float_dtype=dtype)

        loaded_model = load_model(file=path)
        loaded_df = loaded_model.get_joined_dataframes()

        assert loaded_df.shape == model_df.shape

        assert np.allclose(loaded_df.to_numpy(), model_df.to_numpy())

    file_sizes = [path.stat().st_size for path in paths]
    assert sorted(file_sizes) == file_sizes


def test_model_precision_overflow(tmp_path: Path):
    model = Cloud()
    with pytest.raises(OverflowError):
        model.save_model(file=tmp_path / "model_fp16.f5", float_dtype="fp16")
