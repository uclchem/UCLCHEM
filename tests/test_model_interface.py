import numpy as np
import pytest

from uclchem.model import Cloud


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
