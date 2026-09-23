import numpy as np
import pytest

from uclchem.constants import n_species
from uclchem.model import Cloud


def test_out_species_on_oo_model_raises():
    with pytest.raises(TypeError):
        Cloud(out_species=["CO"], run_type="external")


def test_get_final_abundances_for_species():
    model = Cloud(param_dict={"finalTime": 1e0})
    species = ["CO", "H2O", "#CH3"]
    final_abundances = model.get_final_abundances_of_species(species)

    assert len(final_abundances) == len(species)

    phys_df, chem_df = model.get_dataframes(joined=False)
    final_abundances_from_df = []
    for _index, spec in enumerate(species):
        final_abundances_from_df.append(chem_df[spec].iloc[-1])
    assert np.all(final_abundances_from_df == final_abundances)


def test_starting_chemistry_array():
    starting_chemistry = np.random.random(n_species)
    model = Cloud(
        param_dict={"finalTime": 1e0},
        starting_chemistry=starting_chemistry,
    )

    _phys_df, chem_df = model.get_dataframes(joined=False)

    assert np.all(chem_df.iloc[0].to_numpy()[1:] == starting_chemistry)


def test_different_run_type_give_same_result():
    model_managed = Cloud(param_dict={"finalTime": 1e0}, run_type="managed")
    df_managed = model_managed.get_joined_dataframes()

    model_external = Cloud(param_dict={"finalTime": 1e0}, run_type="external")
    model_external.run()
    df_external = model_external.get_joined_dataframes()

    assert df_managed.equals(df_external)


def test_different_run_type_give_same_result_starting_chem():
    starting_chemistry = np.random.random(n_species)

    model_managed = Cloud(
        param_dict={"finalTime": 1e0},
        starting_chemistry=starting_chemistry,
        run_type="managed",
    )
    df_managed = model_managed.get_joined_dataframes()

    model_external = Cloud(
        param_dict={"finalTime": 1e0},
        starting_chemistry=starting_chemistry,
        run_type="external",
    )
    model_external.run()
    df_external = model_external.get_joined_dataframes()

    assert df_managed.equals(df_external)
