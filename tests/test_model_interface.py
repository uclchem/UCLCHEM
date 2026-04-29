import numpy as np

from uclchem.model import Cloud


def test_get_final_abundances_for_species():
    model = Cloud()
    species = ["CO", "H2O", "#CH3"]
    final_abundances = model.get_final_abundances_of_species(species)

    assert len(final_abundances) == len(species)

    phys_df, chem_df = model.get_dataframes(joined=False)
    final_abundances_from_df = []
    for index, spec in enumerate(species):
        final_abundances_from_df.append(chem_df[spec].iloc[-1])
    assert np.all(final_abundances_from_df == final_abundances)
