import numpy as np
import pandas as pd
import pytest
import uclchemwrap

import uclchem.constants
from uclchem.advanced.runtime_network import RuntimeNetwork
from uclchem.analysis import total_element_abundance


def test_get_ionlist():
    network = RuntimeNetwork()
    species_list = network.get_species_list()
    ion_list, n_ions = uclchemwrap.chemistry.get_ion_indices()

    assert np.count_nonzero(ion_list) == n_ions
    for ion_index in ion_list:
        if ion_index != 0:
            assert species_list[ion_index - 1].is_ion()


def test_get_total_elemental_abundance():
    rng = np.random.default_rng()
    abundances = rng.random(size=uclchem.constants.n_species)
    columns = [name.decode("utf-8").strip() for name in uclchemwrap.network.specname]
    abundances_df = pd.DataFrame([abundances], columns=columns)

    elemental_abundances_fortran = uclchemwrap.chemistry.calculate_elemental_abundances(
        abundances
    )
    for idx, element in enumerate(uclchemwrap.network.elem_names):
        element = element.decode("utf-8").strip()

        assert total_element_abundance(
            element, abundances_df
        ).to_numpy() == pytest.approx(elemental_abundances_fortran[idx])
