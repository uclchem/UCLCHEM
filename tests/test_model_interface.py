import re

import numpy as np
import pytest

from uclchem.model import (
    Cloud,
    UnknownModelError,
    UnknownParameterError,
    verify_model_type,
)

valid_model_types = [
    ("Cloud"),
    ("JShock"),
]
invalid_model_types = [
    ("cloud"),
    ("NON_cloud"),
    ("cluod"),
    (""),
]


@pytest.mark.parametrize("model_type", valid_model_types)
def test_verify_model_type(model_type):
    verify_model_type(model_type)


@pytest.mark.parametrize("model_type", invalid_model_types)
def test_verify_model_type_raises(model_type):
    with pytest.raises(UnknownModelError):
        verify_model_type(model_type)


def test_out_species_on_oo_model_raises():
    with pytest.raises(TypeError):
        Cloud(out_species=["CO"], run_type="external")


def test_get_final_abundances_for_species():
    model = Cloud()
    species = ["CO", "H2O", "#CH3"]
    final_abundances = model.get_final_abundances_of_species(species)

    assert len(final_abundances) == len(species)

    _phys_df, chem_df = model.get_dataframes(joined=False)
    final_abundances_from_df = []
    for _index, spec in enumerate(species):
        final_abundances_from_df.append(chem_df[spec].iloc[-1])
    assert np.all(final_abundances_from_df == final_abundances)


def test_unknown_parameter_raises():
    with pytest.raises(UnknownParameterError):
        Cloud({"unknown_parameter": 0})


def test_unknown_parameter_raises_temp():
    with pytest.raises(
        UnknownParameterError,
        match=re.compile(
            "Unknown parameter 'initialtemperature', did you mean 'initialtemp'?",
            flags=re.IGNORECASE,
        ),
    ):
        Cloud({"initialTemperature": 10})
