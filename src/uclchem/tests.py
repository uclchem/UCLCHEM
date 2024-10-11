import os

import numpy as np
from pandas import DataFrame

import uclchem

from .analysis import total_element_abundance
from .uclchemwrap import uclchemwrap as wrap

_ROOT = os.path.dirname(os.path.abspath(__file__))


def test_ode_conservation(element_list=["H", "N", "C", "O"]):
    """Test whether the ODEs conserve elements. Useful to run each time you change network.
    Integrator errors may still cause elements not to be conserved but they cannot be conserved
    if the ODEs are not correct.

    Args:
        element_list (list, optional): A list of elements for which to check the conservation. Defaults to ["H", "N", "C", "O"].

    Returns:
        dict: A dictionary of the elements in element list with values representing the total rate of change of each element.
    """
    species_list = np.loadtxt(
        os.path.join(_ROOT, "species.csv"),
        usecols=[0],
        dtype=str,
        skiprows=1,
        unpack=True,
        delimiter=",",
        comments="%",
    )
    species_list = list(species_list)
    param_dict = {
        "endatfinaldensity": False,
        "freefall": True,
        "initialdens": 1e4,
        "initialtemp": 10.0,
        "finaldens": 1e5,
        "finaltime": 1.0e3,
        "outspecies": len(species_list),
    }
    _, _, abundances, specname, success_flag = wrap.cloud(
        dictionary=param_dict,
        outspeciesin=" ".join(species_list),
        timepoints=1,
        gridpoints=1,
        physicsarray=np.zeros(
            shape=(2, 1, uclchem.constants.N_PHYSICAL_PARAMETERS),
            dtype=np.float64,
            order="F",
        ),
        chemicalabunarray=np.zeros(
            shape=(2, 1, uclchem.constants.n_species),
            dtype=np.float64,
            order="F",
        ),
        returnarray=False,
        givestartabund=False,
    )
    abundances = abundances[: param_dict["outspecies"]]
    param_dict.pop("outspecies")
    input_abund = np.zeros(uclchem.constants.n_species, dtype=np.float64, order="F")
    input_abund[: len(abundances)] = abundances
    rates = wrap.get_odes(param_dict, input_abund)
    # Explicitely clean off the last element, for some reason the shape
    # of the ODE is n_species + 1, and we don't need the last one.
    rates = rates[: len(species_list)]
    df = DataFrame(rates.reshape(1, -1), columns=species_list)
    result = {}
    for element in element_list:
        discrep = total_element_abundance(element, df).values
        result[element] = discrep[0]
    return result
