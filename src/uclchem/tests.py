from .analysis import total_element_abundance
from .uclchemwrap import uclchemwrap as wrap
from pandas import DataFrame
from numpy import zeros,loadtxt
import os
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
    species_list = loadtxt(
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
    abundances, success_flag = wrap.cloud(param_dict," ".join(species_list))
    abundances=abundances[:param_dict["outspecies"]]
    param_dict.pop("outspecies")
    input_abund = zeros(500)
    input_abund[: len(abundances)] = abundances
    rates = wrap.get_odes(param_dict, input_abund)
    df = DataFrame(columns=species_list)
    df.loc[len(df)] = rates[: len(species_list)]
    result = {}
    for element in element_list:
        discrep = total_element_abundance(element, df).values
        result[element] = discrep[0]
    return result
