from .analysis import total_element_abundance
from .uclchemwrap import uclchemwrap as wrap
from pandas import DataFrame
from numpy import zeros

def test_ode_conservation(species_list, element_list=["H", "N", "C", "O"]):
    """Test function which checks whether the ODEs conserve elementsry_

    :param species_list (list): list of each species in the network

    :return: (dict) Dictionary containing total rate of change of important elements
    """
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
