from .analysis import total_element_abundance
from .wrap import cloud,get_odes
from pandas import DataFrame
from numpy import zeros

def test_ode_conservation(species_list, element_list=["H", "N", "C", "O"]):
    """Test function which checks whether the ODEs conserve elementsry_

    :param species_list (list): list of each species in the network

    :return: (dict) Dictionary containing total rate of change of important elements
    """
    param_dict = {
        "phase": 1,
        "switch": 0,
        "collapse": 1,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 1.0e3,
        "outSpecies": " ".join(species_list),
    }
    abundances = cloud(param_dict)
    param_dict.pop("outSpecies")
    input_abund = zeros(500)
    input_abund[: len(abundances)] = abundances
    rates = get_odes(param_dict, input_abund)
    df = DataFrame(columns=species_list)
    df.loc[len(df)] = rates[: len(species_list)]
    result = {}
    for element in element_list:
        discrep = total_element_abundance(element, df).values
        result[element] = discrep[0]
    return result
