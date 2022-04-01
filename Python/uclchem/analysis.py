import numpy as np
elementList = [
    "H",
    "D",
    "HE",
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "CL",
    "LI",
    "NA",
    "MG",
    "SI",
    "PAH",
    "15N",
    "13C",
    "18O",
    "SURFACE",
    "BULK",
]

def get_rates_of_change(rates, reactions, speciesList, species, row):
    """
    Calculate the terms in the rate of equation of a particular species using rates calculated using
    get_species_rates() and a row from the full output of UCLCHEM. See `analysis.py` for intended use.

    :param rates: (ndarray) Array of all reaction rates for a species, as obtained from `get_species_rates`
    :param reactions: (ndarray) Array containing reactions as lists of reactants and products
    :param speciesList: (ndarray) Array of all species names in network
    :param species: (str) Name of species for which rates of change are calculated
    :param row: (panda series) The UCLCHEM output row from the time step at which you want the rate of change
    """
    changes = []
    reactionList = []
    three_phase = "@" in "".join(speciesList)
    for i, reaction in enumerate(reactions):
        change = rates[i]
        reactants = reaction[0:3]
        products = reaction[3:]

        reactant_count = 0
        for reactant in reactants:
            if reactant in speciesList:
                change = change * row[reactant]
                reactant_count += 1
            elif reactant in ["DESOH2", "FREEZE", "LH", "LHDES"]:
                reactant_count += 1

            if reactant in ["DEUVCR", "DESCR", "DESOH2"]:
                change = change / np.max([1.0e-30, row["SURFACE"]])
            if (not three_phase) and (reactant in ["THERM"]):
                change = change * row[reaction[0]] / np.max([1.0e-30, row["SURFACE"]])

        change = change * (row["Density"] ** (reactant_count - 1))

        if species in reactants:
            changes.append(-change)
            reactionList.append(reaction)
        if species in products:
            changes.append(change)
            reactionList.append(reaction)
    A = zip(changes, reactionList)
    A = sorted(A, key=lambda x: np.abs(x[0]), reverse=True)
    changes, reactionList = zip(*A)
    changes = np.asarray(changes)
    return reactionList, changes


def remove_slow_reactions(changes, change_reacs,rate_threshold=0.99):
    """Iterates through a list of reactions adding the fastest reactions to a list until some threshold fraction of the total
    rate of change is reached. This list is returned so that you have the list of reactions that cause rate_threshold of the
    total destruction and formation of a species.

    Args:
        changes (list): List of rates of change due to each reaction a species is involved in
        change_reacs (list): List of corresponding rates of change
        rate_threshold (float, optional): Percentage of overall rate of change to consider before ignoring less important reactions. Defaults to 0.999.

    Returns:
        Total production and destruction rates as a well as list of reactions and rates of change for top rate_threshold reactiosn_
    """
    totalDestruct = sum(changes[np.where(changes < 0)])
    totalProd = sum(changes[np.where(changes > 0)])

    key_reactions = []
    key_changes = []
    form = 0.0
    destruct = 0.0

    for i, reaction in enumerate(change_reacs):
        if (changes[i] > 0) and (form < rate_threshold * totalProd):
            form = form + changes[i]
            key_reactions.append(reaction)
            key_changes.append(changes[i])
        elif (changes[i] < 0) and (abs(destruct) < rate_threshold * abs(totalDestruct)):
            destruct = destruct + changes[i]
            key_reactions.append(reaction)
            key_changes.append(changes[i])

    return totalProd, totalDestruct, key_reactions, key_changes


def write_analysis(output_file,time, total_production, total_destruction, key_reactions, key_changes):
    """Prints key reactions

    Args:
        time (float): Simulation time at which analysis is performed
        total_production (float): Total positive rate of change
        total_destruction (float): Total negative rate of change
        key_reactions (list): A list of all reactions that contribute to the total rate of change
        key_changes (list): A list of rates of change contributing to total
    """
    output_file.write("\n\n***************************\nNew Important Reactions At: {0:.2e} years\n".format(time))
    output_file.write("Formation = {0:.2e} from:".format(total_production))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] > 0:
            outString = f"\n{reaction} : {float(key_changes[k] / total_production):.2%}"
            output_file.write(outString)

    output_file.write("\n\nDestruction = {0:.2e} from:".format(total_destruction))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] < 0:
            outString = f"\n{reaction} : {float(key_changes[k] / total_destruction):.2%}"
            output_file.write(outString)


def format_reactions(reactions):
    """Turn a row of the reaction file into a string.

    Args:
        reactions (list): list of lists, each reaction read from reaction.csv is a list.

    Returns:
        list: list of string, each reaction in readable string form
    """
    formatted_reactions = []
    for reaction in reactions:
        outString = "{x[0]} + {x[1]} + {x[2]} -> {x[3]} + {x[4]} + {x[5]}".format(x=reaction)
        outString = outString.replace(" + NAN", "")
        formatted_reactions.append(outString)
    return formatted_reactions

def count_element(species_list, element):
    """
    Count the number of atoms of an element that appear in each of a list of species,
    return the array of counts

    :param  species_list: (iterable, str), list of species names
    :param element: (str), element

    :return: sums (ndarray) array where each element represents the number of atoms of the chemical element in the corresponding element of species_list
    """
    species_list = pd.Series(species_list)
    # confuse list contains elements whose symbols contain the target eg CL for C
    # We count both sets of species and remove the confuse list counts.
    confuse_list = [x for x in elementList if element in x]
    confuse_list = sorted(confuse_list, key=lambda x: len(x), reverse=True)
    confuse_list.remove(element)
    sums = species_list.str.count(element)
    for i in range(2, 10):
        sums += np.where(species_list.str.contains(element + f"{i:.0f}"), i - 1, 0)
    for spec in confuse_list:
        sums += np.where(species_list.str.contains(spec), -1, 0)
    return sums


def total_element_abundance(element, df):
    """
    Calculates that the total elemental abundance of a species as a function of time. Allows you to check conservation.

    :param element: (str) Element symbol. eg "C"
    :param df: (pandas dataframe) UCLCHEM output in format from `read_output_file`

    :return: Series containing the total abundance of an element at every time step of your output
    """
    sums = count_element(df.columns, element)
    for variable in ["Time", "Density", "gasTemp", "av", "point", "SURFACE", "BULK"]:
        sums = np.where(df.columns == variable, 0, sums)
    return df.mul(sums, axis=1).sum(axis=1)


def check_element_conservation(df, element_list=["H", "N", "C", "O"]):
    """
    Check the conservation of major element by comparing total abundance at start and end of model

    :param	df: (pandas dataframe): UCLCHEM output in format from `read_output_file`

    :return: (dict) Dictionary containing the change in the total abundance of each element as a fraction of initial value
    """
    result = {}
    for element in element_list:
        discrep = total_element_abundance(element, df).values
        discrep = np.abs(discrep[0] - discrep[-1]) / discrep[0]
        result[element] = f"{discrep:.3%}"
    return result