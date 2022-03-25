import numpy as np


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


def remove_slow_reactions(changes, change_reacs,rate_threshold=0.999):
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


def print_analysis(time, total_production, total_destruction, key_reactions, key_changes):
    """Prints key reactions

    Args:
        time (float): Simulation time at which analysis is performed
        total_production (float): Total positive rate of change
        total_destruction (float): Total negative rate of change
        key_reactions (list): A list of all reactions that contribute to the total rate of change
        key_changes (list): A list of rates of change contributing to total
    """
    print("\n***************************\nNew Important Reactions At: {0:.2e} years\n".format(time))
    print("Formation = {0:.2e} from:".format(total_production))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] > 0:
            outString = f"{reaction} : {float(key_changes[k] / total_production):.2%}"
            print(outString)

    print("\nDestruction = {0:.2e} from:".format(total_destruction))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] < 0:
            outString = f"{reaction} : {float(key_changes[k] / total_destruction):.2%}"
            print(outString)


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
