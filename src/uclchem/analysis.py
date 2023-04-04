try:
    from .uclchemwrap import uclchemwrap as wrap
except:
    pass
import numpy as np
from pandas import Series, read_csv
from seaborn import color_palette
import matplotlib.pyplot as plt
import os

_ROOT = os.path.dirname(os.path.abspath(__file__))

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


def read_output_file(output_file):
    """Read the output of a UCLCHEM run created with the outputFile parameter into a pandas DataFrame

    Args:
        output_file (str): path to file containing a full UCLCHEM output

    Returns:
        pandas.DataFrame: A dataframe containing the abundances and physical parameters of the model at every time step.
    """
    f = open(output_file)
    f.readline()
    bits = f.readline().split()
    radfield = float(bits[1])
    data = read_csv(f)
    data["radfield"] = radfield
    data.columns = data.columns.str.strip()
    return data


def create_abundance_plot(df, species, figsize=(16, 9), plot_file=None):
    """Create a plot of the abundance of a list of species through time.

    Args:
        df (pd.DataFrame): Pandas dataframe containing the UCLCHEM output, see `read_output_file`
        species (list): list of strings containing species names. Using a $ instead of # or @ will plot the sum of surface and bulk abundances.
        figsize (tuple, optional): Size of figure, width by height in inches. Defaults to (16, 9).
        plot_file (str, optional): Path to file where figure will be saved. If None, figure is not saved. Defaults to None.

    Returns:
        fig,ax: matplotlib figure and axis objects
    """
    fig, ax = plt.subplots(figsize=figsize, tight_layout=True)

    ax = plot_species(ax, df, species)
    ax.legend(loc=4, fontsize="small")

    ax.set_xlabel("Time / years")
    ax.set_ylabel("X$_{Species}$")

    ax.set_yscale("log")
    if plot_file is not None:
        fig.savefig(plot_file)
    return fig, ax


def plot_species(ax, df, species, **plot_kwargs):
    """Plot the abundance of a list of species through time directly onto an axis.

    Args:
        ax (pyplot.axis): An axis object to plot on
        df (pd.DataFrame): A dataframe created by `read_output_file`
        species (str): A list of species names to be plotted. If species name starts with "$" instead of # or @, plots the sum of surface and bulk abundances

    Returns:
        pyplot.axis: Modified input axis is returned
    """

    color_palette(n_colors=len(species))
    for specIndx, specName in enumerate(species):
        linestyle = "solid"
        if specName[0] == "$":
            abundances = df[specName.replace("$", "#")]
            linestyle = "dashed"
            if specName.replace("$", "@") in df.columns:
                abundances = abundances + df[specName.replace("$", "@")]
        else:
            abundances = df[specName]
        ax.plot(
            df["Time"],
            abundances,
            label=specName,
            lw=2,
            linestyle=linestyle,
            **plot_kwargs,
        )
        ax.set(yscale="log")
        ax.legend()
    return ax


def analysis(species_name, result_file, output_file, rate_threshold=0.99):
    """A function which loops over every time step in an output file and finds the rate of change of a species at that time due to each of the reactions it is involved in.
    From this, the most important reactions are identified and printed to file. This can be used to understand the chemical reason behind a species' behaviour.

    Args:
        species_name (str): Name of species to be analysed
        result_file (str): The path to the file containing the UCLCHEM output
        output_file (str): The path to the file where the analysis output will be written
        rate_threshold (float,optional): Analysis output will contain the only the most efficient reactions that are responsible for rate_threshold of the total production and destruction rate. Defaults to 0.99.
    """
    result_df = read_output_file(result_file)
    species = np.loadtxt(
        os.path.join(_ROOT, "species.csv"),
        usecols=[0],
        dtype=str,
        skiprows=1,
        unpack=True,
        delimiter=",",
        comments="%",
    )
    species = list(species)
    reactions = np.loadtxt(
        os.path.join(_ROOT, "reactions.csv"),
        dtype=str,
        skiprows=1,
        delimiter=",",
        usecols=[0, 1, 2, 3, 4, 5, 6],
        comments="%",
    )

    fortran_reac_indxs = [
        i + 1 for i, reaction in enumerate(reactions) if species_name in reaction
    ]
    reac_indxs = [i for i, reaction in enumerate(reactions) if species_name in reaction]
    species_index = species.index(species_name) + 1  # fortran index of species
    old_key_reactions = []
    old_total_destruct = 0.0
    old_total_form = 0.0
    formatted_reacs = _format_reactions(reactions[reac_indxs])

    with open(output_file, "w") as f:
        f.write("All Reactions\n************************\n")
        for reaction in formatted_reacs:
            f.write(reaction + "\n")
        for i, row in result_df.iterrows():
            # recreate the parameter dictionary needed to get accurate rates
            param_dict = _param_dict_from_output(row)

            # get the rate of all reactions from UCLCHEM along with a few other necessary values
            rates, transfer, swap, bulk_layers = _get_species_rates(
                param_dict, row[species], species_index, fortran_reac_indxs
            )

            # convert reaction rates to total rates of change, this needs manually updating when you add new reaction types!
            change_reacs, changes = _get_rates_of_change(
                rates,
                reactions[reac_indxs],
                species,
                species_name,
                row,
                swap,
                bulk_layers,
            )

            change_reacs = _format_reactions(change_reacs)

            # This whole block adds the transfer of material from surface to bulk as surface grows (or vice versa)
            # it's not a reaction in the network so won't get picked up any other way. We manually add it.
            if species_name[0] == "@":
                if transfer >= 0:
                    change_reacs.append(
                        f"#{species_name[1:]} + SURFACE_TRANSFER -> {species_name}"
                    )
                else:
                    change_reacs.append(
                        f"{species_name} + SURFACE_TRANSFER -> #{species_name[1:]}"
                    )
                changes = np.append(changes, transfer)
            elif species_name[0] == "#":
                if transfer >= 0:
                    change_reacs.append(
                        f"@{species_name[1:]} + SURFACE_TRANSFER -> {species_name}"
                    )
                else:
                    change_reacs.append(
                        f"{species_name} + SURFACE_TRANSFER -> @{species_name[1:]}"
                    )
                changes = np.append(changes, transfer)

            # Then we remove the reactions that are not important enough to be printed by finding
            # which of the top reactions we need to reach rate_threshold*total_rate
            (
                total_formation,
                total_destruct,
                key_reactions,
                key_changes,
            ) = _remove_slow_reactions(
                changes, change_reacs, rate_threshold=rate_threshold
            )

            # only update if list of reactions change or rates change by factor of 10
            if (
                (old_key_reactions != key_reactions)
                or (
                    np.abs(
                        np.log10(np.abs(old_total_destruct))
                        - np.log10(np.abs(total_destruct))
                    )
                    > 1
                )
                or (np.abs(np.log10(old_total_form) - np.log10(total_formation)) > 1)
            ):
                old_key_reactions = key_reactions[:]
                old_total_form = total_formation
                old_total_destruct = total_destruct
                _write_analysis(
                    f,
                    row["Time"],
                    total_formation,
                    total_destruct,
                    key_reactions,
                    key_changes,
                )


def _param_dict_from_output(output_line):
    """
    Generate a parameter dictionary from a UCLCHEM timestep that contains enough of
    the physical variables to recreate the parameter dictionary used to run UCLCHEM.

    :param output_line: (pandas series) any row from the relevant UCLCHEM output
    """
    param_dict = {
        "initialdens": output_line["Density"],
        "initialtemp": output_line["gasTemp"],
        "zeta": output_line["zeta"],
        "radfield": output_line["radfield"],
        "baseav": 0.0,
        "rout": output_line["av"] * (1.6e21) / output_line["Density"],
    }
    return param_dict


def _get_species_rates(param_dict, input_abundances, species_index, reac_indxs):
    """
    Get the rate of up to 500 reactions from UCLCHEM for a given set of parameters and abundances.
    Intended for use within the analysis script.
    :param param_dict:  A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
    :param input_abundances: Abundance of every species in network
    :param reac_indxs: Index of reactions of interest in the network's reaction list.

    :returns: (ndarray) Array containing the rate of every reaction specified by reac_indxs
    """
    input_abund = np.zeros(500)
    input_abund[: len(input_abundances)] = input_abundances
    rate_indxs = np.ones(500)
    rate_indxs[: len(reac_indxs)] = reac_indxs
    rates, success_flag, transfer, swap, bulk_layers = wrap.get_rates(
        param_dict, input_abund, species_index, rate_indxs
    )
    if success_flag < 0:
        raise RuntimeError("UCLCHEM failed to return rates for these parameters")
    return rates[: len(reac_indxs)], transfer, swap, bulk_layers


def _get_rates_of_change(
    rates, reactions, speciesList, species, row, swap, bulk_layers
):
    """Calculate the terms in the rate of equation of a particular species using rates calculated using
    get_species_rates() and a row from the full output of UCLCHEM. See `analysis.py` for intended use.

    Args:
        rates (float, array): Rates of all reactions the species is involved in
        reactions (array): List of all reactions the species is involved in as a list of strings
        speciesList (array): List of species names from network
        species (string): name of species to be analyseds
        row (pd.Series): row from output dataframe
        swap (float): Total swap rate for individual swapping between bulk and surface
        bulk_layers (float): Number of layers in the bulk for individual swapping calc.

    Returns:
        _type_: _description_
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
            elif reactant in ["DESOH2", "FREEZE", "LH", "LHDES", "EXSOLID"]:
                reactant_count += 1
            if reactant in ["DEUVCR", "DESCR", "DESOH2", "SURFSWAP"]:
                change = change / np.max([1.0e-30, row["SURFACE"]])
            if reactant == "SURFSWAP":
                change = change * swap
            if reactant == "BULKSWAP":
                change = change * bulk_layers

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


def _remove_slow_reactions(changes, change_reacs, rate_threshold=0.99):
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


def _write_analysis(
    output_file, time, total_production, total_destruction, key_reactions, key_changes
):
    """Prints key reactions to file

    Args:
        time (float): Simulation time at which analysis is performed
        total_production (float): Total positive rate of change
        total_destruction (float): Total negative rate of change
        key_reactions (list): A list of all reactions that contribute to the total rate of change
        key_changes (list): A list of rates of change contributing to total
    """
    output_file.write(
        "\n\n***************************\nNew Important Reactions At: {0:.2e} years\n".format(
            time
        )
    )    
    # Formation and destruction writing is disabled since the absolute numbers do not appear to be correct.
    # output_file.write("Formation = {0:.2e} from:".format(total_production))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] > 0:
            outString = f"\n{reaction} : {float(key_changes[k] / total_production):.2%}"
            output_file.write(outString)

    # output_file.write("\n\nDestruction = {0:.2e} from:".format(total_destruction))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] < 0:
            outString = (
                f"\n{reaction} : {float(key_changes[k] / total_destruction):.2%}"
            )
            output_file.write(outString)


def _format_reactions(reactions):
    """Turn a row of the reaction file into a string.

    Args:
        reactions (list): list of lists, each reaction read from reaction.csv is a list.

    Returns:
        list: list of string, each reaction in readable string form
    """
    formatted_reactions = []
    for reaction in reactions:
        outString = "{x[0]} + {x[1]} + {x[2]} -> {x[3]} + {x[4]} + {x[5]}".format(
            x=reaction
        )
        outString = outString.replace(" + NAN", "")
        formatted_reactions.append(outString)
    return formatted_reactions


def _count_element(species_list, element):
    """
    Count the number of atoms of an element that appear in each of a list of species,
    return the array of counts

    :param  species_list: (iterable, str), list of species names
    :param element: (str), element

    :return: sums (ndarray) array where each element represents the number of atoms of the chemical element in the corresponding element of species_list
    """
    species_list = Series(species_list)
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
    """Calculates that the total elemental abundance of a species as a function of time. Allows you to check conservation.

    Args:
        element (str): Name of element
        df (pandas.DataFrame): DataFrame from `read_output_file()`

    Returns:
        pandas.Series: Total abundance of element for all time steps in df.
    """
    sums = _count_element(df.columns, element)
    for variable in ["Time", "Density", "gasTemp", "av", "point", "SURFACE", "BULK"]:
        sums = np.where(df.columns == variable, 0, sums)
    return df.mul(sums, axis=1).sum(axis=1)


def check_element_conservation(df, element_list=["H", "N", "C", "O"], percent=True):
    """Check the conservation of major element by comparing total abundance at start and end of model

    Args:
        df (pandas.DataFrame): UCLCHEM output in format from `read_output_file`
        element_list (list, optional): List of elements to check. Defaults to ["H", "N", "C", "O"].

    Returns:
        dict: Dictionary containing the change in the total abundance of each element as a fraction of initial value
    """
    result = {}
    for element in element_list:
        discrep = total_element_abundance(element, df).values
        if percent:
            discrep = np.abs(discrep[0] - discrep[-1]) / discrep[0]
            result[element] = f"{discrep:.3%}"
        else:
            discrep = discrep[0] - discrep[-1]
            result[element] = f"{discrep:.2e}"
    return result
