try:
    from uclchemwrap import uclchemwrap as wrap
except ImportError as E:
    E.add_note("Failed to import wrap.f90 from uclchemwrap, did the installation with f2py succeed?")
    raise
try:     
    from uclchemwrap import surfacereactions
except ImportError as E:
    E.add_note("Failed to import surfacereactions.f90 from uclchemwrap, did the installation with f2py succeed?")
    raise
import os
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import Series, read_csv
from seaborn import color_palette

from uclchem.constants import n_reactions, n_species
from uclchem.makerates import Reaction
from uclchem.makerates.network import Network
from uclchem.makerates.species import Species

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
    data = read_csv(f)
    data.columns = data.columns.str.strip()
    return data


def read_rate_file(rate_file):
    """Read the output of a UCLCHEM run created with the rateFile parameter into a pandas DataFrame

    Args:
        rate_file (str): path to file containing the UCLCHEM reaction rates.

    Returns:
        pandas.DataFrame: A dataframe containing the physical parameters, and reaction rates (s-1) at each timestep.
    """
    f = open(rate_file)
    data = read_csv(f)
    data.columns = data.columns.str.strip()
    return data

def _reactant_count(species: str, reaction_string: str) -> int:
    """Count how many times a species is consumed in a reaction

    Args:
        species (str): species which is maybe consumed in reaction
        reaction_string (str): reaction which maybe consumes species

    Returns:
        int: amount of times species is consumed in reaction
    """
    split_str = reaction_string.split("->")[0].strip()
    if " " not in split_str:
        return species == split_str
    return split_str.split().count(species)


def _product_count(species: str, reaction_string: str) -> int:
    """Count how many times a species is produced in a reaction

    Args:
        species (str): species which is maybe produced by reaction
        reaction_string (str): reaction which maybe produces species

    Returns:
        int: amount of times species is produced by reaction
    """
    split_str = reaction_string.split("->")[1].strip()
    if " " not in split_str:
        return species == split_str
    return split_str.split().count(species)


def _get_rates_change(rate_df: pd.DataFrame, species: str) -> pd.DataFrame:
    phys_param_columns = []
    for i, column in enumerate(rate_df.columns):
        if not "->" in column:  # It is not a reaction, but a physical parameter
            phys_param_columns.append(i)
        if "->" in column:  # Assume reactions come after all the physical parameters
            break
    change_df = rate_df.iloc[:, phys_param_columns]
    for i, column in enumerate(rate_df.columns):
        if not "->" in column:  # It is not a reaction, but a physical parameter
            continue
        rcount = _reactant_count(species, column)
        pcount = _product_count(species, column)
        if rcount == 0 and pcount == 0:
            # Species is unaffected by reaction
            continue
        change_df = pd.concat([change_df, rate_df[column] * (pcount - rcount)], axis=1)
    return change_df


def get_change_df(
    rate_df: pd.DataFrame, species: str, on_grain: bool = False
) -> pd.DataFrame:
    """From a dataframe containing all the reaction rates, get the change of a species over time, due to each reaction.

    Args:
        rate_df (pd.DataFrame): dataframe containing physical parameters and reaction rates over time
        species (str): species to get the change over time
        on_grain (bool): whether to analyse the ice phase of this species

    Returns:
        change_df (pd.DataFrame): change of species over time due to each reaction the species is involved in
    """
    if "#" in species or "@" in species:
        msg = "WARNING: get_change_df IS ONLY FOR ANALYSING ALL OF THE GAS PHASE AND ALL OF THE ICE. "
        msg += "USE on_grain PARAMETER TO INDICATE THIS. IF YOU WANT TO ANALYSE ONLY SURFACE OR ONLY BULK, "
        msg += "USE FUNCTION _get_rates_change WITH SPECIES CONTAINING # OR @ TO INDICATE SURFACE OF BULK."
        raise ValueError(msg)
    if not on_grain:
        return _get_rates_change(rate_df, species)
    df_surf = _get_rates_change(rate_df, "#" + species)
    df_bulk = _get_rates_change(rate_df, "@" + species)
    surf_columns = df_surf.columns
    bulk_columns = df_bulk.columns
    for column in surf_columns:
        if not "->" in column:
            df_bulk.drop(
                columns=column, inplace=True
            )  # Drop the physical parameters from bulk column so we do not have them twice in the final df
            continue
        if column in bulk_columns:
            # If the same reaction is in both dfs, that means that both surf and bulk version of the species is involved in the reaction
            # which means that either surface is lost and bulk forms, or the other way, and so we drop this column
            # from both dfs since it has no effect on the ice overall.
            df_surf.drop(columns=column, inplace=True)
            df_bulk.drop(columns=column, inplace=True)
    # Maybe TODO:
    # Make it such that the columns of the same reactions (but surf and bulk versions) are added
    # such that we have a single reaction rate in the ice, and not seperate surf and bulk reaction rates.
    return pd.concat([df_surf, df_bulk], axis=1)


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


def plot_species(ax, df, species, legend=True, **plot_kwargs):
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
        plot_kwargs["linestyle"] = linestyle
        plot_kwargs["label"] = specName
        # Support legacy code that use either "age" or "Time" as the time variable
        if "age" in df.columns:
            timecolumn = "age"
        elif "Time" in df.columns:
            timecolumn = "Time"
        else:
            raise ValueError("No time variable in dataframe")
        ax.plot(
            df[timecolumn],
            abundances,
            lw=2,
            **plot_kwargs,
        )
        ax.set(yscale="log")
        if legend:
            ax.legend()
    return ax


def read_analysis(filepath, species):
    with open(filepath, "r") as file:
        lines = file.readlines()
    for i, line in enumerate(lines):
        if "All Reactions" in line:
            first_newline = lines.index("\n", i)
            all_reactions = lines[i + 2 : first_newline]
            all_reactions = [reaction.strip("\n") for reaction in all_reactions]
            break

    counts = [all_reactions.count(reac) for reac in all_reactions]
    for i, count in enumerate(counts):
        if count > 1:
            print(
                f"{all_reactions[i]} was found to be in the reaction list {count} times."
            )
            print(f"Renaming this time to '{all_reactions[i]} ({i})'")
            all_reactions[i] = f"{all_reactions[i]} ({i})"

    n_cols = len(all_reactions) + 1
    columns = ["Time"]
    columns.extend(all_reactions)
    df = pd.DataFrame(columns=columns)

    segments_min = [
        i for i, line in enumerate(lines) if "New Important Reactions At:" in line
    ]
    segments_max = [i - 1 for i in segments_min[1:]]
    segments_max.append(len(lines) - 1)
    segments = [
        lines[segments_min[i] : segments_max[i]] for i in range(len(segments_min))
    ]

    for segment_lines in segments:
        new_row = [0] * n_cols
        time = float(segment_lines[0].split()[-2])
        new_row[0] = time
        for i, reaction in enumerate(all_reactions):
            if not any(reaction in line for line in segment_lines):
                rate = 0
            else:
                reac_indx = [
                    j for j, line in enumerate(segment_lines) if reaction in line
                ][0]
                rate = float(segment_lines[reac_indx].split()[-3])
            # Here we use the fact that the sign of the formation/destruction rate is already
            # correct, so we only need to scale it by 1 (in e.g. H+OH->H2O) or 2 (H+H->H2)
            new_row[i + 1] = rate * (
                _product_count(species, reaction) + _reactant_count(species, reaction)
            )

        new_row_dict = dict(zip(columns, new_row))
        new_row_df = pd.DataFrame(new_row_dict, index=[0])
        df = pd.concat(
            [df, new_row_df],
            axis=0,
            ignore_index=True,
        )
    return df, all_reactions


def analysis(species_name, result_file, output_file, rate_threshold=0.99):
    """A function which loops over every time step in an output file and finds the rate of change of a species at that time due to each of the reactions it is involved in.
    From this, the most important reactions are identified and printed to file. This can be used to understand the chemical reason behind a species' behaviour.

    Args:
        species_name (str): Name of species to be analysed
        result_file (str): The path to the file containing the UCLCHEM output
        output_file (str): The path to the file where the analysis output will be written
        rate_threshold (float,optional): Analysis output will contain the only the most efficient reactions that are responsible for rate_threshold of the total production and destruction rate. Defaults to 0.99.
    """
    print(
        "WARNING: THIS VERSION OF THE ANALYSIS TOOL IS CURRENTLY INCORRECT FOR SOLID (SURFACE/BULK) SPECIES,"
    )
    print(
        "ONLY USE IT FOR GAS-PHASE. IF YOU WANT INFORMATION ABOUT FORMATION/DESTRUCTION RATES OF ICE SPECIES, RUN THE MODEL WITH THE 'rateOutput' OPTION SET TO A FILEPATH"
    )
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

    if species_name[0] == "#":
        surftransfer_reacs = [
            f"@{species_name[1:]} + SURFACETRANSFER -> {species_name}",
            f"{species_name} + SURFACETRANSFER -> @{species_name[1:]}",
        ]
        formatted_reacs.extend(surftransfer_reacs)
    if species_name[0] == "@":
        surftransfer_reacs = [
            f"{species_name} + SURFACETRANSFER -> #{species_name[1:]}",
            f"#{species_name[1:]} + SURFACETRANSFER -> {species_name}",
        ]
        formatted_reacs.extend(surftransfer_reacs)

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
            if transfer <= 0:
                if species_name[0] == "#":
                    change_reacs.append(surftransfer_reacs[0])
                    changes = np.append(changes, transfer)
                elif species_name[0] == "@":
                    change_reacs.append(surftransfer_reacs[1])
                    changes = np.append(changes, transfer)
            else:
                if species_name[0] == "#":
                    change_reacs.append(surftransfer_reacs[1])
                    changes = np.append(changes, -transfer)
                elif species_name[0] == "@":
                    change_reacs.append(surftransfer_reacs[0])
                    changes = np.append(changes, transfer)

            # Then we remove the reactions that are not important enough to be printed by finding
            # which of the top reactions we need to reach rate_threshold*total_rate
            # (
            #     total_formation,
            #     total_destruct,
            #     key_reactions,
            #     key_changes,
            # ) = _remove_slow_reactions(
            #     changes, change_reacs, rate_threshold=rate_threshold
            # )

            # only update if list of reactions change or rates change by factor of 10
            # if (
            #     (old_key_reactions != key_reactions)
            #     or (
            #         np.abs(
            #             np.log10(np.abs(old_total_destruct))
            #             - np.log10(np.abs(total_destruct))
            #         )
            #         > 0
            #     )
            #     or (np.abs(np.log10(old_total_form) - np.log10(total_formation)) > 0)
            # ):
            old_key_reactions = key_reactions[:]
            old_total_form = total_formation
            old_total_destruct = total_destruct
            _write_analysis(
                f,
                row["Time"],
                total_formation,
                total_destruct,
                change_reacs,
                changes,
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
        "baseav": output_line["baseAv"],
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
    input_abund = np.zeros(n_species)
    input_abund[: len(input_abundances)] = input_abundances
    rate_indxs = np.ones(n_reactions)
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
    safeMantle = np.max([1.0e-30, row["SURFACE"]])
    for i, reaction in enumerate(reactions):
        reaction_instance = Reaction([*reaction, 0, 0, 0, 0, 0])
        change = rates[i]
        reactants = reaction[0:3]
        products = reaction[3:]

        # Counting the same as Reaction.body_count
        reactant_count = reaction_instance.body_count

        change = change * (row["Density"] ** (reactant_count))
        for reactant in reactants:
            if reactant in speciesList:
                change = change * row[reactant]

            elif reactant == "BULKSWAP":
                change = change * bulk_layers
            elif reactant == "SURFSWAP":
                change = change * swap / safeMantle
            elif reactant in ["DEUVCR", "DESCR", "DESOH2", "ER", "ERDES"]:
                change = change / safeMantle
                if reactant == "DESOH2":
                    change = change * row["H"]
            elif (not three_phase) and (reactant in ["THERM"]):
                change = change * row["Density"] / safeMantle

            if "H2FORM" in reactants:
                # only 1 factor of H abundance in Cazaux & Tielens 2004 H2 formation so stop looping after first iteration
                break

        if "LH" in reactants[2]:
            if "@" in reactants[0]:
                change = change * bulk_layers

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
    output_file.write("Formation = {0:.8e} from:".format(total_production))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] > 0:
            outString = f"\n{reaction} : {float(key_changes[k])} = {float(key_changes[k] / total_production):.2%}"
            output_file.write(outString)

    output_file.write("\n\nDestruction = {0:.8e} from:".format(total_destruction))
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] < 0:
            outString = f"\n{reaction} : {float(key_changes[k])} = {float(key_changes[k] / total_destruction):.2%}"
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


def get_total_swap(rates: pd.DataFrame, abundances: pd.DataFrame, reactions: List[Reaction]) -> np.ndarray:
    """ Obtain the amount of 'random' swapping per timestep

    Args:
        rates (pd.DataFrame): The rates obtained from running an UCLCHEM model
        abundances (pd.DataFrame): The abundances obtained from running an UCLCHEM model
        reactions (List[Reaction]): The reactions used in UCLCHEM

    Returns:
        np.ndarray: The total swap per timestep
    """
    assert len(rates) == len(abundances), "Rates and abundances must be the same length"
    assert rates.shape[1] == len(reactions), "The number of rates and reactions must be equal"
    totalSwap = np.zeros(abundances.shape[0])
    for idx, reac in enumerate(reactions):
        if reac.get_reaction_type() == "BULKSWAP":
            totalSwap += rates.iloc[:, idx] * abundances[reac.get_pure_reactants()[0]]
    return totalSwap


def construct_incidence(species: List[Species], reactions: List[Reaction]) -> np.ndarray:
    """ Construct the incidence matrix, a matrix that describes the in and out degree
    for each of the reactions; useful to matrix multiply by the indvidual fluxes per reaction
    to obtain a flux (dy) per species.

    Args:
        species (List[Species]): A list of species S
        reactions (List[Reaction]): The list of reactions S

    Returns:
        np.ndarray: A RxS incidence matrix
    """
    incidence = np.zeros(
        dtype=np.int8,
        shape=(len(reactions), len(species)),
    )
    for idx, reaction in enumerate(reactions):
        products = reaction.get_pure_products()
        for prod in products:
            incidence[idx, species.index(prod)] += 1
        reactants = reaction.get_pure_reactants()
        for reac in reactants:
            incidence[idx, species.index(reac)] -= 1
    return incidence


def postprocess_rates_to_dy(
    physics: pd.DataFrame,
    abundances: pd.DataFrame,
    rates: pd.DataFrame,
    network: Network = None,
    species: list[Species] = None,
    reactions: list[Reaction] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply postprocessing to obtain the equivalent of GETYDOT from the fortran
    side; It returns two dataframes:
        - ydot: the RHS that is solved in UCLCHEM at every output timestep
        - flux_by_reaction: the individual terms that result in ydot when multiplied by the incidence matrix

    Args:
        physics (pd.DataFrame): The physics output from running a model
        abundances (pd.DataFrame): The abundances output from running a model
        rates (pd.DataFrame): The rates output from running a model
        network (Network, optional): The reaction network used to postprocess the rates. Defaults to None.
        species (list[Species], optional): The species used to postprocess the rates . Defaults to None.
        reactions (list[Reaction], optional): The reactions used to postprocess the rates. Defaults to None.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: dy, flux_by_reaction.
    """
    assert bool(species) == bool(reactions), (
        "If species is specified, reactions also must be and vice ver"
    )
    assert not (network and (species or reactions)), (
        "Choose between providing a network OR (species AND reactions)"
    )
    if network:
        species = network.get_species_list()
        reactions = network.get_reaction_list()
    # Import all of the constants directly from UCLCHEMWRAP to avoid discrepancies
    GAS_DUST_DENSITY_RATIO = surfacereactions.gas_dust_density_ratio
    NUM_SITES_PER_GRAIN = surfacereactions.num_sites_per_grain
    # Compute dynamic quantities that can be precomputed
    bulkLayersReciprocal = (
        NUM_SITES_PER_GRAIN / (GAS_DUST_DENSITY_RATIO * abundances["BULK"])
    ).apply(lambda x: min(1.0, x))
    totalSwap = bulkLayersReciprocal * get_total_swap(rates, abundances, reactions)

    # Create the incidence matrix, we use this to evaluate the rates and
    incidence = construct_incidence(species, reactions)

    # Create a copy so we don't overwrite
    flux_by_reaction = rates.copy()

    # Iterate over each of the rates and compute the contribution to dy for that reaction.
    for idx, reaction in enumerate(network.get_reaction_list()):
        density_multiplier_factor = reaction.body_count
        rate = flux_by_reaction.iloc[:, idx]
        # Multiply by density the right number of times:
        for iiii in range(int(density_multiplier_factor)):
            rate *= physics["Density"]
        # Multiply by the abundances:
        for reactant in reaction.get_sorted_reactants():
            if reactant in list(abundances.columns):
                rate *= abundances[reactant]
        match reaction.get_reaction_type():
            case x if x in ["LH", "LHDES", "BULKSWAP"]:
                if reaction.is_bulk_reaction(include_products=False):
                    rate *= bulkLayersReciprocal
            case "SURFSWAP":
                rate *= totalSwap / abundances["SURFACE"]
            case x if x in ["DESCR", "DEUVCR", "ER", "ERDES"]:
                rate /= abundances["SURFACE"]
            case "DESOH2":
                rate *= abundances["H"] / abundances["SURFACE"]
            case "H2FORM":
                # For some reason, H2form only uses the hydrogen density once
                rate /= abundances["H"]
        flux_by_reaction.iloc[:, idx] = rate

    # Compute the flux at each timestep, adding the appropriate header
    dy = flux_by_reaction @ incidence
    dy.columns = [str(s) for s in species]
    # Compute the SURFACE and BULK:
    dy.loc[:, "SURFACE"] = dy.loc[:, dy.columns.str.startswith("#")].sum(axis=1)
    dy.loc[:, "BULK"] = dy.loc[:, dy.columns.str.startswith("@")].sum(axis=1)

    # Compute the corrections to account for the transport between SURFACE and BULK
    swap_flux_correction = pd.DataFrame()
    # Then correct for the addition or subtraction of material to/from the bulk:
    surfswap_reactions = [r for r in reactions if r.get_reaction_type() == "SURFSWAP"]
    bulkswap_reactions = [r for r in reactions if r.get_reaction_type() == "BULKSWAP"]
    # Sort them in the correct order, s.t. it matches the saving to disk format.
    surfswap_reactions = sorted(
        surfswap_reactions,
        key=lambda x: species.index(x.get_reactants()[0]),
    )
    bulkswap_reactions = sorted(
        bulkswap_reactions,
        key=lambda x: species.index(x.get_reactants()[0]),
    )
    for (idx_j, rate_row), (idx_i, abunds_row) in zip(
        rates.iterrows(), abundances.iterrows()
    ):
        # Walk through the bulkswap and reactionswap pathways:
        _sswap_rates = {}
        _bswap_rates = {}
        # TODO: vectorize this, because this is slower than it has to be.
        for r_bswap, r_sswap in zip(bulkswap_reactions, surfswap_reactions):
            surfaceCoverage = min(1.0, abunds_row["BULK"] / abunds_row["SURFACE"])
            if dy.iloc[idx_j]["SURFACE"] < 0.0:
                surfaceCoverage = min(1.0, abunds_row["BULK"] / abunds_row["SURFACE"])
                # SURFACE is shrinking, so bulk must be growing
                bswap = (
                    dy.iloc[idx_j]["SURFACE"]
                    * surfaceCoverage
                    * abunds_row[r_bswap.get_reactants()[0]]
                    / abunds_row["BULK"]
                )
                _bswap_rates[str(r_bswap).replace("SWAP", "SWAP_TRANSPORT")] = bswap
                _sswap_rates[str(r_sswap).replace("SWAP", "SWAP_TRANSPORT")] = 0.0
                # Immedidiately correct dy:
                dy.loc[idx_j, r_bswap.get_reactants()[0]] -= bswap
                dy.loc[idx_j, r_bswap.get_products()[0]] += bswap
            else:
                surfaceCoverage = 0.5 * GAS_DUST_DENSITY_RATIO / NUM_SITES_PER_GRAIN
                sswap = (
                    dy.iloc[idx_j]["SURFACE"]
                    * surfaceCoverage
                    * abunds_row[r_sswap.get_reactants()[0]]
                )
                _bswap_rates[str(r_bswap).replace("SWAP", "SWAP_TRANSPORT")] = 0.0
                _sswap_rates[str(r_sswap).replace("SWAP", "SWAP_TRANSPORT")] = sswap
                # Immedidiately correct dy:
                dy.loc[idx_j, r_sswap.get_products()[0]] -= sswap
                dy.loc[idx_j, r_sswap.get_reactants()[0]] += sswap
        swap_flux_correction = pd.concat(
            (
                swap_flux_correction,
                (
                    pd.DataFrame.from_dict(
                        _bswap_rates | _sswap_rates,
                        orient="index",
                    ).T
                ),
            )
        )
    swap_flux_correction = swap_flux_correction.reset_index(drop=True)
    flux_by_reaction = pd.concat((flux_by_reaction, swap_flux_correction), axis=1)
    # Correct the change in surface and bulk by summing the constituents:
    dy.loc[:, "SURFACE"] = dy.loc[:, dy.columns.str.startswith("#")].sum(axis=1)
    dy.loc[:, "BULK"] = dy.loc[:, dy.columns.str.startswith("@")].sum(axis=1)
    # Apply the new fluxes to the ydots and compute a new ydot for surface and bulk:
    return (
        dy,
        flux_by_reaction,
    )