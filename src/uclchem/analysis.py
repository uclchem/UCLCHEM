"""UCLCHEM Analysis Module.

Tools for analyzing chemical model outputs and reaction pathways.

This module provides functions to:
- Read and parse UCLCHEM output files
- Analyze chemical reaction pathways for specific species
- Check element conservation in model results
- Create abundance plots and visualizations
- Compare model results across different runs

**Key Functions:**

- :func:`read_output_file` - Read UCLCHEM output files into DataFrames
- :func:`analysis` - Analyze production/destruction pathways for a species
- :func:`check_element_conservation` - Verify element conservation

**Example Usage:**

    >>> import uclchem
    >>>
    >>> model = uclchem.model.Cloud({})
    >>> model.check_error()
    Model ran successfully
    >>>
    >>> physics_df, chemistry_df, rate_constants_df = model.get_dataframes(
    ...     with_rate_constants=True,
    ... )
    >>>
    >>> uclchem.analysis.check_element_conservation(chemistry_df)
    {'H': ..., 'N': ..., 'C': ..., 'O': ...}
    >>>
    >>> # Making a plot of the main formation and destruction reactions
    >>> # at a specific timepoint
    >>> network = uclchem.makerates.network.Network.from_csv()
    >>> dy, reaction_rates = uclchem.analysis.rate_constants_to_dy_and_rates(
    ...     physics_df,
    ...     chemistry_df,
    ...     rate_constants_df,
    ...     network=network,
    ... )
    >>> production_df, destruction_df = uclchem.analysis.get_production_and_destruction(
    ...     "H2O",
    ...     reaction_rates,
    ... )
    >>>

**See Also:**

- :mod:`uclchem.plot` - Dedicated plotting utilities
- :mod:`uclchem.model` - Run chemical models
"""

try:
    from uclchemwrap import uclchemwrap as wrap
except ImportError as E:
    E.add_note(
        "Failed to import wrap.f90 from uclchemwrap, did the installation with f2py succeed?"
    )
    raise
try:
    from uclchemwrap import surfacereactions
except ImportError as E:
    E.add_note(
        "Failed to import surfacereactions.f90 from uclchemwrap, did the installation with f2py succeed?"
    )
    raise
import warnings
from pathlib import Path
from typing import Any, TextIO

import numpy as np
import pandas as pd

from uclchem.constants import default_elements_to_check, n_reactions, n_species
from uclchem.makerates import Reaction
from uclchem.makerates.network import Network
from uclchem.makerates.species import Species, element_list
from uclchem.utils import UCLCHEM_ROOT_DIR, ArrayLike


def read_output_file(output_file: str | Path) -> pd.DataFrame:
    """Read the output of a UCLCHEM run created with the outputFile parameter
    into a pandas DataFrame.

    Args:
        output_file (str | Path): path to file containing a full UCLCHEM output

    Returns:
        data (pd.DataFrame): A dataframe containing the abundances and
            physical parameters of the model at every time step.

    """
    with Path(output_file).open() as f:
        data = pd.read_csv(f)
    data.columns = data.columns.str.strip()
    return data


def read_rate_file(rate_file: str | Path) -> pd.DataFrame:
    """Read the output of a UCLCHEM run created with the rateConstantFile
    parameter into a pandas DataFrame.

    Args:
        rate_file (str | Path): path to file containing the UCLCHEM reaction rates.

    Returns:
        data (pd.DataFrame): A dataframe containing the physical parameters,
            and reaction rates (s-1) at each timestep.

    """
    with Path(rate_file).open() as f:
        data = pd.read_csv(f)
    data.columns = data.columns.str.strip()
    return data


def _reactant_count(species: str, reaction_string: str) -> int:
    """Count how many times a species is consumed in a reaction.

    Args:
        species (str): species which is maybe consumed in reaction
        reaction_string (str): reaction which maybe consumes species

    Returns:
        int: number of times species is consumed in reaction

    """
    # TODO: Put in Reaction
    split_str = reaction_string.split("->")[0].strip()
    if " " not in split_str:
        return species == split_str
    return split_str.split().count(species)


def _product_count(species: str, reaction_string: str) -> int:
    """Count how many times a species is produced in a reaction.

    Args:
        species (str): species which is maybe produced by reaction
        reaction_string (str): reaction which maybe produces species

    Returns:
        int: amount of times species is produced by reaction

    """
    # TODO: Put in Reaction
    split_str = reaction_string.split("->")[1].strip()
    if " " not in split_str:
        return species == split_str
    return split_str.split().count(species)


def _get_rates_change(rate_df: pd.DataFrame, species: str) -> pd.DataFrame:
    phys_param_columns = []
    for i, column in enumerate(rate_df.columns):
        if "->" not in column:  # It is not a reaction, but a physical parameter
            phys_param_columns.append(i)
        if "->" in column:  # Assume reactions come after all the physical parameters
            break
    change_df = rate_df.iloc[:, phys_param_columns]
    for i, column in enumerate(rate_df.columns):
        if "->" not in column:  # It is not a reaction, but a physical parameter
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
    """From a dataframe containing all the reaction rates, get the change of a species over time,
    due to each reaction.

    Args:
        rate_df (pd.DataFrame): dataframe containing physical parameters and
            reaction rate constants over time
        species (str): species to get the change over time
        on_grain (bool): whether to analyse the ice phase of this species

    Returns:
        change_df (pd.DataFrame): change of species over time due to each reaction
            the species is involved in

    Raises:
        DeprecationWarning: Deprecated in UCLCHEM 4.0
        ValueError: If "#" or "@" is in ``species``.

    """
    msg = "This function will be deprecated in UCLCHEM 4.0 and is no longer actively maintained"
    raise DeprecationWarning(msg)

    if "#" in species or "@" in species:
        msg = "WARNING: get_change_df IS ONLY FOR ANALYZING ALL OF THE GAS PHASE AND ALL OF THE ICE. "
        msg += "USE on_grain PARAMETER TO INDICATE THIS. IF YOU WANT TO ANALYZE ONLY SURFACE OR ONLY BULK, "
        msg += "USE FUNCTION _get_rates_change WITH SPECIES CONTAINING # OR @ TO INDICATE SURFACE OF BULK."
        raise ValueError(msg)
    if not on_grain:
        return _get_rates_change(rate_df, species)
    df_surf = _get_rates_change(rate_df, "#" + species)
    df_bulk = _get_rates_change(rate_df, "@" + species)
    surf_columns = df_surf.columns
    bulk_columns = df_bulk.columns
    for column in surf_columns:
        if "->" not in column:
            df_bulk.drop(
                columns=column, inplace=True
            )  # Drop the physical parameters from bulk column so we do not have them twice in the final df
            continue
        if column in bulk_columns:
            # If the same reaction is in both dfs, that means that both surf and bulk version
            # of the species is involved in the reaction, which means that either surface is lost
            # and bulk forms, or the other way, and so we drop this column
            # from both dfs since it has no effect on the ice overall.
            df_surf.drop(columns=column, inplace=True)
            df_bulk.drop(columns=column, inplace=True)
    # Maybe TODO:
    # Make it such that the columns of the same reactions (but surf and bulk versions)
    # are added such that we have a single reaction rate in the ice,
    # and not separate surf and bulk reaction rates.
    return pd.concat([df_surf, df_bulk], axis=1)


def read_analysis(filepath: str | Path, species: str) -> tuple[pd.DataFrame, list[str]]:
    """Read the analysis output.

    Args:
        filepath (str | Path): path to analysis output.
        species (str): Species of interest.

    Returns:
        df (pd.DataFrame): dataframe with rates and time.
        all_reactions (list[str]): list of all reactions that the species is involved in.

    Raises:
        DeprecationWarning: Deprecated in UCLCHEM 4.0

    """
    msg = "This function will be deprecated in UCLCHEM 4.0 and is no longer actively maintained"
    raise DeprecationWarning(msg)
    with Path(filepath).open() as file:
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


def analysis(
    species_name: str,
    output_file: str | Path,
    analysis_file: str | Path,
    rate_threshold: float = 0.99,
) -> None:
    """Loop over every time step in an output file and finds the rate of change
    of a species at that time due to each of the reactions it is involved in.
    From this, the most important reactions are identified and printed to file.
    This can be used to understand the chemical reason behind a species' behavior.

    DEPRECATED

    Args:
        species_name (str): Name of species to be analysed
        output_file (str | Path): The path to the file where the analysis output will be written
        analysis_file (str): The path to the file containing the UCLCHEM output
        rate_threshold (float): Analysis output will contain the only the most efficient
            reactions that are responsible for rate_threshold of the total
            production and destruction rate. Default = 0.99.

    """
    result_df = read_output_file(output_file)
    species_array = np.loadtxt(
        UCLCHEM_ROOT_DIR / "species.csv",
        usecols=[0],
        dtype=str,
        skiprows=1,
        unpack=True,
        delimiter=",",
        comments="%",
    )
    species = list(species_array)
    reactions = np.loadtxt(
        UCLCHEM_ROOT_DIR / "reactions.csv",
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
    old_key_reactions: list[str] = []
    old_total_destruct = 0.0
    old_total_form = 0.0
    formatted_reacs = _format_reactions(list(reactions[reac_indxs]))

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

    with Path(analysis_file).open("w") as f:
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

            # convert reaction rates to total rates of change.
            # this needs manually updating when you add new reaction types!
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

            # This whole block adds transfer of material from surface to bulk as surface grows
            # (or vice versa). It's not a reaction in the network so won't get picked up otherwise.
            # We manually add it.
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
                    > 0
                )
                or (np.abs(np.log10(old_total_form) - np.log10(total_formation)) > 0)
            ):
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


def _param_dict_from_output(
    output_line: dict[str, float] | pd.Series,
) -> dict[str, float]:
    """Generate a parameter dictionary from a UCLCHEM timestep that contains enough of
    the physical variables to recreate the parameter dictionary used to run UCLCHEM.

    Args:
        output_line (dict[str, float] | pd.Series): any row from the relevant UCLCHEM output

    Returns:
        dict[str, float]: dictionary with physical conditions from output.

    """
    param_dict = {
        "initialdens": output_line["Density"],
        "initialtemp": output_line["gasTemp"],
        "zeta": output_line["zeta"],
        "radfield": output_line["radfield"],
        "baseav": output_line["baseAv"],
    }
    return param_dict


def _get_species_rates(
    param_dict: dict[str, Any],
    input_abundances: ArrayLike,
    species_index: int,
    reac_indxs: list[int],
) -> tuple[np.ndarray, float, float, float]:
    """Get the rate of up to 500 reactions from UCLCHEM for a given
    set of parameters and abundances.
    Intended for use within the analysis script.

    Args:
        param_dict (dict[str, Any]): A dictionary of parameters where keys are
            any of the variables in defaultparameters.f90 and values are value for current run.
        input_abundances (ArrayLike): Abundance of every species in network
        species_index (int): Index of the species of interest.
        reac_indxs (list[int]): Indices of reactions of interest in the network's reaction list.

    Returns:
        rates (np.ndarray): Array containing the rate of every reaction specified by reac_indxs
        transfer (float):
        swap (float):
        bulk_layers (float): number of monolayers of bulk ice

    Raises:
        DeprecationWarning: Deprecated in UCLCHEM 4.0
        RuntimeError: If UCLCHEM failed to return the rates for these parameters

    """
    msg = "This function will be deprecated in UCLCHEM 4.0 and is no longer actively maintained"
    raise DeprecationWarning(msg)
    input_abund = np.zeros(n_species)
    input_abund[: len(input_abundances)] = input_abundances
    rate_indxs = np.ones(n_reactions)
    rate_indxs[: len(reac_indxs)] = reac_indxs
    rates, success_flag, transfer, swap, bulk_layers = wrap.get_rates(
        param_dict, input_abund, species_index, rate_indxs
    )
    if success_flag < 0:
        msg = "UCLCHEM failed to return rates for these parameters"
        raise RuntimeError(msg)
    return rates[: len(reac_indxs)], transfer, swap, bulk_layers


def _get_rates_of_change(
    rates: np.ndarray,
    reactions: list[list] | np.ndarray,
    species_list: list[str],
    species: str,
    row: pd.Series,
    swap: float,
    bulk_layers: float,
):
    """Calculate the terms in the rate of equation of a particular species using rates
    calculated using ``get_species_rates()`` and a row from the full output of UCLCHEM.
    See ``analysis.py`` for intended use.

    Args:
        rates (np.ndarray): Rates of all reactions the species is involved in
        reactions (list[list] | np.ndarray): List of all reactions the species is involved in.
        species_list (array): List of species names from network
        species (string): name of species to be analyseds
        row (pd.Series): row from output dataframe
        swap (float): Total swap rate for individual swapping between bulk and surface
        bulk_layers (float): Number of layers in the bulk for individual swapping calc.

    Returns:
        _type_: _description_

    Raises:
        DeprecationWarning: Deprecated in UCLCHEM 4.0

    """
    msg = "This function will be deprecated in UCLCHEM 4.0 and is no longer actively maintained"
    raise DeprecationWarning(msg)
    changes = []
    reactionList = []
    three_phase = "@" in "".join(species_list)
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
            if reactant in species_list:
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
                # only 1 factor of H abundance in Cazaux & Tielens 2004 H2 formation
                # so stop looping after first iteration
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


def _remove_slow_reactions(
    changes: np.ndarray, change_reacs: list[str], rate_threshold: float = 0.99
) -> tuple[float, float, list[str], list[float]]:
    """Iterate through a list of reactions adding the fastest reactions to a list until some
    threshold fraction of the total rate of change is reached. This list is returned so that
    you have the list of reactions that cause rate_threshold of the total destruction and
    formation of a species.

    Args:
        changes (np.ndarray): List of rates of change due to each reaction
            a species is involved in.
        change_reacs (list[str]): List of corresponding rates of change
        rate_threshold (float): Percentage of overall rate of change to consider before ignoring
            less important reactions. Default = 0.99.

    Returns:
        totalProd (float): Total production rate
        totalDestruct (float): Total destruction rate
        key_reactions (list[str]): List of key reactions
        key_changes (list[float]): List of reaction rates of key reactions

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
    output_file: TextIO,
    time: float,
    total_production: float,
    total_destruction: float,
    key_reactions: list[str],
    key_changes: list[float],
) -> None:
    """Print key reactions to a file.

    Args:
        output_file (TextIO): Open file object to write to
        time (float): Simulation time at which analysis is performed
        total_production (float): Total positive rate of change
        total_destruction (float): Total negative rate of change
        key_reactions (list[str]): A list of all reactions that contribute
            to the total rate of change
        key_changes (list[float]): A list of rates of change contributing to total

    """
    output_file.write(
        f"\n\n***************************\nNew Important Reactions At: {time:.2e} years\n"
    )
    # Formation and destruction writing is disabled since the absolute numbers
    # do not appear to be correct.
    output_file.write(f"Formation = {total_production:.8e} from:")
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] > 0:
            outString = f"\n{reaction} : {float(key_changes[k])} = {float(key_changes[k] / total_production):.2%}"
            output_file.write(outString)

    output_file.write(f"\n\nDestruction = {total_destruction:.8e} from:")
    for k, reaction in enumerate(key_reactions):
        if key_changes[k] < 0:
            outString = f"\n{reaction} : {float(key_changes[k])} = {float(key_changes[k] / total_destruction):.2%}"
            output_file.write(outString)


def _format_reactions(reactions: list[list[str]]) -> list[str]:
    """Turn a row of the reaction file into a string.

    Args:
        reactions (list[list[str]]): list of lists, each reaction read from reaction.csv is a list.

    Returns:
        list[str]: list of string, each reaction in readable string form

    """
    # TODO: Replace with str(Reaction()).
    formatted_reactions = []
    for reaction in reactions:
        outString = f"{reaction[0]} + {reaction[1]} + {reaction[2]} -> {reaction[3]} + {reaction[4]} + {reaction[5]}"
        outString = outString.replace(" + NAN", "")
        formatted_reactions.append(outString)
    return formatted_reactions


def _count_element(species_list: list[str], element: str) -> pd.Series:
    """Count the number of atoms of an element that appear in each species of a list of species.

    Args:
        species_list (list[str]): list of species names
        element (str): element

    Returns:
        sums (pd.Series): array where each element represents the number of atoms
            of the chemical element in the corresponding element of species_list

    """
    species = pd.Series(species_list)
    # confuse list contains elements whose symbols contain the target eg CL for C
    # We count both sets of species and remove the confuse list counts.
    confuse_list = [x for x in element_list if element in x]
    confuse_list = sorted(confuse_list, key=lambda x: len(x), reverse=True)
    confuse_list.remove(element)
    sums = species.str.count(element)
    for i in range(2, 10):
        sums += np.where(species.str.contains(element + f"{i:.0f}"), i - 1, 0)
    for spec in confuse_list:
        sums += np.where(species.str.contains(spec), -1, 0)
    return sums


def total_element_abundance(element: str, df: pd.DataFrame) -> pd.Series:
    """Calculate the total elemental abundance of a species as a function of time.
    Allows you to check conservation.

    Args:
        element (str): Name of element
        df (pd.DataFrame): DataFrame from ``read_output_file()``

    Returns:
        pd.Series: Total abundance of element for all time steps in df.

    """
    sums_array = _count_element(list(df.columns), element).to_numpy()
    for variable in ["Time", "Density", "gasTemp", "av", "point", "SURFACE", "BULK"]:
        sums_array = np.where(df.columns == variable, 0, sums_array)
    return df.mul(sums_array, axis=1).sum(axis=1)


def check_element_conservation(
    df: pd.DataFrame, element_list: list[str] | None = None, percent: bool = True
) -> dict[str, str]:
    """Check the conservation of elements by comparing their total
    abundance at start and end of model.

    Args:
        df (pd.DataFrame): UCLCHEM output in format from ``read_output_file``
        element_list (list[str] | None): List of elements to check. If None,
            defaults to ``uclchem.constants.default_elements_to_check``.
        percent (bool): Whether to return the change formatted as a percentage. Default = False.

    Returns:
        dict[str, str]: Dictionary containing the change in the total abundance of each element
            as a fraction of initial value

    """
    if element_list is None:
        element_list = default_elements_to_check
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


def get_total_swap(
    rates: pd.DataFrame, abundances: pd.DataFrame, reactions: list[Reaction]
) -> np.ndarray:
    """Obtain the amount of 'random' swapping per timestep.

    Args:
        rates (pd.DataFrame): The rates obtained from running an UCLCHEM model
        abundances (pd.DataFrame): The abundances obtained from running an UCLCHEM model
        reactions (list[Reaction]): The reactions used in UCLCHEM

    Returns:
        totalSwap (np.ndarray): The total swap per timestep

    Raises:
        ValueError: If ``rates`` and ``abundances`` do not have the same number of timepoints (rows),
            or ``rates`` and ``reactions`` do not have the same number of reactions (columns).
    """
    if len(rates) != len(abundances):
        msg = "Rates and abundances must be the same length"
        raise ValueError(msg)
    if rates.shape[1] != len(reactions):
        msg = "The number of rates and reactions must be equal"
        raise ValueError(msg)
    totalSwap = np.zeros(abundances.shape[0])
    for idx, reac in enumerate(reactions):
        if reac.get_reaction_type() == "BULKSWAP":
            totalSwap += rates.iloc[:, idx] * abundances[reac.get_pure_reactants()[0]]
    return totalSwap


def construct_incidence(
    species: list[str] | list[Species], reactions: list[Reaction]
) -> np.ndarray:
    """Construct the incidence matrix, a matrix that describes the in and out degree
    for each of the reactions; useful to matrix multiply by the individual rates per reaction
    to obtain a rates (dy) per species.

    Args:
        species (list[str]): A list of S species names
        reactions (list[Reaction]): The list of R reactions

    Returns:
        incidence (np.ndarray): An RxS incidence matrix

    Raises:
        ValueError: If a reactant or product of a reaction is not in the ``species``.

    """
    incidence = np.zeros(
        dtype=int,
        shape=(len(reactions), len(species)),
    )
    for reaction_idx, reaction in enumerate(reactions):
        products = reaction.get_pure_products()
        for prod in products:
            if prod not in species:
                msg = f"Product {prod} of reaction {reaction} not in species list"
                raise ValueError(msg)
            incidence[reaction_idx, species.index(prod)] += 1  # type: ignore
        reactants = reaction.get_pure_reactants()
        for reac in reactants:
            if reac not in species:
                msg = f"Reactant {prod} of reaction {reaction} not in species list"
                raise ValueError(msg)
            incidence[reaction_idx, species.index(reac)] -= 1  # type: ignore
    return incidence


def rate_constants_to_dy_and_rates(
    physics: pd.DataFrame,
    abundances: pd.DataFrame,
    rate_constants: pd.DataFrame,
    network: Network | None = None,
    species: list[Species] | None = None,
    reactions: list[Reaction] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply postprocessing to obtain the equivalent of GETYDOT from the fortran side
    and the reaction rates at each timestep.

    Args:
        physics (pd.DataFrame): The physics output from running a model
        abundances (pd.DataFrame): The abundances output from running a model
        rate_constants (pd.DataFrame): The rate constants output from running a model
        network (Network | None): The reaction network used to postprocess the rate constants.
            Defaults to None.
        species (list[Species]): The species used to postprocess the rate constants.
            Defaults to None.
        reactions (list[Reaction] | None): The reactions used to postprocess the rate constants.
            Defaults to None.

    Returns:
        ydot (pd.DataFrame): the RHS that is solved in UCLCHEM at every output timestep
        rate_by_reaction (pd.DataFrame): the individual terms that result in ydot when multiplied
            by the incidence matrix.

    Raises:
        ValueError: If ``species`` is specified, but ``reactions`` is not, or vice versa
        ValueError: If ``species``, ``reactions`` and ``network`` are all specified,
            or all not specified.
        ValueError: If there are any reaction types not processed.
        ValueError: If there is a mismatch between ``SURFSWAP`` and ``BULKSWAP`` reactions, i.e.
            the reactions with the same indices do not have their products as the others reactants.
        NotImplementedError: If no bulk species are detected in the network.
        NotImplementedError: If UCLCHEM was compiled with
            ``surfacereactions.usegarrod2011transfer = .FALSE.``. Only Garrod 2011 geometric transfer
            is currently implemented in this function.
        RuntimeError: If there is no reaction ``@H2 + BULKSWAP -> #H2`` in the network.

    """
    if (species is None) != (reactions is None):
        msg = "If species is specified, reactions also must be and vice versa"
        raise ValueError(msg)

    if (network is None) == (species is None and reactions is None):
        msg = "Choose between providing a network OR (species AND reactions). A network can be obtained using ``uclchem.makerates.network.Network.from_csv()``"
        raise ValueError(msg)

    species_list: list[Species]
    reactions_list: list[Reaction]
    if network is not None:
        species_list = network.get_species_list()
        reactions_list = network.get_reaction_list()
    else:
        species_list = species  # type: ignore[assignment]
        reactions_list = reactions  # type: ignore[assignment]

    if "Point" in rate_constants.columns:
        warnings.warn(
            "Found column `Point` in columns of `rate_constants`. uclchem.analysis.rate_constants_to_dy_and_rates is not designed for multiple points, dropping the column."
        )
        rate_constants = rate_constants.drop(columns=["Point"])

    if "@" not in "".join(spec.get_name() for spec in species_list):
        msg = "uclchem.analysis.rate_constants_to_dy_and_rates is only implemented for three-phase networks,"
        msg += " but no bulk species were found in the network."
        raise NotImplementedError(msg)

    # Import all of the constants directly from UCLCHEMWRAP to avoid discrepancies
    # ruff: noqa: N806
    GAS_DUST_DENSITY_RATIO = surfacereactions.gas_dust_density_ratio
    NUM_SITES_PER_GRAIN = surfacereactions.num_sites_per_grain
    NUM_MONOLAYERS_IS_SURFACE = surfacereactions.num_monolayers_is_surface

    # Calculate safe values to avoid division by zero
    safeBulk = abundances["BULK"].apply(lambda x: max(1.0e-30, x)).copy()
    safeMantle = abundances["SURFACE"].apply(lambda x: max(1.0e-30, x)).copy()
    ratioSurfaceToBulk = (safeMantle / safeBulk).apply(lambda x: min(1.0, x)).copy()

    # Compute dynamic quantities that can be precomputed
    bulkLayersReciprocal = (
        NUM_SITES_PER_GRAIN / (GAS_DUST_DENSITY_RATIO * safeBulk)
    ).apply(lambda x: min(1.0, x))

    totalSwap = ratioSurfaceToBulk * get_total_swap(
        rate_constants, abundances, reactions_list
    )
    # ruff: noqa: N806

    # Create the incidence matrix, we use this to evaluate the rates and
    incidence = construct_incidence(species_list, reactions_list)

    # Create a copy so we don't overwrite
    rate_by_reaction = rate_constants.copy()

    # Keep track of missing reaction types:
    missing_reactions = set()

    # Iterate over each of the rates and compute the contribution to dy for that reaction.
    for reaction_idx, reaction in enumerate(reactions_list):
        rate = rate_by_reaction.iloc[:, reaction_idx]
        # Multiply by density the right number of times:
        rate *= physics["Density"] ** reaction.body_count
        # Multiply by the abundances:
        for reactant in reaction.get_sorted_reactants():
            if reactant in list(abundances.columns):
                rate *= abundances[reactant]

        reaction_type = reaction.get_reaction_type()

        # GAR needs an additional factor of density
        if reaction_type in ["GAR"]:
            rate *= physics["Density"]
        # BULKSWAP reactions use ratioSurfaceToBulk
        elif reaction_type == "BULKSWAP":
            rate *= ratioSurfaceToBulk
        # LH/LHDES bulk reactions
        elif reaction_type in ["LH", "LHDES"]:
            if "@" in reaction.get_sorted_reactants()[0]:
                rate *= bulkLayersReciprocal
        # ED reactions multiply by #H2 abundance
        elif reaction_type == "ED":
            rate *= abundances["#H2"]
        elif reaction_type == "SURFSWAP":
            rate *= totalSwap / safeMantle
        elif reaction_type in ["DESCR", "DEUVCR", "ER", "ERDES"]:
            rate /= safeMantle
        elif reaction_type == "DESOH2":
            rate *= abundances["H"] / safeMantle
        elif reaction_type == "H2FORM":
            # H2FORM only uses 1 factor of H abundance (Cazaux & Tielens 2004)
            rate /= abundances["H"]
        elif reaction_type not in [
            # Standard reactions that require no additional prefactors
            "PHOTON",
            "CRP",
            "CRPHOT",
            "FREEZE",
            "THERM",
            "TWOBODY",
            "IONOPOL1",
            "IONOPOL2",
            "CRS",
            "EXSOLID",
            "EXRELAX",
        ]:
            missing_reactions.add(reaction_type)

        rate_by_reaction.iloc[:, reaction_idx] = rate

    if missing_reactions:
        msg = f"Missing reaction types in rate processing: {missing_reactions}"
        raise ValueError(msg)

    # Compute the rate at each timestep, adding the appropriate header
    dy = rate_by_reaction @ incidence
    dy.columns = [s.get_name() for s in species_list]
    # Compute the SURFACE and BULK:
    dy.loc[:, "SURFACE"] = dy.loc[:, dy.columns.str.startswith("#")].sum(axis=1)
    dy.loc[:, "BULK"] = dy.loc[:, dy.columns.str.startswith("@")].sum(axis=1)

    # Compute the corrections to account for the transport between SURFACE and BULK
    swap_rate_correction = pd.DataFrame()
    # Then correct for the addition or subtraction of material to/from the bulk:
    surfswap_reactions = [
        r for r in reactions_list if r.get_reaction_type() == "SURFSWAP"
    ]
    bulkswap_reactions = [
        r for r in reactions_list if r.get_reaction_type() == "BULKSWAP"
    ]
    # Sort them in the correct order, s.t. it matches the saving to disk format.
    surfswap_reactions = sorted(
        surfswap_reactions,
        key=lambda x: species_list.index(x.get_reactants()[0]),  # type: ignore
    )

    bulkswap_reactions = sorted(
        bulkswap_reactions,
        key=lambda x: species_list.index(x.get_reactants()[0]),  # type: ignore
    )

    H2_bulkswap_index = None
    for i, reaction in enumerate(bulkswap_reactions):
        if reaction.get_reactants()[0] == "@H2":
            H2_bulkswap_index = i
            break
    if H2_bulkswap_index is None:
        msg = "Did not find a BULKSWAP reaction of @H2 + BULKSWAP -> #H2"
        raise RuntimeError(msg)

    surfswap_reactions.insert(
        H2_bulkswap_index,
        Reaction(["#H2", "SURFSWAP", "NAN", "@H2", "NAN", "NAN", "NAN"] + [0] * 5),
    )

    for surfswap_reaction, bulkswap_reaction in zip(
        surfswap_reactions, bulkswap_reactions
    ):
        if surfswap_reaction.get_reactants()[0] != bulkswap_reaction.get_products()[0]:
            msg = "Mismatch for bulkswap and surfswap reactions.\n"
            msg += "\tSurfswap: {surfswap_reaction}\n"
            msg += "\tBulkswap: {bulkswap_reaction}"
            raise ValueError(msg)

    if not surfacereactions.usegarrod2011transfer:
        msg = "Can only calculate transfer reactions for Garrod 2011 transfer."
        msg += " HH transfer is not implemented yet in uclchem.analysis.rate_constants_to_dy_and_rates"
        raise NotImplementedError(msg)

    surfGrowthUncorrected = dy["SURFACE"].copy().values
    for time_idx, abunds_row in abundances.iterrows():
        time_idx = int(time_idx)  # type: ignore[call-overload] # just to make mypy happy
        if surfGrowthUncorrected[time_idx] < 0.0:
            surface_coverage = (
                min(1.0, safeBulk[time_idx] / safeMantle[time_idx])
            ) / safeBulk[time_idx]
        else:
            surface_coverage = GAS_DUST_DENSITY_RATIO / (
                NUM_SITES_PER_GRAIN * NUM_MONOLAYERS_IS_SURFACE
            )
            surface_coverage = (
                min(1.0, surface_coverage * safeMantle[time_idx]) / safeMantle[time_idx]
            )

        # Walk through the bulkswap and reactionswap pathways
        _sswap_rates = {}
        _bswap_rates = {}
        # TODO: vectorize this, because this is slower than it has to be.
        for r_bswap, r_sswap in zip(bulkswap_reactions, surfswap_reactions, strict=True):
            surface_species = r_sswap.get_reactants()[0]
            bulk_species = r_bswap.get_reactants()[0]
            if surfGrowthUncorrected[time_idx] < 0.0:
                # SURFACE is shrinking, so bulk must be shrinking too
                bswap = (
                    -surfGrowthUncorrected[time_idx]
                    * surface_coverage
                    * abunds_row[bulk_species]
                )
                sswap = 0.0
            else:
                # SURFACE is growing, so bulk must grow as well
                sswap = (
                    surfGrowthUncorrected[time_idx]
                    * surface_coverage
                    * abunds_row[surface_species]
                )
                bswap = 0.0

            # Call it SWAP_GEOMETRIC since it is due to the swap induced by the effect
            # of the surface layers growing, this is a geometric bookkeeping thing,
            # So the name geometric makes the most sense.
            _bswap_rates[str(r_bswap).replace("SWAP", "SWAP_GEOMETRIC")] = bswap
            _sswap_rates[str(r_sswap).replace("SWAP", "SWAP_GEOMETRIC")] = sswap

            # Immediately correct dy
            dy.loc[time_idx, bulk_species] += sswap - bswap
            dy.loc[time_idx, surface_species] += bswap - sswap

        swap_rate_correction = pd.concat(
            (
                swap_rate_correction,
                (
                    pd.DataFrame.from_dict(
                        _bswap_rates | _sswap_rates,
                        orient="index",
                    ).T
                ),
            )
        )
    swap_rate_correction = swap_rate_correction.reset_index(drop=True)
    rate_by_reaction = pd.concat((rate_by_reaction, swap_rate_correction), axis=1)

    # Correct the change in surface and bulk by summing the constituents:
    dy.loc[:, "SURFACE"] = dy.loc[:, dy.columns.str.startswith("#")].sum(axis=1)
    dy.loc[:, "BULK"] = dy.loc[:, dy.columns.str.startswith("@")].sum(axis=1)

    return (
        dy,
        rate_by_reaction,
    )


def compute_heating_per_reaction(
    rates: pd.DataFrame,
    network: Network | None = None,
    reactions: list[Reaction] | None = None,
) -> pd.DataFrame:
    """Compute heating/cooling per reaction by multiplying rates by exothermicity.

    Args:
        rates (pd.DataFrame): (time x n_reactions) of reaction rates
        network (Network | None): Network object with exothermicity data. Default = None.
        reactions (list[Reaction] | None): List of Reaction objects (alternative to network)
            Default = None.

    Returns:
        pd.DataFrame: (time x n_reactions) of heating rates in erg/s

    Raises:
        ValueError: If ``network`` and ``reactions`` are both None or both not None.
        ValueError: If the number of reactions and rates are not the same.

    """
    if (network is None) == (reactions is None):
        msg = "Choose between passing either network OR reactions."
        raise ValueError(msg)

    reactions_list: list[Reaction]
    if network is not None:
        reactions_list = network.get_reaction_list()
    else:
        reactions_list = reactions  # type: ignore[assignment]

    if len(reactions_list) != rates.shape[1]:
        msg = "Number of reactions and rates must be equal"
        raise ValueError(msg)

    exothermicities = np.array([r.get_exothermicity() for r in reactions_list])
    return rates * exothermicities


def get_production_and_destruction(
    species: str, dataframe: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split the reaction rates into production and destruction parts for a given species.

    Args:
        species (str): Name of species to split rates for
        dataframe (pd.DataFrame): DataFrame of reaction rates

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: production rates, destruction rates

    Raises:
        RuntimeError: If no production or destruction reactions of ``species`` are found.

    """
    if not dataframe.columns.is_unique:
        # Duplicate column names, can happen with for example UMIST
        # reactions with two temperature ranges.
        dataframe = dataframe.T.groupby(by=dataframe.columns).sum().T

    reactions = [r.strip() for r in list(dataframe.columns)]

    destruction = [r for r in reactions if species in r.split(" -> ")[0].split(" + ")]
    production = [r for r in reactions if species in r.split(" -> ")[-1].split(" + ")]

    if not destruction and not production:
        msg = f"No production or destruction reactions found for species {species}"
        raise RuntimeError(msg)

    production_df = dataframe.loc[:, production]
    destruction_df = dataframe.loc[:, destruction]

    # Take into account that some reactions destroy or form a species twice as fast as
    # the rate of that reaction itself. For example: H + H -> H2, where H is destroyed
    # at double the rate. If we do not take this into account, the rates are wrong,
    # and do not sum to the correct gradient in abundances (dy).
    destruction_count = [
        r.split(" -> ")[0].split(" + ").count(species) for r in destruction
    ]
    production_count = [
        r.split(" -> ")[-1].split(" + ").count(species) for r in production
    ]

    return production_df.mul(production_count), destruction_df.mul(destruction_count)


def derive_phase_from_name(name: str) -> str:
    """Derive the phase of the species from its name.

    Args:
        name (str): name of the species.

    Returns:
        str: Phase. One of `["gas", "surface", "bulk", "ion"]`.
    """
    if name.startswith("@"):
        return "bulk"
    elif name.startswith("#"):
        return "surface"
    elif name.endswith("+"):
        return "ion"
    else:
        return "gas"


def analyze_element_per_phase(element: str, df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the total elemental abundance of a species as a function of time
    within each phase (gas, surface, bulk and ion). Allows you to check conservation
    of elements.

    Args:
        element (str): Name of element
        df (pd.DataFrame): DataFrame from ``read_output_file()``

    Returns:
        content (pd.DataFrame): Total abundance of element for all time steps in df.

    """
    content = pd.DataFrame()
    # Split the columns into ionized, gas, surface and bulk:
    for phase in ["bulk", "surface", "ion", "gas"]:
        species_to_select = [
            s for s in list(df.columns) if derive_phase_from_name(s) == phase
        ]
        _df = df.loc[:, species_to_select]
        sums = _count_element(species_to_select, element)
        content.loc[:, element + "_" + phase] = _df.mul(sums.values, axis=1).sum(axis=1)
    return content
