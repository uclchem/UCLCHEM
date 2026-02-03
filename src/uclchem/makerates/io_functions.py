##########################################################################################
"""
Functions to read in the species and reaction files and write output files
"""

import csv
import fileinput
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict

import numpy as np
import yaml

from uclchem.constants import PHYSICAL_PARAMETERS

from .network import Network
from .reaction import Reaction, reaction_types
from .species import Species


def get_default_coolants() -> list[dict]:
    """Returns the default coolant configuration for UCLCHEM.

    Returns:
        list[dict]: List of coolant dictionaries with 'file' and 'name' keys.
    """
    return [
        {"file": "ly-a.dat", "name": "H"},
        {"file": "12c+_nometa.dat", "name": "C+"},
        {"file": "16o.dat", "name": "O"},
        {"file": "12c.dat", "name": "C"},
        {"file": "12co.dat", "name": "CO"},
        {"file": "p-h2.dat", "name": "p-H2"},
        {"file": "o-h2.dat", "name": "o-H2"},
    ]


def read_species_file(file_name: Path) -> list[Species]:
    """Reads in a Makerates species file

    Args:
        fileName (str): path to file containing the species list

    Returns:
        list: List of Species objects
    """
    species_list = []
    # list to hold user defined bulk species (for adjusting binding energy)
    user_defined_bulk = []
    with open(file_name, "r") as f:
        reader = csv.reader(f, delimiter=",", quotechar="|")
        for idx, row in enumerate(reader):
            try:
                if row[0] != "NAME" and "!" not in row[0]:
                    if "@" in row[0]:
                        user_defined_bulk.append(Species(row))
                    else:
                        species_list.append(Species(row))
            except IndexError as exc:
                print(f"Error reading species file {file_name} at line {idx}")
                raise exc

    return species_list, user_defined_bulk


def read_reaction_file(
    file_name: Path, species_list: list[Species], ftype: str
) -> tuple[list[Reaction], list[Reaction]]:
    """Reads in a reaction file of any kind (user, UMIST, KIDA)
    produces a list of reactions for the network, filtered by species_list

    Args:
        file_name (str): A file name for the reaction file to read.
        species_list (list[Species]): A list of chemical species to be used in the reading.
        ftype (str): 'UMIST','UCL', or 'KIDA' to describe format of file_name

    Returns:
        list,list: Lists of kept and dropped reactions.
    """
    reactions = []
    dropped_reactions = []

    # Every reactant/product of a reaction must be in keep_list to not be dropped
    keep_list = ["", "NAN", "#", "*", "E-", "e-", "ELECTR", "@"]
    keep_list.extend(reaction_types)
    for species in species_list:
        keep_list.append(species.get_name())

    if ftype == "UMIST":
        with open(file_name, "r") as f:
            reader = csv.reader(f, delimiter=":", quotechar="|")
            for row in reader:
                if row[0].startswith("#") or row[0].startswith("!"):
                    continue
                reaction_row = row[2:4] + [""] + row[4:8] + row[9:14] + [""]
                if check_reaction(reaction_row, keep_list):
                    reactions.append(Reaction(reaction_row, reaction_source="UMIST"))
    elif ftype == "UCL":
        with open(file_name, "r") as f:
            reader = csv.reader(f, delimiter=",", quotechar="|")
            for row in reader:
                if (len(row) > 1) and (row[0][0] != "!"):
                    if check_reaction(row, keep_list):
                        reactions.append(Reaction(row, reaction_source="UCL"))
                    else:
                        dropped_reactions.append(row)

    elif ftype == "KIDA":
        for row in kida_parser(file_name):
            if check_reaction(row, keep_list):
                reactions.append(Reaction(row, reaction_source="KIDA"))

    else:
        raise ValueError("Reaction file type must be one of 'UMIST', 'UCL' or 'KIDA'")
    return reactions, dropped_reactions


def check_reaction(reaction_row, keep_list) -> bool:
    """Checks a row parsed from a reaction file and checks it only contains acceptable things.
    It checks if all species in the reaction are present, and adds the temperature range is none is specified.

    Args:
        reaction_row (list): List parsed from a reaction file and formatted to be able to called Reaction(reaction_row)
        keep_list (list): list of elements that are acceptable in the reactant or product bits of row

    Returns:
        bool: Whether the row contains acceptable entries.
    """
    if all(x.upper() in keep_list for x in reaction_row[0:7]):
        if reaction_row[10] == "":
            reaction_row[10] = 0.0
            reaction_row[11] = 10000.0
        if reaction_row[12] == "":
            reaction_row[12] = 0.0
        return True
    else:
        if reaction_row[1] in ["DESORB", "FREEZE"]:
            reac_error = "Desorb or freeze reaction in custom input contains species not in species list"
            reac_error += f"\nReaction was {reaction_row}"
            raise ValueError(reac_error)
        return False


def kida_parser(kida_file):
    """
    KIDA used a fixed format file so we read each line in the chunks they specify
    and use python built in classes to convert to the necessary types.
    NOTE KIDA defines some of the same reaction types to UMIST but with different names
    and coefficients. We fix that by converting them here.
    """

    def str_parse(x):
        return str(x).strip().upper()

    kida_contents = [
        [3, {str_parse: 11}],
        [1, {"skip": 1}],
        [5, {str_parse: 11}],
        [1, {"skip": 1}],
        [3, {float: 10, "skip": 1}],
        [1, {"skip": 27}],
        [2, {int: 6, "skip": 1}],
        [1, {int: 2}],
        [1, {"skip": 11}],
    ]
    rows = []
    with open(kida_file, "r") as f:
        f.readline()  # throw away header
        for line in f:  # then iterate over file
            if line.startswith("!"):
                continue
            row = []
            for item in kida_contents:
                for i in range(item[0]):
                    for func, count in item[1].items():
                        if func != "skip":
                            a = line[:count]
                            row.append(func(a))
                        line = line[count:]

            # Some reformatting required
            # KIDA gives CRP reactions in different units to UMIST
            if row[-1] == 1:
                # Amazingly both UMIST and KIDA use CRP but differently.
                # Translate KIDA names to UMIST
                if row[1] == "CRP":
                    row[1] = "CRPHOT"
                    # with beta=0 and gamma=1, the KIDA formulation of
                    # CRPHOT reactions becomes the UMIST one
                    row[10] = 1.0
                elif row[1] == "CR":
                    row[1] = "CRP"
                # UMIST alpha includes zeta_0 but KIDA doesn't. Since UCLCHEM
                # rate calculation follows UMIST, we convert.
                row[8] = row[8] * 1.36e-17
                rows.append(row[:7] + row[8:-1])
            elif row[-1] in [2, 3]:
                rows.append(row[:7] + row[8:-1])
            elif row[-1] == 4:
                row[2] = "IONOPOL1"
                rows.append(row[:7] + row[8:-1])
            elif row[-1] == 5:
                row[2] = "IONOPOL2"
                rows.append(row[:7] + row[8:-1])
    return rows


def read_grain_assisted_recombination_file(file_name: Path) -> dict:
    with open(file_name, "r") as fh:
        gar_parameters = yaml.safe_load(fh)
    return gar_parameters


def output_drops(
    dropped_reactions: list[Reaction], output_dir: str = None, write_files: bool = True
):
    """Writes the reactions that are dropped to disk/logs

    Args:
        dropped_reactions (list[Reaction]): The reactions that were dropped
        output_dir (str): The directory that dropped_reactions.csv will be written to.
        write_files (bool, optional): Whether or not to write the file. Defaults to True.
    """
    if output_dir is None:
        output_dir = "."
    outputFile = Path(output_dir) / "dropped_reactions.csv"
    # Print dropped reactions from grain file or write if many
    if write_files and dropped_reactions:
        logging.info(f"\nReactions dropped from grain file written to {outputFile}\n")
        with open(outputFile, "w") as f:
            writer = csv.writer(
                f,
                delimiter=",",
                quotechar="|",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator="\n",
            )
            for reaction in dropped_reactions:
                writer.writerow(reaction)
    else:
        logging.info("Reactions dropped from grain file:\n")
        for reaction in dropped_reactions:
            logging.info(reaction)


def write_outputs(
    network: Network,
    output_dir: str = None,
    enable_rates_storage: bool = False,
    gar_database: dict[str, np.array] = None,
    coolants: list[dict] = None,
) -> None:
    """Write the ODE and Network fortran source files to the fortran source.

    Args:
        network (network): The makerates Network class
        output_dir (bool): The directory to write to.
        coolants (list[dict]): List of coolant dictionaries with 'file' and 'name' keys.
                               If None, uses default coolants.
    """
    if output_dir is None:
        output_dir = Path("../src/uclchem")
        fortran_src_dir = Path("../src/fortran_src")
    else:
        output_dir = Path(output_dir)
        fortran_src_dir = Path(output_dir)

    # Use default coolants if none provided
    if coolants is None:
        coolants = get_default_coolants()

    # Create the species file
    filename = output_dir / "species.csv"
    write_species(filename, network.get_species_list())

    filename = output_dir / "reactions.csv"
    write_reactions(filename, network.get_reaction_list())

    # Write the ODEs in the appropriate language format
    filename = fortran_src_dir / "odes.f90"
    write_odes_f90(
        filename,
        network.get_species_list(),
        network.get_reaction_list(),
        enable_rates_storage=enable_rates_storage,
    )

    # Write the network files
    filename = fortran_src_dir / "network.f90"
    write_network_file(
        filename,
        network,
        enable_rates_storage=enable_rates_storage,
        gar_database=gar_database,
    )
    # write the constants needed for wrap.f90

    filename = fortran_src_dir / "f2py_constants.f90"

    f2py_constants = {
        "n_species": len(network.get_species_list()),
        "n_reactions": len(network.get_reaction_list()),
        "n_physical_parameters": len(PHYSICAL_PARAMETERS),
        "n_coolants": len(coolants),
        "coolant_files": [c["file"] for c in coolants],
        "coolant_names": [c["name"] for c in coolants],
    }
    write_f90_constants(f2py_constants, filename)
    # Write some meta information that can be used to read back in the reactions into Python
    # TODO: we can update this s.t. we can directly access the f90 constants from the f2py wrapped Fortran codes
    write_python_constants(f2py_constants, "../src/uclchem/constants.py")


def write_f90_constants(
    replace_dict: Dict[str, int],
    output_file_name: Path,
    template_file_path: Path = "fortran_templates",
) -> None:
    """Write the physical reactions to the f2py_constants.f90 file after every
    run of makerates, this ensures the Fortran and Python bits are compatible.

    Args:
        replace_dict (Dict[str, int]): The dictionary with keys to replace
        output_file_name (Path): The path to target f2py_constants.f90 file
        template_file_path (Path, optional): The file to use as the template.
    """
    _ROOT = Path(__file__).parent
    template_file_path = _ROOT / template_file_path
    with open(template_file_path / "f2py_constants.f90", "r") as fh:
        constants = fh.read()

    # Handle string arrays separately for coolants
    if "coolant_files" in replace_dict and "coolant_names" in replace_dict:
        # Format coolant files
        coolant_files = replace_dict.pop("coolant_files")
        max_file_len = max(len(f) for f in coolant_files)
        coolant_files_str = ",".join(
            f'"{f.ljust(max_file_len)}"' for f in coolant_files
        )
        replace_dict["coolant_file_len"] = max_file_len
        replace_dict["coolant_files"] = "/" + coolant_files_str + "/"

        # Format coolant names
        coolant_names = replace_dict.pop("coolant_names")
        max_name_len = max(len(n) for n in coolant_names)
        coolant_names_str = ",".join(
            f'"{n.ljust(max_name_len)}"' for n in coolant_names
        )
        replace_dict["coolant_name_len"] = max_name_len
        replace_dict["coolant_names"] = "/" + coolant_names_str + "/"

    constants = constants.format(**replace_dict)
    with open(output_file_name, "w") as fh:
        fh.writelines(constants)


def write_python_constants(
    replace_dict: Dict[str, int], python_constants_file: Path
) -> None:
    """Function to write the python constants to the constants.py file after every run,
    this ensure the Python and Fortran bits are compatible with one another.

    Args:
        replace_dict (Dict[str, int]]): Dict with keys to replace and their values
        python_constants_file (Path): Path to the target constant files.
    """
    with fileinput.input(python_constants_file, inplace=True, backup=".bak") as file:
        for line in file:
            # Add a timestamp to the file before the old one:
            if fileinput.isfirstline():
                print(
                    "# This file was machine generated with Makerates on",
                    datetime.now(),
                    end="\n",
                )
                # Don't copy the old timestamp into the new file.
                if line.startswith(
                    "# This file was machine generated with Makerates on"
                ):
                    continue
            # For every line, try to find constants, if we find them, replace them,
            # if not, just print the line.
            hits = {
                constant: line.strip().startswith(constant) for constant in replace_dict
            }
            if any(hits.values()):
                # Filter, also we can only get one hit at a time
                variable = list(filter(hits.get, hits))[0]
                print(f"{variable} = {replace_dict[variable]}")
            else:
                print(line, end="")


def write_species(file_name: Path, species_list: list[Species]) -> None:
    """Write the human readable species file. Note UCLCHEM doesn't use this file.

    Args:
        fileName (str): path to output file
        species_list (list): List of species objects for network
    """
    species_columns = [
        "NAME",
        "MASS",
        "BINDING_ENERGY",
        "SOLID_FRACTION",
        "MONO_FRACTION",
        "VOLCANO_FRACTION",
        "ENTHALPY",
    ]
    with open(file_name, "w") as f:
        writer = csv.writer(
            f,
            delimiter=",",
            quotechar="|",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )
        writer.writerow(species_columns)
        for species in species_list:
            writer.writerow(
                [
                    species.get_name(),
                    species.get_mass(),
                    species.get_binding_energy(),
                    species.get_solid_fraction(),
                    species.get_mono_fraction(),
                    species.get_volcano_fraction(),
                    species.get_enthalpy(),
                ]
            )


# Write the reaction file in the desired format
def write_reactions(fileName, reaction_list) -> None:
    """Write the human readable reaction file.

    Args:
        fileName (str): path to output file
        reaction_list (list): List of reaction objects for network
    """
    reaction_columns = [
        "Reactant 1",
        "Reactant 2",
        "Reactant 3",
        "Product 1",
        "Product 2",
        "Product 3",
        "Product 4",
        "Alpha",
        "Beta",
        "Gamma",
        "T_min",
        "T_max",
        "reduced_mass",
        "extrapolate",
        "exothermicity",
    ]
    with open(fileName, "w") as f:
        writer = csv.writer(
            f,
            delimiter=",",
            quotechar="|",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )
        writer.writerow(reaction_columns)
        for reaction in reaction_list:
            writer.writerow(
                reaction.get_reactants()
                + reaction.get_products()
                + [
                    reaction.get_alpha(),
                    reaction.get_beta(),
                    reaction.get_gamma(),
                    reaction.get_templow(),
                    reaction.get_temphigh(),
                    reaction.get_reduced_mass(),
                    reaction.get_extrapolation(),
                    reaction.get_exothermicity(),
                ]
            )


def write_odes_f90(
    file_name: Path,
    species_list: list[Species],
    reaction_list: list[Reaction],
    enable_rates_storage: bool = False,
) -> None:
    """Write the ODEs in Modern Fortran. This is an actual code file.

    Args:
        file_name (str): Path to file where code will be written
        species_list (list): List of species describing network
        reaction_list (list): List of reactions describing network
    """
    # First generate ODE contributions for all reactions
    species_names = [spec.get_name() for spec in species_list]

    [
        logging.debug(f"{species_names.index(specie) + 1}:{specie}")
        for specie in species_list
    ]

    for i, reaction in enumerate(reaction_list):
        logging.debug(f"RATE({i + 1}):{reaction}")
        reaction.generate_ode_bit(i, species_names)

    # then create ODE code and write to file.
    with open(file_name, mode="w") as output:
        # go through every species and build two strings, one with eq for all destruction routes and one for all formation
        ydotString = build_ode_string(species_list, reaction_list, enable_rates_storage)
        output.write(ydotString)


def write_jacobian(file_name: Path, species_list: list[Species]) -> None:
    """Write jacobian in Modern Fortran. This has never improved UCLCHEM's speed
    and so is not used in the code as it stands.
    Current only works for three phase model.

    Args:
        file_name (str): Path to jacobian file
        species_list (species_list): List of species AFTER being processed by build_ode_string
    """
    output = open(file_name, "w")
    species_names = ""
    for i, species in enumerate(species_list):
        species_names += species.get_name()
        losses = species.losses.split("+")
        gains = species.gains.split("+")
        for j in range(1, len(species_list) + 1):
            if species.get_name() == "SURFACE":
                di_dj = f"J({i + 1},{j})=SUM(J(surfaceList,{j}))\n"
                output.write(di_dj)
            elif species.get_name() == "BULK":
                if species_names.count("@") > 0:
                    di_dj = f"J({i + 1},{j})=SUM(J(bulkList,{j}))\n"
                    output.write(di_dj)
            else:
                # every time an ode bit has our species in it, we remove it (dy/dx=a for y=ax)
                di_dj = [
                    f"-{x}".replace(f"*Y({j})", "", 1)
                    for x in losses
                    if f"*Y({j})" in x
                ]
                di_dj += [
                    f"+{x}".replace(f"*Y({j})", "", 1) for x in gains if f"*Y({j})" in x
                ]
                # of course there might be y=a*x*x so we only replace first instance and if there's still an instance
                # we put a factor of two in since dy/dx=2ax for y=a*x*x
                di_dj = [x + "*2" if f"*Y({j})" in x else x for x in di_dj]

                # safeMantle is a stand in for the surface so do it manually here
                # since it's divided by safemantle, derivative is negative so sign flips and we get another factor of 1/safeMantle
                if species_list[j - 1].get_name() == "SURFACE":
                    di_dj = [f"+{x}/safeMantle" for x in losses if "/safeMantle" in x]
                    di_dj += [f"-{x}/safeMantle" for x in gains if "/safeMantle" in x]
                if len(di_dj) > 0:
                    di_dj = f"J({i + 1},{j})=" + "".join(di_dj) + "\n"
                    output.write(di_dj)

        # tackle density separately.handled
        j = j + 1
        if species.get_name() == "SURFACE":
            di_dj = f"J({i + 1},{j})=SUM(J(surfaceList,{j}))\n"
            output.write(di_dj)
        elif species.get_name() == "BULK":
            if species_names.count("@") > 0:
                di_dj = f"J({i + 1},{j})=SUM(J(bulkList,{j}))\n"
                output.write(di_dj)
        else:
            di_dj = [f"-{x}".replace("*D", "", 1) for x in losses if "*D" in x]
            di_dj += [f"+{x}".replace("*D", "", 1) for x in gains if "*D" in x]
            di_dj = [x + "*2" if "*D" in x else x for x in di_dj]
            if len(di_dj) > 0:
                di_dj = f"J({i + 1},{j})=" + ("".join(di_dj)) + "\n"
                output.write(di_dj)
    i = i + 2
    di_dj = f"J({i},{i})=ddensdensdot(D)\n"
    output.write(di_dj)

    output.close()


def build_ode_string(
    species_list: list[Species],
    reaction_list: list[Reaction],
    enable_rates_storage: bool = False,
) -> str:
    """A long, complex function that does the messy work of creating the actual ODE
    code to calculate the rate of change of each species. Test any change to this code
    thoroughly because ODE mistakes are very hard to spot.

    Args:
        species_list (list): List of species in network
        reaction_list (list): List of reactions in network
        enable_rates_storage (bool): Enable the writing of the rates to the disk.

    Returns:
        str: One long string containing the entire ODE fortran code.
    """

    # We create a string of losses and gains for each species so initialize them all as ""
    species_names = []
    for i, species in enumerate(species_list):
        species_names.append(species.get_name())
        species.losses = ""
        species.gains = ""

    bulk_index = species_names.index("BULK")
    surface_index = species_names.index("SURFACE")
    total_swap = ""

    for i, reaction in enumerate(reaction_list):
        for species in reaction.get_reactants():
            if species in species_names:
                # Eley-Rideal reactions take a share of total freeze out rate which is already accounted for
                # so we add as a loss term to the frozen version of the species rather than the gas version
                if (reaction.get_reaction_type() == "ER") and (
                    not species_list[species_names.index(species)].is_surface_species()
                ):
                    species_list[
                        species_names.index("#" + species)
                    ].losses += reaction.ode_bit
                else:
                    species_list[
                        species_names.index(species)
                    ].losses += reaction.ode_bit
                if reaction.get_reaction_type() == "BULKSWAP":
                    total_swap += reaction.ode_bit
        for species in reaction.get_products():
            if species in species_names:
                species_list[species_names.index(species)].gains += reaction.ode_bit

    ode_string = """MODULE ODES
USE constants
USE network
IMPLICIT NONE
CONTAINS
SUBROUTINE GETYDOT(RATE, Y, ratioSurfaceToBulk, surfaceCoverage, safeMantle, safebulk, D, YDOT)
REAL(dp), INTENT(IN) :: RATE(:), Y(:), ratioSurfaceToBulk, safeMantle, safebulk, D
REAL(dp), INTENT(INOUT) :: YDOT(:), surfaceCoverage
REAL(dp) :: totalSwap, LOSS, PROD
    """
    # Add a logical to determine whether we can write the reaction rates in realtime
    ode_string += truncate_line(f"totalSwap={total_swap[1:]}\n\n")
    # First get total rate of change of bulk and surface by adding ydots
    for n, species in enumerate(species_list):
        if species.get_name()[0] == "@":
            species_list[bulk_index].gains += f"+YDOT({n + 1})"
        elif species.get_name()[0] == "#":
            species_list[surface_index].gains += f"+YDOT({n + 1})"
    if enable_rates_storage:
        for n, reaction in enumerate(reaction_list):
            ode_string += truncate_line(f"REACTIONRATE({n + 1})={reaction.ode_bit}\n")

    for n, species in enumerate(species_list):
        ydot_string = species_ode_string(n, species)
        ode_string += ydot_string

    ode_string += f"    SURFGROWTHUNCORRECTED = YDOT({surface_index + 1})\n"

    # now add bulk transfer to rate of change of surface species after they've already been calculated
    ode_string += "!Update surface species for bulk growth\n"

    ode_string += f"IF (YDOT({surface_index + 1}) .lt. 0) THEN\n"
    ode_string += "    ! Since ydot(surface_index) is negative, bulk is lost and surface forms\n"
    ode_string += "    IF (useGarrod2011Transfer) THEN\n"
    ode_string += "        ! Three-phase treatment of Garrod & Pauly 2011\n"
    ode_string += "        ! Replace surfaceCoverage with alpha_des\n"
    ode_string += "        ! Real value of alpha_des: alpha_des = MIN(1.0D0, safeBulk / safeMantle).\n"
    ode_string += "        ! However, the YDOTs calculated below need to be multiplied with Y(bulkspec)/safeBulk,\n"
    ode_string += "        ! so we divide by safeBulk here to save time\n"
    ode_string += (
        "        surfaceCoverage = MIN(1.0D0, safeBulk/safeMantle)/safeBulk\n"
    )
    ode_string += "    ELSE\n        ! Hasegawa & Herbst 1993\n"
    ode_string += "        surfaceCoverage = MIN(1.0D0, surfaceCoverage*safeMantle)/safeBulk\n"
    ode_string += f"    END IF\n"

    surf_species = [
        i
        for i in species_list
        if i.get_name() not in ["SURFACE", "BULK"] and i.is_surface_species()
    ]
    i = len(reaction_list)
    j = len(reaction_list) + len(surf_species)
    for n, species in enumerate(species_list):
        if species.get_name()[0] == "#":
            i += 1
            j += 1
            bulk_partner = species_names.index(species.get_name().replace("#", "@"))
            if enable_rates_storage:
                ode_string += f"    REACTIONRATE({i}) = -YDOT({surface_index + 1})*surfaceCoverage*Y({bulk_partner + 1})\n"
                ode_string += f"    REACTIONRATE({j}) = 0.0D0\n"
            if not species_list[bulk_partner].is_refractory:
                ode_string += f"    YDOT({n + 1})=YDOT({n + 1})-YDOT({surface_index + 1})*surfaceCoverage*Y({bulk_partner + 1})\n"
        if species.get_name()[0] == "@":
            if not species.is_refractory:
                ode_string += f"    YDOT({n + 1})=YDOT({n + 1})+YDOT({surface_index + 1})*surfaceCoverage*Y({n + 1})\n"
    ode_string += "ELSE\n"
    ode_string += "    ! surfaceCoverage = fractional surface coverage\n"
    ode_string += "    ! Real value of surfaceCoverage: surfaceCoverage = safeMantle / NUM_MONOLAYERS_IS_SURFACE * GAS_DUST_DENSITY_RATIO / NUM_SITES_PER_GRAIN\n"
    ode_string += "    ! However, the YDOTs calculated below need to be multiplied with Y(surfspec)/safeMantle, so we divide by safeMantle here to save time\n"
    ode_string += "    ! In chemistry.f90: surfaceCoverage = 1/NUM_MONOLAYERS_IS_SURFACE * GAS_DUST_DENSITY_RATIO / NUM_SITES_PER_GRAIN\n"
    ode_string += (
        "    surfaceCoverage = MIN(1.0D0, surfaceCoverage*safeMantle)/safeMantle\n"
    )
    i = len(reaction_list)
    j = len(reaction_list) + len(surf_species)
    for n, species in enumerate(species_list):
        if species.get_name() in [
            "#H2",
            "@H2",
        ]:  # Do not allow H2 to transfer from surface to bulk
            if species.get_name() == "@H2":
                i += 1
                j += 1
            continue
        if species.get_name()[0] == "@":
            i += 1
            j += 1
            surface_version = species_names.index(species.get_name().replace("@", "#"))
            if enable_rates_storage:
                ode_string += f"    REACTIONRATE({i}) = 0.0D0\n"
                ode_string += f"    REACTIONRATE({j}) = -YDOT({surface_index + 1})*surfaceCoverage*Y({surface_version + 1})\n"
            ode_string += f"    YDOT({n + 1})=YDOT({n + 1})+YDOT({surface_index + 1})*surfaceCoverage*Y({surface_version + 1})\n"
        if species.get_name()[0] == "#":
            ode_string += f"    YDOT({n + 1})=YDOT({n + 1})-YDOT({surface_index + 1})*surfaceCoverage*Y({n + 1})\n"
    ode_string += "ENDIF\n"

    # once bulk transfer has been added, odes for bulk and surface must be updated to account for it
    ode_string += "!Update total rate of change of bulk and surface for bulk growth\n"
    ode_string += species_ode_string(bulk_index, species_list[bulk_index])
    ode_string += species_ode_string(surface_index, species_list[surface_index])
    ode_string += """    END SUBROUTINE GETYDOT
END MODULE ODES"""
    return ode_string


def species_ode_string(n: int, species: Species) -> str:
    """Build the string of Fortran code for a species once it's loss and gains
    strings have been produced.

    Args:
        n (int): Index of species in python format
        species (Species): species object

    Returns:
        str: the fortran code for the rate of change of the species
    """
    ydot_string = ""
    if species.losses != "":
        loss_string = "    LOSS = " + species.losses[1:] + "\n"
        ydot_string += loss_string
    if species.gains != "":
        prod_string = "    PROD = " + species.gains[1:] + "\n"
        ydot_string += prod_string

    if ydot_string != "":
        ydot_string += f"    YDOT({n + 1}) = "
        # start with empty string and add production and loss terms if they exists
        if species.gains != "":
            ydot_string += "PROD"
        if species.losses != "":
            ydot_string += "-LOSS"
        ydot_string += "\n"
    else:
        ydot_string = f"    YDOT({n + 1}) = {0.0}\n"
    ydot_string = truncate_line(ydot_string)
    return ydot_string


def write_evap_lists(network_file, species_list: list[Species]) -> int:
    """Two phase networks mimic episodic thermal desorption seen in lab (see Viti et al. 2004)
    by desorbing fixed fractions of material at specific temperatures. Three phase networks just
    use binding energy and that fact we set binding energies in bulk to water by default.
    This function writes all necessary arrays to the network file so these processes work.

    Args:
        network_file (file): Open file object to which the network code is being written
        species_list (list[Species]): List of species in network
    """
    gasIceList = []
    surfacelist = []
    solidList = []
    monoList = []
    volcList = []
    binding_energyList = []
    enthalpyList = []
    inertiaProductsList = []
    moleculeIsLinearList = []
    bulkList = []
    iceList = []
    refractoryList = []
    species_names = [spec.get_name() for spec in species_list]
    for i, species in enumerate(species_list):
        if species.get_name()[0] == "#":
            # find gas phase version of grain species. For #CO it looks for first species in list with just CO and then finds the index of that
            try:
                j = species_names.index(species.get_desorb_products()[0])
            except ValueError:
                error = f"{species.get_name()} desorbs as {species.get_desorb_products()[0]}"
                error += "which is not in species list. This desorption is likely user defined.\n"
                error += "Please amend the desorption route in your reaction file and re-run Makerates"
                raise NameError(error)

            # plus ones as fortran and python label arrays differently
            surfacelist.append(i + 1)
            gasIceList.append(j + 1)
            solidList.append(species.get_solid_fraction())
            monoList.append(species.get_mono_fraction())
            volcList.append(species.get_volcano_fraction())
            iceList.append(i + 1)
            binding_energyList.append(species.get_binding_energy())
            enthalpyList.append(species.get_enthalpy())
        elif species.get_name()[0] == "@":
            j = species_names.index(species.get_desorb_products()[0])
            gasIceList.append(j + 1)
            bulkList.append(i + 1)
            iceList.append(i + 1)
            binding_energyList.append(species.get_binding_energy())
            enthalpyList.append(species.get_enthalpy())
            inertiaProductsList.append(species.calculate_rotational_partition_factor())
            moleculeIsLinearList.append(species.is_linear())
            if species.is_refractory:
                refractoryList.append(i + 1)

    # dummy index that will be caught by UCLCHEM
    if len(refractoryList) == 0:
        refractoryList = [-999]

    network_file.write(array_to_string("surfaceList", surfacelist, type="int"))
    if len(bulkList) > 0:
        network_file.write(array_to_string("bulkList", bulkList, type="int"))
    network_file.write(array_to_string("iceList", iceList, type="int"))
    network_file.write(array_to_string("gasIceList", gasIceList, type="int"))
    network_file.write(array_to_string("solidFractions", solidList, type="float"))
    network_file.write(array_to_string("monoFractions", monoList, type="float"))
    network_file.write(array_to_string("volcanicFractions", volcList, type="float"))
    network_file.write(
        array_to_string(
            "bindingEnergy", binding_energyList, type="float", parameter=False
        )
    )
    network_file.write(array_to_string("formationEnthalpy", enthalpyList, type="float"))
    network_file.write(
        array_to_string("inertiaProducts", inertiaProductsList, type="float")
    )
    network_file.write(
        array_to_string("moleculeIsLinear", moleculeIsLinearList, type="logical")
    )
    network_file.write(array_to_string("refractoryList", refractoryList, type="int"))
    return len(iceList)


def truncate_line(input_string: str, lineLength: int = 72) -> str:
    """Take a string and adds line endings at regular intervals
    keeps us from overshooting fortran's line limits and, frankly,
    makes for nicer ode.f90 even if human readability isn't very important

    Args:
        input_string (str): Line of code to be truncated
        lineLength (int, optional): rough line length. Defaults to 72.

    Returns:
        str: Code string with line endings at regular intervals
    """
    result = ""
    i = 0
    j = 0
    # we only want to split at operators to make it look nice
    splits = ["*", ")", "+", ","]
    while len(input_string[i:]) > lineLength:
        j = i + lineLength
        if "\n" in input_string[i:j]:
            j = input_string[i:j].index("\n") + i + 1
            result += input_string[i:j]
        else:
            while input_string[j] not in splits:
                j = j - 1
            result += input_string[i:j] + "&\n    &"
        i = j
    result += input_string[i:]
    return result


def write_network_file(
    file_name: Path,
    network: Network,
    enable_rates_storage: bool = False,
    gar_database=None,
):
    """Write the Fortran code file that contains all network information for UCLCHEM.
    This includes lists of reactants, products, binding energies, formationEnthalpies
    and so on.

    Args:
        file_name (str): The file name where the code will be written.
        network (Network): A Network object built from lists of species and reactions.
    """
    species_list = network.get_species_list()
    reaction_list = network.get_reaction_list()
    openFile = open(file_name, "w")
    openFile.write("MODULE network\nUSE constants\nIMPLICIT NONE\n")

    # write arrays of all species stuff
    names = []
    atoms = []
    masses = []
    for species in species_list:
        names.append(species.get_name())
        masses.append(float(species.mass))
        atoms.append(species.n_atoms)

    speciesIndices = ""
    for name, species_index in network.species_indices.items():
        speciesIndices += "{0}={1},".format(name, species_index)
    if len(speciesIndices) > 72:
        speciesIndices = truncate_line(speciesIndices)
    speciesIndices = speciesIndices[:-1] + "\n"
    openFile.write("    INTEGER(dp), PARAMETER ::" + speciesIndices)
    openFile.write("    LOGICAL, PARAMETER :: THREE_PHASE = .TRUE.\n")
    openFile.write("    REAL(dp) :: SURFGROWTHUNCORRECTED\n")
    openFile.write(array_to_string("    specname", names, type="string"))
    openFile.write(array_to_string("    mass", masses, type="float"))
    openFile.write(array_to_string("    atomCounts", atoms, type="int"))

    # then write evaporation stuff
    n_ice_species = write_evap_lists(openFile, species_list)

    # finally all reactions
    reactant1 = []
    reactant2 = []
    reactant3 = []
    prod1 = []
    prod2 = []
    prod3 = []
    prod4 = []
    alpha = []
    beta = []
    gama = []
    reacTypes = []
    # duplicates = []
    tmins = []
    tmaxs = []
    reduced_masses = []
    extrapolations = []
    exothermicity = []

    # store important reactions
    reaction_indices = ""
    for reaction, index in network.important_reactions.items():
        reaction_indices += reaction + f"={index},"
    reaction_indices = truncate_line(reaction_indices[:-1]) + "\n"
    openFile.write("    INTEGER(dp), PARAMETER ::" + reaction_indices)

    for i, reaction in enumerate(reaction_list):
        reactant1.append(find_reactant(names, reaction.get_reactants()[0]))
        reactant2.append(find_reactant(names, reaction.get_reactants()[1]))
        reactant3.append(find_reactant(names, reaction.get_reactants()[2]))
        prod1.append(find_reactant(names, reaction.get_products()[0]))
        prod2.append(find_reactant(names, reaction.get_products()[1]))
        prod3.append(find_reactant(names, reaction.get_products()[2]))
        prod4.append(find_reactant(names, reaction.get_products()[3]))
        alpha.append(reaction.get_alpha())
        beta.append(reaction.get_beta())
        gama.append(reaction.get_gamma())
        # if reaction.duplicate:
        #     duplicates.append(i + 1)
        tmaxs.append(reaction.get_temphigh())
        tmins.append(reaction.get_templow())
        reduced_masses.append(reaction.get_reduced_mass())
        reacTypes.append(reaction.get_reaction_type())
        extrapolations.append(reaction.get_extrapolation())
        exothermicity.append(reaction.get_exothermicity())
    # if len(duplicates) == 0:
    #     duplicates = [9999]
    #     tmaxs = [0]
    #     tmins = [0]

    reaction_names = []
    for n, reaction in enumerate(reaction_list):
        reaction_names.append(str(reaction))
    for n, species in enumerate(species_list):
        if species.is_surface_species() and species.get_name() not in [
            "SURFACE",
            "BULK",
        ]:
            reaction_name = (
                f"{species.get_name()} + SURFACETRANSFER -> @{species.get_name()[1:]}"
            )
            reaction_names.append(reaction_name)
    for n, species in enumerate(species_list):
        if species.is_surface_species() and species.get_name() not in [
            "SURFACE",
            "BULK",
        ]:
            reaction_name = (
                f"@{species.get_name()[1:]} + SURFACETRANSFER -> {species.get_name()}"
            )
            reaction_names.append(reaction_name)

    # Save some memory by only allocating things we actually want to use:
    if enable_rates_storage:
        openFile.write(
            f"    REAL(dp) :: REACTIONRATE({len(reactant1) + n_ice_species})\n"
        )
        openFile.write("     LOGICAL :: storeRatesComputation=.true.\n")
    else:
        openFile.write("    REAL(dp) :: REACTIONRATE(1)\n")
        openFile.write("    LOGICAL :: storeRatesComputation=.false.\n")
    if any([exo != 0.0 for exo in exothermicity]):
        assert enable_rates_storage, "Chemical heating can only be enabled if rates are being computed and stored in memory. Enable `enable_rates_storage` in the user_settings."
        openFile.write(
            array_to_string(
                "\texothermicities", exothermicity, type="float", parameter=True
            )
        )
        openFile.write("    LOGICAL, PARAMETER :: enableChemicalHeating = .TRUE.\n")
    else:
        openFile.write(
            "    REAL(dp) :: \texothermicities(" + str(len(exothermicity)) + ")\n"
        )
        openFile.write("    LOGICAL, PARAMETER :: enableChemicalHeating = .FALSE.\n")

    openFile.write(array_to_string("\tre1", reactant1, type="int"))
    openFile.write(array_to_string("\tre2", reactant2, type="int"))
    openFile.write(array_to_string("\tre3", reactant3, type="int"))
    openFile.write(array_to_string("\tp1", prod1, type="int"))
    openFile.write(array_to_string("\tp2", prod2, type="int"))
    openFile.write(array_to_string("\tp3", prod3, type="int"))
    openFile.write(array_to_string("\tp4", prod4, type="int"))
    openFile.write(array_to_string("\talpha", alpha, type="float", parameter=False))
    openFile.write(array_to_string("\tbeta", beta, type="float", parameter=False))
    openFile.write(array_to_string("\tgama", gama, type="float", parameter=False))
    # openFile.write(array_to_string("\tduplicates", duplicates, type="int", parameter=True))
    openFile.write(array_to_string("\tminTemps", tmins, type="float", parameter=True))
    openFile.write(array_to_string("\tmaxTemps", tmaxs, type="float", parameter=True))
    openFile.write(
        array_to_string("\treducedMasses", reduced_masses, type="float", parameter=True)
    )
    openFile.write(
        array_to_string(
            "\tExtrapolateRates", extrapolations, type="logical", parameter=True
        )
    )
    reacTypes = np.asarray(reacTypes)

    partners = get_desorption_freeze_partners(reaction_list)
    openFile.write(
        array_to_string("\tfreezePartners", partners, type="int", parameter=True)
    )

    openFile.write(
        array_to_string(
            "\t garParams",
            np.array(list(gar_database.values())) if gar_database else np.zeros((1, 7)),
            type="float",
            parameter=True,
        )
    )

    for reaction_type in reaction_types:
        list_name = reaction_type.lower() + "Reacs"
        indices = np.where(reacTypes == reaction_type)[0]
        if len(indices > 1):
            indices = [indices[0] + 1, indices[-1] + 1]
        else:
            # We still want a dummy array if the reaction type isn't in network
            indices = [99999, 99999]
        openFile.write(
            array_to_string("\t" + list_name, indices, type="int", parameter=True)
        )

    # Write LHDES and ERDES mapping arrays (Feature 3: LH/ER-DES mapping)
    # These arrays map chemical reactive desorption reactions to their parent reactions
    LHDEScorrespondingLHreacs = []
    for reaction in reaction_list:
        if reaction.get_reaction_type() == "LHDES":
            if hasattr(reaction, 'get_partner') and reaction.get_partner() is not None:
                partner = reaction.get_partner()
                reacIndex = reaction_list.index(partner) + 1
                LHDEScorrespondingLHreacs.append(reacIndex)
            else:
                # If no partner set, use dummy index
                LHDEScorrespondingLHreacs.append(99999)

    # Write array (use dummy if empty for backward compatibility)
    if len(LHDEScorrespondingLHreacs) == 0:
        LHDEScorrespondingLHreacs = [99999]
    openFile.write(
        array_to_string(
            "\tLHDEScorrespondingLHreacs",
            LHDEScorrespondingLHreacs,
            type="int",
            parameter=True,
        )
    )

    ERDEScorrespondingERreacs = []
    for reaction in reaction_list:
        if reaction.get_reaction_type() == "ERDES":
            if hasattr(reaction, 'get_partner') and reaction.get_partner() is not None:
                partner = reaction.get_partner()
                reacIndex = reaction_list.index(partner) + 1
                ERDEScorrespondingERreacs.append(reacIndex)
            else:
                # If no partner set, use dummy index
                ERDEScorrespondingERreacs.append(99999)

    # Write array (use dummy if empty for backward compatibility)
    if len(ERDEScorrespondingERreacs) == 0:
        ERDEScorrespondingERreacs = [99999]
    elif len(ERDEScorrespondingERreacs) == 1:
        # Fortran needs at least 2 elements for array
        ERDEScorrespondingERreacs.append(ERDEScorrespondingERreacs[0])
    openFile.write(
        array_to_string(
            "\tERDEScorrespondingERreacs",
            ERDEScorrespondingERreacs,
            type="int",
            parameter=True,
        )
    )
    openFile.write("END MODULE network")
    openFile.close()


def find_reactant(species_list: list[str], reactant: str) -> int:
    """Try to find a reactant in the species list

    Args:
        species_list (list[str]): A list of species in the network
        reactant (str): The reactant to be indexed

    Returns:
        int: The index of the reactant, if it is not found, 9999
    """
    try:
        return species_list.index(reactant) + 1
    except ValueError:
        return 9999


def get_desorption_freeze_partners(reaction_list: list[Reaction]) -> list[Reaction]:
    """Every desorption has a corresponding freeze out eg desorption of #CO and freeze of CO.
    This find the corresponding freeze out for every desorb so that when desorb>>freeze
    we can turn off freeze out in UCLCHEM.

    Args:
        reaction_list (list): Reactions in network

    Returns:
        list: list of indices of freeze out reactions matching order of desorptions.
    """
    freeze_species = [
        x.get_products()[0] for x in reaction_list if x.get_reactants()[1] == "DESCR"
    ]
    partners = []
    for spec in freeze_species:
        for i, reaction in enumerate(reaction_list):
            if reaction.get_reaction_type() == "FREEZE":
                if reaction.get_reactants()[0] == spec:
                    partners.append(i + 1)
                    break
    return partners


def array_to_string(
    name: str, array: np.array, type: str = "int", parameter: bool = True
) -> str:
    """Write an array to fortran source code

    Args:
        name (str): Variable name of array in Fortran
        array (iterable): List of values of array
        type (str, optional): The array's type. Must be one of "int","float", or "string".Defaults to "int".
        parameter (bool, optional): Whether the array is a Fortran PARAMETER (constant). Defaults to True.

    Raises:
        ValueError: Raises an error if type isn't "int","float", or "string"

    Returns:
        str: String containing the Fortran code to declare this array.
    """
    # Check for 2D array
    arr = np.array(array)
    if arr.ndim == 2:
        shape = arr.shape
        flat = arr.flatten(order="F")
        if type == "int":
            dtype = "INTEGER(dp)"
            values = ",".join(str(int(v)) for v in flat)
        elif type == "float":
            dtype = "REAL(dp)"
            values = ",".join("{0:.4e}".format(float(v)) for v in flat)
        elif type == "string":
            strLength = len(max(flat, key=len))
            dtype = f"CHARACTER(Len={strLength})"
            values = ",".join('"' + str(v).ljust(strLength) + '"' for v in flat)
        elif type == "logical":
            dtype = "LOGICAL(dp)"
            values = ",".join(f".{str(v).upper()}." for v in flat)
        else:
            raise ValueError("Not a valid type for array to string")
        param_str = ", PARAMETER" if parameter else ""
        outString = f"{dtype}{param_str} :: {name}({','.join(str(s) for s in shape)}) = RESHAPE((/ {values} /), (/ {', '.join(str(s) for s in shape)} /))\n"
        outString = truncate_line(outString)
        return outString
    else:
        if parameter:
            outString = ", PARAMETER :: " + name + " ({0})=(/".format(len(arr))
        else:
            outString = " :: " + name + " ({0})=(/".format(len(arr))
        if type == "int":
            outString = "INTEGER(dp)" + outString
            for value in arr:
                outString += "{0},".format(value)
        elif type == "float":
            outString = "REAL(dp)" + outString
            for value in arr:
                outString += "{0:.4e},".format(value)
        elif type == "string":
            strLength = len(max(arr, key=len))
            outString = "CHARACTER(Len={0:.0f})".format(strLength) + outString
            for value in arr:
                outString += '"' + value.ljust(strLength) + '",'
        elif type == "logical":
            outString = "LOGICAL(dp)" + outString
            for value in arr:
                outString += ".{0}.,".format(value)
        else:
            raise ValueError("Not a valid type for array to string")
        outString = outString[:-1] + "/)\n"
        outString = truncate_line(outString)
        return outString
