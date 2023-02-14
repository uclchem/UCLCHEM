##########################################################################################
"""
Functions to read in the species and reaction files and write output files
"""
##########################################################################################

import csv
import logging
import numpy as np

from uclchem.makerates.network import Network
from .species import Species, elementList
from .reaction import Reaction, reaction_types
from os.path import join
from datetime import datetime
import yaml


def read_species_file(file_name) -> list[Species]:
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
        for row in reader:
            if row[0] != "NAME" and "!" not in row[0]:
                if "@" in row[0]:
                    user_defined_bulk.append(Species(row))
                else:
                    species_list.append(Species(row))
    return species_list, user_defined_bulk


def read_reaction_file(
    file_name: str, species_list: list[Species], ftype: str
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
        keep_list.append(species.name)

    if ftype == "UMIST":
        with open(file_name, "r") as f:
            reader = csv.reader(f, delimiter=":", quotechar="|")
            for row in reader:
                reaction_row = row[2:4] + [""] + row[4:8] + row[9:]
                if check_reaction(reaction_row, keep_list):
                    reactions.append(Reaction(reaction_row))
    elif ftype == "UCL":
        with open(file_name, "r") as f:
            reader = csv.reader(f, delimiter=",", quotechar="|")
            for row in reader:
                if (len(row) > 1) and (row[0][0] != "!"):
                    if check_reaction(row, keep_list):
                        reactions.append(Reaction(row))
                    else:
                        dropped_reactions.append(row)

    elif ftype == "KIDA":
        for row in kida_parser(file_name):
            if check_reaction(row, keep_list):
                reactions.append(Reaction(row))

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
    str_parse = lambda x: str(x).strip().upper()

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


def output_drops(dropped_reactions: list[Reaction], output_dir: str, write_files=True):
    """Writes the reactions that are dropped to disk/logs

    Args:
        dropped_reactions (list[Reaction]): The reactions that were dropped
        output_dir (str): The directory that dropped_reactions.csv will be written to.
        write_files (bool, optional): Whether or not to write the file. Defaults to True.
    """
    if output_dir is None:
        output_dir = ""
    outputFile = join(output_dir, "dropped_reactions.csv")
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


def write_outputs(network: Network, output_dir: str = None) -> None:
    """Write the ODE and Network fortran source files to the fortran source.

    Args:
        network (network): The makerates Network class
        output_dir (bool): The directory to write to.
    """
    if output_dir is None:
        output_dir = "../src/"
        fortran_src_dir = "../src/fortran_src"
    else:
        fortran_src_dir = output_dir

    # Create the species file
    filename = join(output_dir, "species.csv")
    write_species(filename, network.get_species_list())

    filename = join(output_dir, "reactions.csv")
    write_reactions(filename, network.get_reaction_list())

    # Write the ODEs in the appropriate language format
    filename = join(fortran_src_dir, "odes.f90")
    write_odes_f90(
        filename,
        network.get_species_list(),
        network.get_reaction_list(),
        network.three_phase,
    )

    filename = join(fortran_src_dir, "network.f90")
    write_network_file(filename, network)

    # Write some meta information that can be used to read back in the reactions into Python


def write_species(file_name: str, species_list: list[Species]) -> None:
    """Write the human readable species file. Note UCLCHEM doesn't use this file.

    Args:
        fileName (str): path to output file
        species_list (list): List of species objects for network
    """
    species_columns = [
        "NAME",
        "MASS",
        "BINDING ENERGY",
        "SOLID FRACTION",
        "MONO FRACTION",
        "VOLCANO FRACTION",
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
                    species.name,
                    species.mass,
                    species.binding_energy,
                    species.solidFraction,
                    species.monoFraction,
                    species.volcFraction,
                    species.enthalpy,
                ]
            )


# Write the reaction file in the desired format
def write_reactions(fileName, reaction_list) -> None:
    """Write the human readable reaction file. Note UCLCHEM doesn't use this file.

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
                ]
            )


def write_odes_f90(file_name, species_list, reaction_list, three_phase) -> None:
    """Write the ODEs in Modern Fortran. This is an actual code file.

    Args:
        file_name (str): Path to file where code will be written
        species_list (list): List of species describing network
        reaction_list (list): List of reactions describing network
        three_phase (bool): Flag for whether this is a 3 phase network
    """
    # First generate ODE contributions for all reactions
    species_names = [spec.name for spec in species_list]

    [
        logging.debug(f"{species_names.index(specie)+1}:{specie}")
        for specie in species_list
    ]

    for i, reaction in enumerate(reaction_list):
        logging.debug(f"RATE({i+1}):{reaction}")
        reaction.generate_ode_bit(i, species_names, three_phase)

    # then create ODE code and write to file.
    with open(file_name, mode="w") as output:
        # go through every species and build two strings, one with eq for all destruction routes and one for all formation
        ydotString = build_ode_string(species_list, reaction_list, three_phase)
        output.write(ydotString)


def write_jacobian(file_name, species_list) -> None:
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
        species_names += species.name
        losses = species.losses.split("+")
        gains = species.gains.split("+")
        for j in range(1, len(species_list) + 1):
            if species.name == "SURFACE":
                di_dj = f"J({i+1},{j})=SUM(J(surfaceList,{j}))\n"
                output.write(di_dj)
            elif species.name == "BULK":
                if species_names.count("@") > 0:
                    di_dj = f"J({i+1},{j})=SUM(J(bulkList,{j}))\n"
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
                if species_list[j - 1].name == "SURFACE":
                    di_dj = [f"+{x}/safeMantle" for x in losses if f"/safeMantle" in x]
                    di_dj += [f"-{x}/safeMantle" for x in gains if f"/safeMantle" in x]
                if len(di_dj) > 0:
                    di_dj = f"J({i+1},{j})=" + "".join(di_dj) + "\n"
                    output.write(di_dj)

        # tackle density separately.handled
        j = j + 1
        if species.name == "SURFACE":
            di_dj = f"J({i+1},{j})=SUM(J(surfaceList,{j}))\n"
            output.write(di_dj)
        elif species.name == "BULK":
            if species_names.count("@") > 0:
                di_dj = f"J({i+1},{j})=SUM(J(bulkList,{j}))\n"
                output.write(di_dj)
        else:
            di_dj = [f"-{x}".replace(f"*D", "", 1) for x in losses if f"*D" in x]
            di_dj += [f"+{x}".replace(f"*D", "", 1) for x in gains if f"*D" in x]
            di_dj = [x + "*2" if f"*D" in x else x for x in di_dj]
            if len(di_dj) > 0:
                di_dj = f"J({i+1},{j})=" + ("".join(di_dj)) + "\n"
                output.write(di_dj)
    i = i + 2
    di_dj = f"J({i},{i})=ddensdensdot(D)\n"
    output.write(di_dj)

    output.close()


def build_ode_string(species_list, reaction_list, three_phase) -> str:
    """A long, complex function that does the messy work of creating the actual ODE
    code to calculate the rate of change of each species. Test any change to this code
    thoroughly because ODE mistakes are very hard to spot.

    Args:
        species_list (list): List of species in network
        reaction_list (list): List of reactions in network
        three_phase (bool): Bool denoting if this is a three phase network

    Returns:
        str: One long string containing the entire ODE fortran code.
    """

    # We create a string of losses and gains for each species so initialize them all as ""
    species_names = []
    for i, species in enumerate(species_list):
        species_names.append(species.name)
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
                if ("ER" in reaction.get_reactants()) and (
                    not species_list[species_names.index(species)].is_surface_species()
                ):
                    species_list[
                        species_names.index("#" + species)
                    ].losses += reaction.ode_bit
                else:
                    species_list[
                        species_names.index(species)
                    ].losses += reaction.ode_bit
                if reaction.get_reactants()[1] == "BULKSWAP":
                    total_swap += reaction.ode_bit
        for species in reaction.get_products():
            if species in species_names:
                species_list[species_names.index(species)].gains += reaction.ode_bit

    ode_string = ""
    if three_phase:
        ode_string += truncate_line(f"totalSwap={total_swap[1:]}\n\n")
    # First get total rate of change of bulk and surface by adding ydots
    for n, species in enumerate(species_list):
        if species.name[0] == "@":
            species_list[bulk_index].gains += f"+YDOT({n+1})"
        elif species.name[0] == "#":
            species_list[surface_index].gains += f"+YDOT({n+1})"

    for n, species in enumerate(species_list):
        ydot_string = species_ode_string(n, species)
        ode_string += ydot_string

    # now add bulk transfer to rate of change of surface species after they've already been calculated
    if three_phase:
        ode_string += "!Update surface species for bulk growth, replace surfaceCoverage with alpha_des\n"
        ode_string += (
            "!Since ydot(surface_index) is negative, bulk is lost and surface forms\n"
        )

        ode_string += f"IF (YDOT({surface_index+1}) .lt. 0) THEN\n    surfaceCoverage = MIN(1.0,safeBulk/safeMantle)\n"

        for n, species in enumerate(species_list):
            if species.name[0] == "#":
                bulk_partner = species_names.index(species.name.replace("#", "@"))
                if not species_list[bulk_partner].is_refractory:
                    ode_string += f"    YDOT({n+1})=YDOT({n+1})-YDOT({surface_index+1})*surfaceCoverage*Y({bulk_partner+1})/safeBulk\n"
            if species.name[0] == "@":
                if not species.is_refractory:
                    ode_string += f"    YDOT({n+1})=YDOT({n+1})+YDOT({surface_index+1})*surfaceCoverage*Y({n+1})/safeBulk\n"
        ode_string += "ELSE\n"
        for n, species in enumerate(species_list):
            if species.name[0] == "@":
                surface_version = species_names.index(species.name.replace("@", "#"))
                ode_string += f"    YDOT({n+1})=YDOT({n+1})+YDOT({surface_index+1})*surfaceCoverage*Y({surface_version+1})\n"
            if species.name[0] == "#":
                ode_string += f"    YDOT({n+1})=YDOT({n+1})-YDOT({surface_index+1})*surfaceCoverage*Y({n+1})\n"
        ode_string += "ENDIF\n"

        # once bulk transfer has been added, odes for bulk and surface must be updated to account for it
        ode_string += (
            "!Update total rate of change of bulk and surface for bulk growth\n"
        )
        ode_string += species_ode_string(bulk_index, species_list[bulk_index])
        ode_string += species_ode_string(surface_index, species_list[surface_index])

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
        ydot_string += f"    YDOT({n+1}) = "
        # start with empty string and add production and loss terms if they exists
        if species.gains != "":
            ydot_string += "PROD"
        if species.losses != "":
            ydot_string += "-LOSS"
        ydot_string += "\n"
    else:
        ydot_string = f"    YDOT({n+1}) = {0.0}\n"
    ydot_string = truncate_line(ydot_string)
    return ydot_string


def write_evap_lists(network_file, species_list: list[Species]) -> None:
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
    bulkList = []
    iceList = []
    refractoryList = []
    species_names = [spec.name for spec in species_list]
    for i, species in enumerate(species_list):
        if species.name[0] == "#":
            # find gas phase version of grain species. For #CO it looks for first species in list with just CO and then finds the index of that
            try:
                j = species_names.index(species.get_desorb_products()[0])
            except:
                error = f"{species.name} desorbs as {species.get_desorb_products()[0]}"
                error += "which is not in species list. This desorption is likely user defined.\n"
                error += "Please amend the desorption route in your reaction file and re-run Makerates"
                raise NameError(error)

            # plus ones as fortran and python label arrays differently
            surfacelist.append(i + 1)
            gasIceList.append(j + 1)
            solidList.append(species.solidFraction)
            monoList.append(species.monoFraction)
            volcList.append(species.volcFraction)
            iceList.append(i + 1)
            binding_energyList.append(species.binding_energy)
            enthalpyList.append(species.enthalpy)
        elif species.name[0] == "@":
            j = species_names.index(species.get_desorb_products()[0])
            gasIceList.append(j + 1)
            bulkList.append(i + 1)
            iceList.append(i + 1)
            binding_energyList.append(species.binding_energy)
            enthalpyList.append(species.enthalpy)
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
    network_file.write(array_to_string("refractoryList", refractoryList, type="int"))


def truncate_line(input_string, lineLength=72) -> str:
    """Take a string and adds line endings at regular intervals
    keeps us from overshooting fortran's line limits and, frankly,
    makes for nicer ode.f90 even if human readability isn't very impotant

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


def write_network_file(file_name, network):
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
    openFile.write(
        "    INTEGER, PARAMETER :: nSpec={0}, nReac={1}\n".format(
            len(species_list), len(reaction_list)
        )
    )

    # write arrays of all species stuff
    names = []
    atoms = []
    masses = []
    for species in species_list:
        names.append(species.name)
        masses.append(float(species.mass))
        atoms.append(species.n_atoms)

    speciesIndices = ""
    for name, species_index in network.species_indices.items():
        speciesIndices += "{0}={1},".format(name, species_index)
    if len(speciesIndices) > 72:
        speciesIndices = truncate_line(speciesIndices)
    speciesIndices = speciesIndices[:-1] + "\n"
    openFile.write("    INTEGER, PARAMETER ::" + speciesIndices)
    if network.three_phase:
        openFile.write("    LOGICAL, PARAMETER :: THREE_PHASE = .TRUE.\n")
    else:
        openFile.write("    LOGICAL, PARAMETER :: THREE_PHASE = .FALSE.\n")
    openFile.write(array_to_string("    specname", names, type="string"))
    openFile.write(array_to_string("    mass", masses, type="float"))
    openFile.write(array_to_string("    atomCounts", atoms, type="int"))

    # then write evaporation stuff
    write_evap_lists(openFile, species_list)

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
    duplicates = []
    tmins = []
    tmaxs = []
    # store important reactions
    reaction_indices = ""
    for reaction, index in network.important_reactions.items():
        reaction_indices += reaction + f"={index},"
    reaction_indices = truncate_line(reaction_indices[:-1]) + "\n"
    openFile.write("    INTEGER, PARAMETER ::" + reaction_indices)

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
        reacTypes.append(reaction.get_reaction_type())
    # if len(duplicates) == 0:
    #     duplicates = [9999]
    #     tmaxs = [0]
    #     tmins = [0]

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

    reacTypes = np.asarray(reacTypes)

    partners = get_desorption_freeze_partners(reaction_list)
    openFile.write(
        array_to_string("\tfreezePartners", partners, type="int", parameter=True)
    )

    for reaction_type in reaction_types + ["TWOBODY"]:
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
    except:
        return 9999


def get_desorption_freeze_partners(reaction_list) -> list[Reaction]:
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


def array_to_string(name, array, type="int", parameter=True) -> str:
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
    if parameter:
        outString = ", PARAMETER :: " + name + " ({0})=(/".format(len(array))
    else:
        outString = " :: " + name + " ({0})=(/".format(len(array))
    if type == "int":
        outString = "INTEGER" + outString
        for value in array:
            outString += "{0},".format(value)
    elif type == "float":
        outString = "REAL(dp)" + outString
        for value in array:
            outString += "{0:.4e},".format(value)
    elif type == "string":
        strLength = len(max(array, key=len))
        outString = "CHARACTER(Len={0:.0f})".format(strLength) + outString
        for value in array:
            outString += '"' + value.ljust(strLength) + '",'
    else:
        raise ValueError("Not a valid type for array to string")
    outString = outString[:-1] + "/)\n"
    outString = truncate_line(outString)
    return outString
