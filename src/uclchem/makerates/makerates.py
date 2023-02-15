import uclchem.makerates.io_functions as io
from uclchem.makerates import Network
import os
import yaml
import logging
from logging import Logger
from os import PathLike
from typing import Union

from uclchem.makerates.network import LoadedNetwork


param_list = [
    "species_file",
    "database_reaction_file",
    "database_reaction_type",
    "custom_reaction_file",
    "custom_reaction_type",
    "three_phase",
]


def run_makerates(
    configuration_file: str = "user_settings.yaml", write_files: bool = True
) -> Network:
    """The main run wrapper for makerates, it loads a configuration, parses it in Network
    and then returns the Network. It by default writes to the uclchem fortran directory, but
    this can be skipped.

    Args:
        configuration_file (str, optional): A UCLCHEM Makerates configuration file. Defaults to "user_settings.yaml".
        write_files (bool, optional): Whether to write the fortran files to the src/fortran_src. Defaults to True.

    Raises:
        KeyError: The configuration cannot be found

    Returns:
        Network: A chemical network instance.
    """

    with open(configuration_file, "r") as f:
        user_params = yaml.safe_load(f)

    for param in param_list:
        try:
            # Check if we can access the needed parameter:
            user_params[param]
            logging.info(f"{param} : {user_params[param]}")
        except:
            raise KeyError(f"{param} not found in user_settings.yaml")

    try:
        user_output_dir = user_params["output_directory"]
        if not os.path.exists(user_output_dir):
            os.makedirs(user_output_dir)
    except:
        user_output_dir = None

    # load everything from the configuration file
    reaction_files = [
        user_params["database_reaction_file"],
        user_params["custom_reaction_file"],
    ]
    reaction_types = [
        user_params["database_reaction_type"],
        user_params["custom_reaction_type"],
    ]
    species_file = user_params["species_file"]
    three_phase = user_params["three_phase"]

    # retrieve the network and the dropped reactions
    network, dropped_reactions = _get_network_from_files(
        reaction_files=reaction_files,
        reaction_types=reaction_types,
        species_file=species_file,
        three_phase=three_phase,
    )

    if write_files:
        logging.info(
            "\n################################################\n"
            + "Checks complete, writing output files\n"
            + "################################################\n"
        )
        # Write or output the written files
        io.output_drops(
            dropped_reactions=dropped_reactions,
            output_dir=user_output_dir,
            write_files=write_files,
        )
        logging.info(f"There are {len(dropped_reactions)} droppped reactions")
        io.write_outputs(network, user_output_dir)

    ngrain = len([x for x in network.get_species_list() if x.is_surface_species()])
    logging.info(f"Total number of species = {len(network.get_species_list())}")
    logging.info(f"Number of surface species = {ngrain}")
    logging.info(f"Number of reactions = {len(network.get_reaction_list())}")
    # return the network such that the object can be reused in code/notebooks
    return network


def get_network(
    path_to_input_file: Union[str, bytes, os.PathLike] = None,
    path_to_species_file: Union[str, bytes, os.PathLike] = None,
    path_to_reaction_file: Union[str, bytes, os.PathLike] = None,
    verbosity=None,
):
    """In memory equivalent of Makerates, can either be used on the original input files
    for makerates, or on the output files that makerates generates. So either specify:

    `path_to_input_file ` exclusive OR (`path_to_species_file` and `path_to_reaction_file`)

    The latter scenario allows you to reload a reaction network from a network already written by Makerates.


    Args:
        path_to_input_file (Union[str, bytes, os.PathLike], optional): Path to input file. Defaults to None.
        path_to_species_file (Union[str, bytes, os.PathLike], optional): Path to a species.csv in/from the src directory. Defaults to None.
        path_to_reaction_file (Union[str, bytes, os.PathLike], optional): Path to a reactions.csv in/from the src directory. Defaults to None.
        verbosity (LEVEL, optional): The verbosity level as specified in logging. Defaults to None.

    Raises:
        ValueError: You cannot specify both an input configuration and species+reaction.

    Returns:
        Network: A chemical reaction network.
    """
    if verbosity:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=verbosity)

    if bool(path_to_input_file) and bool(path_to_species_file or path_to_reaction_file):
        raise ValueError(
            "Cannot have both an input Makerates config file and explicit paths to species + reaction files"
        )

    if path_to_input_file:
        return run_makerates(path_to_input_file, write_files=False)
    else:
        # If we load the species/reactions directly from UCLCHEM we can skip the checks
        species_list, _ = io.read_species_file(path_to_species_file)
        reactions_list, _ = io.read_reaction_file(
            path_to_reaction_file, species_list, "UCL"
        )
        return LoadedNetwork(species_list, reactions_list)


def _get_network_from_files(
    species_file: Union[str, bytes, os.PathLike],
    reaction_files: list[Union[str, bytes, os.PathLike]],
    reaction_types: list[str],
    three_phase: bool,
):
    species_list, user_defined_bulk = io.read_species_file(species_file)
    # Check if reaction and type files are lists, if not, make them lists
    if not isinstance(reaction_files, list):
        reaction_files = [reaction_files]
    if not isinstance(reaction_types, list):
        reaction_types = [reaction_types]
    reactions = []
    dropped_reactions = []
    # Support an arbitrary amount of different reaction files and append then in the end.
    for reaction_file, reaction_type in zip(reaction_files, reaction_types):
        temp_reactions, temp_dropped_reactions = io.read_reaction_file(
            reaction_file, species_list, reaction_type
        )
        reactions += temp_reactions
        dropped_reactions += temp_dropped_reactions
    # Create Network
    network = Network(
        species=species_list,
        reactions=reactions,
        three_phase=three_phase,
        user_defined_bulk=user_defined_bulk,
    )

    #################################################################################################

    logging.info(
        "\n################################################\n"
        + "Reading and checking input\n"
        + "################################################\n"
    )

    # check network to see if there are potential problems, in the get wrapper because checking should always happen!
    logging.info("Checking Network")
    network.check_network()
    return network, dropped_reactions
