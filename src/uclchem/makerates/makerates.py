import logging
import os
from typing import Union

from . import io_functions as io
from .config import MakeratesConfig
from .network import Network
from .reaction import Reaction

# Optional parameters that don't raise errors if missing
optional_params = [
    "grain_assisted_recombination_file",
    "output_directory",
    "three_phase",
    "gas_phase_extrapolation",
]


from typing import Dict, List, Optional, Union


def run_makerates(
    configuration_file: str = "user_settings.yaml",
    write_files: bool = True,
    output_directory: Union[str, os.PathLike] | None = None,
) -> Network:
    """
    Main run wrapper for makerates. Loads and validates configuration,
    generates chemical network, and optionally writes output files.

    Args:
        configuration_file: Path to YAML configuration file.
            Defaults to "user_settings.yaml".
        write_files: Whether to write fortran files to src/fortran_src.
            Defaults to True.
        output_directory: Optional override for the output directory where files
            should be written. If None, uses the 'output_directory' from the config
            (if present) or the package defaults.

    Returns:
        Network: A validated chemical network instance.

    Raises:
        ValidationError: If configuration is invalid
        FileNotFoundError: If required files are missing
    """
    # Load and validate configuration using Pydantic
    config = MakeratesConfig.from_yaml(configuration_file)

    # Log the configuration
    config.log_configuration()

    # Prepare output directory (allow caller to override)
    if output_directory is not None:
        user_output_dir = output_directory
        if not os.path.exists(user_output_dir):
            os.makedirs(user_output_dir)
    elif config.output_directory:
        user_output_dir = config.resolve_path(config.output_directory)
        if not os.path.exists(user_output_dir):
            os.makedirs(user_output_dir)
    else:
        user_output_dir = None

    # Get all reaction files and types
    reaction_files = config.get_all_reaction_files()
    reaction_types = config.get_all_reaction_types()

    # Resolve species file path
    species_file = config.resolve_path(config.species_file)

    # Resolve GAR file if present
    gar_file = None
    if config.grain_assisted_recombination_file:
        gar_file = config.resolve_path(config.grain_assisted_recombination_file)

    # Resolve exothermicity files if present
    database_reaction_exothermicity = None
    if config.database_reaction_exothermicity:
        database_reaction_exothermicity = [
            config.resolve_path(ef) for ef in config.database_reaction_exothermicity
        ]

    # Retrieve the network and the dropped reactions
    network, dropped_reactions = _get_network_from_files(
        reaction_files=reaction_files,
        reaction_types=reaction_types,
        species_file=species_file,
        gas_phase_extrapolation=config.gas_phase_extrapolation,
        add_crp_photo_to_grain=config.add_crp_photo_to_grain,
        derive_reaction_exothermicity=config.derive_reaction_exothermicity,
        database_reaction_exothermicity=database_reaction_exothermicity,
    )

    if write_files:
        logging.info(
            "\n################################################\n"
            + "Checks complete, writing output files\n"
            + "################################################\n"
        )
        # Write dropped reactions
        io.output_drops(
            dropped_reactions=dropped_reactions,
            output_dir=user_output_dir,
            write_files=write_files,
        )
        logging.info(f"There are {len(dropped_reactions)} dropped reactions")

        # Check for GAR reactions and validate parameters
        gar_reactions = network.get_reactions_by_types("GAR")
        gar_parameters = None
        if len(gar_reactions) > 0:
            if gar_file is None:
                raise ValueError(
                    "You have GAR reactions in your network, but you did "
                    "not specify a grain_assisted_recombination_file in "
                    "your configuration. Refer to makerates documentation."
                )
            # Get all the individual ions that can recombine
            gar_ions = [gar.get_reactants()[0] for gar in gar_reactions]
            _gar_parameters = io.read_grain_assisted_recombination_file(gar_file)
            if not set(gar_ions).issubset(set(_gar_parameters.keys())):
                missing_ions = set(gar_ions) - set(_gar_parameters.keys())
                raise ValueError(
                    f"You have GAR reactions for ions {missing_ions} but "
                    f"they are not defined in your gar_file {gar_file}"
                )
            # Save the gar parameters in the correct order
            gar_parameters = {ion: _gar_parameters[ion] for ion in gar_ions}
        # Determine which coolants to write. Precedence (highest -> lowest):
        # 1) inline `coolants` in config, 2) `coolants_file` referenced in config,
        # 3) defaults used by write_outputs.
        coolants_to_write = None
        if config.coolants is not None:
            logging.info(f"Using {len(config.coolants)} inline coolants from config")
            coolants_to_write = config.coolants
        elif config.coolants_file:
            coolants_path = config.resolve_path(config.coolants_file)
            try:
                _coolants = io.read_coolants_file(coolants_path)
                logging.info(f"Loaded {len(_coolants)} coolants from {coolants_path}")
                coolants_to_write = _coolants
            except Exception as exc:
                raise ValueError(f"Error reading coolants_file {coolants_path}: {exc}")

        # Pass resolved coolants and coolant_data_dir through to write_outputs
        io.write_outputs(
            network,
            user_output_dir,
            config.enable_rates_storage,
            gar_parameters,
            coolants=coolants_to_write,
            coolant_data_dir=config.coolant_data_dir,
        )

    ngrain = len([x for x in network.get_species_list() if x.is_surface_species()])
    logging.info(f"Total number of species = {len(network.get_species_list())}")
    logging.info(f"Number of surface species = {ngrain}")
    logging.info(f"Number of reactions = {len(network.get_reaction_list())}")

    # Return the network for reuse in code/notebooks
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
        return Network.from_csv(path_to_species_file, path_to_reaction_file)


def _get_network_from_files(
    species_file: Union[str, bytes, os.PathLike],
    reaction_files: list[Union[str, bytes, os.PathLike]],
    reaction_types: list[str],
    gas_phase_extrapolation: bool,
    add_crp_photo_to_grain: bool,
    derive_reaction_exothermicity: Union[bool, str, list[str]],
    database_reaction_exothermicity: list[Union[str, bytes, os.PathLike]] = None,
) -> tuple[Network, list[Reaction]]:
    logging.info(
        f"_get_network_from_files called with database_reaction_exothermicity={database_reaction_exothermicity}"
    )
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

    # Create Network using the build() factory method
    network = Network.build(
        species=species_list,
        reactions=reactions,
        user_defined_bulk=user_defined_bulk,
        gas_phase_extrapolation=gas_phase_extrapolation,
        add_crp_photo_to_grain=add_crp_photo_to_grain,
        derive_reaction_exothermicity=derive_reaction_exothermicity,
        database_reaction_exothermicity=database_reaction_exothermicity,
    )

    #################################################################################################

    logging.info(
        "\n################################################\n"
        + "Reading and checking input\n"
        + "################################################\n"
    )

    # Network checking is now done automatically during build in NetworkBuilder._check_network()
    return network, dropped_reactions
