"""UCLCHEM MakeRates."""

import logging
from collections.abc import Sequence
from pathlib import Path

from uclchem.makerates import io_functions as io
from uclchem.makerates._output_resolver import resolve_output_dirs
from uclchem.makerates.config import MakeratesConfig, ReactionFileTypes
from uclchem.makerates.network import Network

logger = logging.getLogger(__name__)

# Optional parameters that don't raise errors if missing
optional_params = [
    "grain_assisted_recombination_file",
    "output_directory",
    "three_phase",
    "gas_phase_extrapolation",
]


def run_makerates(
    configuration: str | Path | MakeratesConfig = "user_settings.yaml",
    write_files: bool = True,
    output_directory: str | Path | None = None,
) -> Network:
    """Run makerates.

    Main run wrapper for makerates. Loads and validates configuration,
    generates chemical network, and optionally writes output files.

    Args:
        configuration (str | Path | MakeratesConfig): Path to YAML configuration file,
            or ``MakeratesConfig`` instance. Defaults to "user_settings.yaml".
        write_files (bool): Whether to write fortran files to src/fortran_src.
            Defaults to True.
        output_directory (str | Path | None): Optional override for the output directory
            where files should be written. If None, uses the 'output_directory'
            from the config (if present) or the package defaults. Default = None.

    Returns:
        network (Network): A validated chemical network instance.

    Raises:
        ValueError: If ``coolants_file`` is a directory, and not a path to a file.

    """
    if not isinstance(configuration, MakeratesConfig):
        # Load and validate configuration using Pydantic
        config = MakeratesConfig.from_yaml(configuration)
    else:
        config = configuration

    # Log the configuration
    config.log_configuration()

    # Resolve output directories using tiered priority:
    # 1) explicit kwarg, 2) config field, 3) stored project root, 4) legacy relative paths
    explicit_dir = output_directory or config.output_directory
    if explicit_dir:
        # Ensure the directory exists
        explicit_dir = Path(explicit_dir)
        if not explicit_dir.is_dir():
            explicit_dir.mkdir(parents=True)

    output_dir, fortran_src_dir = resolve_output_dirs(
        explicit_dir,
        use_legacy_relative=True,
    )

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

    # Determine which coolants to write. Precedence (highest -> lowest):
    # 1) inline `coolants` in config, 2) `coolants_file` referenced in config,
    # 3) defaults used by write_outputs.
    coolants_to_write = None
    if config.coolants is not None:
        logger.info(f"Using {len(config.coolants)} inline coolants from config")
        coolants_to_write = config.coolants
    elif config.coolants_file:
        coolants_path = config.resolve_path(config.coolants_file)
        # Defensive check: don't try to read a directory as a YAML file
        if coolants_path.is_dir():
            msg = (
                f"coolants_file {coolants_path} resolves to a directory; expected a YAML file listing coolants. "
                "If you intended to set the collisional rate data directory, use 'coolant_data_dir' in your config."
            )
            raise ValueError(msg)
        try:
            _coolants = io.read_coolants_file(coolants_path)
            logger.info(f"Loaded {len(_coolants)} coolants from {coolants_path}")
            coolants_to_write = _coolants
        except Exception as exc:
            msg = f"Error reading coolants_file {coolants_path}"
            raise ValueError(msg) from exc

    if write_files:
        logger.info(
            "\n################################################\n"
            + "Checks complete, writing output files\n"
            + "################################################\n"
        )
        # Write dropped reactions
        io.output_drops(
            dropped_reactions=dropped_reactions,
            output_dir=output_dir,
            write_files=write_files,
        )
        logger.info(f"There are {len(dropped_reactions)} dropped reactions")

        # Check for GAR reactions and validate parameters
        gar_reactions = network.get_reactions_by_types("GAR")
        gar_parameters = None
        if len(gar_reactions) > 0:
            if gar_file is None:
                msg = (
                    "You have GAR reactions in your network, but you did "
                    "not specify a grain_assisted_recombination_file in "
                    "your configuration. Refer to makerates documentation."
                )
                raise ValueError(msg)
            # Get all the individual ions that can recombine
            gar_ions = [gar.get_reactants()[0] for gar in gar_reactions]
            _gar_parameters = io.read_grain_assisted_recombination_file(gar_file)
            if not set(gar_ions).issubset(set(_gar_parameters.keys())):
                missing_ions = set(gar_ions) - set(_gar_parameters.keys())
                msg = (
                    f"You have GAR reactions for ions {missing_ions} but "
                    f"they are not defined in your gar_file {gar_file}"
                )
                raise ValueError(msg)
            # Save the gar parameters in the correct order
            gar_parameters = {ion: _gar_parameters[ion] for ion in gar_ions}

        # Pass resolved output directories and other parameters to write_outputs
        io.write_outputs(
            network,
            output_dir,
            fortran_src_dir,
            enable_rates_storage=config.enable_rates_storage,
            gar_database=gar_parameters,
            coolants=coolants_to_write,
            coolant_data_dir=config.coolant_data_dir,
        )

        # Copy coolant data files to package data directory for installation
        # Only pass coolant_data_dir if it's explicitly set and valid
        source_dir = (
            config.coolant_data_dir
            if config.coolant_data_dir is not None and str(config.coolant_data_dir) != "."
            else None
        )
        io.copy_coolant_files(source_dir=source_dir)

    ngrain = len([x for x in network.get_species_list() if x.is_surface_species()])
    logger.info(f"Total number of species = {len(network.get_species_list())}")
    logger.info(f"Number of surface species = {ngrain}")
    logger.info(f"Number of reactions = {len(network.get_reaction_list())}")

    # Return the network for reuse in code/notebooks
    return network


def get_network(
    path_to_input_file: str | Path | None = None,
    path_to_species_file: str | Path | None = None,
    path_to_reaction_file: str | Path | None = None,
    verbosity: str | None = None,
) -> Network:
    """Get a network into memory.

    In memory equivalent of Makerates, can either be used on the original input files
    for makerates, or on the output files that makerates generates. So either specify:

    `path_to_input_file ` exclusive OR (`path_to_species_file` and `path_to_reaction_file`)

    The latter scenario allows you to reload a reaction network from
    a network already written by Makerates.

    Args:
        path_to_input_file (str | Path | None): Path to input file. Defaults to None.
        path_to_species_file (str | Path | None): Path to a species.csv
            in/from the src directory. Defaults to None.
        path_to_reaction_file (str | Path | None): Path to a reactions.csv in/from
            the src directory. Defaults to None.
        verbosity (str | None): The verbosity level as specified in logger.
            Defaults to None.

    Returns:
        Network: A chemical reaction network.

    Raises:
        ValueError: You cannot specify both an input configuration and species+reaction.


    """
    if verbosity is not None:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=verbosity.upper())

    if bool(path_to_input_file) and bool(path_to_species_file or path_to_reaction_file):
        msg = "Cannot have both an input Makerates config file and explicit paths to species + reaction files"
        raise ValueError(msg)

    if path_to_input_file:
        return run_makerates(path_to_input_file, write_files=False)
    else:
        # If we load the species/reactions directly from UCLCHEM we can skip the checks
        return Network.from_csv(path_to_species_file, path_to_reaction_file)


def _get_network_from_files(
    species_file: str | Path,
    reaction_files: str | Path | Sequence[str | Path],
    reaction_types: ReactionFileTypes | list[ReactionFileTypes],
    gas_phase_extrapolation: bool,
    add_crp_photo_to_grain: bool,
    derive_reaction_exothermicity: bool | str | list[str],
    database_reaction_exothermicity: Sequence[str | Path] | None = None,
) -> tuple[Network, list[list]]:
    """Get a network from files.

    Args:
        species_file (str | Path): path to a species file
        reaction_files (str | Path | Sequence[str | Path]): path to reaction files
        reaction_types (ReactionFileTypes | list[ReactionFileTypes]): list of reaction
            file types corresponding to each file in ``reaction_files``.
        gas_phase_extrapolation (bool): Extrapolate gas-phase temperature.
        add_crp_photo_to_grain (bool): Add CRP/PHOTON to grain.
        derive_reaction_exothermicity (bool | str | list[str]): Reaction types to calculate
            exothermicity for.
        database_reaction_exothermicity (Sequence[str | Path] | None): Custom exothermicity database
            files. Default = None.

    Returns:
        network (Network): Instantiated Network.
        dropped_reactions (list[list]): list of dropped reactions.

    """
    logger.info(
        f"_get_network_from_files called with database_reaction_exothermicity={database_reaction_exothermicity}"
    )
    species_list, user_defined_bulk = io.read_species_file(species_file)
    # Check if reaction and type files are lists, if not, make them lists
    if not isinstance(reaction_files, list):
        reaction_files = [reaction_files]  # type: ignore
    if not isinstance(reaction_types, list):
        reaction_types = [reaction_types]
    reactions = []
    dropped_reactions = []
    # Support an arbitrary amount of different reaction files and append then in the end.
    for reaction_file, reaction_type in zip(reaction_files, reaction_types, strict=True):
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

    logger.info(
        "\n################################################\n"
        + "Reading and checking input\n"
        + "################################################\n"
    )

    # Network checking is now done automatically during build
    # in NetworkBuilder._check_network()
    return network, dropped_reactions
