##########################################################################################
"""Functions to read in the species and reaction files and write output files."""

import csv
import fileinput
import functools
import logging
import shutil
import warnings
from datetime import datetime
from pathlib import Path
from tempfile import mkstemp
from textwrap import dedent
from typing import IO, Any, cast

import numpy as np
import yaml

from uclchem.constants import PHYSICAL_PARAMETERS, ZETA_0
from uclchem.makerates.config import ReactionFileTypes
from uclchem.makerates.network import Network
from uclchem.makerates.reaction import (
    REACTION_TYPES,
    CoupledReaction,
    Reaction,
    reaction_header,
)
from uclchem.makerates.species import (
    Species,
    get_element_counts_per_species,
    species_header,
)
from uclchem.makerates.utils import (
    array_to_string,
    check_reaction,
    get_default_coolants,
    pad_to_max_length,
    replace_value_with_name,
    separate_common_terms,
    strip_comments_from_row,
    truncate_line,
)
from uclchem.utils import (
    MISSING_VALUE_FLOAT,
    MISSING_VALUE_INTEGER,
    NO_REACTANT_OR_PRODUCT,
    UCLCHEM_ROOT_DIR,
)

logger = logging.getLogger(__name__)


_safe_load = functools.partial(yaml.load, Loader=yaml.CSafeLoader)


def get_default_coolant_directory(user_specified: str | Path = "") -> str:
    """Get the collisional rates directory path for use during makerates.

    Parameters
    ----------
    user_specified : str | Path
        User-specified directory from config.
        If empty, searches standard locations relative to CWD. (Default value = '')

    Returns
    -------
    str
        Absolute directory path to collisional rate data files.

    """
    if user_specified:
        return str(user_specified)
    # Search standard locations (makerates runs from Makerates/ directory)
    candidates = [
        Path.cwd() / "data" / "collisional_rates",  # Makerates/data/collisional_rates
        Path.cwd().parent
        / "data"
        / "collisional_rates",  # project_root/data/collisional_rates
        Path.cwd() / "Makerates" / "data" / "collisional_rates",  # from project root
        Path.cwd().parent
        / "Makerates"
        / "data"
        / "collisional_rates",  # from subdirectory
    ]
    for candidate in candidates:
        # Check if the directory exists and contains coolant data files
        if candidate.is_dir():
            # Look for at least one expected file to confirm it's the right directory
            if list(candidate.glob("*.dat")):  # Has individual .dat files
                return str(candidate.resolve())
    return ""


def read_species_file(
    file_name: str | Path,
) -> tuple[list[Species], list[Species]]:
    """Read in a Makerates species file.

    Parameters
    ----------
    file_name : str | Path
        path to file containing the species list

    Returns
    -------
    species_list : list[Species]
        List of Species objects
    user_defined_bulk : list[Species]
        List of user defined bulk species

    Raises
    ------
    IndexError
        If there is an error parsing a line in the file.

    """
    species_list = []
    # list to hold user defined bulk species (for adjusting binding energy)
    user_defined_bulk = []
    with Path(file_name).open() as f:
        reader = csv.reader(f, delimiter=",", quotechar="|")
        for idx, row in enumerate(reader):
            try:
                if row[0] != "NAME" and "!" not in row[0]:
                    row = strip_comments_from_row(row)
                    if "@" in row[0]:
                        user_defined_bulk.append(Species(cast("list[str | float]", row)))
                    else:
                        species_list.append(Species(cast("list[str | float]", row)))
            except IndexError as exc:
                print(f"Error reading species file {file_name} at line {idx}")
                raise exc

    return species_list, user_defined_bulk


def read_reaction_file(
    file_name: str | Path,
    species_list: list[Species],
    ftype: ReactionFileTypes,
) -> tuple[list[Reaction], list[list[str]]]:
    """Read in a reaction file of any kind (UCL, UMIST, KIDA), and.

    produces a list of reactions for the network, filtered by species_list.

    Parameters
    ----------
    file_name : str | Path
        A file name for the reaction file to read.
    species_list : list[Species]
        A list of chemical species to be used in the reading.
    ftype : ReactionFileTypes
        'UMIST','UCL', or 'KIDA' to describe format of file_name

    Returns
    -------
    reactions : list[Reaction]
        List of kept reactions.
    dropped_reactions : list[list[str]]
        List of raw CSV rows for dropped reactions.

    Raises
    ------
    ValueError
        If reaction file type is not one of ["UMIST", "UCL", "KIDA"].

    """
    reactions = []
    dropped_reactions = []

    # Every reactant/product of a reaction must be in keep_list to not be dropped
    keep_list = ["", "NAN", "#", "*", "E-", "e-", "ELECTR", "@"]
    keep_list.extend(REACTION_TYPES)
    keep_list.extend(species.get_name() for species in species_list)

    if ftype == "UMIST":
        with Path(file_name).open() as f:
            reader = csv.reader(f, delimiter=":", quotechar="|")
            for row in reader:
                if row[0].startswith("#") or row[0].startswith("!"):
                    continue
                row = strip_comments_from_row(row)
                reaction_row = [*row[2:4], "", *row[4:8], *row[9:14], ""]
                if check_reaction(reaction_row, keep_list):
                    reactions.append(
                        Reaction(
                            cast("list[str | float]", reaction_row),
                            reaction_source="UMIST",
                        )
                    )
    elif ftype == "UCL":
        with Path(file_name).open() as f:
            reader = csv.reader(f, delimiter=",", quotechar="|")
            for row in reader:
                if (len(row) > 1) and (row[0][0] != "!"):
                    row = strip_comments_from_row(row)
                    if check_reaction(row, keep_list):
                        reactions.append(
                            Reaction(
                                cast("list[str | float]", row),
                                reaction_source="UCL",
                            )
                        )
                    else:
                        dropped_reactions.append(row)

    elif ftype == "KIDA":
        reactions.extend(
            Reaction(row, reaction_source="KIDA")
            for row in kida_parser(file_name)
            if check_reaction(row, keep_list)
        )

    else:
        msg = "Reaction file type must be one of 'UMIST', 'UCL' or 'KIDA'"
        raise ValueError(msg)
    return reactions, dropped_reactions


def kida_parser(kida_file: str | Path) -> list[list[str | int | float]]:
    """Parse rows of a KIDA file.

    KIDA used a fixed format file so we read each line in the chunks they specify
    and use python built in classes to convert to the necessary types.
    NOTE KIDA defines some of the same reaction types to UMIST but with different names
    and coefficients. We fix that by converting them here.

    Parameters
    ----------
    kida_file : str | Path
        path to KIDA file

    Returns
    -------
    list[list[str | int | float]]
        rows (list[list[Any]])

    """

    def str_parse(x: Any) -> str:
        return str(x).strip().upper()

    kida_contents: list[tuple[int, dict[Any, int]]] = [
        (3, {str_parse: 11}),
        (1, {"skip": 1}),
        (5, {str_parse: 11}),
        (1, {"skip": 1}),
        (3, {float: 10, "skip": 1}),
        (1, {"skip": 27}),
        (2, {int: 6, "skip": 1}),
        (1, {int: 2}),
        (1, {"skip": 11}),
    ]
    rows = []
    with Path(kida_file).open() as f:
        f.readline()  # throw away header
        for line in f:  # then iterate over file
            if line.startswith("!"):
                continue
            row: list[str | int | float] = []
            for item in kida_contents:
                for _ in range(item[0]):
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
                row[8] = cast("float", row[8]) * ZETA_0
                rows.append(row[:7] + row[8:-1])
            elif row[-1] in {2, 3}:
                rows.append(row[:7] + row[8:-1])
            elif row[-1] == 4:  # ruff: ignore[magic-value-comparison]
                row[2] = "IONOPOL1"
                rows.append(row[:7] + row[8:-1])
            elif row[-1] == 5:  # ruff: ignore[magic-value-comparison]
                row[2] = "IONOPOL2"
                rows.append(row[:7] + row[8:-1])
    return rows


def read_grain_assisted_recombination_file(
    file_name: str | Path,
) -> dict[str, np.ndarray]:
    """Read a grain assisted recombination file.

    Parameters
    ----------
    file_name : str | Path
        file path of grain assisted recombination data

    Returns
    -------
    gar_parameters : dict[str, np.ndarray]
        Database for grain-activated recombination
        reactions.

    """
    with Path(file_name).open() as fh:
        gar_parameters = _safe_load(fh)
    return gar_parameters


def read_coolants_file(file_name: str | Path) -> list[dict]:
    """Read a YAML file specifying coolants.

    The file should contain either a single mapping or a list of mappings where each
    mapping contains 'file' and 'name' keys. 'file' must be a bare filename (no path).

    Parameters
    ----------
    file_name : str | Path
        Path to coolants file.

    Returns
    -------
    list[dict]
        Normalized list of coolant dicts.

    Raises
    ------
    ValueError
        If the yaml-parsed data is not a dictionary, or list of dictionaries.
    ValueError
        If the "file" entries in coolants_file are not bare filenames.

    """
    with Path(file_name).open() as fh:
        data = _safe_load(fh)

    if data is None:
        return []
    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list):
        msg = "Coolants file must contain a mapping or list of mappings"
        raise ValueError(msg)

    normalized = []

    for item in data:
        if not isinstance(item, dict):
            msg = "Each coolant entry must be a mapping with 'file' and 'name' keys"
            raise ValueError(msg)
        if "file" not in item or "name" not in item:
            msg = "Each coolant mapping must contain 'file' and 'name' keys"
            raise ValueError(msg)
        file_val = str(item["file"])
        if Path(file_val).name != file_val or Path(file_val).parent != Path():
            msg = "Coolant 'file' entries in coolants_file must be bare filenames (no directories)"
            raise ValueError(msg)
        entry: dict[str, Any] = {"file": file_val, "name": str(item["name"])}
        if "parent_species" in item:
            entry["parent_species"] = str(item["parent_species"])
        if "conversion_factor" in item:
            entry["conversion_factor"] = float(item["conversion_factor"])
        normalized.append(entry)
    return normalized


def output_drops(
    dropped_reactions: list[list[str]],
    output_dir: str | Path | None = None,
    *,
    write_files: bool = True,
) -> None:
    """Write the reactions that are dropped to disk/logs.

    Parameters
    ----------
    dropped_reactions : list[list[str]]
        The reactions that were dropped
    output_dir : str | Path | None
        The directory that dropped_reactions.csv will be written to.
        Default = None (writes to current directory).
    write_files : bool
        Whether or not to write the file. Defaults to True.

    """
    if output_dir is None:
        output_dir = "."
    outputFile = Path(output_dir) / "dropped_reactions.csv"
    # Print dropped reactions from grain file or write if many
    if write_files and dropped_reactions:
        logger.info(f"\nReactions dropped from grain file written to {outputFile}\n")
        with Path(outputFile).open("w") as f:
            writer = csv.writer(
                f,
                delimiter=",",
                quotechar="|",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator="\n",
            )
            writer.writerows(dropped_reactions)
    else:
        logger.info("Reactions dropped from grain file:\n")
        for reaction in dropped_reactions:
            logger.info(reaction)


def write_outputs(
    network: Network,
    python_src_dir: Path,
    fortran_src_dir: Path,
    gar_database: dict[str, np.ndarray] | None = None,
    coolants: list[dict] | None = None,
    coolant_data_dir: str | Path | None = "",
    *,
    enable_rates_storage: bool = False,
) -> None:
    """Write the ODE and Network fortran source files to the fortran source.

    Parameters
    ----------
    network : Network
        The makerates Network class
    python_src_dir : Path
        Directory to write Python source files
        (species.csv, reactions.csv).
    fortran_src_dir : Path
        Directory to write Fortran source files
        (odes.f90, network.f90, f2py_constants.f90).
    gar_database : dict[str, np.ndarray] | None
        Database for grain-activated recombination
        reactions. Default = None.
    coolants : list[dict] | None
        List of coolants or None. If None,
        use default list of coolants. See `get_default_coolants().`
    coolant_data_dir : str | Path | None
        User-specified directory from config.
        If empty, searches standard locations relative to CWD. (Default value = '')
    enable_rates_storage : bool
        Enable storage of writing rates to files.
        Default = False.

    Raises
    ------
    ValueError
        If coolants entries do not have a key "file" or the file names are
        not bare file names.
    ValueError
        If coolants start with "o-" or "p-" for ortho and para, but no
        `conversion_factor` is given in the dictionary.
    ValueError
        If coolants have parent species specified that are not in the network.

    """
    # Use default coolants if none provided
    if coolants is None:
        coolants = get_default_coolants()

    # Validate that coolant 'file' entries are bare filenames (not paths)
    for c in coolants:
        f = c.get("file")
        if f is None:
            msg = "Each coolant dict must contain a 'file' key"
            raise ValueError(msg)
        if Path(f).name != f or Path(f).parent != Path():
            msg = (
                "Coolant file names must be bare filenames (no directories). "
                "Set the coolant directory at runtime via coolantDataDir."
            )
            raise ValueError(msg)

    # Compute energy level counts from coolant data files
    from uclchem._coolant_utils import (  # ruff: ignore[import-outside-top-level] heavy-extension
        get_energy_levels_info,
        validate_coolant_frequencies,
    )

    coolant_data_directory = get_default_coolant_directory(
        coolant_data_dir if coolant_data_dir is not None else ""
    )
    n_total_levels, n_se_stats_per_coolant = get_energy_levels_info(
        coolant_names=[c["name"] for c in coolants],
        coolant_files=[c["file"] for c in coolants],
        data_dir=coolant_data_directory,
    )

    # Validate frequency consistency and compute suggested tolerance
    freq_deviations = validate_coolant_frequencies(
        coolant_names=[c["name"] for c in coolants],
        coolant_files=[c["file"] for c in coolants],
        data_dir=coolant_data_directory,
    )
    DEFAULT_SUGGESTED_FREQ_TOL = 0.01
    if freq_deviations:
        max_deviation = max(freq_deviations.values())
        # Add 10% margin above the largest observed deviation
        suggested_freq_rel_tol = max_deviation * 1.1
        # Enforce a minimum tolerance of 0.01 (1%)
        suggested_freq_rel_tol = max(suggested_freq_rel_tol, DEFAULT_SUGGESTED_FREQ_TOL)

        # Log per-coolant frequency deviations (sorted largest first)
        sorted_devs = sorted(freq_deviations.items(), key=lambda x: -x[1])
        name_width = max(len(name) for name, _ in sorted_devs)
        logger.info("Coolant frequency deviations (|E_i-E_j|/h vs LAMDA):")
        logger.info(f"  {'Coolant':<{name_width}}  {'Deviation':>10}")
        logger.info(f"  {'-' * name_width}  {'-' * 10}")
        for name, dev in sorted_devs:
            marker = " <<<" if dev > DEFAULT_SUGGESTED_FREQ_TOL else ""
            logger.info(f"  {name:<{name_width}}  {dev * 100:9.4f}%{marker}")
        logger.info(
            f"  Auto-setting freq_rel_tol = {suggested_freq_rel_tol:.4f} "
            f"(max deviation {max_deviation * 100:.2f}% + 10% margin)"
        )
    else:
        suggested_freq_rel_tol = DEFAULT_SUGGESTED_FREQ_TOL
        max_deviation = 0.0

    # Compute coolant conversion arrays
    parent_names = []
    conversion_factors = []
    conversion_modes = []
    for c in coolants:
        name = c["name"]
        explicit_parent = c.get("parent_species")
        explicit_factor = c.get("conversion_factor")

        if explicit_parent is not None and explicit_factor is not None:
            # Fully specified: fixed factor mode
            parent_names.append(explicit_parent)
            conversion_factors.append(explicit_factor)
            conversion_modes.append(0)
        elif name == "p-H2":
            parent_names.append(explicit_parent or "H2")
            conversion_factors.append(
                explicit_factor if explicit_factor is not None else 0.0
            )
            conversion_modes.append(0 if explicit_factor is not None else 1)
        elif name == "o-H2":
            parent_names.append(explicit_parent or "H2")
            conversion_factors.append(
                explicit_factor if explicit_factor is not None else 0.0
            )
            conversion_modes.append(0 if explicit_factor is not None else 2)
        elif name.startswith(("o-", "p-")):
            # Other ortho/para species — require explicit conversion_factor
            if explicit_factor is None:
                msg = (
                    f"Coolant '{name}' appears to be an ortho/para species but has no "
                    f"'conversion_factor'. Please specify 'conversion_factor' and optionally "
                    f"'parent_species' in the coolant configuration."
                )
                raise ValueError(msg)
            parent = explicit_parent or name[2:]
            parent_names.append(parent)
            conversion_factors.append(explicit_factor)
            conversion_modes.append(0)
        else:
            # Normal species
            parent_names.append(explicit_parent or name)
            conversion_factors.append(
                explicit_factor if explicit_factor is not None else 1.0
            )
            conversion_modes.append(0)

    # Validate that all parent species exist in the network
    species_dict = network.get_species_dict()
    missing_parents = []
    for i, (coolant, parent) in enumerate(zip(coolants, parent_names, strict=False)):
        if parent not in species_dict:
            missing_parents.append((i, coolant["name"], parent))

    if missing_parents:
        error_msg = "ERROR: The following coolants reference parent species that don't exist in the network:\n"
        for idx, coolant_name, parent_name in missing_parents:
            error_msg += f"  - Coolant #{idx + 1} '{coolant_name}' → parent species '{parent_name}' not found\n"
        error_msg += (
            "\nAvailable species: "
            + ", ".join(sorted(species_dict.keys())[:20])
            + ", ...\n"
        )
        error_msg += "Fix by adding 'parent_species: <existing_species>' in the coolant YAML config."
        raise ValueError(error_msg)

    f2py_constants: dict[str, Any] = {
        "n_species": len(network.get_species_list()),
        "n_reactions": len(network.get_reaction_list()),
        "n_physical_parameters": len(PHYSICAL_PARAMETERS),
        "n_dvode_stats": 19,
        "n_coolants": len(coolants),
        "max_coolants": len(coolants),
        "n_total_levels": n_total_levels,
        "n_se_stats_per_coolant": n_se_stats_per_coolant,
        "coolant_files": [c["file"] for c in coolants],
        "coolant_names": [c["name"] for c in coolants],
        "parent_names": parent_names,
        "conversion_factors": conversion_factors,
        "conversion_modes": conversion_modes,
        "coolant_data_dir": coolant_data_dir or "",
        "suggested_freq_rel_tol": suggested_freq_rel_tol,
        "missing_value_integer": MISSING_VALUE_INTEGER,
        "missing_value_float": MISSING_VALUE_FLOAT,
        "no_reactant_or_product": NO_REACTANT_OR_PRODUCT,
    }
    # Write all outputs to temporary files first; only replace finals if all succeed.
    tmp_paths = []
    try:
        _, tmp_species = mkstemp(dir=python_src_dir, prefix="species_", suffix=".csv.tmp")
        tmp_paths.append(Path(tmp_species))
        _, tmp_reactions = mkstemp(
            dir=python_src_dir, prefix="reactions_", suffix=".csv.tmp"
        )
        tmp_paths.append(Path(tmp_reactions))
        _, tmp_odes = mkstemp(dir=fortran_src_dir, prefix="odes_", suffix=".f90.tmp")
        tmp_paths.append(Path(tmp_odes))
        _, tmp_network = mkstemp(
            dir=fortran_src_dir, prefix="network_", suffix=".f90.tmp"
        )
        tmp_paths.append(Path(tmp_network))
        _, tmp_constants = mkstemp(
            dir=fortran_src_dir, prefix="f2py_constants_", suffix=".f90.tmp"
        )
        tmp_paths.append(Path(tmp_constants))

        write_species(Path(tmp_species), network.get_species_list())
        write_reactions(Path(tmp_reactions), network.get_reaction_list())
        write_odes_f90(
            Path(tmp_odes),
            network.get_species_list(),
            network.get_reaction_list(),
            enable_rates_storage=enable_rates_storage,
        )
        write_network_file(
            Path(tmp_network),
            network,
            enable_rates_storage=enable_rates_storage,
            gar_database=gar_database,
        )
        write_f90_constants(f2py_constants, Path(tmp_constants))

        # All writes succeeded — atomically replace the final files.
        shutil.move(tmp_species, python_src_dir / "species.csv")
        shutil.move(tmp_reactions, python_src_dir / "reactions.csv")
        shutil.move(tmp_odes, fortran_src_dir / "odes.f90")
        shutil.move(tmp_network, fortran_src_dir / "network.f90")
        shutil.move(tmp_constants, fortran_src_dir / "f2py_constants.f90")
        tmp_paths.clear()
    finally:
        for p in tmp_paths:
            if p.exists():
                p.unlink()
    # Note: constants.py now reads directly from f2py_constants module,
    # so we no longer need to write it during MakeRates.
    # After running MakeRates, just reinstall to update the Python constants.


def write_f90_constants(
    replace_dict: dict[str, Any],
    output_file_name: Path,
    template_file_path: str | Path = "fortran_templates",
) -> None:
    """Write the physical reactions to the f2py_constants.f90 file after every.

    run of makerates, this ensures the Fortran and Python bits are compatible.

    Parameters
    ----------
    replace_dict : dict[str, Any]
        The dictionary with keys to replace
    output_file_name : Path
        The path to target f2py_constants.f90 file
    template_file_path : str | Path
        The file to use as the template. (Default value = 'fortran_templates')

    """
    template_file_path = UCLCHEM_ROOT_DIR / "makerates" / template_file_path
    constants = Path(template_file_path / "f2py_constants.f90").read_text()

    # Handle string arrays separately for coolants
    if "coolant_files" in replace_dict and "coolant_names" in replace_dict:
        # Format coolant files
        coolant_files = replace_dict.pop("coolant_files")
        max_file_len = max(len(f) for f in coolant_files)
        coolant_files_str = truncate_line(
            ",".join(f'"{f.ljust(max_file_len)}"' for f in coolant_files), line_length=60
        )
        replace_dict["coolant_file_len"] = max_file_len
        replace_dict["coolant_files"] = truncate_line("/" + coolant_files_str + "/")

        # Format coolant names
        coolant_names = replace_dict.pop("coolant_names")
        max_name_len = max(len(n) for n in coolant_names)
        coolant_names_str = ",".join(f'"{n.ljust(max_name_len)}"' for n in coolant_names)
        replace_dict["coolant_name_len"] = max_name_len
        replace_dict["coolant_names"] = truncate_line(
            "/" + coolant_names_str + "/", line_length=60
        )

        # Format parent names (same pattern as coolant names)
        if "parent_names" in replace_dict:
            parent_names = replace_dict.pop("parent_names")
            max_parent_len = max(len(n) for n in parent_names)
            parent_names_str = ",".join(
                f'"{n.ljust(max_parent_len)}"' for n in parent_names
            )
            replace_dict["parent_name_len"] = max_parent_len
            replace_dict["parent_names"] = truncate_line(
                "/" + parent_names_str + "/", line_length=60
            )

    # Extract numeric arrays to be written via array_to_string (handles line limits)
    conversion_factors = replace_dict.pop("conversion_factors", None)
    conversion_modes = replace_dict.pop("conversion_modes", None)
    suggested_freq_rel_tol = replace_dict.pop("suggested_freq_rel_tol", None)

    constants = constants.format(**replace_dict)

    # Insert array declarations before end module using array_to_string
    extra_lines = ""
    if conversion_factors is not None:
        extra_lines += "    ! Conversion factors for abundance scaling (used when coolantConversionMode=0)\n"
        extra_lines += "    " + array_to_string(
            "coolantConversionFactors",
            np.array(conversion_factors),
            value_type="float",
            length_name="MAX_COOLANTS",
        )
    if conversion_modes is not None:
        extra_lines += "    ! Conversion mode: 0=fixed factor, 1=thermal OPR para, 2=thermal OPR ortho\n"
        extra_lines += "    " + array_to_string(
            "coolantConversionMode",
            np.array(conversion_modes),
            value_type="int",
            length_name="MAX_COOLANTS",
        )
    # Generate coolant_active array: all coolants enabled by default
    if "n_coolants" in replace_dict or conversion_factors is not None:
        n_coolants = (
            len(conversion_factors)
            if conversion_factors is not None
            else replace_dict.get("n_coolants", 0)
        )
        if n_coolants > 0:
            coolant_active_defaults = np.ones(n_coolants, dtype=bool)
            extra_lines += "    ! Per-coolant on/off toggle (can be changed at runtime via HeatingSettings)\n"
            extra_lines += "    " + array_to_string(
                "coolant_active",
                coolant_active_defaults,
                value_type="logical",
                parameter=False,
                length_name="MAX_COOLANTS",
            )
    if suggested_freq_rel_tol is not None:
        extra_lines += "    ! Suggested freq_rel_tol based on max observed deviation (with 10% margin)\n"
        extra_lines += f"    real(dp), parameter :: suggested_freq_rel_tol = {suggested_freq_rel_tol:.6e}_dp\n"

    if extra_lines:
        constants = constants.replace(
            "end module F2PY_CONSTANTS", extra_lines + "end module F2PY_CONSTANTS"
        )

    with Path(output_file_name).open("w") as fh:
        fh.writelines(constants)


def write_python_constants(
    replace_dict: dict[str, int], python_constants_file: Path
) -> None:
    """Write the python constants to the constants.py file (deprecated).

    As of the latest version, constants.py reads directly from the f2py_constants
    module, so this function is no longer needed. It's kept for backward compatibility
    but does nothing.

    Parameters
    ----------
    replace_dict : dict[str, int]
        Dict with keys to replace and their values (ignored)
    python_constants_file : Path
        Path to the target constant files (ignored)

    """
    warnings.warn(
        "write_python_constants() is deprecated. "
        "constants.py now reads directly from f2py_constants module.",
        DeprecationWarning,
        stacklevel=2,
    )
    # Do nothing - constants.py is now self-updating
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
                if line.startswith("# This file was machine generated with Makerates on"):
                    continue
            # For every line, try to find constants, if we find them, replace them,
            # if not, just print the line.
            hits = {
                constant: line.strip().startswith(constant) for constant in replace_dict
            }
            if any(hits.values()):
                # Filter, also we can only get one hit at a time
                variable = next(filter(hits.get, hits))
                print(f"{variable} = {replace_dict[variable]}")
            else:
                print(line, end="")


def write_species(file_name: str | Path, species_list: list[Species]) -> None:
    """Write the human readable species file. Note UCLCHEM doesn't use this file.

    Parameters
    ----------
    file_name : str | Path
        path to output file
    species_list : list[Species]
        List of species objects for network

    """
    with Path(file_name).open("w") as f:
        writer = csv.writer(
            f,
            delimiter=",",
            quotechar="|",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )
        writer.writerow(species_header)

        # Order is the same as in uclchem.species.species_header
        writer.writerows(
            [
                species.get_name(),
                species.get_mass(),
                species.get_binding_energy(),
                species.get_solid_fraction(),
                species.get_mono_fraction(),
                species.get_volcano_fraction(),
                species.get_enthalpy(),
                species.get_vdes(),
                species.get_diffusion_barrier(),
                species.get_vdiff(),
                species.get_Ix(),
                species.get_Iy(),
                species.get_Iz(),
                species.get_symmetry_factor(),
            ]
            for species in species_list
        )


# Write the reaction file in the desired format
def write_reactions(file_name: Path, reaction_list: list[Reaction]) -> None:
    """Write the human readable reaction file.

    Parameters
    ----------
    file_name : Path
        path to output file
    reaction_list : list[Reaction]
        List of reaction objects for network

    """
    with Path(file_name).open("w") as f:
        writer = csv.writer(
            f,
            delimiter=",",
            quotechar="|",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )
        writer.writerow(reaction_header)
        writer.writerows(
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
            for reaction in reaction_list
        )


def write_odes_f90(
    file_name: Path,
    species_list: list[Species],
    reaction_list: list[Reaction],
    *,
    enable_rates_storage: bool = False,
) -> None:
    """Write the ODEs in Modern Fortran. This is an actual code file.

    Parameters
    ----------
    file_name : Path
        Path to file where code will be written
    species_list : list[Species]
        List of species describing network
    reaction_list : list[Reaction]
        List of reactions describing network
    enable_rates_storage : bool
        Enable storage of writing rates to files.
        Default = False.

    """
    # First generate ODE contributions for all reactions
    species_names = [spec.get_name() for spec in species_list]

    for specie in species_list:
        logger.debug(f"{species_names.index(specie.get_name()) + 1}:{specie}")

    for i, reaction in enumerate(reaction_list):
        logger.debug(f"RATE_CONSTANTS({i + 1}):{reaction}")
        reaction.generate_ode_bit(i, species_names)

    ydotString = build_ode_string(
        species_list, reaction_list, enable_rates_storage=enable_rates_storage
    )
    Path(file_name).write_text(ydotString)


def write_jacobian(file_name: Path, species_list: list[Species]) -> None:
    """Write jacobian in Modern Fortran.

    This has never improved UCLCHEM's speed, and so is not used in the code as it stands.
    Current only works for three phase model.

    Parameters
    ----------
    file_name : Path
        Path to jacobian file
    species_list : list[Species]
        List of species AFTER being processed by build_ode_string

    """
    species_names = ""
    with Path(file_name).open("w") as output:
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
                    di_dj_parts: list[str] = [
                        f"-{x}".replace(f"*Y({j})", "", 1)
                        for x in losses
                        if f"*Y({j})" in x
                    ]
                    di_dj_parts += [
                        f"+{x}".replace(f"*Y({j})", "", 1)
                        for x in gains
                        if f"*Y({j})" in x
                    ]
                    # of course there might be y=a*x*x so we only replace first instance and if
                    # there's still an instance we put a factor of two in since
                    # dy/dx=2ax for y=a*x*x
                    di_dj_parts = [
                        x + "*2" if f"*Y({j})" in x else x for x in di_dj_parts
                    ]

                    # safeMantle is a stand in for the surface so do it manually here
                    # since it's divided by safemantle, derivative is negative so sign flips
                    # and we get another factor of 1/safeMantle
                    if species_list[j - 1].get_name() == "SURFACE":
                        di_dj_parts = [
                            f"+{x}/safeMantle" for x in losses if "/safeMantle" in x
                        ]
                        di_dj_parts += [
                            f"-{x}/safeMantle" for x in gains if "/safeMantle" in x
                        ]
                    if len(di_dj_parts) > 0:
                        output.write(f"J({i + 1},{j})=" + "".join(di_dj_parts) + "\n")

            # tackle density separately handled
            j += 1
            if species.get_name() == "SURFACE":
                di_dj = f"J({i + 1},{j})=SUM(J(surfaceList,{j}))\n"
                output.write(di_dj)
            elif species.get_name() == "BULK":
                if species_names.count("@") > 0:
                    di_dj = f"J({i + 1},{j})=SUM(J(bulkList,{j}))\n"
                    output.write(di_dj)
            else:
                di_dj_d_parts: list[str] = [
                    f"-{x}".replace("*D", "", 1) for x in losses if "*D" in x
                ]
                di_dj_d_parts += [
                    f"+{x}".replace("*D", "", 1) for x in gains if "*D" in x
                ]
                di_dj_d_parts = [x + "*2" if "*D" in x else x for x in di_dj_d_parts]
                if len(di_dj_d_parts) > 0:
                    output.write(f"J({i + 1},{j})=" + "".join(di_dj_d_parts) + "\n")
        i += 2
        di_dj = f"J({i},{i})=ddensdensdot(D)\n"
        output.write(di_dj)


def build_ode_string(
    species_list: list[Species],
    reaction_list: list[Reaction],
    *,
    enable_rates_storage: bool = False,
) -> str:
    """Build the ODE string.

    A long, complex function that does the messy work of creating the actual ODE
    code to calculate the rate of change of each species. Test any change to this code
    thoroughly because ODE mistakes are very hard to spot.

    Parameters
    ----------
    species_list : list[Species]
        List of species in network
    reaction_list : list[Reaction]
        List of reactions in network
    enable_rates_storage : bool
        Enable the writing of the rates to the disk.
        Default = False.

    Returns
    -------
    ode_string : str
        One long string containing the entire ODE fortran code.

    """
    # We create a string of losses and gains for each species so initialize them all as ""
    species_names = []
    for species in species_list:
        species_names.append(species.get_name())
        species.losses = ""
        species.gains = ""

    bulk_index = species_names.index("BULK")
    surface_index = species_names.index("SURFACE")
    total_swap = ""

    for reaction in reaction_list:
        for species_name in reaction.get_reactants():
            if species_name in species_names:
                # Eley-Rideal reactions take a share of total freeze out rate
                # which is already accounted for so we add as a loss term to the
                # frozen version of the species rather than the gas version
                if (reaction.get_reaction_type() == "ER") and (
                    not species_list[
                        species_names.index(species_name)
                    ].is_surface_species()
                ):
                    species_list[
                        species_names.index("#" + species_name)
                    ].losses += reaction.ode_bit
                else:
                    species_list[
                        species_names.index(species_name)
                    ].losses += reaction.ode_bit
                if reaction.get_reaction_type() == "BULKSWAP":
                    total_swap += reaction.ode_bit
        for species_name in reaction.get_products():
            if species_name in species_names:
                species_list[species_names.index(species_name)].gains += reaction.ode_bit

    total_swap = separate_common_terms(total_swap[1:], "ratioSurfaceToBulk")

    # 10 August 2026, Tobias Dijkhuis:
    #    safeMantle = MAX(MIN_ABUND, sum(Y(surfaceList))) # ruff: ignore[commented-out-code]
    # could be replaced with
    #    safeMantle = MAX(MIN_ABUND, Y(nSurface)) # ruff: ignore[commented-out-code]
    # and same for bulk, but if I did that, I suddenly got conservation errors for hotcore models.
    # Probably because then safeMantle is the same as some surface species abunds because it is
    # clamped twice (also in Y_safe, which is passed here in chemistry.f90).
    ode_string = dedent("""    module ODES
        use constants, only: dp, MIN_ABUND
        use f2py_constants, only: nReac, nSpec
        use network, only: refractoryList, bulkList, surfaceList, REACTIONRATE, &
            nSurface, nBulk
        use surfacereactions, only: useGarrod2011Transfer, NUM_SITES_PER_GRAIN, GAS_DUST_DENSITY_RATIO

        implicit none

        public
    contains
        pure subroutine GETYDOT(RATE_CONSTANTS, Y, surfaceCoverage, D, YDOT, surfGrowthUncorrected)
            real(dp), intent(in) :: RATE_CONSTANTS(nReac), Y(nSpec+2)
            real(dp), intent(in) :: D
            real(dp), intent(in) :: surfaceCoverage
            real(dp), intent(out) :: YDOT(nSpec+2)
            real(dp), intent(out) :: surfGrowthUncorrected

            real(dp) :: totalSwap, LOSS, PROD, alpha
            real(dp) :: safeMantle, safeBulk, ratioSurfaceToBulk, bulkLayersReciprocal

            safeMantle = MAX(MIN_ABUND, sum(Y(surfaceList)))
            safeBulk   = MAX(MIN_ABUND, sum(Y(bulkList)))
            if (refractoryList(1) > 0) then
                safeBulk = MAX(MIN_ABUND, safeBulk - SUM(Y(refractoryList)))
            end if
            ratioSurfaceToBulk   = MIN(1.0_dp, safeMantle/safeBulk)
            bulkLayersReciprocal = MIN(1.0_dp, NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
    """)
    ode_string += truncate_line(f"        totalSwap={total_swap}\n\n")

    if enable_rates_storage:
        for n, reaction in enumerate(reaction_list):
            ode_string += truncate_line(
                f"        REACTIONRATE({n + 1})={reaction.ode_bit}\n"
            )

    # Get total rate of change of bulk and surface by adding ydots
    species_list[bulk_index].gains = "+sum(YDOT(bulkList))"
    species_list[surface_index].gains = "+sum(YDOT(surfaceList))"

    for n, species in enumerate(species_list):
        ydot_string = species_ode_string(n, species)
        ode_string += ydot_string

    ode_string += f"\n        surfGrowthUncorrected = YDOT({surface_index + 1})\n"

    # now add bulk transfer to rate of change of surface species
    # after they've already been calculated
    ode_string += dedent("""

    ! Update surface species for bulk growth
    if (YDOT(nSurface) < 0) then
        ! Since ydot(surface_index) is negative, bulk is lost and surface forms
        if (useGarrod2011Transfer) then
            ! Three-phase treatment of Garrod & Pauly 2011
            ! Real value of alpha_des: alpha_des = MIN(1.0_dp, safeBulk / safeMantle).
            ! However, the YDOTs calculated below need to be multiplied with Y(bulkspec)/safeBulk,
            ! so we divide by safeBulk here to save time
            alpha = MIN(1.0_dp, safeBulk/safeMantle)/safeBulk
        else
            ! Hasegawa & Herbst 1993
            alpha = MIN(1.0_dp, surfaceCoverage*safeMantle)/safeBulk
        end if
    """)

    surf_species = [
        i
        for i in species_list
        if i.get_name() not in {"SURFACE", "BULK"} and i.is_surface_species()
    ]
    i = len(reaction_list)
    j = len(reaction_list) + len(surf_species)
    for n, species in enumerate(species_list):
        if species.get_name()[0] == "#":
            i += 1
            j += 1
            bulk_partner = species_names.index(species.get_name().replace("#", "@"))
            if enable_rates_storage:
                ode_string += f"        REACTIONRATE({i}) = -YDOT(nSurface)*alpha*Y({bulk_partner + 1})\n"
                ode_string += f"        REACTIONRATE({j}) = 0.0_dp\n"
            if not species_list[bulk_partner].is_refractory:
                ode_string += f"        YDOT({n + 1})=YDOT({n + 1})-YDOT(nSurface)*alpha*Y({bulk_partner + 1})\n"
        if species.get_name()[0] == "@" and not species.is_refractory:
            ode_string += (
                f"        YDOT({n + 1})=YDOT({n + 1})+YDOT(nSurface)*alpha*Y({n + 1})\n"
            )
    ode_string += "    else\n"
    ode_string += "        ! surfaceCoverage = fractional surface coverage\n"
    ode_string += "        ! Real value of surfaceCoverage: surfaceCoverage = safeMantle / NUM_MONOLAYERS_IS_SURFACE * GAS_DUST_DENSITY_RATIO / NUM_SITES_PER_GRAIN\n"
    ode_string += "        ! However, the YDOTs calculated below need to be multiplied with Y(surfspec)/safeMantle, so we divide by safeMantle here to save time\n"
    ode_string += "        ! In chemistry.f90: surfaceCoverage = 1/NUM_MONOLAYERS_IS_SURFACE * GAS_DUST_DENSITY_RATIO / NUM_SITES_PER_GRAIN\n"
    ode_string += "        alpha = MIN(1.0_dp, surfaceCoverage*safeMantle)/safeMantle\n"
    i = len(reaction_list)
    j = len(reaction_list) + len(surf_species)
    for n, species in enumerate(species_list):
        if species.get_name()[0] == "@":
            i += 1
            j += 1
            surface_version = species_names.index(species.get_name().replace("@", "#"))
            if enable_rates_storage:
                ode_string += f"        REACTIONRATE({i}) = 0.0_dp\n"
                ode_string += f"        REACTIONRATE({j}) = -YDOT(nSurface)*alpha*Y({surface_version + 1})\n"
            ode_string += f"        YDOT({n + 1})=YDOT({n + 1})+YDOT(nSurface)*alpha*Y({surface_version + 1})\n"
        if species.get_name()[0] == "#":
            ode_string += (
                f"        YDOT({n + 1})=YDOT({n + 1})-YDOT(nSurface)*alpha*Y({n + 1})\n"
            )
    ode_string += "    end if\n"

    # once bulk transfer has been added, odes for bulk and surface must account for it
    ode_string += (
        "    ! Update total rate of change of bulk and surface for bulk growth\n"
    )
    ode_string += species_ode_string(bulk_index, species_list[bulk_index])
    ode_string += species_ode_string(surface_index, species_list[surface_index])
    ode_string += """    end subroutine GETYDOT
end module ODES"""
    return ode_string


def species_ode_string(n: int, species: Species) -> str:
    """Build the string of Fortran code for a species once it's loss and gains.

    strings have been produced.

    Parameters
    ----------
    n : int
        Index of species in python format
    species : Species
        species object

    Returns
    -------
    str
        the fortran code for the rate of change of the species

    """
    ydot_string = ""
    if species.losses:
        # 10 August 2026, Tobias Dijkhuis:
        # Could also separate common term f"Y({n+1})", here, but if I did that
        # I got much worse runtime performance. Probably the compiler cannot
        # optimize it as much?
        loss_string = "        LOSS = " + species.losses[1:] + "\n"
        ydot_string += loss_string
    if species.gains:
        prod_string = "        PROD = " + species.gains[1:] + "\n"
        ydot_string += prod_string

    if ydot_string:
        ydot_string += f"        YDOT({n + 1}) = "
        # start with empty string and add production and loss terms if they exists
        if species.gains:
            ydot_string += "PROD"
        if species.losses:
            ydot_string += "-LOSS"
        ydot_string += "\n"
    else:
        ydot_string = f"        YDOT({n + 1}) = 0.0_dp\n"
    ydot_string = truncate_line(ydot_string)
    return ydot_string


def write_evap_lists(network_file: IO[str], species_list: list[Species]) -> None:
    """Write evaporation list to network file.

    Two phase networks mimic episodic thermal desorption seen in lab (see Viti et al. 2004)
    by desorbing fixed fractions of material at specific temperatures.
    Three phase networks just use binding energy and that fact we set binding energies
    in bulk to water by default. This function writes all necessary arrays to
    the network file so these processes work.

    Parameters
    ----------
    network_file : IO[str]
        Open file handle to which the network code is being written
    species_list : list[Species]
        List of species in network

    Raises
    ------
    NameError
        If a species desorbs as another species that is not in the species
        list.

    """
    gasIceList = []
    surfacelist = []
    solidList = []
    monoList = []
    volcList = []
    binding_energyList = []
    customVdesList = []
    diffusion_barriersList = []
    customVdiffList = []
    inertiaProducts = []
    isLinears = []

    enthalpyList = []
    bulkList = []
    iceList = []
    refractoryList = []
    species_names = [spec.get_name() for spec in species_list]
    for i, species in enumerate(species_list):
        if species.get_name()[0] == "#":
            # find gas phase version of grain species.
            # For #CO it looks for first species in list with just CO
            # and then finds the index of that
            try:
                j = species_names.index(species.get_standard_desorb_products()[0])
            except ValueError:
                # Standard gas counterpart not in species list (e.g. isomer-only networks).
                # Fall back to the user-defined DESORB product if one was supplied.
                desorb_fallback = species.get_desorb_products()[0]
                if (
                    desorb_fallback not in {"NAN", ""}
                    and desorb_fallback in species_names
                ):
                    j = species_names.index(desorb_fallback)
                else:
                    error = (
                        f"{species.get_name()} standard desorb product is "
                        f"{species.get_standard_desorb_products()[0]}"
                    )
                    error += " which is not in species list.\n"
                    error += "If this species desorbs to a non-standard gas product, add a single-product DESORB\n"
                    error += "reaction in your reaction file to specify the gasIceList entry, then re-run Makerates."
                    raise NameError(error) from None

            # plus ones as fortran and python label arrays differently
            surfacelist.append(i + 1)
            gasIceList.append(j + 1)
            solidList.append(species.get_solid_fraction())
            monoList.append(species.get_mono_fraction())
            volcList.append(species.get_volcano_fraction())
            iceList.append(i + 1)

            binding_energyList.append(species.get_binding_energy())
            customVdesList.append(species.get_vdes())
            diffusion_barriersList.append(species.get_diffusion_barrier())
            customVdiffList.append(species.get_vdiff())

            isLinears.append(species.is_linear())
            inertiaProducts.append(species.calculate_rotational_partition_factor())

            enthalpyList.append(species.get_enthalpy())
        elif species.get_name()[0] == "@":
            try:
                j = species_names.index(species.get_standard_desorb_products()[0])
            except ValueError:
                desorb_fallback = species.get_desorb_products()[0]
                if (
                    desorb_fallback not in {"NAN", ""}
                    and desorb_fallback in species_names
                ):
                    j = species_names.index(desorb_fallback)
                else:
                    error = (
                        f"{species.get_name()} standard desorb product is "
                        f"{species.get_standard_desorb_products()[0]}"
                    )
                    error += " which is not in species list.\n"
                    error += "If this species desorbs to a non-standard gas product, add a single-product DESORB\n"
                    error += "reaction in your reaction file to specify the gasIceList entry, then re-run Makerates."
                    raise NameError(error) from None
            gasIceList.append(j + 1)
            bulkList.append(i + 1)
            iceList.append(i + 1)

            binding_energyList.append(species.get_binding_energy())
            customVdesList.append(species.get_vdes())
            diffusion_barriersList.append(species.get_diffusion_barrier())
            customVdiffList.append(species.get_vdiff())

            isLinears.append(species.is_linear())
            inertiaProducts.append(species.calculate_rotational_partition_factor())

            enthalpyList.append(species.get_enthalpy())
            if species.is_refractory:
                refractoryList.append(i + 1)

    # dummy index that will be caught by UCLCHEM
    if len(refractoryList) == 0:
        refractoryList = [-999]

    network_file.write(f"integer, parameter :: N_ICE_SPECIES = {len(iceList)}\n")
    network_file.write(f"integer, parameter :: N_SURFACE_SPECIES = {len(surfacelist)}\n")
    network_file.write(f"integer, parameter :: N_BULK_SPECIES = {len(bulkList)}\n")

    network_file.write(
        array_to_string(
            "surfaceList", surfacelist, value_type="int", length_name="N_SURFACE_SPECIES"
        )
    )
    if len(bulkList) > 0:
        network_file.write(
            array_to_string(
                "bulkList", bulkList, value_type="int", length_name="N_BULK_SPECIES"
            )
        )
    network_file.write(
        array_to_string("iceList", iceList, value_type="int", length_name="N_ICE_SPECIES")
    )
    network_file.write(
        array_to_string(
            "gasIceList", gasIceList, value_type="int", length_name="N_ICE_SPECIES"
        )
    )
    network_file.write(
        array_to_string(
            "solidFractions",
            solidList,
            value_type="float",
            length_name="N_SURFACE_SPECIES",
        )
    )
    network_file.write(
        array_to_string(
            "monoFractions", monoList, value_type="float", length_name="N_SURFACE_SPECIES"
        )
    )
    network_file.write(
        array_to_string(
            "volcanicFractions",
            volcList,
            value_type="float",
            length_name="N_SURFACE_SPECIES",
        )
    )
    network_file.write(
        array_to_string(
            "bindingEnergy",
            binding_energyList,
            value_type="float",
            parameter=False,
            length_name="N_ICE_SPECIES",
        )
    )
    network_file.write(
        replace_value_with_name(
            array_to_string(
                "customVdes",
                customVdesList,
                value_type="float",
                length_name="N_ICE_SPECIES",
            ),
            MISSING_VALUE_FLOAT,
            "MISSING_VALUE_FLOAT",
        )
    )
    network_file.write(
        replace_value_with_name(
            array_to_string(
                "diffusionBarrier",
                diffusion_barriersList,
                value_type="float",
                parameter=False,
                length_name="N_ICE_SPECIES",
            ),
            MISSING_VALUE_FLOAT,
            "MISSING_VALUE_FLOAT",
        )
    )
    network_file.write(
        replace_value_with_name(
            array_to_string(
                "customVdiff",
                customVdiffList,
                value_type="float",
                length_name="N_ICE_SPECIES",
            ),
            MISSING_VALUE_FLOAT,
            "MISSING_VALUE_FLOAT",
        )
    )

    network_file.write(
        array_to_string(
            "moleculeIsLinear",
            isLinears,
            value_type="logical",
            parameter=False,
            length_name="N_ICE_SPECIES",
        )
    )
    network_file.write(
        replace_value_with_name(
            array_to_string(
                "inertiaProducts",
                inertiaProducts,
                value_type="float",
                parameter=False,
                length_name="N_ICE_SPECIES",
            ),
            MISSING_VALUE_FLOAT,
            "MISSING_VALUE_FLOAT",
        )
    )
    network_file.write(
        replace_value_with_name(
            array_to_string(
                "formationEnthalpy",
                enthalpyList,
                value_type="float",
                parameter=False,
                length_name="N_ICE_SPECIES",
            ),
            MISSING_VALUE_FLOAT,
            "MISSING_VALUE_FLOAT",
        )
    )
    network_file.write(
        array_to_string("refractoryList", refractoryList, value_type="int")
    )


def write_network_file(
    file_name: Path,
    network: Network,
    gar_database: dict[str, np.ndarray] | None = None,
    *,
    enable_rates_storage: bool = False,
) -> None:
    """Write the Fortran code file that contains all network information for UCLCHEM.

    This includes lists of reactants, products, binding energies, formationEnthalpies
    and so on.

    Parameters
    ----------
    file_name : Path
        The file name where the code will be written.
    network : Network
        A Network object built from lists of species and reactions.
    gar_database : dict[str, np.ndarray] | None
        Database for grain-activated recombination
        reactions. Default = None.
    enable_rates_storage : bool
        Enable storage of writing rates to files.
        Default = False.

    Raises
    ------
    AssertionError
        If exothermicity is non-zero but ``enable_rates_storage`` is ``False``.
    RuntimeError
        If a :class:`~uclchem.makerates.reaction.CoupledReaction` with partner None is found.

    """
    species_list = network.get_species_list()
    reaction_list = network.get_reaction_list()

    with file_name.open("w") as openFile:
        openFile.write(
            dedent("""        module network
            use constants, only: dp, REAC_NOT_PRESENT
            use f2py_constants, only: nSpec, nReac, MISSING_VALUE_FLOAT, MISSING_VALUE_INTEGER, &
                NO_REACTANT_OR_PRODUCT

            implicit none

            public

        """)
        )

        # write arrays of all species stuff
        n_species = len(species_list)
        names = [""] * n_species
        atoms = [0] * n_species
        masses = [0] * n_species
        for i, species in enumerate(species_list):
            names[i] = species.get_name()
            masses[i] = species.get_mass()
            atoms[i] = species.get_n_atoms()

        speciesIndices = ",".join(
            f"{name}={species_index}"
            for name, species_index in network.species_indices.items()
        )
        openFile.write(truncate_line("integer, parameter :: " + speciesIndices + "\n"))
        openFile.write("logical, parameter :: THREE_PHASE = .true.\n")
        openFile.write(
            array_to_string("specname", names, value_type="string", length_name="nSpec")
        )
        openFile.write(
            array_to_string("mass", masses, value_type="float", length_name="nSpec")
        )
        openFile.write(
            array_to_string("atomCounts", atoms, value_type="int", length_name="nSpec")
        )

        unique_elements, elem_count_2d = get_element_counts_per_species(
            network.get_species_list()
        )

        padded_elems = pad_to_max_length(unique_elements)

        openFile.write(f"integer, parameter :: n_elem_tracked = {len(unique_elements)}\n")
        openFile.write(
            array_to_string(
                "elem_names",
                padded_elems,
                value_type="string",
                length_name="n_elem_tracked",
            )
        )
        openFile.write(
            array_to_string(
                "elem_count",
                elem_count_2d,
                value_type="int",
                length_name="nSpec, n_elem_tracked",
            )
        )

        # then write evaporation stuff
        write_evap_lists(openFile, species_list)

        # finally all reactions
        n_reactions = len(reaction_list)
        reactant1 = [0] * n_reactions
        reactant2 = [0] * n_reactions
        reactant3 = [0] * n_reactions
        prod1 = [0] * n_reactions
        prod2 = [0] * n_reactions
        prod3 = [0] * n_reactions
        prod4 = [0] * n_reactions
        alpha = [0.0] * n_reactions
        beta = [0.0] * n_reactions
        gama = [0.0] * n_reactions
        reacTypes = [""] * n_reactions
        tmins = [0.0] * n_reactions
        tmaxs = [0.0] * n_reactions
        reduced_masses = [0.0] * n_reactions
        extrapolations = [False] * n_reactions
        exothermicity = [0.0] * n_reactions
        partner_indices = [0] * n_reactions

        # store important reactions
        reaction_indices = ""
        for reaction_name, reaction_idx in network.important_reactions.items():
            reaction_indices += reaction_name + f"={reaction_idx},"
        reaction_indices = truncate_line(reaction_indices[:-1]) + "\n"
        openFile.write("integer, parameter :: " + reaction_indices)

        species_dict = {
            species.get_name(): idx + 1 for idx, species in enumerate(species_list)
        }
        for idx, reaction in enumerate(reaction_list):
            reactants = reaction.get_reactants()
            products = reaction.get_products()
            reactant1[idx] = species_dict.get(reactants[0], NO_REACTANT_OR_PRODUCT)
            reactant2[idx] = species_dict.get(reactants[1], NO_REACTANT_OR_PRODUCT)
            reactant3[idx] = species_dict.get(reactants[2], NO_REACTANT_OR_PRODUCT)
            prod1[idx] = species_dict.get(products[0], NO_REACTANT_OR_PRODUCT)
            prod2[idx] = species_dict.get(products[1], NO_REACTANT_OR_PRODUCT)
            prod3[idx] = species_dict.get(products[2], NO_REACTANT_OR_PRODUCT)
            prod4[idx] = species_dict.get(products[3], NO_REACTANT_OR_PRODUCT)
            alpha[idx] = reaction.get_alpha()
            beta[idx] = reaction.get_beta()
            gama[idx] = reaction.get_gamma()
            tmaxs[idx] = reaction.get_temphigh()
            tmins[idx] = reaction.get_templow()
            reduced_masses[idx] = reaction.get_reduced_mass()
            reacTypes[idx] = reaction.get_reaction_type()
            extrapolations[idx] = reaction.get_extrapolation()
            exothermicity[idx] = reaction.get_exothermicity()

            if isinstance(reaction, CoupledReaction):
                partner = reaction.get_partner()
                if partner is None:
                    msg = f"Found CoupledReaction '{reaction}' with partner None."
                    raise RuntimeError(msg)
                partner_indices[idx] = reaction_list.index(partner) + 1
            else:
                partner_indices[idx] = MISSING_VALUE_INTEGER
        reacTypes_array = np.array(reacTypes)

        reaction_names = [str(reaction) for reaction in reaction_list]
        for species in species_list:
            if species.is_surface_species() and species.get_name() not in {
                "SURFACE",
                "BULK",
            }:
                reaction_name = (
                    f"{species.get_name()} + SURFACETRANSFER -> @{species.get_name()[1:]}"
                )
                reaction_names.append(reaction_name)
        for species in species_list:
            if species.is_surface_species() and species.get_name() not in {
                "SURFACE",
                "BULK",
            }:
                reaction_name = (
                    f"@{species.get_name()[1:]} + SURFACETRANSFER -> {species.get_name()}"
                )
                reaction_names.append(reaction_name)

        # Save some memory by only allocating things we actually want to use:
        if enable_rates_storage:
            openFile.write("real(dp) :: REACTIONRATE(nReac + N_ICE_SPECIES)\n")
            openFile.write("logical, parameter :: storeRatesComputation=.true.\n")
        else:
            openFile.write("real(dp) :: REACTIONRATE(1)\n")
            openFile.write("logical, parameter :: storeRatesComputation=.false.\n")
        if any(exo != 0 for exo in exothermicity):
            if not enable_rates_storage:
                msg = "Chemical heating can only be enabled if rates are being computed and stored in memory. Enable `enable_rates_storage` in the user_settings."
                raise AssertionError(msg)
            openFile.write(
                array_to_string(
                    "exothermicities",
                    exothermicity,
                    value_type="float",
                    parameter=True,
                    length_name="nReac",
                )
            )
            openFile.write("logical, parameter :: enableChemicalHeating = .true.\n")
        else:
            openFile.write("real(dp) :: exothermicities(nReac)\n")
            openFile.write("logical, parameter :: enableChemicalHeating = .false.\n")

        openFile.write(
            replace_value_with_name(
                array_to_string("re1", reactant1, value_type="int", length_name="nReac"),
                NO_REACTANT_OR_PRODUCT,
                "NO_REACTANT_OR_PRODUCT",
            )
        )
        openFile.write(
            replace_value_with_name(
                array_to_string("re2", reactant2, value_type="int", length_name="nReac"),
                NO_REACTANT_OR_PRODUCT,
                "NO_REACTANT_OR_PRODUCT",
            )
        )
        openFile.write(
            replace_value_with_name(
                array_to_string("re3", reactant3, value_type="int", length_name="nReac"),
                NO_REACTANT_OR_PRODUCT,
                "NO_REACTANT_OR_PRODUCT",
            )
        )
        openFile.write(
            replace_value_with_name(
                array_to_string("p1", prod1, value_type="int", length_name="nReac"),
                NO_REACTANT_OR_PRODUCT,
                "NO_REACTANT_OR_PRODUCT",
            )
        )
        openFile.write(
            replace_value_with_name(
                array_to_string("p2", prod2, value_type="int", length_name="nReac"),
                NO_REACTANT_OR_PRODUCT,
                "NO_REACTANT_OR_PRODUCT",
            )
        )
        openFile.write(
            replace_value_with_name(
                array_to_string("p3", prod3, value_type="int", length_name="nReac"),
                NO_REACTANT_OR_PRODUCT,
                "NO_REACTANT_OR_PRODUCT",
            )
        )
        openFile.write(
            replace_value_with_name(
                array_to_string("p4", prod4, value_type="int", length_name="nReac"),
                NO_REACTANT_OR_PRODUCT,
                "NO_REACTANT_OR_PRODUCT",
            )
        )
        openFile.write(
            array_to_string(
                "alpha", alpha, value_type="float", parameter=False, length_name="nReac"
            )
        )
        openFile.write(
            array_to_string(
                "beta", beta, value_type="float", parameter=False, length_name="nReac"
            )
        )
        openFile.write(
            array_to_string(
                "gama", gama, value_type="float", parameter=False, length_name="nReac"
            )
        )
        openFile.write(
            array_to_string(
                "minTemps", tmins, value_type="float", parameter=True, length_name="nReac"
            )
        )
        openFile.write(
            array_to_string(
                "maxTemps", tmaxs, value_type="float", parameter=True, length_name="nReac"
            )
        )
        openFile.write(
            array_to_string(
                "reducedMasses",
                reduced_masses,
                value_type="float",
                parameter=True,
                length_name="nReac",
            )
        )
        openFile.write(
            array_to_string(
                "extrapolateRateConstants",
                extrapolations,
                value_type="logical",
                parameter=True,
                length_name="nReac",
            )
        )
        openFile.write(
            replace_value_with_name(
                array_to_string(
                    "partnerIndices",
                    partner_indices,
                    value_type="int",
                    parameter=True,
                    length_name="nReac",
                ),
                MISSING_VALUE_INTEGER,
                "MISSING_VALUE_INTEGER",
            )
        )

        partners = get_desorption_freeze_partners(reaction_list)
        openFile.write(f"integer, parameter :: N_FREEZE_PARTNERS = {len(partners)}\n")
        openFile.write(
            array_to_string(
                "freezePartners",
                partners,
                value_type="int",
                parameter=True,
                length_name="N_FREEZE_PARTNERS",
            )
        )

        openFile.write(
            f"integer, parameter :: N_GAR_SPECIES = {len(gar_database) if gar_database else 1}\n"
        )
        openFile.write("integer, parameter :: N_GAR_PARAMS = 7\n")
        openFile.write(
            array_to_string(
                "garParams",
                np.array(list(gar_database.values()))
                if gar_database
                else np.zeros((1, 7)),
                value_type="float",
                parameter=True,
                length_name="N_GAR_SPECIES, N_GAR_PARAMS",
            )
        )

        for reaction_type in sorted(REACTION_TYPES):
            list_name = reaction_type.lower() + "Reacs"
            reac_indices = np.where(reacTypes_array == reaction_type)[0]
            indices: list[int]
            if len(reac_indices > 1):
                indices = [int(reac_indices[0]) + 1, int(reac_indices[-1]) + 1]
            else:
                # We still want a dummy array if the reaction type isn't in network
                indices = [99999, 99999]
            openFile.write(
                array_to_string(
                    list_name, indices, value_type="int", parameter=True
                ).replace("99999", "REAC_NOT_PRESENT")
            )

        # Write LHDES and ERDES mapping arrays (Feature 3: LH/ER-DES mapping)
        # These arrays map chemical reactive desorption reactions to their parent reactions
        LHDEScorrespondingLHreacs = []
        for reaction in reaction_list:
            if reaction.get_reaction_type() == "LHDES":
                if (
                    hasattr(reaction, "get_partner")
                    and reaction.get_partner() is not None
                ):
                    partner = reaction.get_partner()
                    reacIndex = reaction_list.index(partner) + 1
                    LHDEScorrespondingLHreacs.append(reacIndex)
                else:
                    # If no partner set, use dummy index
                    LHDEScorrespondingLHreacs.append(99999)

        openFile.write(
            f"integer, parameter :: N_LHDES_REACTIONS = {len(LHDEScorrespondingLHreacs)}\n"
        )

        # Write array (use dummy if empty for backward compatibility)
        if len(LHDEScorrespondingLHreacs) == 0:
            LHDEScorrespondingLHreacs = [99999]
        openFile.write(
            array_to_string(
                "LHDEScorrespondingLHreacs",
                LHDEScorrespondingLHreacs,
                value_type="int",
                parameter=True,
                length_name="N_LHDES_REACTIONS"
                if len(LHDEScorrespondingLHreacs) > 2  # ruff: ignore[magic-value-comparison]
                else None,
            ).replace("99999", "REAC_NOT_PRESENT")
        )

        ERDEScorrespondingERreacs = []
        for reaction in reaction_list:
            if reaction.get_reaction_type() == "ERDES":
                if (
                    hasattr(reaction, "get_partner")
                    and reaction.get_partner() is not None
                ):
                    partner = reaction.get_partner()
                    reacIndex = reaction_list.index(partner) + 1
                    ERDEScorrespondingERreacs.append(reacIndex)
                else:
                    # If no partner set, use dummy index
                    ERDEScorrespondingERreacs.append(99999)

        openFile.write(
            f"integer, parameter :: N_ERDES_REACTIONS = {len(ERDEScorrespondingERreacs)}\n"
        )
        # Write array (use dummy if empty for backward compatibility)
        if len(ERDEScorrespondingERreacs) == 0:
            ERDEScorrespondingERreacs = [99999]
        if len(ERDEScorrespondingERreacs) == 1:
            # Fortran needs at least 2 elements for array
            ERDEScorrespondingERreacs.append(ERDEScorrespondingERreacs[0])
        openFile.write(
            array_to_string(
                "ERDEScorrespondingERreacs",
                ERDEScorrespondingERreacs,
                value_type="int",
                parameter=True,
                length_name="N_ERDES_REACTIONS"
                if len(ERDEScorrespondingERreacs) > 2  # ruff: ignore[magic-value-comparison]
                else None,
            ).replace("99999", "REAC_NOT_PRESENT")
        )
        openFile.write("end module network")


def get_desorption_freeze_partners(reaction_list: list[Reaction]) -> list[int]:
    """Every desorption has a corresponding freeze out eg desorption of #CO and freeze of CO.

    This find the corresponding freeze out for every desorb so that when desorb>>freeze
    we can turn off freeze out in UCLCHEM.

    Parameters
    ----------
    reaction_list : list[Reaction]
        Reactions in network

    Returns
    -------
    partners : list[int]
        list of indices of freeze out reactions
        matching order of desorptions.

    """
    freeze_species = [
        x.get_products()[0] for x in reaction_list if x.get_reactants()[1] == "DESCR"
    ]
    partners = []
    for spec in freeze_species:
        for i, reaction in enumerate(reaction_list):
            if (
                reaction.get_reaction_type() == "FREEZE"
                and reaction.get_reactants()[0] == spec
            ):
                partners.append(i + 1)
                break
    return partners


def copy_coolant_files(source_dir: str | Path | None = None) -> None:
    """Copy coolant data files to the package data directory for installation.

    This function copies .dat files from the source coolant directory
    (typically Makerates/data/collisional_rates/) to src/uclchem/data/collisional_rates/
    so they can be installed with the package via meson.

    Parameters
    ----------
    source_dir : str | Path | None
        Optional source directory.
        If None, uses get_default_coolant_directory(). Default = None.

    Raises
    ------
    FileNotFoundError
        If source directory doesn't exist or contains no .dat files.

    """
    # Determine source directory
    if source_dir is None:
        source_dir = get_default_coolant_directory()

    source_path = Path(source_dir)
    if not source_path.is_dir():
        msg = (
            f"Source coolant directory not found: {source_path}\n"
            f"Expected to find .dat files for coolant data."
        )
        raise FileNotFoundError(msg)

    # Find .dat files in source directory
    dat_files = list(source_path.glob("*.dat"))
    if not dat_files:
        msg = (
            f"No .dat files found in source directory: {source_path}\n"
            f"Cannot copy coolant data files."
        )
        raise FileNotFoundError(msg)

    # Target directory in package structure
    target_path = UCLCHEM_ROOT_DIR / "data" / "collisional_rates"
    target_path.mkdir(parents=True, exist_ok=True)

    # Copy each .dat file
    logger.info(f"Copying {len(dat_files)} coolant data files to {target_path}")
    for dat_file in dat_files:
        target_file = target_path / dat_file.name
        shutil.copy2(dat_file, target_file)
        logger.debug(f"  Copied {dat_file.name}")

    logger.info("Successfully copied coolant data files for package installation")
