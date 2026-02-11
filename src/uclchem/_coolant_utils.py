"""Utilities for loading and parsing coolant data files."""

import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# Number of statistical evolution statistics per coolant
# (converged flag, iterations, max_rel_change)
N_SE_STATS_PER_COOLANT = 3


def _normalize_for_comparison(text: str) -> str:
    """Normalize text for robust comparison.

    Strips whitespace, converts to lowercase, removes special characters.

    Args:
        text: Text to normalize

    Returns:
        Normalized text with only alphanumeric characters
    """
    return re.sub(r"[^a-z0-9]", "", text.strip().lower())


def get_energy_levels_info(
    coolant_names: List[str],
    coolant_files: List[str],
    data_dir: str,
) -> Tuple[int, int]:
    """Compute energy level information from coolant data files.

    This is the core function that works both at makerates time (without uclchemwrap)
    and at runtime (via wrapper function). No fail-safes - raises errors if anything
    goes wrong.

    Args:
        coolant_names: List of coolant species names (e.g., ['H', 'C+', 'O', ...])
        coolant_files: List of coolant data file names (e.g., ['ly-a.dat', ...])
        data_dir: Directory containing the coolant data files

    Returns:
        Tuple of (n_total_levels, n_se_stats_per_coolant)
        where n_total_levels is the sum of all energy levels across all coolants
        and n_se_stats_per_coolant is always 3 (converged, iterations, max_rel_change)

    Raises:
        FileNotFoundError: If data directory or any coolant file doesn't exist
        ValueError: If coolant file format is invalid or energy levels not found
        RuntimeError: If any parsing error occurs
    """
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileNotFoundError(f"Coolant data directory not found: {data_dir}")

    n_total_levels = 0
    target_marker = _normalize_for_comparison("NUMBER OF ENERGY LEVELS")

    for coolant_name, coolant_file in zip(coolant_names, coolant_files):
        filepath = data_path / coolant_file

        if not filepath.exists():
            raise FileNotFoundError(f"Coolant file not found: {filepath}")

        with open(filepath, "r") as f:
            lines = f.readlines()

        # Find NUMBER OF ENERGY LEVELS line (robust search)
        found_levels = False
        for j, line in enumerate(lines):
            if target_marker in _normalize_for_comparison(line):
                found_levels = True
                nlevel = int(lines[j + 1].strip())
                n_total_levels += nlevel
                logger.debug(f"Coolant {coolant_name}: {nlevel} levels")
                break

        if not found_levels:
            raise ValueError(f"No 'NUMBER OF ENERGY LEVELS' marker found in {filepath}")

    logger.info(
        f"Total energy levels across {len(coolant_names)} coolants: {n_total_levels}"
    )
    return n_total_levels, N_SE_STATS_PER_COOLANT


def get_energy_levels_info_from_runtime() -> Tuple[int, int]:
    """Runtime wrapper that fetches parameters from uclchemwrap.

    Returns:
        Tuple of (n_total_levels, n_se_stats_per_coolant)

    Raises:
        ImportError: If uclchemwrap is not available
        Various exceptions from get_energy_levels_info if data is invalid
    """
    from uclchemwrap import f2py_constants

    from uclchem.advanced import HeatingSettings

    coolant_names = [str(name.decode()).strip() for name in f2py_constants.coolantnames]
    coolant_files = [str(fname.decode()).strip() for fname in f2py_constants.coolantfiles]

    # Get data directory
    heating_settings = HeatingSettings()
    data_dir = heating_settings.get_coolant_directory()

    return get_energy_levels_info(coolant_names, coolant_files, data_dir)


def load_coolant_level_names() -> Dict[int, List[str]]:
    """Load coolant level information from disk for meaningful column names.

    Returns:
        Dict mapping coolant index to list of level names
        (e.g., {0: ['H_2_S_1/2', 'H_2_P_1/2'], 1: ['C+_2_P_1/2', ...]})

    Raises:
        FileNotFoundError: If data directory or coolant files don't exist
        ValueError: If coolant file format is invalid
        RuntimeError: If parsing fails
    """
    from uclchemwrap import f2py_constants

    from uclchem.advanced import HeatingSettings

    coolant_names = [str(name.decode()).strip() for name in f2py_constants.coolantnames]
    coolant_files = [str(fname.decode()).strip() for fname in f2py_constants.coolantfiles]

    # Get data directory
    heating_settings = HeatingSettings()
    data_dir = heating_settings.get_coolant_directory()

    # Verify directory exists
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileNotFoundError(f"Coolant data directory not found: {data_dir}")

    level_names = {}
    target_marker = _normalize_for_comparison("NUMBER OF ENERGY LEVELS")

    for i, (coolant_name, coolant_file) in enumerate(zip(coolant_names, coolant_files)):
        filepath = data_path / coolant_file

        if not filepath.exists():
            raise FileNotFoundError(f"Coolant file not found: {filepath}")

        level_names[i] = []
        with open(filepath, "r") as f:
            lines = f.readlines()

        # Find NUMBER OF ENERGY LEVELS line (robust search)
        found_levels = False
        for j, line in enumerate(lines):
            if target_marker in _normalize_for_comparison(line):
                found_levels = True
                nlevel = int(lines[j + 1].strip())
                # Read level info - parse quantum numbers
                level_count = 0
                for k in range(j + 2, len(lines)):
                    level_line = lines[k].strip()
                    if not level_line or level_line.startswith("!"):
                        if level_count >= nlevel:
                            break
                        continue
                    # Format: level_num energy weight qnum
                    parts = level_line.split()
                    if len(parts) >= 4:
                        qnum = parts[3]  # quantum number
                        level_names[i].append(f"{coolant_name}_{qnum}")
                        level_count += 1
                        if level_count >= nlevel:
                            break
                break

        if not found_levels:
            raise ValueError(f"No energy levels found in {filepath}")

    # Validate all coolants were parsed
    if any(not names for names in level_names.values()):
        raise RuntimeError(
            "Some coolants have no parsed levels. Check coolant file format."
        )

    logger.debug(f"Loaded level names for {len(level_names)} coolants")
    return level_names
