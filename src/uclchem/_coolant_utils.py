"""Utilities for loading and parsing coolant data files."""

import logging
import re
from pathlib import Path

from uclchem.constants import PLANCK_CONSTANT_CGS, SPEED_OF_LIGHT_CGS

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
    coolant_names: list[str],
    coolant_files: list[str],
    data_dir: str,
) -> tuple[int, int]:
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

    """
    data_path = Path(data_dir)
    if not data_path.exists():
        msg = f"Coolant data directory not found: {data_dir}"
        raise FileNotFoundError(msg)

    n_total_levels = 0
    target_marker = _normalize_for_comparison("NUMBER OF ENERGY LEVELS")

    for coolant_name, coolant_file in zip(coolant_names, coolant_files):
        filepath = data_path / coolant_file

        if not filepath.exists():
            msg = f"Coolant file not found: {filepath}"
            raise FileNotFoundError(msg)

        with filepath.open() as f:
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
            msg = f"No 'NUMBER OF ENERGY LEVELS' marker found in {filepath}"
            raise ValueError(msg)

    logger.info(
        f"Total energy levels across {len(coolant_names)} coolants: {n_total_levels}"
    )
    return n_total_levels, N_SE_STATS_PER_COOLANT


def get_energy_levels_info_from_runtime() -> tuple[int, int]:
    """Runtime wrapper that fetches parameters from uclchemwrap.

    Returns:
        Tuple of (n_total_levels, n_se_stats_per_coolant)

    """
    from uclchemwrap import f2py_constants

    from uclchem.advanced import HeatingSettings

    coolant_names = [str(name.decode()).strip() for name in f2py_constants.coolantnames]
    coolant_files = [str(fname.decode()).strip() for fname in f2py_constants.coolantfiles]

    # Get data directory
    heating_settings = HeatingSettings()
    data_dir = heating_settings.get_coolant_directory()

    return get_energy_levels_info(coolant_names, coolant_files, data_dir)


def validate_coolant_frequencies(
    coolant_names: list[str],
    coolant_files: list[str | Path],
    data_dir: str | Path,
) -> dict[str, float]:
    """Validate frequency consistency in LAMDA files at makerates time.

    For each coolant, computes freq = |E_i - E_j| / h from energy levels and
    compares to the frequency stored in the LAMDA file. Returns per-coolant
    maximum relative deviation.

    Args:
        coolant_names (list[str]): List of coolant species names
        coolant_files (list[str | Path]): List of coolant data file names
        data_dir (str | Path): Directory containing the coolant data files

    Returns:
        max_deviations (dict[str, float]): Dict mapping coolant name to max relative frequency deviation

    """
    data_path = Path(data_dir)
    level_marker = _normalize_for_comparison("NUMBER OF ENERGY LEVELS")
    trans_marker = _normalize_for_comparison("NUMBER OF RADIATIVE TRANSITIONS")

    max_deviations = {}

    for coolant_name, coolant_file in zip(coolant_names, coolant_files):
        filepath = data_path / coolant_file
        if not filepath.exists():
            continue

        with filepath.open() as f:
            lines = f.readlines()

        # Parse energy levels
        energies = {}  # level_index -> energy in erg
        nlevel = 0
        level_start = -1
        for j, line in enumerate(lines):
            if level_marker in _normalize_for_comparison(line):
                nlevel = int(lines[j + 1].strip())
                level_start = j + 2  # skip the comment line after count
                break

        if nlevel == 0:
            continue

        # Read energy levels (format: index energy_cm-1 weight [qnum])
        count = 0
        for k in range(level_start, len(lines)):
            line = lines[k].strip()
            if not line or line.startswith("!"):
                if count >= nlevel:
                    break
                continue
            parts = line.split()
            if len(parts) >= 3:
                idx = int(parts[0])
                energy_cm = float(parts[1])
                energies[idx] = (
                    energy_cm * SPEED_OF_LIGHT_CGS * PLANCK_CONSTANT_CGS
                )  # cm^-1 -> erg
                count += 1
                if count >= nlevel:
                    break

        # Parse radiative transitions
        max_dev = 0.0
        worst_transition = ""
        for j, line in enumerate(lines):
            if trans_marker in _normalize_for_comparison(line):
                ntrans = int(lines[j + 1].strip())
                trans_start = j + 2  # skip comment line
                count = 0
                for k in range(trans_start, len(lines)):
                    line = lines[k].strip()
                    if not line or line.startswith("!"):
                        if count >= ntrans:
                            break
                        continue
                    parts = line.split()
                    if len(parts) >= 5:
                        up = int(parts[1])
                        low = int(parts[2])
                        freq_ghz = float(parts[4])
                        freq_file = freq_ghz * 1.0e9  # GHz -> Hz

                        if up in energies and low in energies and freq_file > 0:
                            freq_calc = (
                                abs(energies[up] - energies[low]) / PLANCK_CONSTANT_CGS
                            )
                            rel_dev = abs(freq_calc - freq_file) / freq_file
                            if rel_dev > max_dev:
                                max_dev = rel_dev
                                worst_transition = f"{up}->{low}"
                        count += 1
                        if count >= ntrans:
                            break
                break

        max_deviations[coolant_name] = max_dev

        if max_dev > 0.01:  # Warn for > 1% deviation
            logger.warning(
                f"Coolant '{coolant_name}' ({coolant_file}): max frequency deviation "
                f"{max_dev:.4f} ({max_dev * 100:.2f}%) at transition {worst_transition}"
            )

    return max_deviations


def load_coolant_level_names() -> dict[int, list[str]]:
    """Load coolant level information from disk for meaningful column names.

    Returns:
        level_names (dict[int, list[str]]): Dict mapping coolant index to list of level names
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
        msg = f"Coolant data directory not found: {data_dir}"
        raise FileNotFoundError(msg)

    level_names = {}
    target_marker = _normalize_for_comparison("NUMBER OF ENERGY LEVELS")

    for i, (coolant_name, coolant_file) in enumerate(zip(coolant_names, coolant_files)):
        filepath = data_path / coolant_file

        if not filepath.exists():
            msg = f"Coolant file not found: {filepath}"
            raise FileNotFoundError(msg)

        level_names[i] = []
        with filepath.open() as f:
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
            msg = f"No energy levels found in {filepath}"
            raise ValueError(msg)

    # Validate all coolants were parsed
    if any(not names for names in level_names.values()):
        msg = "Some coolants have no parsed levels. Check coolant file format."
        raise RuntimeError(msg)

    logger.debug(f"Loaded level names for {len(level_names)} coolants")
    return level_names
