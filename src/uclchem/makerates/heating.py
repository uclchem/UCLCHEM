"""
Heating and cooling calculations for UCLCHEM reactions.

Provides functions to set reaction exothermicities from thermochemical
databases or custom CSV files with various units.
"""

import logging
import re
from typing import List, Optional, Tuple

import pandas as pd

from .reaction import Reaction

# Physical constants (2019 SI definitions)
AVOGADRO_NUMBER = 6.02214076e23  # mol^-1
EV_TO_JOULE = 1.602176634e-19  # J/eV
CALORIE_TO_JOULE = 4.184  # J/cal (thermochemical)
ERG_TO_JOULE = 1.0e-7  # J/erg

# Base conversion factors to erg per reaction
EV_TO_ERG = EV_TO_JOULE / ERG_TO_JOULE
JOULE_TO_ERG = 1.0 / ERG_TO_JOULE
KCAL_TO_ERG = CALORIE_TO_JOULE * 1000.0 / ERG_TO_JOULE
ERG_TO_ERG = 1.0

# Cache for parsed unit conversions
_UNIT_CACHE = {}

# Base unit mappings
_BASE_UNITS = {
    "ev": EV_TO_ERG,
    "joule": JOULE_TO_ERG,
    "j": JOULE_TO_ERG,
    "kj": JOULE_TO_ERG * 1e3,
    "cal": CALORIE_TO_JOULE,
    "kcal": KCAL_TO_ERG,
    "erg": ERG_TO_ERG,
}

# Denominator mappings
_DENOMINATORS = {
    "reaction": 1.0,
    "mol": 1.0 / AVOGADRO_NUMBER,
}


def parse_species_from_row(row: pd.Series, prefix: str) -> List[str]:
    """Parse species list from CSV row.

    Args:
        row: DataFrame row
        prefix: 'reactant' or 'product'

    Returns:
        List of species names (uppercase, NAN for missing)
    """
    species = []
    idx = 1
    while f"{prefix}{idx}" in row.index:
        val = row[f"{prefix}{idx}"]
        if pd.isna(val) or val == "" or str(val).upper().strip() == "NAN":
            species.append("NAN")
        else:
            species.append(str(val).strip().upper())
        idx += 1
    return species if species else ["NAN"]


def _parse_unit(unit: str) -> float:
    """Parse unit string and return conversion factor to erg per reaction.

    Parses units like: ev, ev_per_reaction, ev/mol, joule_per_mol, etc.

    Args:
        unit: Unit string (case-insensitive)

    Returns:
        Conversion factor to erg per reaction
    """
    unit_lower = unit.strip().lower()

    # Check cache first
    if unit_lower in _UNIT_CACHE:
        return _UNIT_CACHE[unit_lower]

    # Split by separators (/ or _per_)
    # Replace _per_ with / for uniform handling
    normalized = re.sub(r"_per_", "/", unit_lower)
    parts = normalized.split("/")
    if len(parts) == 1:
        # Just a base unit (e.g., "ev", "joule")
        base = parts[0].strip()
        if base not in _BASE_UNITS:
            raise ValueError(
                f"Unknown unit '{unit}'. Available: {list(_BASE_UNITS.keys())}"
            )
        factor = _BASE_UNITS[base]
        # Default to per reaction
        factor *= _DENOMINATORS["reaction"]
    elif len(parts) == 2:
        # Unit with denominator (e.g., "ev/reaction", "kcal/mol")
        base = parts[0].strip()
        denom = parts[1].strip()
        if base not in _BASE_UNITS:
            raise ValueError(
                f"Unknown base unit '{base}' in '{unit}'. "
                f"Available: {list(_BASE_UNITS.keys())}"
            )
        if denom not in _DENOMINATORS:
            raise ValueError(
                f"Unknown denominator '{denom}' in '{unit}'. "
                f"Available: {list(_DENOMINATORS.keys())}"
            )
        factor = _BASE_UNITS[base] * _DENOMINATORS[denom]
    else:
        raise ValueError(
            f"Cannot parse unit '{unit}'. Format: <unit> or <unit>/<denominator> or <unit>_per_<denominator>"
        )
    # Cache the result
    _UNIT_CACHE[unit_lower] = factor
    return factor


def convert_to_erg(value: float, unit: str) -> float:
    """Convert exothermicity to erg per reaction.
    Args:
        value: Exothermicity value
        unit: Unit string (case-insensitive)

    Returns:
        Value in erg per reaction
    """
    factor = _parse_unit(unit)
    return value * factor


def match_reaction(
    reactants: List[str], products: List[str], reactions: List[Reaction]
) -> Optional[Reaction]:
    """Find matching reaction in list.

    Args:
        reactants: List of reactant names
        products: List of product names
        reactions: List to search

    Returns:
        Matching Reaction or None
    """
    sorted_r = sorted(reactants)
    sorted_p = sorted(products)

    for reaction in reactions:
        if (
            sorted(reaction.get_reactants()) == sorted_r
            and sorted(reaction.get_products()) == sorted_p
        ):
            return reaction
    return None


def load_custom_exothermicities(csv_path: str) -> pd.DataFrame:
    """Load custom exothermicities from CSV.

    Expected columns: reactant1-3, product1-4, exothermicity, unit

    Args:
        csv_path: Path to CSV file

    Returns:
        DataFrame with custom exothermicities
    """
    df = pd.read_csv(csv_path, comment="#")

    required = ["exothermicity", "unit"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")

    has_reactants = any(c.startswith("reactant") for c in df.columns)
    has_products = any(c.startswith("product") for c in df.columns)

    if not has_reactants or not has_products:
        raise ValueError("CSV must have reactant and product columns")

    return df


def set_custom_exothermicities(
    reactions: List[Reaction], csv_path: str, overwrite: bool = True
) -> Tuple[int, int]:
    """Set reaction exothermicities from custom CSV.

    Args:
        reactions: List of Reaction objects to modify
        csv_path: Path to CSV with custom exothermicities
        overwrite: If False, only set reactions with zero exothermicity

    Returns:
        tuple: (num_matched, num_unmatched)
    """
    df = load_custom_exothermicities(csv_path)
    matched = 0
    unmatched = 0

    for _, row in df.iterrows():
        reactants = parse_species_from_row(row, "reactant")
        products = parse_species_from_row(row, "product")

        try:
            exo_erg = convert_to_erg(row["exothermicity"], row["unit"])
        except ValueError as e:
            logging.warning(f"Skipping row: {e}")
            continue

        reaction = match_reaction(reactants, products, reactions)

        if reaction:
            if overwrite or reaction.get_exothermicity() == 0.0:
                reaction.set_exothermicity(exo_erg)
                matched += 1
        else:
            unmatched += 1
            logging.warning(
                f"No match: {' + '.join(r for r in reactants if r != 'NAN')} "
                f"-> {' + '.join(p for p in products if p != 'NAN')}"
            )

    logging.info(f"Custom exothermicities: {matched} matched, {unmatched} unmatched")
    return matched, unmatched
