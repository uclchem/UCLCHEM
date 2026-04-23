#!/usr/bin/env python3
"""ATCT-UCLCHEM Interactive Species Matcher.

Interactive tool for matching UCLCHEM species with ATCT thermochemical data.
Handles exact matches, isomer detection, and user selection for ambiguous cases.
Supports save/resume functionality for interrupted sessions.
Outputs clean CSV files by default with YAML backup for compatibility.

Usage:
    from atct_uclchem_user_interactive import SpeciesMatcher

    matcher = SpeciesMatcher("atct_cleaned_v1.220.csv")
    mapping = matcher.match_species(uclchem_species_list)
    matcher.save_mapping(mapping, "species_mapping_v1.220.csv")  # Default CSV output
    matcher.save_mapping_yaml(mapping, "species_mapping_v1.220.yaml")  # Optional YAML
"""

import argparse
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def clean_numeric_value(value: Any) -> float | None:
    """Convert numpy scalars and other numeric types to clean Python floats.

    Args:
        value (Any): Numeric value that might be a numpy scalar, pandas Series, etc.

    Returns:
        float | None: Clean Python float or None if the value is NaN/missing

    """
    if pd.isna(value):
        return None

    # Handle numpy scalars and arrays
    if hasattr(value, "item"):
        clean_val = float(value.item())
    else:
        clean_val = float(value)

    # Return None for NaN values
    return None if np.isnan(clean_val) else clean_val


class FormulaParser:
    """Chemical formula parsing utilities for UCLCHEM and ATCT formats."""

    @staticmethod
    def parse_uclchem_formula(formula: str) -> dict[str, int]:
        """Parse UCLCHEM all-caps formulas (e.g., 'SIC2+', 'HE+').

        Args:
            formula (str): UCLCHEM-type formula

        Returns:
            elements (dict[str, int]): Dictionary with counter of elements

        """
        if not isinstance(formula, str):
            return {}

        clean_formula = formula.rstrip("+-")
        charge = 1 if formula.endswith("+") else -1 if formula.endswith("-") else 0

        # Known two-letter elements in UCLCHEM format
        two_letter_elements = ["HE", "SI", "CL", "MG", "NA"]

        elements = {"charge": charge}
        i = 0

        while i < len(clean_formula):
            if clean_formula[i].isalpha():
                # Check for two-letter element first
                found_element = None
                if i + 1 < len(clean_formula):
                    two_char = clean_formula[i : i + 2]
                    if two_char in two_letter_elements:
                        found_element = two_char
                        i += 2

                # If no two-letter match, take single character
                if found_element is None:
                    found_element = clean_formula[i]
                    i += 1

                # Parse count
                count_str = ""
                while i < len(clean_formula) and clean_formula[i].isdigit():
                    count_str += clean_formula[i]
                    i += 1

                count = int(count_str) if count_str else 1

                # Normalize element name to proper case
                element_map = {
                    "HE": "He",
                    "SI": "Si",
                    "CL": "Cl",
                    "MG": "Mg",
                    "NA": "Na",
                }
                normalized = element_map.get(found_element, found_element)

                if normalized in elements:
                    elements[normalized] += count
                else:
                    elements[normalized] = count
            else:
                i += 1

        return elements

    @staticmethod
    def parse_atct_formula(formula: str) -> dict[str, int]:
        """Parse ATCT proper-case formulas (e.g., 'SiC2+', 'He+').

        Args:
            formula (str): ATCT formula

        Returns:
            element_counts (dict[str, int]): Count of number of elements

        """
        if not isinstance(formula, str):
            return {}

        clean_formula = formula.rstrip("+-")
        charge = 1 if formula.endswith("+") else -1 if formula.endswith("-") else 0

        # Regular expression to match elements and their counts
        pattern = r"([A-Z][a-z]*)(\d*)"
        matches = re.findall(pattern, clean_formula)

        element_counts = {"charge": charge}

        for element, count_str in matches:
            count = int(count_str) if count_str else 1
            if element in element_counts:
                element_counts[element] += count
            else:
                element_counts[element] = count

        return element_counts

    @classmethod
    def formulas_match(cls, uclchem_formula: str, atct_formula: str) -> bool:
        """Check if UCLCHEM and ATCT formulas represent the same species.

        Args:
            uclchem_formula (str): UCLCHEM formula
            atct_formula (str): ATCT formula

        Returns:
            bool: whether the elemental composition of the two formulas is the same.

        """
        uclchem_parsed = cls.parse_uclchem_formula(uclchem_formula)
        atct_parsed = cls.parse_atct_formula(atct_formula)
        return uclchem_parsed == atct_parsed


class SpeciesMatcher:
    """Interactive species matching between UCLCHEM and ATCT databases."""

    def __init__(self, atct_csv_path: str):
        """Initialize matcher with ATCT data.

        Args:
            atct_csv_path (str): Path to cleaned ATCT CSV file

        """
        self.atct_data = pd.read_csv(atct_csv_path)
        self.atct_gas = self.atct_data[
            (self.atct_data["Phase"].str.startswith("g", na=False))
            & (self.atct_data["Phase"] != "graphite")
        ].copy()
        self.parser = FormulaParser()

        print(f"Loaded ATCT: {len(self.atct_gas)} gas-phase species")

    def match_species(
        self, uclchem_species: list[str], resume_file: str | Path | None = None
    ) -> dict[str, dict[str, Any]]:
        """Match UCLCHEM species with ATCT data.

        Args:
            uclchem_species (list[str]): List of UCLCHEM species names
            resume_file (str | Path | None): Optional path to resume from previous session.
                Default = None.

        Returns:
            dict[str, dict[str, Any]]: Dictionary mapping UCLCHEM species to ATCT matches

        """
        # Filter valid species
        target_species = [
            s
            for s in uclchem_species
            if pd.notna(s) and s not in {"NAN", ""} and not s.startswith(("#", "@"))
        ]

        # Resume from previous session if requested
        if resume_file and Path(resume_file).exists():
            return self._resume_matching(resume_file)

        # Stage 1: Find exact matches
        exact_matches = self._find_exact_matches(target_species)

        # Stage 2: Find isomer matches for unmatched species
        unmatched_species = [s for s in target_species if s not in exact_matches]
        isomer_matches = self._find_isomer_matches(unmatched_species)

        print(
            f"Matched {len(exact_matches)} exact, "
            f"{len(isomer_matches)} isomers from {len(target_species)} species"
        )

        # Stage 3: Interactive selection for multi-option species
        canonical_matches = exact_matches.copy()

        # Use default session file if none provided
        if resume_file is None and len(isomer_matches) > 0:
            resume_file = "species_matching_session.yaml"

        canonical_matches.update(self._interactive_selection(isomer_matches, resume_file))

        return canonical_matches

    def _find_exact_matches(self, species_list: list[str]) -> dict[str, dict[str, Any]]:
        """Find exact formula matches between UCLCHEM and ATCT.

        Args:
            species_list (list[str]): list of UCLCHEM species

        Returns:
            exact_matches (dict[str, dict[str, Any]]): mapping from UCLCHEM species
                to exact matches in ATCT.

        """
        exact_matches = {}

        for species in species_list:
            # Try direct formula matching (case-insensitive)
            formula_matches = self.atct_gas[
                self.atct_gas["Species_formula"].str.upper() == species.upper()
            ]

            if len(formula_matches) > 0:
                match = formula_matches.iloc[0]
                exact_matches[species] = {
                    "atct_name": match["Species_Name"],
                    "atct_formula": match["Species_formula"],
                    "enthalpy_298k": clean_numeric_value(match["DfH_298K"]),
                    "match_type": "exact",
                    "selection_reason": "exact_formula_match",
                }

        return exact_matches

    def _find_isomer_matches(
        self, species_list: list[str]
    ) -> dict[str, list[dict[str, Any]]]:
        """Find formula-based isomer matches.

        Args:
            species_list (list[str]): list of species to find all ATcT entries for.

        Returns:
            isomer_matches (dict[str, list[dict[str, Any]]]): mapping from UCLCHEM species
                to list of ATcT matches.

        """
        isomer_matches = {}

        for species in species_list:
            matches = []
            seen_formulas = set()  # Track seen formulas to avoid duplicates

            # Compare formulas, preferring simpler phase descriptions
            sorted_atct = self.atct_gas.sort_values(
                "Phase"
            )  # 'g' comes before 'g, ortho'

            for _, atct_row in sorted_atct.iterrows():
                if self.parser.formulas_match(species, atct_row["Species_formula"]):
                    formula_key = (
                        atct_row["Species_Name"],
                        atct_row["Species_formula"],
                    )

                    # Skip if we already have this exact name+formula combination
                    if formula_key not in seen_formulas:
                        matches.append(
                            {
                                "atct_name": atct_row["Species_Name"],
                                "atct_formula": atct_row["Species_formula"],
                                "enthalpy_298k": clean_numeric_value(
                                    atct_row["DfH_298K"]
                                ),
                            }
                        )
                        seen_formulas.add(formula_key)

            if matches:
                isomer_matches[species] = matches

        return isomer_matches

    def _interactive_selection(
        self,
        isomer_matches: dict[str, list[dict[str, Any]]],
        session_file: str | Path | None = None,
    ) -> dict[str, dict[str, Any]]:
        """Interactive selection for species with multiple matches.

        Args:
            isomer_matches (dict[str, list[dict[str, Any]]]): mapping from UCLCHEM species
                to list of ATcT matches.
            session_file (str | Path | None): path to session file. Default = None.

        Returns:
            canonical_matches (dict[str, dict[str, Any]]): mapping from UCLCHEM species
                to selection of the isomer matches.

        """
        canonical_matches = {}

        # Auto-select species with only one option
        single_option = {s: m for s, m in isomer_matches.items() if len(m) == 1}
        multi_option = {s: m for s, m in isomer_matches.items() if len(m) > 1}

        for species, matches in single_option.items():
            canonical_matches[species] = {
                **matches[0],
                "match_type": "isomer",
                "selection_reason": "auto_selected_only_option",
                "total_options": 1,
            }

        if not multi_option:
            return canonical_matches

        print(f"\nInteractive selection: {len(multi_option)} species")
        print("Enter choice number, 's' to skip, 'q' to quit and save")
        print("=" * 60)

        # Check if user wants to continue from previous session or start fresh
        if session_file is not None and Path(session_file).exists():
            session_file = Path(session_file)
            while True:
                try:
                    continue_choice = (
                        input(
                            f"\nFound existing session file: {session_file}\n"
                            f"Would you like to:\n"
                            f"  'c' - Continue from where you left off\n"
                            f"  'n' - Start fresh (previous work backed up)\n"
                            f"  'q' - Quit without making changes\n"
                            f"Your choice: "
                        )
                        .strip()
                        .lower()
                    )

                    if continue_choice == "q":
                        print("Exiting without changes.")
                        return canonical_matches
                    elif continue_choice == "c":
                        # Load previous session
                        try:
                            with session_file.open() as f:
                                previous_session = yaml.safe_load(f)
                            prev = previous_session.get("completed", {})
                            canonical_matches.update(prev)
                            print(f"Loaded {len(prev)} previous selections")
                            multi_option = {
                                s: m for s, m in multi_option.items() if s not in prev
                            }
                            if not multi_option:
                                print("All selections completed")
                                return canonical_matches
                            break
                        except Exception as e:
                            logger.warning(f"Error loading session: {e}")
                            break
                    elif continue_choice == "n":
                        # Backup existing session instead of deleting it
                        backup_name = self._backup_session_file(session_file)
                        print(f"Session backed up: {backup_name}")
                        break
                    else:
                        print("Invalid choice. Please enter 'c', 'n', or 'q'")
                except KeyboardInterrupt:
                    print("\nExiting...")
                    return canonical_matches

        session_data = {"completed": canonical_matches, "remaining": multi_option}

        for i, (species, matches) in enumerate(multi_option.items(), 1):
            print(f"\n[{i}/{len(multi_option)}] {species}")
            print(f"Found {len(matches)} isomer matches:")

            # Display options
            for j, match in enumerate(matches, 1):
                enthalpy_str = (
                    f"{match['enthalpy_298k']:.1f} kJ/mol"
                    if match["enthalpy_298k"] is not None
                    else "No enthalpy data"
                )
                print(
                    f"   {j}. {match['atct_name']:35} "
                    f"({match['atct_formula']:10}) | "
                    f"ΔfH°: {enthalpy_str:>12}"
                )

            # Get user selection
            while True:
                try:
                    user_input = (
                        input(f"\nChoice (1-{len(matches)}, s, q): ").strip().lower()
                    )

                    if user_input == "q":
                        print("Quitting")
                        if session_file:
                            self._save_session(session_data, session_file)
                            print("Session saved")
                        return canonical_matches

                    if user_input == "s":
                        print(f"Skipped: {species}")
                        break

                    choice = int(user_input)
                    if 1 <= choice <= len(matches):
                        selected = matches[choice - 1]
                        canonical_matches[species] = {
                            **selected,
                            "match_type": "isomer",
                            "selection_reason": "user_selected_option",
                            "total_options": len(matches),
                        }
                        print(f"Selected: {selected['atct_name']}")
                        break
                    else:
                        print(f"Invalid: use 1-{len(matches)}, s, or q")

                except ValueError:
                    print(f"Invalid: use 1-{len(matches)}, s, or q")
                except KeyboardInterrupt:
                    print("\nSession interrupted. Saving progress...")
                    if session_file:
                        self._save_session(session_data, session_file)
                    return canonical_matches

            # Save progress after each selection
            if session_file:
                session_data["completed"] = canonical_matches
                remaining_multi = {
                    s: m for s, m in multi_option.items() if s not in canonical_matches
                }
                session_data["remaining"] = remaining_multi
                self._save_session(session_data, session_file)

        return canonical_matches

    @staticmethod
    def _backup_session_file(session_file: str | Path) -> str:
        """Create a simple backup of existing session file.

        Args:
            session_file (str | Path): path to session file to back up.

        Returns:
            str: String representing backup. If the backup failed,
                is "backup_failed{extension}", otherwise is
                "{basename}_backup{extension}".

        """
        base_name = Path(session_file).stem
        extension = Path(session_file).suffix
        backup_name = f"{base_name}_backup{extension}"

        try:
            shutil.copy2(session_file, backup_name)
            return backup_name
        except Exception as e:
            logger.warning(f"Could not backup session file: {e}")
            return f"backup_failed{extension}"

    @staticmethod
    def _save_session(session_data: dict, session_file: str | Path) -> None:
        """Save current matching session.

        Args:
            session_data (dict): data of the current session
            session_file (str | Path): where to save the current session

        """
        try:
            with Path(session_file).open("w") as f:
                yaml.dump(session_data, f, indent=2, default_flow_style=False)

            completed = len(session_data.get("completed", {}))
            remaining = len(session_data.get("remaining", {}))
            print(f"Progress: {completed} done, {remaining} left")
        except Exception as e:
            logger.warning(f"Warning: Could not save session: {e}")

    def _resume_matching(self, resume_file: str | Path) -> dict[str, dict[str, Any]]:
        """Resume matching from saved session.

        Args:
            resume_file (str | Path): Path to saved session file.

        Returns:
            completed (dict[str, dict[str, Any]]): species mapping dictionary.

        """
        with Path(resume_file).open() as f:
            session_data = yaml.safe_load(f)

        completed = session_data.get("completed", {})
        remaining = session_data.get("remaining", {})

        print(f"Resumed: {len(completed)} done, {len(remaining)} left")

        if remaining:
            completed.update(self._interactive_selection(remaining, resume_file))

        return completed

    @staticmethod
    def save_mapping_yaml(
        mapping: dict[str, dict[str, Any]], output_path: str | Path
    ) -> None:
        """Save species mapping to YAML file.

        Args:
            mapping (dict[str, dict[str, Any]]): Species mapping dictionary
            output_path (str | Path): Output YAML file path

        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with output_path.open("w") as f:
            yaml.dump(mapping, f, indent=2, default_flow_style=False)

        print(f"Saved YAML: {len(mapping)} species → {output_path}")

    def save_mapping(
        self, mapping: dict[str, dict[str, Any]], output_path: str | Path
    ) -> None:
        """Save species mapping to CSV file (default format).

        Args:
            mapping (dict[str, dict[str, Any]]): Species mapping dictionary
            output_path (str | Path): Output CSV file path

        """
        self.export_mapping_csv(mapping, output_path)

    @staticmethod
    def load_mapping(mapping_path: str | Path) -> dict[str, dict[str, Any]]:
        """Load species mapping from YAML file.

        Args:
            mapping_path (str | Path): Path to mapping YAML file

        Returns:
            mapping (dict[str, dict[str, Any]]): Species mapping dictionary

        """
        with Path(mapping_path).open() as f:
            mapping = yaml.safe_load(f)

        print(f"Loaded: {len(mapping)} species from {mapping_path}")
        return mapping

    @staticmethod
    def export_mapping_csv(
        mapping: dict[str, dict[str, Any]], output_path: str | Path
    ) -> None:
        """Export mapping to CSV for easy review.

        Args:
            mapping (dict[str, dict[str, Any]]): Species mapping dictionary
            output_path (str | Path): Output CSV file path

        """
        mapping_data = []
        for species, match_info in mapping.items():
            mapping_data.append(
                {
                    "uclchem_species": species,
                    "atct_name": match_info["atct_name"],
                    "atct_formula": match_info["atct_formula"],
                    "enthalpy_298k": match_info["enthalpy_298k"],
                    "match_type": match_info["match_type"],
                    "selection_reason": match_info["selection_reason"],
                    "has_enthalpy": "YES"
                    if match_info["enthalpy_298k"] is not None
                    else "NO",
                }
            )

        df = pd.DataFrame(mapping_data)
        df = df.sort_values("uclchem_species").reset_index(drop=True)
        df.to_csv(output_path, index=False)

        has_enthalpy = len(df[df["has_enthalpy"] == "YES"])
        print(f"Exported CSV: {len(df)} species ({has_enthalpy} with data)")

    @staticmethod
    def write_back_to_file(
        species_df: pd.DataFrame,
        mapping: dict[str, dict[str, Any]],
        original_csv_path: str | Path,
    ) -> None:
        """Write back the original species DataFrame with enthalpy data.

        Merges ATCT enthalpy data into the existing ENTHALPY column without
        overwriting any existing non-zero values. Creates a backup of the
        original file as NAME_backup.csv.

        Args:
            species_df (pd.DataFrame): Original UCLCHEM species DataFrame
            mapping (dict[str, dict[str, Any]]): Species mapping dictionary from match_species()
            original_csv_path (str | Path): Path to the original species CSV file

        """
        # Create backup file path
        original_path = Path(original_csv_path)
        backup_path = original_path.with_name(
            f"{original_path.stem}_backup{original_path.suffix}"
        )

        # Create backup
        try:
            shutil.copy2(original_csv_path, backup_path)
        except Exception as e:
            logger.warning(f"Backup failed: {e}")

        # Create a copy of the DataFrame to avoid modifying the original
        updated_df = species_df.copy()

        # Ensure ENTHALPY column exists
        if "ENTHALPY" not in updated_df.columns:
            updated_df["ENTHALPY"] = 0.0

        # Track statistics
        existing_nonzero_count = 0
        new_values_added = 0
        skipped_existing = 0

        # Merge enthalpy values into existing ENTHALPY column
        for idx, row in updated_df.iterrows():
            species_name = row["NAME"]
            current_enthalpy = row.get("ENTHALPY", 0.0)

            # Check if current value is non-zero (already has data)
            has_existing_value = (
                pd.notna(current_enthalpy) and abs(current_enthalpy) > 1e-10
            )

            if has_existing_value:
                existing_nonzero_count += 1
                # Never overwrite existing non-zero values
                if (
                    species_name in mapping
                    and mapping[species_name]["enthalpy_298k"] is not None
                ):
                    skipped_existing += 1
                continue

            # Fill in missing/zero values from ATCT mapping
            if species_name in mapping:
                atct_enthalpy = mapping[species_name]["enthalpy_298k"]
                if atct_enthalpy is not None:
                    # Convert kJ/mol to kcal/mol
                    updated_df.at[idx, "ENTHALPY"] = atct_enthalpy / 4.184
                    new_values_added += 1

        # Write back to original file
        try:
            updated_df.to_csv(original_csv_path, index=False)
            total_with_enthalpy = existing_nonzero_count + new_values_added
            print(
                f"Updated {original_csv_path}: "
                f"{new_values_added} added, "
                f"{existing_nonzero_count} preserved, "
                f"{total_with_enthalpy}/{len(updated_df)} total"
            )
        except Exception as e:
            logger.warning(f"Write failed: {e}")
            # Restore backup if write failed
            try:
                shutil.copy2(backup_path, original_csv_path)
                print("Restored from backup")
            except Exception as restore_error:
                logger.warning(f"Restore failed: {restore_error}")


def main() -> None:
    """Command-line interface for species matcher."""
    # Configure logging
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Interactive tool for matching UCLCHEM species with "
        "ATCT thermochemical data."
    )

    parser.add_argument(
        "--atct_csv",
        required=True,
        help=("Path to the cleaned ATCT CSV file containing thermochemical data"),
    )

    parser.add_argument(
        "--uclchem_species_csv",
        required=True,
        help=("Path to CSV file containing UCLCHEM species (must have 'NAME' column)"),
    )

    parser.add_argument(
        "--output_mapping",
        help="Path for full output CSV mapping file.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--overwrite_uclchem_species_csv",
        action="store_true",
        help=(
            "If set, enhance the original UCLCHEM species CSV with "
            "enthalpy data based on the matching."
        ),
    )

    args = parser.parse_args()

    try:
        # Load UCLCHEM species
        species_df = pd.read_csv(args.uclchem_species_csv)
        species_list = species_df["NAME"].tolist()

        # Run matching
        matcher = SpeciesMatcher(args.atct_csv)
        mapping = matcher.match_species(species_list)

        # Save results (default to CSV)
        if args.output_mapping is not None:
            matcher.save_mapping(mapping, args.output_mapping)
        if args.overwrite_uclchem_species_csv:
            matcher.write_back_to_file(species_df, mapping, args.uclchem_species_csv)

        print(f"Complete: {len(mapping)} species mapped")

    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
