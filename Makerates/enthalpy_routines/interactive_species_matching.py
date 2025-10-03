#!/usr/bin/env python3
"""
ATCT-UCLCHEM Interactive Species Matcher

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

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import yaml


def clean_numeric_value(value: Any) -> Optional[float]:
    """Convert numpy scalars and other numeric types to clean Python floats.

    Args:
        value: Numeric value that might be a numpy scalar, pandas Series, etc.

    Returns:
        Clean Python float or None if the value is NaN/missing
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
    def parse_uclchem_formula(formula: str) -> Dict[str, int]:
        """Parse UCLCHEM all-caps formulas (e.g., 'SIC2+', 'HE+')."""
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
    def parse_atct_formula(formula: str) -> Dict[str, int]:
        """Parse ATCT proper-case formulas (e.g., 'SiC2+', 'He+')."""
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
        """Check if UCLCHEM and ATCT formulas represent the same species."""
        uclchem_parsed = cls.parse_uclchem_formula(uclchem_formula)
        atct_parsed = cls.parse_atct_formula(atct_formula)
        return uclchem_parsed == atct_parsed


class SpeciesMatcher:
    """Interactive species matching between UCLCHEM and ATCT databases."""

    def __init__(self, atct_csv_path: str):
        """Initialize matcher with ATCT data.

        Args:
            atct_csv_path: Path to cleaned ATCT CSV file
        """
        self.atct_data = pd.read_csv(atct_csv_path)
        self.atct_gas = self.atct_data[
            (self.atct_data["Phase"].str.startswith("g", na=False))
            & (self.atct_data["Phase"] != "graphite")
        ].copy()
        self.parser = FormulaParser()

        print(f"✓ Loaded ATCT data: {len(self.atct_data)} total species")
        print(f"✓ Gas phase species: {len(self.atct_gas)} species")

    def match_species(
        self, uclchem_species: List[str], resume_file: Optional[str] = None
    ) -> Dict[str, Dict[str, Any]]:
        """Match UCLCHEM species with ATCT data.

        Args:
            uclchem_species: List of UCLCHEM species names
            resume_file: Optional path to resume from previous session

        Returns:
            Dictionary mapping UCLCHEM species to ATCT matches
        """
        # Filter valid species
        target_species = [
            s
            for s in uclchem_species
            if pd.notna(s) and s not in ["NAN", ""] and not s.startswith(("#", "@"))
        ]

        print(f"Matching {len(target_species)} UCLCHEM species...")

        # Resume from previous session if requested
        if resume_file and Path(resume_file).exists():
            return self._resume_matching(target_species, resume_file)

        # Stage 1: Find exact matches
        exact_matches = self._find_exact_matches(target_species)
        print(f"✓ Found {len(exact_matches)} exact matches")

        # Stage 2: Find isomer matches for unmatched species
        unmatched_species = [s for s in target_species if s not in exact_matches]
        isomer_matches = self._find_isomer_matches(unmatched_species)
        print(f"✓ Found isomer matches for {len(isomer_matches)} species")

        # Stage 3: Interactive selection for multi-option species
        canonical_matches = exact_matches.copy()

        # Use default session file if none provided
        if resume_file is None and len(isomer_matches) > 0:
            resume_file = "species_matching_session.yaml"

        canonical_matches.update(
            self._interactive_selection(isomer_matches, resume_file)
        )

        return canonical_matches

    def _find_exact_matches(self, species_list: List[str]) -> Dict[str, Dict[str, Any]]:
        """Find exact formula matches between UCLCHEM and ATCT."""
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
        self, species_list: List[str]
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Find formula-based isomer matches."""
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
        isomer_matches: Dict[str, List[Dict[str, Any]]],
        session_file: Optional[str] = None,
    ) -> Dict[str, Dict[str, Any]]:
        """Interactive selection for species with multiple matches."""
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

        print(f"✓ Auto-selected {len(single_option)} species with single option")

        if not multi_option:
            return canonical_matches

        print(f"\n Interactive selection needed for {len(multi_option)} species...")
        print("=" * 60)
        print("  INSTRUCTIONS:")
        print("   • Progress is automatically saved after each selection")
        print("   • You can quit anytime with 'q' and resume later")
        print("   • Enter numbers (1,2,3...) to select species")
        print("   • Enter 's' to skip a species")
        print("=" * 60)

        # Check if user wants to continue from previous session or start fresh
        if session_file and Path(session_file).exists():
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
                            with open(session_file, "r") as f:
                                previous_session = yaml.safe_load(f)
                            previous_completed = previous_session.get("completed", {})
                            canonical_matches.update(previous_completed)
                            print(
                                f"✓ Loaded {len(previous_completed)} previous selections"
                            )
                            # Remove already completed species from multi_option
                            multi_option = {
                                s: m
                                for s, m in multi_option.items()
                                if s not in previous_completed
                            }
                            if not multi_option:
                                print("✓ All selections already completed!")
                                return canonical_matches
                            break
                        except Exception as e:
                            print(f"Error loading session: {e}")
                            print("Starting fresh...")
                            break
                    elif continue_choice == "n":
                        # Backup existing session instead of deleting it
                        backup_name = self._backup_session_file(session_file)
                        print(f"✓ Previous session backed up as: {backup_name}")
                        print("Starting fresh session...")
                        break
                    else:
                        print("Invalid choice. Please enter 'c', 'n', or 'q'")
                except KeyboardInterrupt:
                    print("\nExiting...")
                    return canonical_matches

        session_data = {"completed": canonical_matches, "remaining": multi_option}

        for i, (species, matches) in enumerate(multi_option.items(), 1):
            print(f"\n[{i}/{len(multi_option)}] UCLCHEM Species: {species}")
            print(f"Chemical formula: {self.parser.parse_uclchem_formula(species)}")
            print(f"Found {len(matches)} isomer matches in ATCT:")

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
                        input(
                            f"\nEnter choice (1-{len(matches)}), 's' to skip, 'q' to quit: "
                        )
                        .strip()
                        .lower()
                    )

                    if user_input == "q":
                        print("Session saved. You can resume later.")
                        if session_file:
                            self._save_session(session_data, session_file)
                        return canonical_matches

                    if user_input == "s":
                        print(f"→ SKIPPED: {species}")
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
                        print(f"→ SELECTED: {selected['atct_name']} (option {choice})")
                        break
                    else:
                        print(f"Invalid choice. Enter 1-{len(matches)}, 's', or 'q'")

                except ValueError:
                    print(f"Invalid input. Enter 1-{len(matches)}, 's', or 'q'")
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

    def _backup_session_file(self, session_file: str) -> str:
        """Create a simple backup of existing session file."""
        import shutil

        base_name = Path(session_file).stem
        extension = Path(session_file).suffix
        backup_name = f"{base_name}_backup{extension}"

        try:
            shutil.copy2(session_file, backup_name)
            return backup_name
        except Exception as e:
            print(f"⚠️  Could not backup session file: {e}")
            return f"backup_failed{extension}"

    def _save_session(self, session_data: Dict, session_file: str) -> None:
        """Save current matching session."""
        try:
            with open(session_file, "w") as f:
                yaml.dump(session_data, f, indent=2, default_flow_style=False)

            completed_count = len(session_data.get("completed", {}))
            remaining_count = len(session_data.get("remaining", {}))
            print(
                f"✓ Progress saved: {completed_count} completed, "
                f"{remaining_count} remaining → {session_file}"
            )
        except Exception as e:
            logging.warning(f"Warning: Could not save session: {e}")

    def _resume_matching(
        self, target_species: List[str], resume_file: str
    ) -> Dict[str, Dict[str, Any]]:
        """Resume matching from saved session."""
        with open(resume_file, "r") as f:
            session_data = yaml.safe_load(f)

        completed = session_data.get("completed", {})
        remaining = session_data.get("remaining", {})

        print(
            f"✓ Resumed session: {len(completed)} completed, {len(remaining)} remaining"
        )

        if remaining:
            completed.update(self._interactive_selection(remaining, resume_file))

        return completed

    def save_mapping_yaml(
        self, mapping: Dict[str, Dict[str, Any]], output_path: str
    ) -> None:
        """Save species mapping to YAML file.

        Args:
            mapping: Species mapping dictionary
            output_path: Output YAML file path
        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, "w") as f:
            yaml.dump(mapping, f, indent=2, default_flow_style=False)

        print(f"✓ Saved YAML mapping for {len(mapping)} species to {output_path}")

    def save_mapping(
        self, mapping: Dict[str, Dict[str, Any]], output_path: str
    ) -> None:
        """Save species mapping to CSV file (default format).

        Args:
            mapping: Species mapping dictionary
            output_path: Output CSV file path
        """
        self.export_mapping_csv(mapping, output_path)

    def load_mapping(self, mapping_path: str) -> Dict[str, Dict[str, Any]]:
        """Load species mapping from YAML file.

        Args:
            mapping_path: Path to mapping YAML file

        Returns:
            Species mapping dictionary
        """
        with open(mapping_path, "r") as f:
            mapping = yaml.safe_load(f)

        print(f"✓ Loaded mapping for {len(mapping)} species from {mapping_path}")
        return mapping

    def export_mapping_csv(
        self, mapping: Dict[str, Dict[str, Any]], output_path: str
    ) -> None:
        """Export mapping to CSV for easy review.

        Args:
            mapping: Species mapping dictionary
            output_path: Output CSV file path
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

        print(f"✓ Exported mapping CSV to {output_path}")
        print(f"  - {len(df)} total mappings")
        print(f"  - {len(df[df['has_enthalpy'] == 'YES'])} with enthalpy data")

    def write_back_to_file(
        self,
        species_df: pd.DataFrame,
        mapping: Dict[str, Dict[str, Any]],
        original_csv_path: str
    ) -> None:
        """Write back the original species DataFrame with enthalpy data.

        Creates a backup of the original file as NAME_backup.csv and writes
        the updated DataFrame with enthalpy column to the original file.

        Args:
            species_df: Original UCLCHEM species DataFrame
            mapping: Species mapping dictionary from match_species()
            original_csv_path: Path to the original species CSV file
        """
        import shutil
        from pathlib import Path

        # Create backup file path
        original_path = Path(original_csv_path)
        backup_path = original_path.with_name(
            f"{original_path.stem}_backup{original_path.suffix}"
        )

        # Create backup
        try:
            shutil.copy2(original_csv_path, backup_path)
            print(f"✓ Created backup: {backup_path}")
        except Exception as e:
            print(f"⚠️  Could not create backup: {e}")

        # Create a copy of the DataFrame to avoid modifying the original
        updated_df = species_df.copy()

        # Check if gas_enthalpy column already exists with non-zero values
        if "GAS_ENTHALPY" in updated_df.columns:
            existing_nonzero = updated_df[
                (updated_df["GAS_ENTHALPY"].notna()) &
                (updated_df["GAS_ENTHALPY"] != 0.0)
            ]
            if len(existing_nonzero) > 0:
                logging.warning(
                    f"Warning: Found {len(existing_nonzero)} existing non-zero "
                    f"GAS_ENTHALPY values in the file. These will be "
                    f"overwritten."
                )

        # Add enthalpy column
        enthalpy_values = []
        for _, row in updated_df.iterrows():
            species_name = row["NAME"]
            if species_name in mapping:
                enthalpy_values.append(mapping[species_name]["enthalpy_298k"])
            else:
                enthalpy_values.append(None)

        updated_df["GAS_ENTHALPY"] = enthalpy_values

        # Convert kJ/mol to kcal/mol
        updated_df["GAS_ENTHALPY"] /= 4.184
          
        # Write back to original file
        try:
            updated_df.to_csv(original_csv_path, index=False)
            mapped_count = len([e for e in enthalpy_values if e is not None])
            print(f"✓ Updated {original_csv_path} with enthalpy data")
            print(f"  - {mapped_count}/{len(updated_df)} species have "
                  f"enthalpy values")
        except Exception as e:
            print(f"❌ Error writing to file: {e}")
            # Restore backup if write failed
            try:
                shutil.copy2(backup_path, original_csv_path)
                print("✓ Restored original file from backup")
            except Exception as restore_error:
                print(f"❌ Could not restore backup: {restore_error}")


def main():
    """Command-line interface for species matcher."""
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Interactive tool for matching UCLCHEM species with "
        "ATCT thermochemical data."
    )

    parser.add_argument(
        "--atct_csv",
        required=True,
        help="Path to the cleaned ATCT CSV file containing thermochemical data",
    )

    parser.add_argument(
        "--uclchem_species_csv",
        required=True,
        help="Path to CSV file containing UCLCHEM species (must have 'NAME' column)",
    )

    parser.add_argument(
        "--output_mapping",
        help="Path for full output CSV mapping file.",
        required=False,
        default=None,
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
        else:
            matcher.write_back_to_file(species_df, mapping, args.uclchem_species_csv)

        print(f"\n✅ Matching complete! Mapped {len(mapping)} species.")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
