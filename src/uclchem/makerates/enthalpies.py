from pathlib import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import Path

from uclchem.makerates.reaction import Reaction, reaction_types


class ReactionEnthalpyCalculator:
    """Calculate reaction enthalpies using ATCT thermochemical data."""

    def __init__(self, uclchem_enthalpy_mapping: str):
        """Initialize calculator with species mapping.

        Args:
            species_mapping_path: Path to CSV file with UCLCHEM-ATCT mapping
        """
        # Ignore the reaction type strings and electrons
        self.ignore_list = reaction_types + [
            "E-",
        ]
        self.atct_entry_by_species = self._load_mapping(uclchem_enthalpy_mapping)

    def _load_mapping(self, mapping_path: str) -> Dict[str, Dict[str, Any]]:
        """Load species mapping from CSV file."""
        mapping_df = pd.read_csv(mapping_path)
        mapping = {}

        # This can all be generalized if we provide a standardized format for BURCAT as well
        for _, row in mapping_df.iterrows():
            species = row["uclchem_species"]
                
            # First load the nice_to_have columns if they exist:
            mapping[species] = {
                "atct_name": row.get("atct_name", None),
                "atct_formula": row.get("atct_formula", None),
                "match_type": row.get("match_type", None),
                "selection_reason": row.get("selection_reason", None),
            }
            # Ensure we load the enthalpy around 298K:
            try:
                mapping[species]["enthalpy_298k"] = float(row["enthalpy_298k"])
            except KeyError as e:
                e.add_note(f"Missing enthalpy_298k entry in mapping file for {species}")
                raise e
            # Assign zero enthalpy to all ignored species
            if species in self.ignore_list:
                logging.debug(f"Ignoring species {species} in enthalpy calculations")
                mapping[species]["enthalpy_298k"] = 0.0
        return mapping

    def process_reactions(self, reactions: list[Reaction]) -> pd.DataFrame:
        """Process reactions CSV and calculate enthalpies.

        Args:
            reactions_csv_path: Path to UCLCHEM reactions.csv file

        Returns:
            DataFrame with reaction enthalpies added
        """

        # Calculate enthalpies and directly assign them
        for reaction in reactions:
            delta_h, missing_species = self._calculate_reaction_enthalpy(reaction)
            if len(missing_species) > 0:
                delta_h = 0.0
            reaction.set_delta_enthalpy(delta_h)

    def _calculate_reaction_enthalpy(
        self, reaction: Reaction
    ) -> Tuple[Optional[float], List[str]]:
        """Calculate enthalpy change for a single reaction.

        Args:
            reaction_row: Row from reactions DataFrame

        Returns:
            Tuple of (enthalpy_change, missing_species_list)
        """
        # Extract reactants and products
        reactants = reaction.get_reactants()
        products = reaction.get_products()

        # Check for missing species
        all_species = reactants + products
        missing = [
            s
            for s in all_species
            if s not in self.species_lookup and s not in self.ignore_list
        ]

        # Disable enthalpy in case more than one species is missing:
        if missing:
            return 0.0, missing

        # Calculate ΔH = sum(products) - sum(reactants)
        products_enthalpy = sum(
            self.species_lookup.get(species, 0.0) * coeff
            for species, coeff in products.items()
        )

        reactants_enthalpy = sum(
            self.species_lookup.get(species, 0.0) * coeff
            for species, coeff in reactants.items()
        )

        delta_h = products_enthalpy - reactants_enthalpy
        return delta_h, missing
