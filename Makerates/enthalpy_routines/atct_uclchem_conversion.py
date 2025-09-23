#!/usr/bin/env python3
"""
ATCT-UCLCHEM Reaction Enthalpy Converter

Calculates reaction enthalpies for UCLCHEM reaction networks using ATCT
thermochemical data. Takes UCLCHEM reactions.csv, species.csv, and ATCT
species mapping to produce reactions_with_enthalpy.csv.

Features:
- Reads CSV species mapping (CSV format)
- Calculates reaction enthalpies with 2 decimal precision
- Clean output without species reports

Usage:
    from atct_uclchem_conversion import ReactionEnthalpyCalculator
    
    calc = ReactionEnthalpyCalculator("species_mapping.csv")
    enhanced_reactions = calc.process_reactions("reactions.csv")
    calc.save_enhanced_reactions(enhanced_reactions, "reactions_with_enthalpy.csv")
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


class ReactionEnthalpyCalculator:
    """Calculate reaction enthalpies using ATCT thermochemical data."""
    
    def __init__(self, species_mapping_path: str):
        """Initialize calculator with species mapping.
        
        Args:
            species_mapping_path: Path to CSV file with UCLCHEM-ATCT mapping
        """
        # Species to ignore (not chemical species) - must be defined first
        self.ignore_list = [
            "PHOTON", "CRP", "CRPHOT", "FREEZE", "DESORB", "THERM",
            "DESOH2", "DESCR", "DEUVCR", "H2FORM", "ER", "ERDES",
            "LH", "LHDES", "BULKSWAP", "SURFSWAP", "IONOPOL1", "IONOPOL2",
            "CRS", "EXSOLID", "EXRELAX", "GAR", "E-"
        ]
        
        self.species_mapping = self._load_mapping(species_mapping_path)
        self.species_lookup = self._create_species_lookup()
        
        print(f"✓ Loaded enthalpy data for {len(self.species_lookup)} species")
    
    def _load_mapping(self, mapping_path: str) -> Dict[str, Dict[str, Any]]:
        """Load species mapping from CSV file."""
        mapping_df = pd.read_csv(mapping_path)
        mapping = {}
        
        for _, row in mapping_df.iterrows():
            species = row['uclchem_species']
            mapping[species] = {
                'atct_name': row['atct_name'],
                'atct_formula': row['atct_formula'],
                'enthalpy_298k': row['enthalpy_298k'],
                'match_type': row['match_type'],
                'selection_reason': row['selection_reason']
            }
        
        return mapping
    
    def _create_species_lookup(self) -> Dict[str, float]:
        """Create lookup table for species enthalpies."""
        lookup = {}
        
        for species, match_info in self.species_mapping.items():
            enthalpy = match_info.get('enthalpy_298k')
            if pd.notna(enthalpy):
                lookup[species] = float(enthalpy)
        
        # Add ignored species with zero enthalpy
        for ignored in self.ignore_list:
            lookup[ignored] = 0.0
        
        return lookup
    
    def process_reactions(self, reactions_csv_path: str) -> pd.DataFrame:
        """Process reactions CSV and calculate enthalpies.
        
        Args:
            reactions_csv_path: Path to UCLCHEM reactions.csv file
            
        Returns:
            DataFrame with reaction enthalpies added
        """
        reactions_df = pd.read_csv(reactions_csv_path)
        print(f"Processing {len(reactions_df)} reactions...")
        
        # Calculate enthalpies
        results = []
        for idx, row in reactions_df.iterrows():
            delta_h, missing_species = self._calculate_reaction_enthalpy(row)
            results.append({
                'reaction_id': idx,
                'delta_h_kj_mol': delta_h,
                'status': 'SUCCESS' if delta_h is not None else 'MISSING_SPECIES',
                'missing_species': missing_species
            })
        
        # Add enthalpy column to original dataframe
        results_df = pd.DataFrame(results)
        reactions_df['delta_h_kj_mol'] = 0.0  # Initialize with zeros
        
        # Fill in calculated values
        for _, row in results_df.iterrows():
            reaction_id = row['reaction_id']
            if pd.notna(row['delta_h_kj_mol']):
                reactions_df.loc[reaction_id, 'delta_h_kj_mol'] = row['delta_h_kj_mol']
        
        # Report statistics
        success_count = len(results_df[results_df['status'] == 'SUCCESS'])
        success_rate = success_count / len(reactions_df) * 100
        
        print(f"✓ Processed {len(reactions_df)} reactions:")
        print(f"  - {success_count} successful ({success_rate:.1f}%)")
        print(f"  - {len(reactions_df) - success_count} with missing species")
        
        # Report most common missing species
        all_missing = []
        for _, row in results_df.iterrows():
            if row['missing_species']:
                all_missing.extend(row['missing_species'])
        
        if all_missing:
            missing_counts = pd.Series(all_missing).value_counts()
            print(f"\nMost common missing species:")
            for species, count in missing_counts.head(10).items():
                print(f"  - {species}: {count} reactions")
        
        return reactions_df
    
    def _calculate_reaction_enthalpy(self, reaction_row: pd.Series) -> Tuple[Optional[float], List[str]]:
        """Calculate enthalpy change for a single reaction.
        
        Args:
            reaction_row: Row from reactions DataFrame
            
        Returns:
            Tuple of (enthalpy_change, missing_species_list)
        """
        # Extract reactants and products
        reactants = self._extract_species(reaction_row, ['Reactant 1', 'Reactant 2', 'Reactant 3'])
        products = self._extract_species(reaction_row, ['Product 1', 'Product 2', 'Product 3', 'Product 4'])
        
        # Check for missing species
        all_species = list(reactants.keys()) + list(products.keys())
        missing = [s for s in all_species if s not in self.species_lookup and s not in self.ignore_list]
        
        if missing:
            return None, missing
        
        # Calculate ΔH = Σ(products) - Σ(reactants)
        try:
            products_enthalpy = sum(
                self.species_lookup.get(species, 0.0) * coeff
                for species, coeff in products.items()
            )
            
            reactants_enthalpy = sum(
                self.species_lookup.get(species, 0.0) * coeff
                for species, coeff in reactants.items()
            )
            
            delta_h = products_enthalpy - reactants_enthalpy
            return delta_h, []
            
        except Exception:
            return None, missing
    
    def _extract_species(self, reaction_row: pd.Series, columns: List[str]) -> Dict[str, float]:
        """Extract species and stoichiometry from reaction row.
        
        Args:
            reaction_row: Row from reactions DataFrame
            columns: Column names to extract species from
            
        Returns:
            Dictionary of {species: stoichiometry}
        """
        species_dict = {}
        
        for col in columns:
            if col in reaction_row:
                species = reaction_row[col]
                if isinstance(species, str) and species.upper() not in ["NAN", ""]:
                    species = species.upper()
                    # Default stoichiometry is 1.0
                    species_dict[species] = 1.0
        
        return species_dict
    
    def save_enhanced_reactions(self, 
                               reactions_df: pd.DataFrame, 
                               output_path: str) -> None:
        """Save enhanced reactions to CSV file.
        
        Args:
            reactions_df: DataFrame with enthalpy data
            output_path: Output CSV file path
        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Format delta_h_kj_mol to 2 decimal places
        reactions_df_formatted = reactions_df.copy()
        reactions_df_formatted['delta_h_kj_mol'] = reactions_df_formatted['delta_h_kj_mol'].round(2)
        
        reactions_df_formatted.to_csv(output_file, index=False)
        
        # Calculate statistics
        total_reactions = len(reactions_df)
        with_enthalpy = len(reactions_df[reactions_df['delta_h_kj_mol'] != 0])
        
        print(f"✓ Saved {total_reactions} reactions to {output_path}")
        print(f"  - {with_enthalpy} reactions with calculated enthalpies")
        print(f"  - {total_reactions - with_enthalpy} reactions with zero enthalpy (missing data)")
    
    def generate_species_report(self, 
                               reactions_df: pd.DataFrame,
                               output_path: Optional[str] = None) -> pd.DataFrame:
        """Generate report on species coverage and missing data.
        
        Args:
            reactions_df: DataFrame with reaction data
            output_path: Optional path to save report CSV
            
        Returns:
            DataFrame with species statistics
        """
        # Collect all species from reactions
        all_species = set()
        
        species_columns = ['Reactant 1', 'Reactant 2', 'Reactant 3', 
                          'Product 1', 'Product 2', 'Product 3', 'Product 4']
        
        for col in species_columns:
            if col in reactions_df.columns:
                species_list = reactions_df[col].dropna()
                species_list = species_list[species_list.str.upper() != "NAN"]
                all_species.update(species_list.str.upper())
        
        # Create species report
        report_data = []
        for species in sorted(all_species):
            has_enthalpy = species in self.species_lookup
            is_ignored = species in self.ignore_list
            
            # Count occurrences in reactions
            occurrence_count = 0
            for col in species_columns:
                if col in reactions_df.columns:
                    occurrence_count += (reactions_df[col].str.upper() == species).sum()
            
            status = "IGNORED" if is_ignored else ("MAPPED" if has_enthalpy else "MISSING")
            enthalpy_value = self.species_lookup.get(species, np.nan)
            
            report_data.append({
                'species': species,
                'status': status,
                'enthalpy_kj_mol': enthalpy_value,
                'occurrence_count': occurrence_count,
                'has_atct_data': has_enthalpy and not is_ignored
            })
        
        report_df = pd.DataFrame(report_data)
        report_df = report_df.sort_values(['status', 'occurrence_count'], 
                                         ascending=[True, False])
        
        # Print summary
        total_species = len(report_df)
        mapped_species = len(report_df[report_df['status'] == 'MAPPED'])
        missing_species = len(report_df[report_df['status'] == 'MISSING'])
        ignored_species = len(report_df[report_df['status'] == 'IGNORED'])
        
        print(f"\n📊 Species Coverage Report:")
        print(f"  - Total unique species: {total_species}")
        print(f"  - Mapped with ATCT data: {mapped_species} ({mapped_species/total_species*100:.1f}%)")
        print(f"  - Missing enthalpy data: {missing_species} ({missing_species/total_species*100:.1f}%)")
        print(f"  - Ignored (non-chemical): {ignored_species}")
        
        if missing_species > 0:
            print(f"\nTop missing species by occurrence:")
            missing_df = report_df[report_df['status'] == 'MISSING']
            for _, row in missing_df.head(10).iterrows():
                print(f"  - {row['species']}: {row['occurrence_count']} reactions")
        
        if output_path:
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            report_df.to_csv(output_file, index=False)
            print(f"✓ Species report saved to {output_path}")
        
        return report_df
    
    def update_species_csv(self, 
                          species_csv_path: str,
                          output_path: str) -> None:
        """Update UCLCHEM species.csv with ATCT enthalpy data.
        
        Args:
            species_csv_path: Path to original UCLCHEM species.csv
            output_path: Path for enhanced species.csv
        """
        species_df = pd.read_csv(species_csv_path)
        print(f"Updating {len(species_df)} species with ATCT enthalpy data...")
        
        updated_count = 0
        for idx, row in species_df.iterrows():
            species_name = row['NAME'].upper()
            if species_name in self.species_lookup:
                atct_enthalpy = self.species_lookup[species_name]
                if atct_enthalpy != 0.0:  # Only update non-zero values
                    species_df.loc[idx, 'ENTHALPY'] = atct_enthalpy
                    updated_count += 1
        
        # Save updated species file
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        species_df.to_csv(output_file, index=False)
        
        print(f"✓ Updated {updated_count}/{len(species_df)} species with ATCT data")
        print(f"✓ Enhanced species file saved to {output_path}")


def main():
    """Command-line interface for reaction enthalpy calculator."""
    import sys
    
    if len(sys.argv) < 4:
        print("Usage: python atct_uclchem_conversion.py "
              "<species_mapping.csv> <reactions.csv> <output_reactions.csv>")
        sys.exit(1)
    
    mapping_file, reactions_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
    
    try:
        # Initialize calculator
        calc = ReactionEnthalpyCalculator(mapping_file)
        
        # Process reactions
        enhanced_reactions = calc.process_reactions(reactions_file)
        
        # Save results
        calc.save_enhanced_reactions(enhanced_reactions, output_file)
        
        print(f"\n✅ Conversion complete!")
        print(f"Enhanced reactions saved to: {output_file}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()