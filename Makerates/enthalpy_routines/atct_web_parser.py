#!/usr/bin/env python3
"""
ATCT Web Database Parser

Converts ATCT HTML database files to clean CSV format for thermochemical
analysis. Extracts species names, formulas, phases, and enthalpy of formation
values.

Usage:
    from atct_web_parser import ATCTParser

    parser = ATCTParser()
    data = parser.parse_html_file("ATCTDatabase_v1.220.html")
    parser.save_to_csv(data, "atct_cleaned_v1.220.csv")
"""

import re
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

try:
    from bs4 import BeautifulSoup
except ImportError:
    raise ImportError(
        "BeautifulSoup4 is required. Install with `pip install beautifulsoup4`."
    )


class ATCTParser:
    """Efficient ATCT HTML database parser with robust error handling."""

    def __init__(self):
        self.columns = [
            "Species_Name",
            "Species_formula",
            "Phase",
            "Image",
            "DfH_0K",
            "DfH_298K",
            "Uncertainty",
            "Units",
            "Relative_Molecular_Mass",
            "ATcT_ID",
        ]

    def parse_html_file(self, html_file_path: str) -> Optional[pd.DataFrame]:
        """Parse ATCT HTML file and return cleaned DataFrame.

        Args:
            html_file_path: Path to ATCT HTML database file

        Returns:
            DataFrame with parsed thermochemical data or None if parsing fails
        """
        html_path = Path(html_file_path)
        if not html_path.exists():
            raise FileNotFoundError(f"ATCT HTML file not found: {html_file_path}")

        try:
            with open(html_path, "r", encoding="utf-8") as f:
                html_content = f.read()

            soup = BeautifulSoup(html_content, "html.parser")

            # Find main data table
            data_table = self._find_data_table(soup)
            if data_table is None:
                raise ValueError("Could not locate main data table in HTML file")

            # Extract and clean data
            raw_data = self._extract_table_data(data_table)
            if not raw_data:
                raise ValueError("No data rows found in table")

            # Create DataFrame and clean
            df = pd.DataFrame(raw_data, columns=self.columns[: len(raw_data[0])])
            return self._clean_dataframe(df)

        except Exception as e:
            raise RuntimeError(f"Failed to parse ATCT HTML file: {e}") from e

    def _find_data_table(self, soup: BeautifulSoup) -> Optional[Any]:
        """Find the main thermochemical data table in HTML."""
        tables = soup.find_all("table")

        for table in tables:
            rows = table.find_all("tr")
            if len(rows) > 100:  # Large data table
                sample_text = table.get_text()[:500].lower()
                if any(
                    keyword in sample_text
                    for keyword in ["species", "formula", "enthalpy"]
                ):
                    return table
        return None

    def _extract_table_data(self, table: Any) -> list:
        """Extract raw data from HTML table."""
        rows = table.find_all("tr")
        data = []

        for row in rows[1:]:  # Skip header
            cells = row.find_all(["td", "th"])
            row_data = []

            for cell in cells:
                text = cell.get_text().strip()
                # Prefer link text if available
                links = cell.find_all("a")
                if links and links[0].get_text().strip():
                    text = links[0].get_text().strip()
                row_data.append(text)

            if len(row_data) >= 3:
                # Split formula and phase: "SiC2+ (g)" -> ["SiC2+", "g"]
                if len(row_data) > 1:
                    formula_phase = row_data[1]
                    match = re.match(r"^(.+?)\s+\((.+)\)$", formula_phase)
                    if match:
                        species_formula, phase = match.groups()
                        row_data = [row_data[0], species_formula, phase] + row_data[2:]

                data.append(row_data)

        return data

    def _clean_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean and standardize the parsed DataFrame."""
        # Clean Unicode issues
        for col in df.columns:
            df[col] = df[col].apply(self._clean_unicode_string)

        # Ensure correct column names
        df.columns = self.columns[: len(df.columns)]

        # Convert numeric columns
        numeric_cols = ["DfH_0K", "DfH_298K", "Uncertainty"]
        for col in numeric_cols:
            if col in df.columns:
                df[col] = df[col].str.extract(r"(-?\d+\.?\d*)", expand=False)
                df[col] = pd.to_numeric(df[col], errors="coerce")

        # Remove invalid rows
        df = df[df["Species_Name"].str.len() > 0]

        return df.reset_index(drop=True)

    def _clean_unicode_string(self, text: Any) -> str:
        """Clean problematic Unicode characters from text."""
        if pd.isna(text) or text is None:
            return ""

        text = str(text)
        cleaned = ""
        for char in text:
            try:
                char.encode("utf-8")
                cleaned += char
            except UnicodeEncodeError:
                continue

        return cleaned

    def save_to_csv(self, data: pd.DataFrame, output_path: str) -> None:
        """Save parsed data to CSV file.

        Args:
            data: Parsed ATCT data
            output_path: Output CSV file path
        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        data.to_csv(output_file, index=False)
        print(f"✓ Saved {len(data)} species to {output_path}")

    def get_summary_stats(self, data: pd.DataFrame) -> Dict[str, Any]:
        """Get summary statistics for parsed data.

        Args:
            data: Parsed ATCT data

        Returns:
            Dictionary with summary statistics
        """
        stats = {
            "total_species": len(data),
            "species_with_298k": data["DfH_298K"].notna().sum(),
            "species_with_0k": data["DfH_0K"].notna().sum(),
            "gas_phase": len(data[data["Phase"].str.startswith("g", na=False)]),
            "unique_phases": data["Phase"].nunique(),
        }
        return stats

    def validate_data(self, data: pd.DataFrame, min_species: int = 3000) -> None:
        """Validate parsed data meets expected criteria.

        Args:
            data: Parsed ATCT data
            min_species: Minimum expected number of species

        Raises:
            ValueError: If data doesn't meet validation criteria
        """
        stats = self.get_summary_stats(data)

        if stats["total_species"] < min_species:
            raise ValueError(
                f"Expected ≥{min_species} species, got {stats['total_species']}"
            )

        if stats["species_with_298k"] < min_species:
            raise ValueError(
                f"Expected ≥{min_species} species with 298K data, got {stats['species_with_298k']}"
            )

        print(
            f"✓ Validation passed: {stats['total_species']} species with {stats['species_with_298k']} having 298K data"
        )


def main():
    """Command-line interface for ATCT parser."""
    import sys

    if len(sys.argv) != 3:
        print("Usage: python atct_web_parser.py <input_html> <output_csv>")
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]

    try:
        parser = ATCTParser()
        print(f"Parsing {input_file}...")

        data = parser.parse_html_file(input_file)
        parser.validate_data(data)

        stats = parser.get_summary_stats(data)
        print(f"✓ Successfully parsed {stats['total_species']} species:")
        print(f"  - {stats['species_with_298k']} with ΔfH°(298K) data")
        print(f"  - {stats['species_with_0k']} with ΔfH°(0K) data")
        print(f"  - {stats['gas_phase']} gas phase species")

        parser.save_to_csv(data, output_file)
        print(f"✓ Output saved to {output_file}")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
