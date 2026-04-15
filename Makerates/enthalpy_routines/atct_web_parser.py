#!/usr/bin/env python3
"""ATCT Web Database Parser.

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
from typing import Any

import pandas as pd

try:
    from bs4 import BeautifulSoup
except ImportError:
    msg = "BeautifulSoup4 is required. Install with `pip install beautifulsoup4`."
    raise ImportError(msg)


class ATCTParser:
    """Efficient ATCT HTML database parser with robust error handling."""

    def __init__(self):
        """Create an ATCTParser."""
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

    def parse_html_file(self, html_file_path: str | Path) -> pd.DataFrame | None:
        """Parse ATCT HTML file and return cleaned DataFrame.

        Args:
            html_file_path (str | Path): Path to ATCT HTML database file

        Returns:
            pd.DataFrame with parsed thermochemical data or None if parsing fails

        Raises:
            FileNotFoundError: If the ATcT HTML file could not be found.
            ValueError: If the main data table could not be located,
                or no data was found in the table.
            RuntimeError: If the parsing of the ATcT failed.

        """
        html_path = Path(html_file_path)
        if not html_path.exists():
            msg = f"ATCT HTML file not found: {html_file_path}"
            raise FileNotFoundError(msg)

        try:
            with html_path.open(encoding="utf-8") as f:
                html_content = f.read()

            soup = BeautifulSoup(html_content, "html.parser")

            # Find main data table
            data_table = self._find_data_table(soup)
            if data_table is None:
                msg = "Could not locate main data table in HTML file"
                raise ValueError(msg)

            # Extract and clean data
            raw_data = self._extract_table_data(data_table)
            if not raw_data:
                msg = "No data rows found in table"
                raise ValueError(msg)

            # Create DataFrame and clean
            df = pd.DataFrame(raw_data, columns=self.columns[: len(raw_data[0])])
            return self._clean_dataframe(df)

        except Exception as e:
            msg = f"Failed to parse ATCT HTML file: {e}"
            raise RuntimeError(msg) from e

    @staticmethod
    def _find_data_table(soup: BeautifulSoup) -> Any | None:
        """Find the main thermochemical data table in HTML.

        Args:
            soup (BeautifulSoup): soup instance to scrape ATcT.

        Returns:
            Any | None: Table if it could be found, else None.

        """
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

    @staticmethod
    def _extract_table_data(table: Any) -> list:
        """Extract raw data from HTML table.

        Args:
            table (Any): HTML table to extract data from

        Returns:
            list: extracted data.

        """
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
        """Clean and standardize the parsed DataFrame.

        Args:
            df (pd.DataFrame): dataframe to clean

        Returns:
            pd.DataFrame: DataFrame with unnecessary columns removed, and cleaned
                column names.
        """
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

    @staticmethod
    def _clean_unicode_string(text: Any) -> str:
        """Clean problematic Unicode characters from text.

        Args:
            text (Any): object

        Returns:
            cleaned (str): string with problematic Unicode characters removed.

        """
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

    @staticmethod
    def save_to_csv(data: pd.DataFrame, output_path: str | Path) -> None:
        """Save parsed data to CSV file.

        Args:
            data (pd.DataFrame): Parsed ATCT data
            output_path (str | Path): Output CSV file path

        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        data.to_csv(output_file, index=False)
        print(f"✓ Saved {len(data)} species to {output_path}")

    @staticmethod
    def get_summary_stats(data: pd.DataFrame) -> dict[str, int]:
        """Get summary statistics for parsed data.

        Args:
            data (pd.DataFrame): Parsed ATCT data

        Returns:
            dict[str, int]: Dictionary with summary statistics

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
            msg = f"Expected ≥{min_species} species, got {stats['total_species']}"
            raise ValueError(msg)

        if stats["species_with_298k"] < min_species:
            msg = f"Expected ≥{min_species} species with 298K data, got {stats['species_with_298k']}"
            raise ValueError(msg)

        print(
            f"✓ Validation passed: {stats['total_species']} species with {stats['species_with_298k']} having 298K data"
        )


def main() -> None:
    """Command-line interface for ATCT parser.

    Raises:
        RuntimeError: If parsing of input file failed.

    """
    import sys

    if len(sys.argv) != 3:
        print("Usage: python atct_web_parser.py <input_html> <output_csv>")
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]

    try:
        parser = ATCTParser()
        print(f"Parsing {input_file}...")

        data = parser.parse_html_file(input_file)
        if data is None:
            msg = f"Parsing of {input_file} failed"
            raise RuntimeError(msg)

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
