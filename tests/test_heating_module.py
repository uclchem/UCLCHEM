"""
Tests for the heating module.

Tests custom reaction exothermicity loading from CSV files with
various unit formats.
"""

import tempfile
from pathlib import Path

import pytest

import uclchem.makerates.heating as heating

# Import directly to avoid package-level import issues
from uclchem.makerates.reaction import Reaction


class TestCustomExothermicities:
    """Test custom exothermicity loading from CSV."""

    @pytest.fixture
    def sample_custom_csv(self):
        """Create a temporary custom exothermicity CSV file."""
        csv_content = (
            "reactant1,reactant2,reactant3,product1,product2,"
            "product3,product4,exothermicity,unit\n"
            "H,H,NAN,H2,NAN,NAN,NAN,-4.52,ev\n"
            "H,H,NAN,H2,NAN,NAN,NAN,-4.52,ev_per_reaction\n"
            "CO,O,NAN,CO2,NAN,NAN,NAN,-5.5,ev/reaction\n"
            "OH,H,NAN,H2O,NAN,NAN,NAN,-100,kcal/mol\n"
            "C,H,H,CH4,H,H,H,-10.0,ev_per_mol\n"
            "O,O,NAN,O2,NAN,NAN,NAN,-50000,joule/mol\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_path = f.name
        yield temp_path
        Path(temp_path).unlink()

    @pytest.fixture
    def sample_reactions_custom(self):
        """Create sample Reaction objects for custom exothermicity testing."""
        # H + H -> H2
        r1 = Reaction(
            [
                "H",
                "H",
                "NAN",
                "H2",
                "NAN",
                "NAN",
                "NAN",
                1e-10,
                0.0,
                0.0,
                10.0,
                300.0,
                1.0,
                False,
                0.0,
            ]
        )

        # CO + O -> CO2
        r2 = Reaction(
            [
                "CO",
                "O",
                "NAN",
                "CO2",
                "NAN",
                "NAN",
                "NAN",
                1e-10,
                0.0,
                0.0,
                10.0,
                300.0,
                1.0,
                False,
                0.0,
            ]
        )

        # OH + H -> H2O
        r3 = Reaction(
            [
                "OH",
                "H",
                "NAN",
                "H2O",
                "NAN",
                "NAN",
                "NAN",
                1e-10,
                0.0,
                0.0,
                10.0,
                300.0,
                1.0,
                False,
                0.0,
            ]
        )

        # O + O -> O2
        r4 = Reaction(
            [
                "O",
                "O",
                "NAN",
                "O2",
                "NAN",
                "NAN",
                "NAN",
                1e-10,
                0.0,
                0.0,
                10.0,
                300.0,
                1.0,
                False,
                0.0,
            ]
        )

        return [r1, r2, r3, r4]

    def test_load_custom_exothermicities(self, sample_custom_csv):
        """Test loading custom exothermicities from CSV."""
        df = heating.load_custom_exothermicities(sample_custom_csv)

        assert len(df) == 6
        assert "exothermicity" in df.columns
        assert "unit" in df.columns
        assert "reactant1" in df.columns
        assert "product1" in df.columns

    def test_convert_to_erg_dynamic_parsing(self):
        """Test dynamic unit parsing with various formats."""
        # Test base units
        erg1 = heating.convert_to_erg(-4.52, "ev")
        expected1 = -4.52 * heating.EV_TO_ERG
        assert erg1 == pytest.approx(expected1)

        # Test with _per_ separator
        erg2 = heating.convert_to_erg(-4.52, "ev_per_reaction")
        assert erg2 == pytest.approx(expected1)

        # Test with / separator
        erg3 = heating.convert_to_erg(-4.52, "ev/reaction")
        assert erg3 == pytest.approx(expected1)

        # Test per mol
        erg4 = heating.convert_to_erg(-100.0, "kcal/mol")
        expected4 = -100.0 * heating.KCAL_TO_ERG / heating.AVOGADRO_NUMBER
        assert erg4 == pytest.approx(expected4)

        # Test joule
        erg5 = heating.convert_to_erg(-50000, "joule/mol")
        expected5 = -50000 * heating.JOULE_TO_ERG / heating.AVOGADRO_NUMBER
        assert erg5 == pytest.approx(expected5)

        # Test case insensitivity
        erg6 = heating.convert_to_erg(1.0, "EV/MOL")
        erg7 = heating.convert_to_erg(1.0, "ev/mol")
        assert erg6 == pytest.approx(erg7)

    def test_unit_caching(self):
        """Test that unit conversions are cached."""
        # First call should populate cache
        heating.convert_to_erg(1.0, "ev/mol")
        assert "ev/mol" in heating._UNIT_CACHE

        # Second call should use cached value
        cached_factor = heating._UNIT_CACHE["ev/mol"]
        erg = heating.convert_to_erg(2.0, "ev/mol")
        expected = 2.0 * cached_factor
        assert erg == pytest.approx(expected)

    def test_invalid_unit_parsing(self):
        """Test error handling for invalid units."""
        with pytest.raises(ValueError, match="Unknown unit"):
            heating.convert_to_erg(1.0, "invalid_unit")

        with pytest.raises(ValueError, match="Unknown base unit"):
            heating.convert_to_erg(1.0, "fake/mol")

        with pytest.raises(ValueError, match="Unknown denominator"):
            heating.convert_to_erg(1.0, "ev/fake")

    def test_parse_species_from_row(self, sample_custom_csv):
        """Test parsing species from CSV rows."""
        df = heating.load_custom_exothermicities(sample_custom_csv)

        # Test first row (H + H -> H2)
        row = df.iloc[0]
        reactants = heating.parse_species_from_row(row, "reactant")
        assert reactants == ["H", "H", "NAN"]

        products = heating.parse_species_from_row(row, "product")
        assert products == ["H2", "NAN", "NAN", "NAN"]

        # Test row with multiple products
        row = df.iloc[4]  # C + H + H -> CH4 + H + H + H
        reactants = heating.parse_species_from_row(row, "reactant")
        assert "C" in reactants
        assert "H" in reactants

        products = heating.parse_species_from_row(row, "product")
        assert "CH4" in products
        assert "H" in products

    def test_match_reaction(self, sample_reactions_custom):
        """Test matching reactions by reactants and products."""
        reactions = sample_reactions_custom

        # Should find H + H -> H2
        found = heating.match_reaction(
            ["H", "H", "NAN"], ["H2", "NAN", "NAN", "NAN"], reactions
        )
        assert found is not None
        assert sorted(found.get_reactants()) == sorted(["H", "H", "NAN"])

        # Should not find non-existent reaction
        found = heating.match_reaction(
            ["X", "Y", "NAN"], ["Z", "NAN", "NAN", "NAN"], reactions
        )
        assert found is None

    def test_set_custom_exothermicities(
        self, sample_custom_csv, sample_reactions_custom
    ):
        """Test setting custom exothermicities on reactions."""
        reactions = sample_reactions_custom

        matched, unmatched = heating.set_custom_exothermicities(
            reactions, sample_custom_csv, overwrite=True
        )

        # Should match: 2x H+H->H2, CO+O->CO2, OH+H->H2O, O+O->O2 = 5 matches
        assert matched == 5
        # C + H + H -> CH4 + H + H + H should not match (wrong stoichiometry)
        assert unmatched == 1

        # Check that exothermicities were set
        h2_reaction = reactions[0]  # H + H -> H2
        exo = h2_reaction.get_exothermicity()
        assert exo != 0.0
        # Should be -4.52 eV in erg
        expected = -4.52 * heating.EV_TO_ERG
        assert exo == pytest.approx(expected, rel=1e-6)

    def test_overwrite_flag(self, sample_custom_csv, sample_reactions_custom):
        """Test that overwrite flag works correctly."""
        reactions = sample_reactions_custom

        # Set initial values
        reactions[0].set_exothermicity(-1e10)  # Pre-existing non-zero value

        # First call with overwrite=False (should skip)
        matched, _ = heating.set_custom_exothermicities(
            reactions, sample_custom_csv, overwrite=False
        )

        # H+H->H2 should not be updated
        assert reactions[0].get_exothermicity() == pytest.approx(-1e10)

        # Second call with overwrite=True (should update)
        matched, _ = heating.set_custom_exothermicities(
            reactions, sample_custom_csv, overwrite=True
        )

        # H+H->H2 should now be updated
        expected = -4.52 * heating.EV_TO_ERG
        assert reactions[0].get_exothermicity() == pytest.approx(expected, rel=1e-6)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
