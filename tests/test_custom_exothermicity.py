"""
Unit tests for heating module.

Tests conversion constants and exothermicity functions.
"""

import tempfile
from pathlib import Path

import pytest

from uclchem.makerates import heating


class TestConversionConstants:
    """Test conversion constants."""

    def test_avogadro_number(self):
        assert heating.AVOGADRO_NUMBER == 6.02214076e23

    def test_electron_volt_to_joule(self):
        assert heating.EV_TO_JOULE == 1.602176634e-19

    def test_derived_conversions(self):
        expected_ev_to_erg = heating.EV_TO_JOULE / heating.ERG_TO_JOULE
        assert abs(heating.EV_TO_ERG - expected_ev_to_erg) < 1e-10

        expected_kcal_to_erg = heating.CALORIE_TO_JOULE * 1000.0 / heating.ERG_TO_JOULE
        assert abs(heating.KCAL_TO_ERG - expected_kcal_to_erg) < 1e-10


class TestHeatingFunctions:
    """Test heating module functions."""

    @pytest.fixture
    def sample_csv(self):
        csv_content = (
            "reactant1,reactant2,reactant3,product1,product2,"
            "product3,product4,exothermicity,unit\n"
            "H,H,NAN,H2,NAN,NAN,NAN,-4.52,ev_per_reaction\n"
            "CO,O,NAN,CO2,NAN,NAN,NAN,-5.5,ev\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_path = f.name
        yield temp_path
        Path(temp_path).unlink()

    def test_load_custom_exothermicities(self, sample_csv):
        df = heating.load_custom_exothermicities(sample_csv)
        assert len(df) == 2
        assert "exothermicity" in df.columns
        assert "unit" in df.columns

    def test_convert_to_erg(self):
        erg_val = heating.convert_to_erg(-4.52, "ev_per_reaction")
        # ev_per_reaction means eV per reaction, which converts to erg using EV_TO_ERG * 1.0
        expected = -4.52 * heating.EV_TO_ERG
        assert abs(erg_val - expected) < 1e-20

        # Test case insensitivity
        erg1 = heating.convert_to_erg(1.0, "EV")
        erg2 = heating.convert_to_erg(1.0, "ev")
        assert erg1 == erg2

    def test_invalid_unit(self):
        with pytest.raises(ValueError, match="Unknown unit"):
            heating.convert_to_erg(1.0, "invalid_unit")

    def test_parse_species_from_row(self, sample_csv):
        df = heating.load_custom_exothermicities(sample_csv)
        row = df.iloc[0]

        reactants = heating.parse_species_from_row(row, "reactant")
        assert reactants == ["H", "H", "NAN"]

        products = heating.parse_species_from_row(row, "product")
        assert products == ["H2", "NAN", "NAN", "NAN"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
