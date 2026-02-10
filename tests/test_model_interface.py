import pytest

from uclchem.model import Cloud


def test_invalid_out_species_name_raises():
    # Use run_type='external' to prevent the model from running during init
    with pytest.raises(ValueError):
        Cloud(out_species=["NOT_A_SPECIES"], run_type="external")


def test_non_string_out_species_entry_raises():
    with pytest.raises(ValueError):
        Cloud(out_species=[123], run_type="external")
