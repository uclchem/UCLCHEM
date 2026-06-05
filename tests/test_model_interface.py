import pytest

from uclchem.model import Cloud


def test_out_species_on_oo_model_raises():
    with pytest.raises(TypeError):
        Cloud(out_species=["CO"], run_type="external")
