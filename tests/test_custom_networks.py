"""
Tests to verify the prefix network and small chemistry network do function.
Checks that small_chemistry and small_chemistry_prefix networks
can be loaded and have expected properties.
"""

from pathlib import Path

import pytest

from uclchem.makerates.makerates import run_makerates

NETWORKS_DIR = Path(__file__).parent / "networks"


@pytest.fixture(scope="module")
def small_chemistry_network():
    """Load the small_chemistry test network."""
    settings = NETWORKS_DIR / "small_chemistry" / "user_settings.yaml"
    return run_makerates(str(settings), write_files=False)


@pytest.fixture(scope="module")
def small_chemistry_prefix_network():
    """Load the small_chemistry_prefix test network."""
    settings = NETWORKS_DIR / "small_chemistry_prefix" / "user_settings.yaml"
    return run_makerates(str(settings), write_files=False)


# ============================================================================
# Small Chemistry Network Tests
# ============================================================================


def test_small_chemistry_network_loads(small_chemistry_network):
    """Test that small_chemistry network loads without errors."""
    assert small_chemistry_network is not None


def test_small_chemistry_has_species(small_chemistry_network):
    """Test that small_chemistry network has species."""
    species_list = small_chemistry_network.get_species_list()
    assert len(species_list) > 0
    species_names = {s.get_name() for s in species_list}
    assert "H" in species_names
    assert "CH4" in species_names
    assert "#CH4" in species_names


def test_small_chemistry_has_reactions(small_chemistry_network):
    """Test that small_chemistry network has reactions."""
    reactions = small_chemistry_network.get_reaction_list()
    assert len(reactions) > 0


def test_small_chemistry_no_prefix_species(small_chemistry_network):
    """Test that small_chemistry does not have prefix species."""
    species_names = {s.get_name() for s in small_chemistry_network.get_species_list()}
    assert "a-CH4" not in species_names
    assert "b-CH4" not in species_names
    assert "a-NH3" not in species_names
    assert "b-NH3" not in species_names


# ============================================================================
# Small Chemistry with Prefixes Network Tests
# ============================================================================


def test_small_chemistry_prefix_network_loads(small_chemistry_prefix_network):
    """Test that small_chemistry_prefix network loads without errors."""
    assert small_chemistry_prefix_network is not None


def test_small_chemistry_prefix_has_species(small_chemistry_prefix_network):
    """Test that small_chemistry_prefix has species."""
    species_list = small_chemistry_prefix_network.get_species_list()
    assert len(species_list) > 0
    species_names = {s.get_name() for s in species_list}
    assert "H" in species_names
    assert "CH4" in species_names


def test_small_chemistry_prefix_has_prefix_species(small_chemistry_prefix_network):
    """Test that small_chemistry_prefix has all expected prefix species."""
    species_names = {s.get_name() for s in small_chemistry_prefix_network.get_species_list()}

    # Check CH4 isomers
    assert "a-CH4" in species_names
    assert "b-CH4" in species_names
    assert "#CH4" in species_names

    # Check NH3 isomers
    assert "a-NH3" in species_names
    assert "b-NH3" in species_names
    assert "#a-NH3" in species_names
    assert "#b-NH3" in species_names


def test_small_chemistry_prefix_species_properties(small_chemistry_prefix_network):
    """Test that prefix species have correct properties."""
    species_dict = {s.get_name(): s for s in small_chemistry_prefix_network.get_species_list()}

    # Check a-CH4 and b-CH4 have correct mass
    assert species_dict["a-CH4"].get_mass() == 16
    assert species_dict["b-CH4"].get_mass() == 16

    # Check a-NH3 and b-NH3 have correct mass
    assert species_dict["a-NH3"].get_mass() == 17
    assert species_dict["b-NH3"].get_mass() == 17


def test_small_chemistry_prefix_freeze_reactions(small_chemistry_prefix_network):
    """Test that freeze reactions exist for prefix species."""
    freeze_rxns = [
        r for r in small_chemistry_prefix_network.get_reaction_list()
        if r.get_reaction_type() == "FREEZE"
    ]
    freeze_reactants = [r.get_reactants()[0] for r in freeze_rxns]

    # a-CH4 and b-CH4 should freeze to #CH4
    assert "a-CH4" in freeze_reactants
    assert "b-CH4" in freeze_reactants

    # a-NH3 and b-NH3 should freeze to their own grain versions
    assert "a-NH3" in freeze_reactants
    assert "b-NH3" in freeze_reactants


def test_small_chemistry_prefix_desorb_reactions(small_chemistry_prefix_network):
    """Test that desorption reactions exist for prefix species."""
    desorb_types = {"DESORB", "DESOH2", "DESCR", "DEUVCR", "THERM"}
    desorb_rxns = [
        r for r in small_chemistry_prefix_network.get_reaction_list()
        if r.get_reaction_type() in desorb_types
    ]

    # #CH4 should have desorption reactions
    ch4_desorbs = [
        r for r in desorb_rxns
        if r.get_reactants()[0] == "#CH4"
    ]
    assert len(ch4_desorbs) > 0

    # Check that both a-CH4 and b-CH4 appear as products
    all_products = [p for r in ch4_desorbs for p in r.get_products() if p != "NAN"]
    assert "a-CH4" in all_products
    assert "b-CH4" in all_products


def test_small_chemistry_prefix_has_reactions(small_chemistry_prefix_network):
    """Test that small_chemistry_prefix network has reactions."""
    reactions = small_chemistry_prefix_network.get_reaction_list()
    assert len(reactions) > 0
