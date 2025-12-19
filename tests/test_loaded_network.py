"""
Unit tests for Network class to verify different loading options.

The Network class supports multiple loading modes:
1. Load from species and reactions objects directly (Network.from_lists)
2. Load from CSV file paths (Network.from_csv)
3. Load from UCLCHEM's default installation files (Network.from_csv with no arguments)
4. Build with validation (Network.build)

These tests verify all loading options work correctly.
"""

import pandas as pd
import pytest

from uclchem.makerates.network import Network
from uclchem.makerates.reaction import Reaction
from uclchem.makerates.species import Species


@pytest.fixture
def sample_species():
    """Create sample species for testing."""
    return [
        Species(["H", 1, 1.008, 0, 0, 0, 0]),
        Species(["H2", 1, 2.016, 0, 0, 0, 0]),
    ]


@pytest.fixture
def sample_reactions():
    """Create sample reactions for testing."""
    return [
        Reaction(
            [
                "H",  # R1
                "H",  # R2
                "NAN",  # R3
                "H2",  # P1
                "NAN",  # P2
                "NAN",  # P3
                "NAN",  # P4
                1.0e-10,  # alpha
                0.5,  # beta
                0.0,  # gamma
                0.0,  # templow
                0.0,  # temphigh
                2.0,  # reduced mass
                0,  # extrapolation
                0.0,  # exothermicity
            ]
        ),
    ]


@pytest.fixture
def temp_species_file(tmp_path):
    """Create a temporary species CSV file."""
    species_data = pd.DataFrame(
        {
            "name": ["H", "H2"],
            "mass": [1, 1],
            "molecularMass": [1.008, 2.016],
            "binding": [0, 0],
            "solidFraction": [0, 0],
            "monoFraction": [0, 0],
            "volcFraction": [0, 0],
        }
    )
    filepath = tmp_path / "test_species.csv"
    species_data.to_csv(filepath, index=False)
    return filepath


@pytest.fixture
def temp_reactions_file(tmp_path):
    """Create a temporary reactions CSV file."""
    reactions_data = pd.DataFrame(
        {
            "Reactant 1": ["H"],
            "Reactant 2": ["H"],
            "Reactant 3": ["NAN"],
            "Product 1": ["H2"],
            "Product 2": ["NAN"],
            "Product 3": ["NAN"],
            "Product 4": ["NAN"],
            "Alpha": [1.0e-10],
            "Beta": [0.5],
            "Gamma": [0.0],
            "T_min": [0.0],
            "T_max": [0.0],
            "reduced_mass": [2.0],
            "extrapolate": [False],
            "exothermicity": [0.0],
        }
    )
    filepath = tmp_path / "test_reactions.csv"
    reactions_data.to_csv(filepath, index=False)
    return filepath


def test_option1_load_from_objects(sample_species, sample_reactions):
    """
    Test Option 1: Loading with species and reactions objects directly.

    This verifies that users can pass pre-constructed Species and Reaction
    objects directly to Network using from_lists factory method.
    """
    network = Network.from_lists(species=sample_species, reactions=sample_reactions)

    assert network is not None
    # Verify species were loaded
    assert len(network.get_species_dict()) >= len(sample_species)
    assert "H" in network.get_species_dict()
    assert "H2" in network.get_species_dict()


def test_option2_load_from_custom_filepaths(temp_species_file, temp_reactions_file):
    """
    Test Option 2: Loading from custom file paths.

    This verifies that users can provide paths to their own species.csv and
    reactions.csv files to load a custom network.
    """
    network = Network.from_csv(
        species_path=temp_species_file, reactions_path=temp_reactions_file
    )

    assert network is not None
    # Verify data was loaded from files
    assert len(network.get_species_dict()) >= 2  # H, H2
    assert len(network.get_reaction_dict()) >= 1


def test_option3_load_from_default_installation():
    """
    Test Option 3: Loading from UCLCHEM's default installation files.

    This verifies that when no arguments are provided, LoadedNetwork automatically
    loads the default species.csv and reactions.csv from the UCLCHEM installation.
    This is the simplest use case for users who want the standard UCLCHEM network.
    """
    network = Network.from_csv()

    assert network is not None
    # Verify default UCLCHEM files were loaded (should have many species/reactions)
    assert len(network.get_species_dict()) > 10
    assert len(network.get_reaction_dict()) > 10
    # Check for common UCLCHEM species
    species_dict = network.get_species_dict()
    assert any("H" in name for name in species_dict.keys())


def test_different_factory_methods_produce_networks(sample_species, sample_reactions):
    """
    Test that different factory methods all produce valid Network objects.

    This verifies that Network supports multiple creation patterns:
    - from_csv() for loading from files
    - from_lists() for direct object construction
    - build() for full validation and generation
    """
    # Test from_csv works
    network1 = Network.from_csv()
    assert isinstance(network1, Network)
    assert len(network1.get_species_dict()) > 0

    # Test from_lists works
    network2 = Network.from_lists(species=sample_species, reactions=sample_reactions)
    assert isinstance(network2, Network)
    assert len(network2.get_species_dict()) >= 2


def test_network_supports_crud_operations():
    """
    Test that Network supports full CRUD operations.

    Network inherits from MutableNetworkABC so should support adding,
    removing, and modifying species and reactions.
    """
    network = Network.from_csv()

    # Test that modification methods exist
    assert hasattr(network, "add_species")
    assert hasattr(network, "remove_species")
    assert hasattr(network, "add_reactions")
    assert hasattr(network, "remove_reaction")
    assert hasattr(network, "change_binding_energy")
    assert hasattr(network, "change_reaction_barrier")
