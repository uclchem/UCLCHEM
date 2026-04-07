"""
Tests to verify the prefix network and small chemistry network do function.
Checks that small_chemistry and small_chemistry_prefix networks
can be loaded and have expected properties.
"""

import subprocess
import sys
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
    species_names = {
        s.get_name() for s in small_chemistry_prefix_network.get_species_list()
    }

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
    species_dict = {
        s.get_name(): s for s in small_chemistry_prefix_network.get_species_list()
    }

    # Check a-CH4 and b-CH4 have correct mass
    assert species_dict["a-CH4"].get_mass() == 16
    assert species_dict["b-CH4"].get_mass() == 16

    # Check a-NH3 and b-NH3 have correct mass
    assert species_dict["a-NH3"].get_mass() == 17
    assert species_dict["b-NH3"].get_mass() == 17


def test_small_chemistry_prefix_freeze_reactions(small_chemistry_prefix_network):
    """Test that freeze reactions exist for prefix species."""
    freeze_rxns = [
        r
        for r in small_chemistry_prefix_network.get_reaction_list()
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
        r
        for r in small_chemistry_prefix_network.get_reaction_list()
        if r.get_reaction_type() in desorb_types
    ]

    # #CH4 should have desorption reactions
    ch4_desorbs = [r for r in desorb_rxns if r.get_reactants()[0] == "#CH4"]
    assert len(ch4_desorbs) > 0

    # Check that both a-CH4 and b-CH4 appear as products
    all_products = [p for r in ch4_desorbs for p in r.get_products() if p != "NAN"]
    assert "a-CH4" in all_products
    assert "b-CH4" in all_products


def test_small_chemistry_prefix_has_reactions(small_chemistry_prefix_network):
    """Test that small_chemistry_prefix network has reactions."""
    reactions = small_chemistry_prefix_network.get_reaction_list()
    assert len(reactions) > 0


# ============================================================================
# Dynamic Network Installation Tests
# ============================================================================
# These tests validate that each custom network can be installed and rebuilt
# using makerates. The workflow is:
# 1. Run makerates (generate rates)
# 2. Load the network
# 3. Run makerates again (verify network installation doesn't break rates)


def get_all_networks():
    """Discover all custom networks in the tests/networks directory."""
    networks_path = NETWORKS_DIR
    if not networks_path.exists():
        return []

    networks = []
    for network_dir in sorted(networks_path.iterdir()):
        if network_dir.is_dir():
            settings_file = network_dir / "user_settings.yaml"
            if settings_file.exists():
                networks.append((network_dir.name, settings_file))

    return networks


@pytest.mark.parametrize("network_name,settings_file", get_all_networks())
def test_network_installable(network_name, settings_file):
    """Test that each custom network can be installed with makerates.

    Workflow:
    1. Run makerates to generate reaction rates
    2. Load the network from the generated rates
    3. Run makerates again to verify the network doesn't break rates generation

    This validates that networks are self-consistent and can be rebuilt.
    """
    # Step 1: Run makerates (generate rates)
    network = run_makerates(str(settings_file), write_files=False)
    assert network is not None, f"Failed to generate rates for {network_name}"

    # Step 2: Verify network loaded correctly
    species = network.get_species_list()
    reactions = network.get_reaction_list()
    assert len(species) > 0, f"Network {network_name} has no species"
    assert len(reactions) > 0, f"Network {network_name} has no reactions"

    # Step 3: Run makerates again (verify network installation consistency)
    # This ensures the network can be rebuilt without errors
    network_2 = run_makerates(str(settings_file), write_files=False)
    assert network_2 is not None, f"Failed to rebuild rates for {network_name}"

    # Verify the rebuilt network has same number of species and reactions
    species_2 = network_2.get_species_list()
    reactions_2 = network_2.get_reaction_list()
    assert len(species_2) == len(species), (
        f"Network {network_name} inconsistent: "
        f"first build {len(species)} species, second build {len(species_2)}"
    )
    assert len(reactions_2) == len(reactions), (
        f"Network {network_name} inconsistent: "
        f"first build {len(reactions)} reactions, second build {len(reactions_2)}"
    )


@pytest.mark.parametrize("network_name,settings_file", get_all_networks())
def test_network_has_minimum_content(network_name, settings_file):
    """Test that each custom network has minimum expected content.

    Validates:
    - Network has species definitions
    - Network has reactions
    - Network can be accessed without errors
    """
    network = run_makerates(str(settings_file), write_files=False)

    species_names = {s.get_name() for s in network.get_species_list()}
    reactions = network.get_reaction_list()

    # Every network must have at least one species and reaction
    assert len(species_names) > 0, f"Network {network_name} has no species"
    assert len(reactions) > 0, f"Network {network_name} has no reactions"

    # Every network should have at least hydrogen (fundamental species)
    assert "H" in species_names, f"Network {network_name} missing fundamental species 'H'"


# ============================================================================
# Integration Tests - Makerates → Install → Cloud Model
# ============================================================================
# These tests run a complete workflow: generate network files with makerates,
# reinstall the package to recompile Fortran, and run a minimal cloud model.
# This validates that networks are not only internally consistent, but also
# work end-to-end with the compiled binary.


@pytest.mark.timeout(600)
@pytest.mark.parametrize("network_name,settings_file", get_all_networks())
def test_network_cloud_model_runs(network_name, settings_file):
    """Full integration test: makerates → pip install → run Cloud model.

    Workflow:
    1. Run makerates writing Fortran source files to src/
    2. Reinstall the package to recompile the extension with the new network
    3. Run a minimal cloud model (1e4 yr) in a fresh subprocess
       (to load the newly compiled binary, not the stale in-process one)

    This validates that networks can be fully installed and used in simulation.
    """
    uclchem_root = Path(__file__).resolve().parent.parent

    # Step 1: Run makerates and write Fortran source files
    run_makerates(str(settings_file), write_files=True)

    # Step 2: Reinstall to recompile the Fortran extension
    result = subprocess.run(
        "pip install .",
        shell=True,
        text=True,
        capture_output=True,
        cwd=str(uclchem_root),
        timeout=300,
    )
    assert result.returncode == 0, (
        f"pip install failed for {network_name}:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # Step 3: Run a minimal cloud model in a fresh subprocess
    # (avoids stale .so binary locked in memory by the test process)
    script = """
import uclchem
cloud = uclchem.model.Cloud(
    param_dict={
        "finalTime": 1e4,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "freefall": False,
        "endAtFinalDensity": False,
    },
    out_species=["H", "H2", "CO"]
)
cloud.check_error()
print("OK")
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, (
        f"Cloud model failed for {network_name}:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )
