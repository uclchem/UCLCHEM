"""
Test building a network with both computed enthalpies and custom exothermicity files.
"""

import pytest
import yaml

from uclchem.makerates.makerates import get_network


@pytest.fixture
def custom_exothermicity_file(tmp_path):
    """Create a custom exothermicity file with dummy incrementing values.

    Uses reactions that actually exist in UMIST22 database:
    - CH+ + e- -> C + H (dissociative recombination)
    - H3+ + e- -> H2 + H (dissociative recombination)
    - C10H+ + e- -> C10 + H (dissociative recombination)
    """
    csv_content = """reactant1,reactant2,reactant3,product1,product2,product3,product4,exothermicity,unit,reference
CH+,E-,NAN,C,H,NAN,NAN,-1.0,eV,Test
H3+,E-,NAN,H2,H,NAN,NAN,-2.0,eV,Test
C10H+,E-,NAN,C10,H,NAN,NAN,-3.0,eV,Test
C2H+,E-,NAN,C2,H,NAN,NAN,-4.0,eV,Test
C3H+,E-,NAN,C3,H,NAN,NAN,-5.0,eV,Test
C4H+,E-,NAN,C4,H,NAN,NAN,-6.0,eV,Test
C5H+,E-,NAN,C5,H,NAN,NAN,-7.0,eV,Test
C6H+,E-,NAN,C6,H,NAN,NAN,-8.0,eV,Test
C7H+,E-,NAN,C7,H,NAN,NAN,-9.0,eV,Test
C8H+,E-,NAN,C8,H,NAN,NAN,-10.0,eV,Test
C9H+,E-,NAN,C9,H,NAN,NAN,-11.0,eV,Test
NH+,E-,NAN,N,H,NAN,NAN,-12.0,eV,Test
OH+,E-,NAN,O,H,NAN,NAN,-13.0,eV,Test
"""

    exo_file = tmp_path / "custom_exothermicity.csv"
    exo_file.write_text(csv_content)
    return exo_file


@pytest.fixture
def config_file_with_exothermicity(tmp_path, custom_exothermicity_file):
    """Create a complete configuration file with exothermicity settings."""
    import os

    workspace_root = os.getcwd()

    config = {
        "species_file": f"{workspace_root}/Makerates/data/default/default_species.csv",
        "database_reaction_file": f"{workspace_root}/Makerates/data/databases/umist22.csv",
        "database_reaction_type": "UMIST",
        "custom_reaction_file": f"{workspace_root}/Makerates/data/default/default_grain_network.csv",
        "custom_reaction_type": "UCL",
        "add_crp_photo_to_grain": False,
        "enable_rates_storage": False,
        "three_phase": True,
        "gas_phase_extrapolation": False,
        "derive_reaction_exothermicity": "GAS",
        "database_reaction_exothermicity": str(custom_exothermicity_file),
        "output_directory": str(tmp_path / "output"),
    }

    config_path = tmp_path / "test_config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    return config_path


@pytest.fixture
def config_file_multiple_exo_files(tmp_path, custom_exothermicity_file):
    """Create a configuration file with multiple exothermicity files."""
    import os

    workspace_root = os.getcwd()

    # Create a second exothermicity file that overrides some values
    exo_csv_2 = tmp_path / "custom_exothermicities_2.csv"
    content = """reactant1,reactant2,reactant3,product1,product2,product3,product4,exothermicity,unit,reference
CH+,E-,NAN,C,H,NAN,NAN,-100.0,ev,Override Test Data
H3+,E-,NAN,H2,H,NAN,NAN,-200.0,ev,Override Test Data
"""
    exo_csv_2.write_text(content)

    config = {
        "species_file": f"{workspace_root}/Makerates/data/default/default_species.csv",
        "database_reaction_file": f"{workspace_root}/Makerates/data/databases/umist22.csv",
        "database_reaction_type": "UMIST",
        "custom_reaction_file": f"{workspace_root}/Makerates/data/default/default_grain_network.csv",
        "custom_reaction_type": "UCL",
        "add_crp_photo_to_grain": False,
        "enable_rates_storage": False,
        "three_phase": True,
        "gas_phase_extrapolation": False,
        "derive_reaction_exothermicity": "GAS",
        "database_reaction_exothermicity": [
            str(custom_exothermicity_file),
            str(exo_csv_2),
        ],
        "output_directory": str(tmp_path / "output"),
    }

    config_path = tmp_path / "test_config_multi.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    return config_path


def test_network_with_custom_exothermicity(config_file_with_exothermicity):
    """Test building a network with custom exothermicity files."""
    # Debug: Print config file content
    import yaml

    with open(config_file_with_exothermicity) as f:
        config = yaml.safe_load(f)
    print(
        f"\n=== Config database_reaction_exothermicity: {config.get('database_reaction_exothermicity')} ==="
    )

    # Build the network
    network = get_network(path_to_input_file=config_file_with_exothermicity)

    assert network is not None
    assert len(network.get_species_list()) > 0
    assert len(network.get_reaction_list()) > 0

    # Check that custom exothermicities were applied
    reactions_with_exo = [
        r
        for r in network.get_reaction_list()
        if r.get_exothermicity() is not None and r.get_exothermicity() != 0.0
    ]

    # We should have reactions with exothermicity set
    assert len(reactions_with_exo) > 0, "No reactions have exothermicity set"

    # Find specific reactions and check their exothermicity values
    # CH+ + E- -> C + H should have -1.0 eV (dummy test data)
    ch_recomb = None
    for reaction in network.get_reaction_list():
        reactants = [r for r in reaction.get_reactants() if r != "NAN"]
        products = [p for p in reaction.get_products() if p != "NAN"]
        if set(reactants) == {"CH+", "E-"} and set(products) == {"C", "H"}:
            ch_recomb = reaction
            break

    assert ch_recomb is not None, "CH+ + E- -> C + H reaction not found"
    assert (
        ch_recomb.get_exothermicity() is not None
    ), "CH+ recombination has no exothermicity"

    # Convert from erg to eV for comparison (1 eV = 1.602176634e-12 erg)
    EV_TO_ERG = 1.602176634e-12
    exo_ev = ch_recomb.get_exothermicity() / EV_TO_ERG

    # Should be approximately -1.0 eV (dummy test data)
    assert abs(exo_ev - (-1.0)) < 0.01, f"Expected -1.0 eV, got {exo_ev} eV"

    print(f"\n✓ CH+ + E- -> C + H has exothermicity: {exo_ev:.2f} eV")


def test_network_with_multiple_database_reaction_exothermicity(
    config_file_multiple_exo_files,
):
    """Test that multiple exothermicity files are processed in order with later files overriding."""
    # Build the network
    network = get_network(path_to_input_file=config_file_multiple_exo_files)

    assert network is not None

    # Find CH+ + E- -> C + H reaction
    ch_recomb = None
    for reaction in network.get_reaction_list():
        reactants = [r for r in reaction.get_reactants() if r != "NAN"]
        products = [p for p in reaction.get_products() if p != "NAN"]
        if set(reactants) == {"CH+", "E-"} and set(products) == {"C", "H"}:
            ch_recomb = reaction
            break

    assert ch_recomb is not None, "CH+ + E- -> C + H reaction not found"
    assert ch_recomb.get_exothermicity() is not None

    # Convert from erg to eV
    EV_TO_ERG = 1.602176634e-12
    exo_ev = ch_recomb.get_exothermicity() / EV_TO_ERG

    # Should be -100.0 eV (from second file, overriding -1.0 from first)
    assert (
        abs(exo_ev - (-100.0)) < 0.01
    ), f"Expected -100.0 eV (override), got {exo_ev} eV"

    print(f"\n✓ CH+ + E- -> C + H has overridden exothermicity: {exo_ev:.2f} eV")

    # Find H3+ + E- -> H2 + H reaction
    h3_recomb = None
    for reaction in network.get_reaction_list():
        reactants = [r for r in reaction.get_reactants() if r != "NAN"]
        products = [p for p in reaction.get_products() if p != "NAN"]
        if set(reactants) == {"H3+", "E-"} and set(products) == {"H2", "H"}:
            h3_recomb = reaction
            break

    assert h3_recomb is not None, "H3+ + E- -> H2 + H reaction not found"
    assert h3_recomb.get_exothermicity() is not None

    exo_ev_h3 = h3_recomb.get_exothermicity() / EV_TO_ERG

    # Should be -200.0 eV (from second file, overriding -2.0 from first)
    assert (
        abs(exo_ev_h3 - (-200.0)) < 0.01
    ), f"Expected -200.0 eV (override), got {exo_ev_h3} eV"

    print(f"✓ H3+ + E- -> H2 + H has overridden exothermicity: {exo_ev_h3:.2f} eV")


def test_network_exothermicity_without_custom_file(tmp_path):
    """Test that network works with derive_reaction_exothermicity but no custom files."""
    import os

    workspace_root = os.getcwd()

    config = {
        "species_file": f"{workspace_root}/Makerates/data/default/default_species.csv",
        "database_reaction_file": f"{workspace_root}/Makerates/data/databases/umist22.csv",
        "database_reaction_type": "UMIST",
        "custom_reaction_file": f"{workspace_root}/Makerates/data/default/default_grain_network.csv",
        "custom_reaction_type": "UCL",
        "add_crp_photo_to_grain": False,
        "enable_rates_storage": False,
        "three_phase": True,
        "gas_phase_extrapolation": False,
        "derive_reaction_exothermicity": "GAS",  # Enable exothermicity but no custom files
        "output_directory": str(tmp_path / "output"),
    }

    config_path = tmp_path / "test_config_no_custom.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    # This should work without errors (computed enthalpies only)
    network = get_network(path_to_input_file=config_path)

    assert network is not None
    assert len(network.get_species_list()) > 0
    assert len(network.get_reaction_list()) > 0

    print(
        "\n✓ Network builds successfully with derive_reaction_exothermicity=GAS and no custom files"
    )


def test_exothermicity_unit_conversion(config_file_with_exothermicity):
    """Test that exothermicity values are correctly converted to erg."""
    network = get_network(path_to_input_file=config_file_with_exothermicity)

    # All exothermicity values should be in erg
    # Check a few reactions with dummy incrementing values
    test_cases = [
        ("CH+", "E-", "C", "H", -1.0),  # eV
        ("H3+", "E-", "H2", "H", -2.0),  # eV
        ("C2H+", "E-", "C2", "H", -4.0),  # eV (note: -3.0 is C10H+ which may not exist)
    ]

    EV_TO_ERG = 1.602176634e-12

    for reactant1, reactant2, product1, product2, expected_ev in test_cases:
        reaction_found = None
        for reaction in network.get_reaction_list():
            reactants = [r for r in reaction.get_reactants() if r != "NAN"]
            products = [p for p in reaction.get_products() if p != "NAN"]
            if (
                reactant1 in reactants
                and reactant2 in reactants
                and product1 in products
                and product2 in products
            ):
                reaction_found = reaction
                break

        assert (
            reaction_found is not None
        ), f"{reactant1} + {reactant2} -> {product1} + {product2} not found"

        if reaction_found.get_exothermicity() is not None:
            exo_ev = reaction_found.get_exothermicity() / EV_TO_ERG
            assert abs(exo_ev - expected_ev) < 0.01, (
                f"{reactant1} + {reactant2} -> {product1} + {product2}: "
                f"Expected {expected_ev} eV, got {exo_ev} eV"
            )
            print(
                f"✓ {reactant1} + {reactant2} -> {product1} + {product2}: {exo_ev:.2f} eV"
            )


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
