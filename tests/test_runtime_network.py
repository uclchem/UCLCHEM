import numpy as np
import pytest
import uclchemwrap

import uclchem
import uclchem.advanced as advanced
from uclchem.makerates.reaction import Reaction


class TestRuntimeNetwork:
    """Test suite for RuntimeNetwork class."""

    @pytest.fixture
    def network(self) -> advanced.RuntimeNetwork:
        """Create a RuntimeNetwork instance for testing."""
        return advanced.RuntimeNetwork()

    def test_get_ice_list_index(self, network: advanced.RuntimeNetwork):
        ice_list_index = network._get_ice_list_index("#H")
        spec_index_from_ice_list = network.fortran_module.icelist[ice_list_index] - 1
        spec_index_from_spec_list = network._get_species_index("#H")
        assert spec_index_from_spec_list == spec_index_from_ice_list

    def test_initialization(self, network: advanced.RuntimeNetwork) -> None:
        """Test that RuntimeNetwork initializes correctly."""
        assert network is not None
        # Should have species and reactions from Fortran
        assert len(network.get_species_list()) > 0
        assert len(network.get_reaction_list()) > 0

    def test_inherits_from_network_abc(self, network: advanced.RuntimeNetwork) -> None:
        """Test that RuntimeNetwork inherits from NetworkABC."""
        from uclchem.makerates.network import NetworkABC

        assert isinstance(network, NetworkABC)

    def test_read_operations_work(self, network: advanced.RuntimeNetwork) -> None:
        """Test that all read operations work correctly."""
        # Test species reads
        species_list = network.get_species_list()
        species_dict = network.get_species_dict()
        assert len(species_list) == len(species_dict)

        # Test getting a specific species
        if len(species_list) > 0:
            first_species_name = species_list[0].get_name()
            species = network.get_specie(first_species_name)
            assert species is not None

        # Test reaction reads
        reactions_list = network.get_reaction_list()
        reactions_dict = network.get_reaction_dict()
        assert len(reactions_list) == len(reactions_dict)

        # Test getting a specific reaction
        if len(reactions_list) > 0:
            reaction = network.get_reaction(0)
            assert reaction is not None

    def test_modify_reaction_parameters(self, network: advanced.RuntimeNetwork) -> None:
        """Test that reaction parameter modification works."""
        if len(network.get_reaction_list()) == 0:
            pytest.skip("No reactions to test")

        # Get original values
        original_alpha = uclchemwrap.network.alpha[0]

        # Modify via RuntimeNetwork
        network.modify_reaction_parameters(0, alpha=999.0)

        # Verify Fortran was modified
        assert np.allclose(float(uclchemwrap.network.alpha[0]), 999.0)

        # Reset
        network.reset_to_initial_state()

        # Verify reset worked
        assert np.allclose(float(uclchemwrap.network.alpha[0]), float(original_alpha))

    def test_modify_reaction_parameters_coupled(
        self, network: advanced.RuntimeNetwork
    ) -> None:

        # Get a reaction that we know is in there for sure
        original_reaction = Reaction(["#H", "#H", "LH", "#H2", *["NAN"] * 8])
        expected_partners = [
            Reaction(["#H", "#H", "LHDES", "H2", *["NAN"] * 8]),
            Reaction(["@H", "@H", "LH", "@H2", *["NAN"] * 8]),
        ]
        partners = network.get_all_coupled_to_reaction(original_reaction)

        for partner in partners:
            assert partner in expected_partners

        reactions = network.get_reaction_list()
        original_idx = reactions.index(original_reaction)

        new_value = 999.0
        network.change_reaction_barrier(original_reaction, new_value)

        original_reaction_modified = network.get_reaction(original_idx)
        assert original_reaction_modified.get_gamma() == new_value
        partners = network.get_all_coupled_to_reaction(original_reaction)
        for partner in partners:
            assert partner.get_gamma() == new_value

        # Reset
        network.reset_to_initial_state()

    def test_disable_reaction(self, network: advanced.RuntimeNetwork) -> None:
        """Test that disabling a reaction sets alpha to zero."""
        if len(network.get_reaction_list()) == 0:
            pytest.skip("No reactions to test")

        # Disable a reaction
        network.disable_reaction(0)

        # Verify alpha is zero
        assert np.allclose(float(uclchemwrap.network.alpha[0]), 0.0)

        # Reset
        network.reset_to_initial_state()

    def test_add_species_not_supported(self, network: advanced.RuntimeNetwork) -> None:
        """Test that adding species raises NotImplementedError."""
        from uclchem.makerates.species import Species

        test_species = Species(["TEST", 10, 0.0, 0, 0, 0, 0])

        with pytest.raises(NotImplementedError, match="Cannot add species"):
            network.add_species(test_species)

    def test_add_reactions_not_supported(self, network: advanced.RuntimeNetwork) -> None:
        """Test that adding reactions raises NotImplementedError."""

        test_reaction = Reaction(
            ["H", "H", "NAN", "H2", "NAN", "NAN", "NAN", 1e-10, 0, 0, 0, 1e6, 0, 0, 0]
        )

        with pytest.raises(NotImplementedError, match="Cannot add reactions"):
            network.add_reactions(test_reaction)

    def test_remove_species_not_supported(self, network: advanced.RuntimeNetwork) -> None:
        """Test that removing species raises NotImplementedError."""
        if len(network.get_species_list()) == 0:
            pytest.skip("No species to test")

        first_species = network.get_species_list()[0].get_name()
        with pytest.raises(NotImplementedError, match="Cannot remove species"):
            network.remove_species(first_species)

    def test_remove_reaction_not_supported(
        self, network: advanced.RuntimeNetwork
    ) -> None:
        """Test that removing reactions raises NotImplementedError."""
        if len(network.get_reaction_list()) == 0:
            pytest.skip("No reactions to test")

        with pytest.raises(NotImplementedError, match="disable_reaction"):
            network.remove_reaction_by_index(0)

    def test_csv_validation(self, network: advanced.RuntimeNetwork) -> None:
        """Test that CSV files are validated against Fortran dimensions."""
        # This should pass during initialization
        # If there's a mismatch, ValueError would be raised in __init__
        assert hasattr(network, "_species_csv")
        assert hasattr(network, "_reactions_csv")
        assert len(network._species_csv) == len(network.get_species_list())
        assert len(network._reactions_csv) == len(network.get_reaction_list())

    def test_change_reaction_barrier_different_chemistry(
        self, network: advanced.RuntimeNetwork
    ) -> None:
        network.reset_to_initial_state()

        param_dict = {"initialTemp": 10, "initialDens": 1e4, "finalTime": 1e3}

        model = uclchem.model.Cloud(param_dict)
        _, chem_df = model.get_split_dataframes()

        reaction = Reaction(["#H", "#OH", "LH", "#H2O", "NAN", "NAN", "NAN", *[0] * 10])
        reaction_idx = network.get_reaction_index(reaction)
        current_reaction = network.get_reaction(reaction_idx)
        current_gamma = current_reaction.get_gamma()
        assert current_gamma < 10000

        network.change_reaction_barrier(reaction, 10000)

        reaction_idx = network.get_reaction_index(reaction)
        current_reaction = network.get_reaction(reaction_idx)
        current_gamma = current_reaction.get_gamma()
        assert current_gamma == 10000

        model_adj = uclchem.model.Cloud(param_dict)
        _, chem_df_adj = model_adj.get_split_dataframes()

        assert not np.allclose(
            chem_df.to_numpy(), chem_df_adj.to_numpy(), rtol=2, atol=1e-20
        )

        network.reset_to_initial_state()

    def test_change_diffusion_barrier_different_chemistry(
        self, network: advanced.RuntimeNetwork
    ) -> None:
        param_dict = {
            "initialTemp": 10,
            "initialDens": 1e4,
            "finalTime": 1e3,
            "useCustomDiffusionBarriers": True,
        }

        model = uclchem.model.Cloud(param_dict)
        _, chem_df = model.get_split_dataframes()

        network.change_diffusion_barrier("#H2", 10000)
        model = uclchem.model.Cloud(param_dict)
        _, chem_df_adj = model.get_split_dataframes()

        assert not np.allclose(
            chem_df.to_numpy(), chem_df_adj.to_numpy(), rtol=2, atol=1e-20
        )

        network.reset_to_initial_state()
