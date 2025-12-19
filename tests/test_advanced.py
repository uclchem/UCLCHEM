"""
Unit tests for advanced module classes.

Tests the runtime interfaces for:
- HeatingSettings: Controlling heating and cooling mechanisms
- NetworkState: Accessing and resetting the chemical network
- GeneralSettings: Accessing and modifying UCLCHEM settings
"""

import numpy as np
import pytest
import uclchemwrap

from uclchem import advanced


@pytest.fixture(autouse=True)
def reset_fortran_state() -> None:
    """Reset Fortran module state after each test to ensure test isolation.

    This prevents state leakage between tests when they modify global Fortran values.
    """
    yield
    # Cleanup after test
    network = advanced.NetworkState()
    network.reset_state()


class TestHeatingSettings:
    """Test suite for HeatingSettings class."""

    @pytest.fixture
    def settings(self) -> advanced.HeatingSettings:
        """Create a HeatingSettings instance for testing."""
        return advanced.HeatingSettings()

    def test_initialization(self, settings: advanced.HeatingSettings) -> None:
        """Test that HeatingSettings initializes correctly on the Python side.

        Python side: Tests Python object attributes and structure.
        """
        assert settings is not None
        assert hasattr(settings, "PHOTOELECTRIC")
        assert hasattr(settings, "H2_FORMATION")
        assert hasattr(settings, "H2_PHOTODISSOCIATION")
        assert hasattr(settings, "H2_FUV_PUMPING")
        assert hasattr(settings, "CARBON_IONIZATION")
        assert hasattr(settings, "COSMIC_RAY")
        assert hasattr(settings, "TURBULENT")
        assert hasattr(settings, "GAS_GRAIN_COLLISIONS")
        assert hasattr(settings, "ATOMIC_LINE_COOLING")
        assert hasattr(settings, "H2_COLLISIONALLY_INDUCED")
        assert hasattr(settings, "COMPTON_COOLING")
        assert hasattr(settings, "CONTINUUM_EMISSION")
        assert hasattr(settings, "MOLECULAR_LINE_COOLING")
        assert hasattr(settings, "DUST_TEMP_HOCUK")
        assert hasattr(settings, "DUST_TEMP_HOLLENBACH")

        # Verify PHOTOELECTRIC is a dict with expected keys
        assert isinstance(settings.PHOTOELECTRIC, dict)
        assert "BAKES" in settings.PHOTOELECTRIC
        assert "WEINGARTNER" in settings.PHOTOELECTRIC

    def test_photoelectric_mutual_exclusion(
        self, settings: advanced.HeatingSettings
    ) -> None:
        """Test that photoelectric mechanisms are mutually exclusive.
        Fortran side: Tests Fortran heating module states via get_heating_modules().
        """
        settings.set_heating_mechanism(settings.PHOTOELECTRIC["WEINGARTNER"], True)

        # Check that Bakes is disabled
        heating_modules = settings.get_heating_modules()
        assert not heating_modules["PhotoelectricBakes"]
        assert heating_modules["PhotoelectricWeingartner"]

        # Enable Bakes method
        settings.set_heating_mechanism(settings.PHOTOELECTRIC["BAKES"], True)

        # Check that Weingartner is now disabled
        heating_modules = settings.get_heating_modules()
        assert heating_modules["PhotoelectricBakes"]
        assert not heating_modules["PhotoelectricWeingartner"]

    def test_toggle_heating_mechanism(self, settings: advanced.HeatingSettings) -> None:
        """Test enabling/disabling heating mechanisms.

        Fortran side: Tests toggling Fortran heating module states.
        """
        heating_modules = settings.get_heating_modules()

        # Disable mechanism
        settings.set_heating_mechanism(settings.H2_FORMATION, False)
        heating_modules = settings.get_heating_modules()
        assert not heating_modules["H2Formation"]

        # Re-enable mechanism
        settings.set_heating_mechanism(settings.H2_FORMATION, True)
        heating_modules = settings.get_heating_modules()
        assert heating_modules["H2Formation"]

    def test_toggle_cooling_mechanism(self, settings: advanced.HeatingSettings) -> None:
        """Test enabling/disabling cooling mechanisms.

        Fortran side: Tests toggling Fortran cooling module states.
        """
        settings.set_cooling_mechanism(settings.COMPTON_COOLING, False)
        cooling_modules = settings.get_cooling_modules()
        assert not cooling_modules["Compton"]

        # Re-enable cooling
        settings.set_cooling_mechanism(settings.COMPTON_COOLING, True)
        cooling_modules = settings.get_cooling_modules()
        assert cooling_modules["Compton"]

    def test_get_heating_modules(self, settings: advanced.HeatingSettings) -> None:
        """Test getting all heating mechanism states.

        Fortran side: Retrieves heating module states from Fortran wrapper.
        """
        heating_modules = settings.get_heating_modules()
        assert isinstance(heating_modules, dict)
        assert len(heating_modules) > 0
        # Check some expected keys
        assert (
            "H2Formation" in heating_modules and "PhotoelectricBakes" in heating_modules
        )

    def test_get_cooling_modules(self, settings: advanced.HeatingSettings) -> None:
        """Test getting all cooling mechanism states.

        Fortran side: Retrieves cooling module states from Fortran wrapper.
        """
        cooling_modules = settings.get_cooling_modules()
        assert isinstance(cooling_modules, dict)
        assert len(cooling_modules) >= 0  # May be empty if no cooling enabled

    def test_pah_abundance_access(self, settings: advanced.HeatingSettings) -> None:
        """Test reading and writing PAH abundance parameter.

        Fortran side: Tests get/set operations on Fortran parameter.
        """
        original = settings.get_pah_abundance()
        # Should be a numeric type (scalar or 0-d array)
        assert isinstance(original, (float, np.floating, np.ndarray))

        # Set new value
        set_value = 1.5e-7
        settings.set_pah_abundance(set_value)
        new_value = settings.get_pah_abundance()
        assert np.allclose(new_value, set_value)

    def test_line_solver_attempts_access(
        self, settings: advanced.HeatingSettings
    ) -> None:
        """Test reading and writing line solver attempts parameter.

        Fortran side: Tests get/set operations on Fortran parameter.
        """
        original = settings.get_line_solver_attempts()
        # Should be a numeric type (scalar or 0-d array)
        assert isinstance(original, (int, np.integer, np.ndarray))

        # Set new value
        settings.set_line_solver_attempts(9)
        new_value = settings.get_line_solver_attempts()
        assert new_value == 9

    def test_invalid_mechanism_index(self, settings: advanced.HeatingSettings) -> None:
        """Test that invalid mechanism indices raise errors.

        Python-Fortran integration: Tests error handling when accessing invalid Fortran array index.
        """
        with pytest.raises((IndexError, ValueError, RuntimeError)):
            settings.set_heating_mechanism(999, True)


class TestNetworkState:
    """Test suite for NetworkState class."""

    @pytest.fixture
    def network(self) -> advanced.NetworkState:
        """Create a NetworkState instance for testing."""
        return advanced.NetworkState()

    def test_initialization(self, network: advanced.NetworkState) -> None:
        """Test that NetworkState initializes correctly.

        Python side: Tests Python object attributes.
        """
        assert network is not None
        assert hasattr(network, "species_list")
        assert hasattr(network, "reaction_list")

    def test_species_list_loaded(self, network: advanced.NetworkState) -> None:
        """Test that species are loaded from CSV.

        Python-Fortran integration: Compares Python list size to Fortran network.
        """
        assert len(network.species_list) > 0
        # Should match compiled network size
        assert len(network.species_list) == len(uclchemwrap.network.specname)

    def test_reaction_list_loaded(self, network: advanced.NetworkState) -> None:
        """Test that reactions are loaded from CSV.

        Python-Fortran integration: Compares Python list size to Fortran network.
        """
        assert len(network.reaction_list) > 0
        # Should match compiled network size
        assert len(network.reaction_list) == len(uclchemwrap.network.alpha)

    def test_species_objects(self, network: advanced.NetworkState) -> None:
        """Test that species_list contains proper Species objects.

        Python side: Tests Python wrapper object attributes.
        """
        species = network.species_list[0]

        assert hasattr(species, "get_name")
        assert hasattr(species, "get_binding_energy")
        assert isinstance(species.get_name(), str)

    def test_reaction_objects(self, network: advanced.NetworkState) -> None:
        """Test that reaction_list contains proper Reaction objects.

        Python side: Tests Python wrapper object attributes.
        """
        reaction = network.reaction_list[0]

        assert hasattr(reaction, "get_alpha")
        assert hasattr(reaction, "get_beta")
        assert hasattr(reaction, "get_gamma")
        assert isinstance(reaction.get_alpha(), (int, float))

    def test_validation(self, network: advanced.NetworkState) -> None:
        """Test that validation passes for unmodified network.

        Python-Fortran integration: Python method validates Fortran network state.
        """
        network.validate()

    def test_reset_network(self, network: advanced.NetworkState) -> None:
        """Test resetting network from CSV.

        Fortran side: Modifies and resets Fortran network arrays.
        """
        import uclchemwrap

        original_alpha = float(uclchemwrap.network.alpha[0])
        uclchemwrap.network.alpha[0] = original_alpha * 10.0

        # Reset to initial state
        network.reset_state()

        # Verify it's reset
        new_alpha = float(uclchemwrap.network.alpha[0])
        assert np.allclose(original_alpha, new_alpha)

    def test_network_consistency(self, network: advanced.NetworkState) -> None:
        """Test that on-disk and in-memory networks are consistent."""
        # Python-Fortran integration: Validates consistency between Python and Fortran
        # This is tested by validate() in initialization - just verify counts match
        import uclchemwrap

        assert len(network.species_list) == len(uclchemwrap.network.specname)
        assert len(network.reaction_list) == len(uclchemwrap.network.alpha)

    def test_reset_preserves_counts(self, network: advanced.NetworkState) -> None:
        """Test that reset doesn't change array sizes.

        Python-Fortran integration: Verifies Python list sizes match Fortran after reset.
        """
        original_species_count = len(network.species_list)
        original_reaction_count = len(network.reaction_list)

        network.reset_state()

        assert len(network.species_list) == original_species_count
        assert len(network.reaction_list) == original_reaction_count

    def test_reset_restores_alpha_beta_gama(
        self, network: advanced.NetworkState
    ) -> None:
        """Test that reset restores alpha, beta, gama correctly.

        Fortran side: Tests Fortran array values are restored correctly.
        """
        import uclchemwrap

        # Store original values
        originals = [
            (
                float(uclchemwrap.network.alpha[i]),
                float(uclchemwrap.network.beta[i]),
                float(uclchemwrap.network.gama[i]),
            )
            for i in range(5)
        ]

        # Modify
        for i in range(5):
            uclchemwrap.network.alpha[i] *= 2.0
            uclchemwrap.network.beta[i] += 1.0
            uclchemwrap.network.gama[i] *= 0.5

        # Reset
        network.reset_state()

        # Verify restoration
        for i in range(5):
            assert np.allclose(float(uclchemwrap.network.alpha[i]), originals[i][0])
            assert np.allclose(float(uclchemwrap.network.beta[i]), originals[i][1])
            assert np.allclose(float(uclchemwrap.network.gama[i]), originals[i][2])

    def test_multiple_resets(self, network: advanced.NetworkState) -> None:
        """Test that multiple resets work correctly.

        Python-Fortran integration: Python method repeatedly resets Fortran state.
        """
        network.reset_state()
        network.validate()

        # Second reset
        network.reset_state()
        network.validate()

        # Should not raise errors


class TestRuntimeNetwork:
    """Test suite for RuntimeNetwork class."""

    @pytest.fixture
    def network(self) -> advanced.RuntimeNetwork:
        """Create a RuntimeNetwork instance for testing."""
        return advanced.RuntimeNetwork()

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
        from uclchem.makerates.reaction import Reaction
        test_reaction = Reaction(["H", "H", "NAN", "H2", "NAN", "NAN", "NAN", 
                                  1e-10, 0, 0, 0, 1e6, 0, 0, 0])
        
        with pytest.raises(NotImplementedError, match="Cannot add reactions"):
            network.add_reactions(test_reaction)

    def test_remove_species_not_supported(self, network: advanced.RuntimeNetwork) -> None:
        """Test that removing species raises NotImplementedError."""
        if len(network.get_species_list()) == 0:
            pytest.skip("No species to test")
        
        first_species = network.get_species_list()[0].get_name()
        with pytest.raises(NotImplementedError, match="Cannot remove species"):
            network.remove_species(first_species)

    def test_remove_reaction_not_supported(self, network: advanced.RuntimeNetwork) -> None:
        """Test that removing reactions raises NotImplementedError."""
        if len(network.get_reaction_list()) == 0:
            pytest.skip("No reactions to test")
        
        with pytest.raises(NotImplementedError, match="disable_reaction"):
            network.remove_reaction_by_index(0)

    def test_csv_validation(self, network: advanced.RuntimeNetwork) -> None:
        """Test that CSV files are validated against Fortran dimensions."""
        # This should pass during initialization
        # If there's a mismatch, ValueError would be raised in __init__
        assert hasattr(network, '_species_csv')
        assert hasattr(network, '_reactions_csv')
        assert len(network._species_csv) == len(network.get_species_list())
        assert len(network._reactions_csv) == len(network.get_reaction_list())


class TestGeneralSettings:
    """Test suite for GeneralSettings class."""

    @pytest.fixture
    def settings(self) -> advanced.GeneralSettings:
        """Create a GeneralSettings instance for testing."""
        return advanced.GeneralSettings()

    def test_initialization(self, settings: advanced.GeneralSettings) -> None:
        """Test that GeneralSettings initializes correctly.

        Python side: Tests Python object attributes.
        """
        assert settings is not None
        assert hasattr(settings, "defaultparameters")
        assert hasattr(settings, "network")
        assert hasattr(settings, "heating")

    def test_setting_get(self, settings: advanced.GeneralSettings) -> None:
        """Test getting setting values.

        Fortran side: Retrieves value from Fortran module.
        """
        initial_dens = settings.defaultparameters.initialdens.get()
        # Should be a numeric type (scalar or 0-d array)
        assert isinstance(initial_dens, (int, float, np.number, np.ndarray))
        assert initial_dens > 0

    def test_setting_set(self, settings: advanced.GeneralSettings) -> None:
        """Test setting new values.

        Fortran side: Sets value in Fortran module.
        """
        # Set new value
        set_value = 500.0
        settings.defaultparameters.initialdens.set(set_value)
        new_value = settings.defaultparameters.initialdens.get()
        assert np.allclose(new_value, set_value)

    def test_file_path_protection(self, settings: advanced.GeneralSettings) -> None:
        """Test that file path parameters cannot be modified via GeneralSettings.

        Python side: Tests Python error handling for file path parameters.
        """
        # Test outputfile
        with pytest.raises(
            RuntimeError, match="file paths should be set via param_dict"
        ):
            settings.defaultparameters.outputfile.set("test.dat")

        # Test abundsavefile
        with pytest.raises(
            RuntimeError, match="file paths should be set via param_dict"
        ):
            settings.defaultparameters.abundsavefile.set("test-abund.dat")

        # Test columnfile
        with pytest.raises(
            RuntimeError, match="file paths should be set via param_dict"
        ):
            settings.defaultparameters.columnfile.set("test-column.dat")

    def test_parameter_protection(self, settings: advanced.GeneralSettings) -> None:
        """Test that PARAMETER arrays cannot be modified.

        Python side: Tests Python error handling for protected Fortran PARAMETERs.
        """
        with pytest.raises(RuntimeError, match="PARAMETER"):
            settings.network.mintemps.set(123.0)

    def test_edited_tracking(self, settings: advanced.GeneralSettings) -> None:
        """Test that is_edited flag tracks modifications.

        Python side: Tests Python tracking flag (not Fortran state).
        """
        setting = settings.defaultparameters.initialdens
        original = setting.get()

        # Should not be edited initially (may be from previous tests)
        # So let's reset first
        setting.reset()
        assert not setting.is_edited

        # Modify and check
        setting.set(original * 2)
        assert setting.is_edited

        # Reset and check
        setting.reset()
        assert not setting.is_edited

    def test_reset_to_default(self, settings: advanced.GeneralSettings) -> None:
        """Test resetting settings to default values.

        Python-Fortran integration: Python method resets Fortran value to default.
        """
        setting = settings.defaultparameters.initialdens
        default = float(setting.default_value)

        # Modify
        setting.set(999.0)
        assert setting.get() != default

        # Reset
        setting.reset()
        assert np.allclose(setting.get(), default)
        assert not setting.is_edited

    def test_setting_attributes(self, settings: advanced.GeneralSettings) -> None:
        """Test Setting object attributes.

        Python side: Tests Python object attributes.
        """
        setting = settings.defaultparameters.initialdens

        assert hasattr(setting, "current_value")
        assert hasattr(setting, "default_value")
        assert hasattr(setting, "is_edited")
        assert hasattr(setting, "dtype")
        assert hasattr(setting, "is_parameter")
        assert hasattr(setting, "is_internal")

    def test_parameter_flag(self, settings: advanced.GeneralSettings) -> None:
        """Test that PARAMETER arrays are flagged correctly.

        Python side: Tests Python classification flag.
        """
        assert settings.network.mintemps.is_parameter

        # Regular variable
        assert not settings.defaultparameters.initialdens.is_parameter

    def test_internal_flag(self, settings: advanced.GeneralSettings) -> None:
        """Test that internal solver parameters are flagged.

        Python side: Tests Python classification flag.
        """
        if hasattr(settings.physicscore, "dstep"):
            assert settings.physicscore.dstep.is_internal

    def test_search_functionality(self, settings: advanced.GeneralSettings) -> None:
        """Test searching for settings across modules.

        Python side: Tests Python search logic.
        """
        results = settings.search("initial")

        assert isinstance(results, dict)
        assert len(results) > 0
        assert any("initial" in key.lower() for key in results.keys())

    def test_module_list_settings(self, settings: advanced.GeneralSettings) -> None:
        """Test listing settings in a module.

        Python side: Tests Python list generation.
        """
        settings_dict = settings.defaultparameters.list_settings()

        assert isinstance(settings_dict, dict)
        assert len(settings_dict) > 0
        assert "initialdens" in settings_dict

    def test_reset_all(self, settings: advanced.GeneralSettings) -> None:
        """Test resetting all settings to defaults.

        Python-Fortran integration: Python method resets all Fortran values.
        """
        settings.defaultparameters.initialdens.set(999.0)
        settings.defaultparameters.initialtemp.set(50.0)

        # Reset all (confirm=False to avoid stdin prompt)
        settings.reset_all(confirm=False)

        # Check they're reset
        assert not settings.defaultparameters.initialdens.is_edited
        assert not settings.defaultparameters.initialtemp.is_edited

    def test_memory_consistency_warning(
        self, settings: advanced.GeneralSettings
    ) -> None:
        """Test that get() warns if memory value changed externally.

        Python-Fortran integration: Tests Python tracking vs Fortran memory state.
        """
        import numpy as np
        import uclchemwrap

        setting = settings.defaultparameters.initialdens
        original = setting.get()
        original_float = float(original)  # Convert to Python float

        # Modify directly in Fortran (bypass Setting interface)
        new_value = original_float * 2
        uclchemwrap.defaultparameters.initialdens = new_value

        # Getting should warn when value changed externally
        with pytest.warns(
            UserWarning, match="has been modified outside of GeneralSettings"
        ):
            memory_value = setting.get(check_memory=True)
        assert np.allclose(float(memory_value), new_value)

        # Restore
        uclchemwrap.defaultparameters.initialdens = original_float

    def test_invalid_module_access(self, settings: advanced.GeneralSettings) -> None:
        """Test that accessing invalid modules raises error.

        Python side: Tests Python error handling.
        """
        with pytest.raises(AttributeError):
            _ = settings.nonexistent_module

    def test_invalid_setting_access(self, settings: advanced.GeneralSettings) -> None:
        """Test that accessing invalid settings raises error.

        Python side: Tests Python error handling.
        """
        with pytest.raises(AttributeError):
            _ = settings.defaultparameters.nonexistent_setting

    def test_temporary_changes_context_manager(
        self, settings: advanced.GeneralSettings
    ) -> None:
        """Test that temporary_changes context manager restores state.

        Python-Fortran integration: Python context manager saves/restores Fortran state.
        """
        settings.defaultparameters.initialdens.set(1000.0)
        settings.defaultparameters.initialtemp.set(20.0)
        initial_dens = settings.defaultparameters.initialdens.get()
        initial_temp = settings.defaultparameters.initialtemp.get()

        # Use context manager
        with settings.temporary_changes():
            settings.defaultparameters.initialdens.set(5000.0)
            settings.defaultparameters.initialtemp.set(50.0)
            assert np.allclose(settings.defaultparameters.initialdens.get(), 5000.0)
            assert np.allclose(settings.defaultparameters.initialtemp.get(), 50.0)

        # Verify restoration
        assert np.allclose(settings.defaultparameters.initialdens.get(), initial_dens)
        assert np.allclose(settings.defaultparameters.initialtemp.get(), initial_temp)

    def test_temporary_changes_exception_handling(
        self, settings: advanced.GeneralSettings
    ) -> None:
        """Test that temporary_changes restores even on exception.

        Python-Fortran integration: Tests context manager restoration on exception.
        """
        settings.defaultparameters.initialdens.set(1000.0)
        initial_dens = settings.defaultparameters.initialdens.get()

        # Use context manager with exception
        try:
            with settings.temporary_changes():
                settings.defaultparameters.initialdens.set(5000.0)
                raise ValueError("Test exception")
        except ValueError:
            pass

        # Verify restoration despite exception
        assert np.allclose(settings.defaultparameters.initialdens.get(), initial_dens)

    def test_temporary_changes_nested(self, settings: advanced.GeneralSettings) -> None:
        """Test that nested temporary_changes work correctly.

        Python-Fortran integration: Tests nested context manager with Fortran state.
        """
        settings.defaultparameters.initialdens.set(1000.0)
        initial_dens = settings.defaultparameters.initialdens.get()

        with settings.temporary_changes():
            settings.defaultparameters.initialdens.set(2000.0)
            mid_dens = settings.defaultparameters.initialdens.get()

            with settings.temporary_changes():
                settings.defaultparameters.initialdens.set(3000.0)
                assert np.allclose(settings.defaultparameters.initialdens.get(), 3000.0)

            # Should restore to mid_dens
            assert np.allclose(settings.defaultparameters.initialdens.get(), mid_dens)

        # Should restore to initial_dens
        assert np.allclose(settings.defaultparameters.initialdens.get(), initial_dens)


def test_network_reset_preserves_all_parameters():
    """Test that network reset restores all alpha, beta, gamma, and binding energies from CSV.

    This test retrieves all reaction parameters before reset, modifies them,
    then resets the network from CSV and verifies all parameters match their
    CSV-stored values exactly. Note: "exact" here means the values from CSV,
    which may differ in the last few digits from the initially compiled values
    due to CSV text format precision limits (~15-17 decimal digits).
    """
    network = advanced.NetworkState()

    # Get expected values from CSV (via reaction_list and species_list)
    # These are what reset should restore to
    csv_alpha = np.array([rxn.get_alpha() for rxn in network.reaction_list])
    csv_beta = np.array([rxn.get_beta() for rxn in network.reaction_list])
    csv_gamma = np.array([rxn.get_gamma() for rxn in network.reaction_list])

    # Get binding energies directly from the network
    n_binding = len(network._network.bindingenergy)
    csv_binding = np.array(
        [float(network._network.bindingenergy[i]) for i in range(n_binding)]
    )

    # Modify Fortran arrays directly
    for i in range(min(10, len(network._network.alpha))):
        network._network.alpha[i] *= 2.0
        network._network.beta[i] += 5.0
        network._network.gama[i] *= 0.5

    for i in range(min(10, n_binding)):
        network._network.bindingenergy[i] *= 1.5

    # Verify modifications took effect
    assert not np.allclose(
        network._network.alpha[:10], csv_alpha[:10]
    ), "Alpha modification didn't work"

    # Reset network to initial state
    network.reset_state()

    # Get values after reset
    reset_alpha = np.array(
        [float(network._network.alpha[i]) for i in range(len(csv_alpha))]
    )
    reset_beta = np.array(
        [float(network._network.beta[i]) for i in range(len(csv_beta))]
    )
    reset_gamma = np.array(
        [float(network._network.gama[i]) for i in range(len(csv_gamma))]
    )
    reset_binding = np.array(
        [float(network._network.bindingenergy[i]) for i in range(n_binding)]
    )

    # Check exact matches - arrays should be identical element-wise
    assert np.allclose(
        csv_alpha, reset_alpha
    ), "Alpha values not exactly matching CSV after reset"
    assert np.allclose(
        csv_beta, reset_beta
    ), "Beta values not exactly matching CSV after reset"
    assert np.allclose(
        csv_gamma, reset_gamma
    ), "Gamma values not exactly matching CSV after reset"
    assert np.allclose(
        csv_binding, reset_binding
    ), "Binding energies not exactly matching CSV after reset"


def test_network_reset_isolated():
    """Test that network reset operation works in isolation.

    This test isolates the network reset operation and verifies it completes
    without errors. It does not run a full model to avoid long test times.
    """
    # Reset network to initial state
    network = advanced.NetworkState()

    # Verify network is loaded
    assert len(network.reaction_list) > 0
    assert len(network.species_list) > 0

    # Reset should complete without raising exceptions
    network.reset_state()

    # Verify network still valid after reset
    assert len(network.reaction_list) > 0
    assert len(network.species_list) > 0

    # Verify some basic network properties are still intact
    assert len(network._network.alpha) > 0
    assert len(network._network.beta) > 0
    assert len(network._network.gama) > 0
    assert len(network._network.bindingenergy) > 0
