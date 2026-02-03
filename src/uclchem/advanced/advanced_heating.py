"""
Heating and cooling mechanism configuration for UCLCHEM.

This module provides class-based interface for accessing and modifying the Fortran
heating module that controls UCLCHEM's heating and cooling mechanisms during execution.

**Thread Safety Warning:**
HeatingSettings modifies global Fortran module state and is **NOT thread-safe**.
Do not use with multiprocessing, multithreading, or concurrent model runs.
Settings should only be modified during initialization, before running models.

Note: Changes made through HeatingSettings affect the global Fortran state and persist
across model runs in the same Python session.
"""

from typing import Dict

import numpy as np
from uclchemwrap import f2py_constants as f2py_constants_module
from uclchemwrap import heating as heating_module


class HeatingSettings:
    """Class-based wrapper for UCLCHEM heating and cooling configuration.

    This class provides access to the heating and cooling mechanism controls,
    solver parameters, and computed values from heating.f90.

    Attributes:
        PHOTOELECTRIC (dict): Photoelectric heating implementations
            - 'BAKES': Bakes & Tielens method (index 1)
            - 'WEINGARTNER': Weingartner method (index 2)
            (Only one can be enabled at a time)
        H2_FORMATION (int): Index for H2 formation heating (3)
        H2_PHOTODISSOCIATION (int): Index for H2 photodissociation heating (4)
        H2_FUV_PUMPING (int): Index for H2 FUV pumping heating (5)
        CARBON_IONIZATION (int): Index for carbon ionization heating (6)
        COSMIC_RAY (int): Index for cosmic ray heating (7)
        TURBULENT (int): Index for turbulent heating (8)
        GAS_GRAIN_COLLISIONS (int): Index for gas-grain collisional heating/cooling (9)

        ATOMIC_LINE_COOLING (int): Index for atomic line cooling
        H2_COLLISIONALLY_INDUCED (int): Index for H2 collisionally induced emission
        COMPTON_COOLING (int): Index for Compton cooling
        CONTINUUM_EMISSION (int): Index for continuum emission cooling
        MOLECULAR_LINE_COOLING (int): Index for molecular line cooling

        DUST_TEMP_HOCUK (int): Hocuk et al. 2017 dust temperature method
        DUST_TEMP_HOLLENBACH (int): Hollenbach 1991 dust temperature method

    Example:
        >>> from uclchem.advanced import HeatingSettings
        >>> settings = HeatingSettings()
        >>> settings.print_configuration()
        >>> # Switch photoelectric method (auto-disables Bakes when enabling Weingartner)
        >>> settings.set_heating_mechanism(settings.PHOTOELECTRIC['WEINGARTNER'], True)
        >>> # Disable H2 formation
        >>> settings.set_heating_mechanism(settings.H2_FORMATION, False)
    """

    # Heating mechanism indices (matching heating.f90)
    # Mechanisms with multiple implementations are dicts; only one can be enabled at a time
    # Fortran Indexed:
    PHOTOELECTRIC = {
        "BAKES": 1,  # Bakes & Tielens photoelectric heating
        "WEINGARTNER": 2,  # Weingartner photoelectric heating
    }
    H2_FORMATION = 3
    H2_PHOTODISSOCIATION = 4
    H2_FUV_PUMPING = 5
    CARBON_IONIZATION = 6
    COSMIC_RAY = 7
    TURBULENT = 8
    GAS_GRAIN_COLLISIONS = 9

    # Cooling mechanism indices (matching heating.f90)
    ATOMIC_LINE_COOLING = 1
    H2_COLLISIONALLY_INDUCED = 2
    COMPTON_COOLING = 3
    CONTINUUM_EMISSION = 4
    MOLECULAR_LINE_COOLING = 5

    # Dust-gas coupling methods
    DUST_TEMP_HOCUK = 1  # Hocuk et al. 2017 parametric formulation
    DUST_TEMP_HOLLENBACH = 2  # Hollenbach 1991 detailed balance method

    def __init__(self):
        """Initialize the HeatingSettings wrapper.

        This connects to the Fortran heating module and provides access to
        its configuration parameters. Captures the current state as the default
        configuration for reset_to_defaults().
        """
        self._heating_module = heating_module
        self._f2py_constants_module = f2py_constants_module

        # Build a mapping of mechanism IDs to their groups (for mutual exclusion)
        self._heating_groups = {}
        for attr_name in dir(self):
            if attr_name.isupper() and not attr_name.startswith("_"):
                attr = getattr(self, attr_name)
                if isinstance(attr, dict):
                    # This is a group with multiple implementations
                    for impl_id in attr.values():
                        self._heating_groups[impl_id] = list(attr.values())

        # Backup initial state for reset_to_defaults()
        self._default_heating_modules = self._heating_module.heating_modules.copy()
        self._default_cooling_modules = self._heating_module.cooling_modules.copy()
        self._default_dust_gas_coupling_method = (
            self._heating_module.dust_gas_coupling_method
        )
        self._default_line_solver_attempts = self._heating_module.line_solver_attempts
        self._default_pahabund = self._heating_module.pahabund
        self._default_coolant_data_dir = self._f2py_constants_module.coolantdatadir

    def set_heating_mechanism(self, mechanism_id: int, enabled: bool = True):
        """Enable or disable a specific heating mechanism.

        For mechanisms with multiple implementations (like PHOTOELECTRIC),
        enabling one automatically disables others in the same group.

        Args:
            mechanism_id: Index of the heating mechanism (1-9). Use class
                constants (e.g., self.PHOTOELECTRIC['BAKES'], self.H2_FORMATION)
            enabled: True to enable, False to disable. Default: True

        Example:
            >>> settings = HeatingSettings()
            >>> # Enable Weingartner photoelectric (auto-disables Bakes)
            >>> settings.set_heating_mechanism(settings.PHOTOELECTRIC['WEINGARTNER'], True)
            >>> # Or disable H2 formation
            >>> settings.set_heating_mechanism(settings.H2_FORMATION, False)
        """
        if not 1 <= mechanism_id <= self._heating_module.nheating:
            raise ValueError(
                f"mechanism_id must be between 1 and {self._heating_module.nheating}"
            )

        # If this mechanism is part of a mutually exclusive group and we're enabling it,
        # disable all others in the group first
        if enabled and mechanism_id in self._heating_groups:
            group = self._heating_groups[mechanism_id]
            for other_id in group:
                if other_id != mechanism_id:
                    self._heating_module.heating_modules[other_id - 1] = False

        # Now set the requested mechanism
        self._heating_module.heating_modules[mechanism_id - 1] = enabled

    def set_cooling_mechanism(self, mechanism_id: int, enabled: bool = True):
        """Enable or disable a specific cooling mechanism.

        Args:
            mechanism_id: Index of the cooling mechanism (1-5). Use class
                constants (e.g., self.ATOMIC_LINE_COOLING)
            enabled: True to enable, False to disable. Default: True

        Example:
            >>> settings = HeatingSettings()
            >>> settings.set_cooling_mechanism(settings.MOLECULAR_LINE_COOLING, False)
        """
        if not 1 <= mechanism_id <= self._heating_module.ncooling:
            raise ValueError(
                f"mechanism_id must be between 1 and {self._heating_module.ncooling}"
            )
        self._heating_module.cooling_modules[mechanism_id - 1] = enabled

    def get_heating_modules(self) -> Dict[str, bool]:
        """Obtain the state (on/off) of all heating mechanisms.

        Returns:
            Dictionary mapping mechanism names to their enabled state

        Example:
            >>> settings = HeatingSettings()
            >>> state = settings.get_heating_modules()
            >>> print(state['H2Formation'])
            True
        """
        labels = [
            str(np.char.decode(label)).strip()
            for label in self._heating_module.heatinglabels
        ]
        return {
            label: bool(self._heating_module.heating_modules[i])
            for i, label in enumerate(labels)
        }

    def get_cooling_modules(self) -> Dict[str, bool]:
        """Obtain the state (on/off) of all cooling mechanisms.

        Returns:
            Dictionary mapping mechanism names to their enabled state

        Example:
            >>> settings = HeatingSettings()
            >>> state = settings.get_cooling_modules()
            >>> print(state['AtomicLineCooling'])
            True
        """
        labels = [
            str(np.char.decode(label)).strip()
            for label in self._heating_module.coolinglabels
        ]
        return {
            label: bool(self._heating_module.cooling_modules[i])
            for i, label in enumerate(labels)
        }

    def set_dust_gas_coupling_method(self, method: int):
        """Set the dust-gas temperature coupling method.

        Args:
            method: Coupling method to use:
                - 1 (DUST_TEMP_HOCUK): Hocuk et al. 2017
                - 2 (DUST_TEMP_HOLLENBACH): Hollenbach 1991

        Example:
            >>> settings = HeatingSettings()
            >>> settings.set_dust_gas_coupling_method(settings.DUST_TEMP_HOLLENBACH)
        """
        if method not in [1, 2]:
            raise ValueError("method must be 1 (Hocuk) or 2 (Hollenbach)")
        self._heating_module.dust_gas_coupling_method = method

    def get_dust_gas_coupling_method(self) -> int:
        """Get the current dust-gas temperature coupling method.

        Returns:
            Current method (1=Hocuk, 2=Hollenbach)
        """
        return self._heating_module.dust_gas_coupling_method

    def set_line_solver_attempts(self, attempts: int):
        """Set the number of line cooling solver iterations.

        The line cooling is computed multiple times and the median is used.

        Args:
            attempts: Number of solver attempts (odd number recommended).
                Default is 5. Recommended range: 3-11.

        Example:
            >>> settings = HeatingSettings()
            >>> settings.set_line_solver_attempts(7)
        """
        if attempts < 1:
            raise ValueError("attempts must be at least 1")
        if attempts % 2 == 0:
            import warnings

            warnings.warn(
                "Even number of attempts may not produce proper median. "
                "Consider using an odd number."
            )
        self._heating_module.line_solver_attempts = attempts

    def get_line_solver_attempts(self) -> int:
        """Get the current number of line cooling solver attempts.

        Returns:
            Number of solver attempts
        """
        return self._heating_module.line_solver_attempts

    def set_pah_abundance(self, abundance: float):
        """Set the PAH (Polycyclic Aromatic Hydrocarbon) abundance.

        Args:
            abundance: PAH abundance relative to hydrogen. Default: 6e-7

        Example:
            >>> settings = HeatingSettings()
            >>> settings.set_pah_abundance(1e-6)
        """
        if abundance < 0:
            raise ValueError("abundance must be non-negative")
        self._heating_module.pahabund = abundance

    def get_pah_abundance(self) -> float:
        """Get the current PAH abundance.

        Returns:
            PAH abundance
        """
        return self._heating_module.pahabund

    def set_coolant_directory(self, directory: str):
        """Set the directory containing collisional rate data files.

        This directory must contain the LAMDA-format collisional rate files
        (e.g., co.dat, o-h2.dat, etc.) used for molecular line cooling calculations.

        Args:
            directory: Path to directory containing LAMDA-format rate files.
                      Must end with '/'. Max 255 characters.

        Raises:
            ValueError: If path is too long or doesn't end with '/'
            FileNotFoundError: If directory doesn't exist

        Example:
            >>> settings = HeatingSettings()
            >>> settings.set_coolant_directory("/custom/rates/")
        """
        from pathlib import Path

        # Validate
        if not directory.endswith("/"):
            directory += "/"

        dir_path = Path(directory)
        if not dir_path.exists():
            raise FileNotFoundError(f"Directory not found: {directory}")
        if not dir_path.is_dir():
            raise ValueError(f"Not a directory: {directory}")
        if len(directory) > 255:
            raise ValueError(f"Path too long (max 255): {directory}")

        # Pad and set
        current = self._f2py_constants_module.coolantdatadir
        max_len = int(current.dtype.itemsize)
        padded = directory.ljust(max_len)
        self._f2py_constants_module.coolantdatadir = padded

    def get_coolant_directory(self) -> str:
        """Get the current collisional rate data directory.

        Returns:
            Directory path as string (stripped of trailing spaces)
        """
        return str(np.char.decode(self._f2py_constants_module.coolantdatadir)).strip()

    def reset_to_defaults(self):
        """Reset all heating and cooling mechanisms to their initial values.

        Restores the configuration that was present when this HeatingSettings
        instance was first initialized. This allows reverting any changes made
        during the session back to the starting state.

        Example:
            >>> settings = HeatingSettings()
            >>> settings.set_heating_mechanism(settings.H2_FORMATION, False)
            >>> settings.reset_to_defaults()  # Restores original state
        """
        # Restore backed-up initial state
        self._heating_module.heating_modules[:] = self._default_heating_modules
        self._heating_module.cooling_modules[:] = self._default_cooling_modules
        self._heating_module.dust_gas_coupling_method = (
            self._default_dust_gas_coupling_method
        )
        self._heating_module.line_solver_attempts = self._default_line_solver_attempts
        self._heating_module.pahabund = self._default_pahabund
        self._f2py_constants_module.coolantdatadir = self._default_coolant_data_dir

    def print_configuration(self):
        """Print the current heating and cooling configuration.

        Example:
            >>> settings = HeatingSettings()
            >>> settings.print_configuration()
        """
        print("=" * 60)
        print("UCLCHEM Heating & Cooling Configuration")
        print("=" * 60)

        print("\nHeating Mechanisms:")
        print("-" * 40)
        heating_state = self.get_heating_modules()
        for name, enabled in heating_state.items():
            status = "ENABLED " if enabled else "DISABLED"
            print(f"  {name:30s} : {status}")

        print("\nCooling Mechanisms:")
        print("-" * 40)
        cooling_state = self.get_cooling_modules()
        for name, enabled in cooling_state.items():
            status = "ENABLED " if enabled else "DISABLED"
            print(f"  {name:30s} : {status}")

        print("\nSolver Parameters:")
        print("-" * 40)
        method = self.get_dust_gas_coupling_method()
        method_name = "Hocuk et al. 2017" if method == 1 else "Hollenbach 1991"
        print(f"  Dust-Gas Coupling Method       : {method} ({method_name})")
        print(f"  Line Solver Attempts           : {self.get_line_solver_attempts()}")
        print(f"  PAH Abundance                  : {self.get_pah_abundance():.2e}")
        print(f"  Coolant Data Directory         : {self.get_coolant_directory()}")
        print("=" * 60)
