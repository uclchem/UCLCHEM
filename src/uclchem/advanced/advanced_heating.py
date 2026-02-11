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
import uclchemwrap
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

        COOLANT_WARM (int): Warm restart mode - initialize LTE, rescale on density change (default)
        COOLANT_FORCE_LTE (int): Always reset to LTE before SE iteration (original behavior)
        COOLANT_FORCE_GROUND (int): Always reset to ground state before SE iteration

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

    # Coolant population restart modes
    COOLANT_WARM = 0  # Initialize to LTE on first call, then rescale on density change (default, most efficient)
    COOLANT_FORCE_LTE = 1  # Always reset to LTE before SE iteration (original behavior)
    COOLANT_FORCE_GROUND = 2  # Always reset to ground state before SE iteration

    def __init__(self):
        """Initialize the HeatingSettings wrapper.

        This connects to the Fortran heating module and provides access to
        its configuration parameters. Captures the current state as the default
        configuration for reset_to_defaults().
        """
        self._heating_module = heating_module
        self._f2py_constants_module = f2py_constants_module
        self._uclchemwrap = uclchemwrap

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

        # Check if coolant restart mode functions are available
        self._coolant_functions_available = hasattr(
            self._uclchemwrap, "get_coolant_restart_mode_wrap"
        )
        if self._coolant_functions_available:
            self._default_coolant_restart_mode = (
                self._uclchemwrap.get_coolant_restart_mode_wrap()
            )
        else:
            # Functions not yet exposed, use default value
            self._default_coolant_restart_mode = 0  # WARM mode

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

    # TODO: refactor once Fortran is exposed
    def set_coolant_restart_mode(self, mode: int):
        """Set the coolant population restart mode.

        Args:
            mode: Restart mode (0=WARM, 1=FORCE_LTE, 2=FORCE_GROUND)

        Example:
            >>> settings = HeatingSettings()
            >>> settings.set_coolant_restart_mode(settings.COOLANT_WARM)
        """
        if mode not in [0, 1, 2]:
            raise ValueError(f"mode must be 0, 1, or 2, got {mode}")
        self._uclchemwrap.uclchemwrap.set_coolant_restart_mode_wrap(mode)
        assert self.get_coolant_restart_mode() == mode, "Failed to set coolant restart mode"

    # TODO: refactor once Fortran is exposed
    def get_coolant_restart_mode(self) -> int:
        """Get the current coolant population restart mode.

        Returns:
            Current restart mode (0=WARM, 1=FORCE_LTE, 2=FORCE_GROUND)
        """
        return self._uclchemwrap.uclchemwrap.get_coolant_restart_mode_wrap()


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
        if self._coolant_functions_available:
            self._uclchemwrap.set_coolant_restart_mode_wrap(
                self._default_coolant_restart_mode
            )

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

        restart_mode = self.get_coolant_restart_mode()
        restart_names = {0: "WARM", 1: "FORCE_LTE", 2: "FORCE_GROUND"}
        restart_name = restart_names.get(restart_mode, "UNKNOWN")
        print(f"  Coolant Restart Mode           : {restart_mode} ({restart_name})")
        print("=" * 60)


def initialize_coolant_directory() -> str:
    """
    Locate and return the collisional rate data directory.

    This function searches for coolant data files in the following order:
    1. UCLCHEM_COOLANT_DATA environment variable (if set)
    2. Installed package data via importlib.resources (normal installation)
    3. Development mode: Makerates/data/collisional_rates/ (relative to project root)

    Returns:
        str: Absolute path to the coolant data directory (with trailing slash)

    Raises:
        RuntimeError: If the Fortran heating module is not available (not compiled)
        FileNotFoundError: If coolant data directory cannot be found in any location

    Example:
        >>> from uclchem.advanced import initialize_coolant_directory
        >>> coolant_dir = initialize_coolant_directory()
        >>> print(f"Coolant data at: {coolant_dir}")
    """
    import logging
    import os
    from pathlib import Path

    # Check if heating module is available
    try:
        import uclchemwrap
        from uclchemwrap import f2py_constants
    except (ImportError, AttributeError) as e:
        raise RuntimeError(
            "UCLCHEM heating module not available. "
            "The Fortran extension may not be compiled. "
            f"Install UCLCHEM with: pip install . (error: {e})"
        )

    # Priority 1: Environment variable
    env_dir = os.environ.get("UCLCHEM_COOLANT_DATA")
    if env_dir:
        env_path = Path(env_dir)
        if env_path.is_dir() and list(env_path.glob("*.dat")):
            coolant_dir = str(env_path.resolve())
            if not coolant_dir.endswith("/"):
                coolant_dir += "/"
            logging.info(f"Using coolant data from UCLCHEM_COOLANT_DATA: {coolant_dir}")
            return coolant_dir
        else:
            logging.warning(
                f"UCLCHEM_COOLANT_DATA set to {env_dir}, but directory not found or empty. "
                "Searching other locations..."
            )

    # Priority 2: Installed package data (importlib.resources for Python 3.9+)
    try:
        # Try new API first (Python 3.9+)
        try:
            from importlib.resources import files
            package_data_path = files("uclchem") / "data" / "collisional_rates"
            # Convert to Path object
            if hasattr(package_data_path, "as_posix"):  # Traversable
                package_data_path = Path(str(package_data_path))
        except (ImportError, TypeError):
            # Fallback to older API (Python 3.7-3.8)
            from importlib.resources import path as resource_path
            with resource_path("uclchem.data", "collisional_rates") as p:
                package_data_path = Path(p)

        if package_data_path.is_dir() and list(package_data_path.glob("*.dat")):
            coolant_dir = str(package_data_path.resolve())
            if not coolant_dir.endswith("/"):
                coolant_dir += "/"
            logging.debug(f"Using installed coolant data: {coolant_dir}")
            return coolant_dir
    except (ImportError, FileNotFoundError, AttributeError) as e:
        logging.debug(f"Installed package data not found: {e}")

    # Priority 3: Development mode - search for Makerates/data/collisional_rates/
    try:
        from uclchem.utils import UCLCHEM_ROOT_DIR

        # Try relative to UCLCHEM_ROOT_DIR (src/uclchem/)
        candidates = [
            UCLCHEM_ROOT_DIR.parent.parent / "Makerates" / "data" / "collisional_rates",  # from src/uclchem to project root
            Path.cwd() / "Makerates" / "data" / "collisional_rates",  # from current working directory
            Path.cwd().parent / "Makerates" / "data" / "collisional_rates",  # one level up
        ]

        for candidate in candidates:
            if candidate.is_dir() and list(candidate.glob("*.dat")):
                coolant_dir = str(candidate.resolve())
                if not coolant_dir.endswith("/"):
                    coolant_dir += "/"
                logging.info(f"Using development mode coolant data: {coolant_dir}")
                return coolant_dir
    except Exception as e:
        logging.debug(f"Development mode search failed: {e}")

    # Not found in any location
    raise FileNotFoundError(
        "Could not locate coolant data files (.dat files for collisional rates). "
        "Searched:\n"
        "  1. UCLCHEM_COOLANT_DATA environment variable\n"
        "  2. Installed package data (uclchem/data/collisional_rates/)\n"
        "  3. Development mode (Makerates/data/collisional_rates/)\n"
        "\n"
        "To fix:\n"
        "  - For installed package: Run 'python makerates.py' then 'pip install .'\n"
        "  - For development: Ensure Makerates/data/collisional_rates/*.dat files exist\n"
        "  - Or set UCLCHEM_COOLANT_DATA=/path/to/coolant/data/"
    )


def auto_initialize_coolant_directory() -> bool:
    """
    Automatically initialize the coolant data directory for the Fortran module.

    This is a convenience wrapper around initialize_coolant_directory() that:
    - Attempts to locate coolant data files
    - Sets the coolant directory in the Fortran module if found
    - Logs warnings instead of raising exceptions if initialization fails

    This function is called automatically when the uclchem module is imported.

    Returns:
        bool: True if initialization succeeded, False if it failed

    Example:
        >>> from uclchem.advanced import auto_initialize_coolant_directory
        >>> if auto_initialize_coolant_directory():
        ...     print("Coolant data initialized successfully")
    """
    import logging

    try:
        coolant_dir = initialize_coolant_directory()
        settings = HeatingSettings()
        settings.set_coolant_directory(coolant_dir)
        logging.debug(f"Auto-initialized coolant directory: {coolant_dir}")
        return True
    except RuntimeError as e:
        # Heating module not available - this is expected for makerates-only builds
        logging.debug(f"Coolant initialization skipped: {e}")
        return False
    except FileNotFoundError as e:
        # Could not find coolant data - warn user
        logging.warning(
            f"Could not auto-initialize coolant data directory: {e}\n"
            "Heating/cooling calculations may fail. "
            "Run 'python makerates.py' and reinstall if needed."
        )
        return False
    except Exception as e:
        # Unexpected error - warn but don't crash
        logging.warning(f"Unexpected error during coolant initialization: {e}" + 
                        "\nEnabling heating and cooling might cause errors at runtime.")
        return False
