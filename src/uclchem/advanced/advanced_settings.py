"""
General settings interface for UCLCHEM Fortran modules.

This module provides class-based interfaces for accessing and modifying runtime
settings across all UCLCHEM Fortran modules:

- Setting: Represents a single setting with metadata and edit tracking
- ModuleSettings: Container for all settings in a Fortran module
- GeneralSettings: Top-level interface to all modules

**Thread Safety Warning:**
All classes in this module modify global Fortran module state and are **NOT thread-safe**.
Do not use with multiprocessing, multithreading, or concurrent model runs.
Settings should only be modified during initialization, before running models.

Note: Changes made through these classes affect the global Fortran state and persist
across model runs in the same Python session.
"""

import warnings
from contextlib import contextmanager
from typing import Any, Dict, Set, Union

import numpy as np
import numpy.typing as npt
import uclchemwrap

# Import parameter classifications from constants module
from .constants import FILE_PATH_PARAMETERS, FORTRAN_PARAMETERS, INTERNAL_PARAMETERS


class Setting:
    """Represents a single runtime setting from a Fortran module.

    Tracks the current value, edit status, default value, and metadata
    for a Fortran module variable.

    Attributes:
        name: Setting name
        module_name: Parent Fortran module name
        current_value: Current value (cached on last read)
        is_edited: Whether value has been modified from default
        default_value: Original default value at initialization
        dtype: NumPy dtype or Python type
        is_parameter: True if this is a compile-time constant (read-only)
        is_internal: True if this is an internal solver parameter
        is_file_path: True if this is a file path (should use param_dict)
        shape: Array shape (None for scalars)
    """

    def __init__(
        self,
        name: str,
        module_name: str,
        fortran_module,
        is_parameter: bool = False,
        is_internal: bool = False,
        is_file_path: bool = False,
    ):
        """Initialize a Setting object.

        Args:
            name: Setting name
            module_name: Parent module name
            fortran_module: Reference to the Fortran module
            is_parameter: Whether this is a PARAMETER (read-only)
            is_internal: Whether this is an internal solver parameter
            is_file_path: Whether this is a file path parameter (should use param_dict)
        """
        self.name = name
        self.module_name = module_name
        self._fortran_module = fortran_module
        self.is_parameter = is_parameter
        self.is_internal = is_internal
        self.is_file_path = is_file_path

        # Read initial state
        value = getattr(fortran_module, name)
        self.default_value = self._copy_value(value)
        self.current_value = self._copy_value(value)
        self.is_edited = False

        # Store metadata
        if isinstance(value, np.ndarray):
            self.dtype = value.dtype
            self.shape = value.shape
        else:
            self.dtype = type(value)
            self.shape = None

    def _copy_value(self, value):
        """Make a copy of a value (handle arrays and scalars)."""
        if isinstance(value, np.ndarray):
            return value.copy()
        else:
            return value

    def get(self, check_memory: bool = True) -> Union[float, int, npt.NDArray]:
        """Get the current value of the setting.

        Args:
            check_memory: If True, compare cached value with actual Fortran memory
                         and warn if they differ

        Returns:
            Current value from Fortran memory
        """
        memory_value = getattr(self._fortran_module, self.name)

        if check_memory:
            # Compare with cached value
            if isinstance(memory_value, np.ndarray):
                if not np.array_equal(memory_value, self.current_value):
                    warnings.warn(
                        f"{self.module_name}.{self.name} has been modified "
                        f"outside of GeneralSettings (cache out of sync)",
                        UserWarning,
                        stacklevel=2,
                    )
            else:
                if memory_value != self.current_value:
                    warnings.warn(
                        f"{self.module_name}.{self.name} has been modified "
                        f"outside of GeneralSettings (cache out of sync)",
                        UserWarning,
                        stacklevel=2,
                    )

        # Update cache
        self.current_value = self._copy_value(memory_value)
        return memory_value

    def set(self, value: Union[float, int, npt.NDArray]) -> None:
        """Set the value of the setting.

        Args:
            value: New value to set

        Raises:
            RuntimeError: If attempting to modify a PARAMETER or file path parameter
        """
        if self.is_parameter:
            raise RuntimeError(
                f"Cannot modify {self.module_name}.{self.name}: "
                f"it is a Fortran PARAMETER (compile-time constant). "
                f"To change this value, you must edit the Fortran source, "
                f"regenerate with MakeRates, and rebuild."
            )

        if self.is_file_path:
            raise RuntimeError(
                f"Cannot modify {self.module_name}.{self.name}: "
                f"file paths should be set via param_dict when calling model functions "
                f"(e.g., uclchem.model.cloud(param_dict={{'outputFile': '...'}})). "
                f"File I/O paths are handled specially by the model wrapper and parsers."
            )

        # Set the value in Fortran memory
        setattr(self._fortran_module, self.name, value)

        # Update cache and edit status
        self.current_value = self._copy_value(value)
        self.is_edited = (
            not np.array_equal(self.current_value, self.default_value)
            if isinstance(self.current_value, np.ndarray)
            else (self.current_value != self.default_value)
        )

    def reset(self):
        """Reset the setting to its default value."""
        if self.is_parameter:
            raise RuntimeError(
                f"Cannot reset {self.module_name}.{self.name}: "
                f"it is a Fortran PARAMETER (compile-time constant)"
            )

        self.set(self.default_value)
        self.is_edited = False

    def __repr__(self) -> str:
        """String representation of the setting."""
        status = []
        if self.is_parameter:
            status.append("PARAMETER")
        if self.is_internal:
            status.append("INTERNAL")
        if self.is_file_path:
            status.append("FILE_PATH")
        if self.is_edited:
            status.append("EDITED")

        status_str = f" [{', '.join(status)}]" if status else ""

        if isinstance(self.current_value, np.ndarray):
            # For small arrays, show values; for large ones, just shape
            if self.current_value.size <= 3:
                value_str = str(self.current_value)
            else:
                value_str = f"<array shape={self.shape}, dtype={self.dtype}>"
        else:
            value_str = str(self.current_value)

        return f"Setting({self.module_name}.{self.name} = {value_str}{status_str})"


class ModuleSettings:
    """Container for all settings from a single Fortran module.

    Provides dict-like access to Setting objects with attribute-style syntax.
    """

    def __init__(
        self,
        module_name: str,
        fortran_module,
        parameter_names: Set[str],
        internal_names: Set[str],
        file_path_names: Set[str],
    ):
        """Initialize settings for a module.

        Args:
            module_name: Name of the Fortran module
            fortran_module: Reference to the actual Fortran module
            parameter_names: Set of names that are PARAMETERs
            internal_names: Set of names that are internal solver parameters
            file_path_names: Set of names that are file paths (should use param_dict)
        """
        self.module_name = module_name
        self._fortran_module = fortran_module
        self._settings = {}

        # Discover all attributes
        attrs = [a for a in dir(fortran_module) if not a.startswith("_")]

        for attr in attrs:
            try:
                value = getattr(fortran_module, attr)
                if callable(value):
                    continue

                # Determine classification
                is_param = attr.lower() in parameter_names
                is_internal = attr.lower() in internal_names
                is_file_path = attr.lower() in file_path_names

                # Create Setting object
                setting = Setting(
                    attr,
                    module_name,
                    fortran_module,
                    is_param,
                    is_internal,
                    is_file_path,
                )
                self._settings[attr.lower()] = setting

            except Exception:
                # Skip attributes that can't be accessed
                pass

    def __getattr__(self, name):
        """Get a Setting object by name."""
        if name.startswith("_"):
            return object.__getattribute__(self, name)

        name_lower = name.lower()
        if name_lower in self._settings:
            return self._settings[name_lower]

        raise AttributeError(
            f"Module '{self.module_name}' has no setting '{name}'. "
            f"Available: {', '.join(list(self._settings.keys())[:10])}..."
        )

    def __setattr__(self, name, value):
        """Set a setting value."""
        if name.startswith("_") or name in ["module_name"]:
            object.__setattr__(self, name, value)
            return

        name_lower = name.lower()
        if name_lower in self._settings:
            self._settings[name_lower].set(value)
        else:
            raise AttributeError(f"Module '{self.module_name}' has no setting '{name}'")

    def __dir__(self):
        """List all settings."""
        return list(self._settings.keys())

    def list_settings(
        self, include_internal: bool = False, include_parameters: bool = False
    ) -> Dict[str, "Setting"]:
        """List settings with their current values.

        Args:
            include_internal: Include internal solver parameters
            include_parameters: Include read-only PARAMETERs

        Returns:
            Dict mapping setting names to Setting objects
        """
        result = {}
        for name, setting in self._settings.items():
            if not include_internal and setting.is_internal:
                continue
            if not include_parameters and setting.is_parameter:
                continue
            result[name] = setting
        return result

    def print_settings(
        self, include_internal: bool = False, include_parameters: bool = False
    ):
        """Print all settings in a readable format."""
        print(f"\n{'=' * 70}")
        print(f"Module: {self.module_name}")
        print(f"{'=' * 70}\n")

        settings = self.list_settings(include_internal, include_parameters)
        for name in sorted(settings.keys()):
            setting = settings[name]
            flags = []
            if setting.is_parameter:
                flags.append("PARAM")
            if setting.is_internal:
                flags.append("INTERNAL")
            if setting.is_file_path:
                flags.append("FILE_PATH")
            if setting.is_edited:
                flags.append("EDITED")
            flag_str = f" [{','.join(flags)}]" if flags else ""

            if isinstance(setting.current_value, np.ndarray):
                if setting.current_value.size <= 5:
                    value_str = str(setting.current_value)
                else:
                    value_str = f"<array shape={setting.shape}>"
            else:
                value_str = str(setting.current_value)

            print(f"  {name:35s} = {value_str}{flag_str}")
        print()


class GeneralSettings:
    """General interface to all UCLCHEM settings across all Fortran modules.

    Provides dynamic access to modifiable parameters in any uclchemwrap module,
    with automatic detection of read-only PARAMETERs and internal solver settings.

    Each setting is represented by a Setting object that tracks:
    - Current value (with caching)
    - Edit status
    - Default value
    - Metadata (dtype, shape, flags)

    **Thread Safety Warning:**
        This class modifies global Fortran module state and is **NOT thread-safe**.
        Do not use with multiprocessing or concurrent model runs. Settings should
        only be modified during initialization, before running models.

    Usage:
        >>> from uclchem.advanced import GeneralSettings
        >>> settings = GeneralSettings()
        >>>
        >>> # Access a setting
        >>> setting = settings.defaultparameters.initialdens
        >>> print(setting.get())
        >>>
        >>> # Modify a setting
        >>> setting.set(500.0)
        >>>
        >>> # Or use attribute-style syntax
        >>> settings.defaultparameters.initialdens = 500.0
        >>>
        >>> # List all settings in a module
        >>> settings.defaultparameters.print_settings()
        >>>
        >>> # Use context manager for temporary changes
        >>> with settings.temporary_changes():
        ...     settings.defaultparameters.initialdens = 1000.0
        ...     # Run model with modified settings
        >>> # Settings automatically restored after context
    """

    def __init__(self):
        """Initialize with all available uclchemwrap modules."""
        self._modules = {}
        self._discover_modules()

    def _discover_modules(self):
        """Discover all available modules in uclchemwrap."""
        # Known module names
        module_names = [
            "defaultparameters",
            "network",
            "heating",
            "physicscore",
            "constants",
            "cloud_mod",
            "collapse_mod",
            "cshock_mod",
            "jshock_mod",
            "hotcore",
            "chemistry",
            "rates",
            "photoreactions",
            "surfacereactions",
            "io",
            "f2py_constants",
            "postprocess_mod",
            "sputtering",
        ]

        for name in module_names:
            if hasattr(uclchemwrap, name):
                module = getattr(uclchemwrap, name)
                self._modules[name] = ModuleSettings(
                    name,
                    module,
                    FORTRAN_PARAMETERS,
                    INTERNAL_PARAMETERS,
                    FILE_PATH_PARAMETERS,
                )

    def __getattr__(self, name):
        """Access modules as attributes."""
        if name.startswith("_"):
            return object.__getattribute__(self, name)

        if name in self._modules:
            return self._modules[name]

        raise AttributeError(
            f"No module '{name}' available.\n"
            f"Available modules: {', '.join(sorted(self._modules.keys()))}"
        )

    def __dir__(self):
        """List all available modules."""
        return list(self._modules.keys())

    def list_modules(self):
        """List all available modules with statistics."""
        print(f"\n{'=' * 70}")
        print("Available UCLCHEM Modules")
        print(f"{'=' * 70}\n")

        for name in sorted(self._modules.keys()):
            mod_settings = self._modules[name]
            all_settings = mod_settings._settings

            n_user = sum(
                1
                for s in all_settings.values()
                if not s.is_parameter and not s.is_internal
            )
            n_param = sum(1 for s in all_settings.values() if s.is_parameter)
            n_internal = sum(1 for s in all_settings.values() if s.is_internal)

            print(
                f"  {name:25s} {n_user:3d} user settings, "
                f"{n_param:3d} parameters, {n_internal:3d} internal"
            )
        print()

    def search(
        self,
        pattern: str,
        include_internal: bool = False,
        include_parameters: bool = False,
    ) -> Dict[str, Setting]:
        """Search for settings matching a pattern across all modules.

        Args:
            pattern: String pattern to search for (case-insensitive)
            include_internal: Include internal solver parameters
            include_parameters: Include read-only PARAMETERs

        Returns:
            Dict mapping "module.setting" to Setting objects
        """
        pattern = pattern.lower()
        results = {}

        for module_name, mod_settings in self._modules.items():
            for setting_name, setting in mod_settings._settings.items():
                if not include_internal and setting.is_internal:
                    continue
                if not include_parameters and setting.is_parameter:
                    continue

                if pattern in setting_name:
                    full_name = f"{module_name}.{setting_name}"
                    results[full_name] = setting

        return results

    def print_all_edited(self):
        """Print all settings that have been modified from defaults."""
        print(f"\n{'=' * 70}")
        print("Modified Settings")
        print(f"{'=' * 70}\n")

        found_any = False
        for module_name in sorted(self._modules.keys()):
            mod_settings = self._modules[module_name]
            edited = {
                name: s for name, s in mod_settings._settings.items() if s.is_edited
            }

            if edited:
                found_any = True
                print(f"\n{module_name}:")
                for name in sorted(edited.keys()):
                    setting = edited[name]
                    print(
                        f"  {name:35s} {setting.default_value} → {setting.current_value}"
                    )

        if not found_any:
            print("  No settings have been modified")
        print()

    def print_all_settings(self):
        """Print all settings across all modules."""
        for module_name in sorted(self._modules.keys()):
            mod_settings = self._modules[module_name]
            mod_settings.print_settings(include_internal=True, include_parameters=True)

    def reset_all(self, confirm: bool = True):
        """Reset all settings to their default values.

        Args:
            confirm: If True, require user confirmation
        """
        if confirm:
            response = input("Reset ALL settings to defaults? (yes/no): ")
            if response.lower() != "yes":
                print("Reset cancelled")
                return

        for mod_settings in self._modules.values():
            for setting in mod_settings._settings.values():
                if not setting.is_parameter and setting.is_edited:
                    try:
                        setting.reset()
                    except Exception as e:
                        print(f"Warning: Could not reset {setting.name}: {e}")

        print("✓ All settings reset to defaults")

    @contextmanager
    def temporary_changes(self):
        """Context manager for temporary setting modifications.

        Saves current state on entry and restores it on exit, even if an
        exception occurs. Useful for running models with temporary parameter
        changes without affecting the global state.

        Example:
            >>> settings = GeneralSettings()
            >>> settings.defaultparameters.initialdens = 1000.0
            >>> print(settings.defaultparameters.initialdens.get())  # 1000.0
            >>>
            >>> with settings.temporary_changes():
            ...     settings.defaultparameters.initialdens = 5000.0
            ...     print(settings.defaultparameters.initialdens.get())  # 5000.0
            ...     # Run model here
            >>>
            >>> print(settings.defaultparameters.initialdens.get())  # 1000.0 (restored)

        Yields:
            self: The GeneralSettings instance for chaining
        """
        # Save all current values (not just edited ones)
        # Skip PARAMETERs (compile-time constants), INTERNAL_PARAMETERS (solver state), and arrays
        saved_states: Dict[str, Dict[str, Any]] = {}
        for module_name, mod_settings in self._modules.items():
            saved_states[module_name] = {}
            for setting_name, setting in mod_settings._settings.items():
                # Skip PARAMETERs (can't be modified), internal parameters (unsafe to restore),
                # and arrays (Fortran memory management issues)
                if (
                    setting.is_parameter
                    or setting.is_internal
                    or setting.shape is not None
                ):
                    continue
                try:
                    saved_states[module_name][setting_name] = {
                        "value": setting._copy_value(setting.get(check_memory=False)),
                        "is_edited": setting.is_edited,
                    }
                except Exception as e:
                    # Warn about settings that can't be safely saved
                    warnings.warn(
                        f"Could not save {module_name}.{setting_name} for temporary_changes: {e}",
                        UserWarning,
                        stacklevel=2,
                    )

        try:
            yield self
        finally:
            # Restore all values
            for module_name, settings_dict in saved_states.items():
                mod_settings = self._modules[module_name]
                for setting_name, state in settings_dict.items():
                    setting = mod_settings._settings[setting_name]
                    try:
                        setting.set(state["value"])
                        setting.is_edited = state["is_edited"]
                    except Exception as e:
                        warnings.warn(
                            f"Could not restore {module_name}.{setting_name}: {e}",
                            UserWarning,
                            stacklevel=2,
                        )
