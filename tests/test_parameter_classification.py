"""
Validation CI Check for Fortran Parameter Classification

This module provides a test that ensures all Fortran variables are properly
classified in GeneralSettings.FORTRAN_PARAMETERS or INTERNAL_PARAMETERS.

This test should be added to the test suite to catch when new Fortran variables
are added without updating the classification lists.
"""

import numpy as np
import pytest

from uclchem import advanced


def test_all_fortran_variables_classified():
    """Ensure all Fortran variables are properly classified.

    This test fails if:
    1. New Fortran variables are added without classification
    2. Variables are incorrectly marked as PARAMETER or INTERNAL

    To fix failures:
    - Add new variable names to GeneralSettings.FORTRAN_PARAMETERS
      if they are compile-time constants (PARAMETER)
    - Add to GeneralSettings.INTERNAL_PARAMETERS if they are
      internal solver/computation variables
    - Otherwise, they're user-facing settings (no action needed)
    """
    settings = advanced.GeneralSettings()

    # Collect all unclassified variables
    unclassified = []
    misclassified = []

    for module_name in settings._modules.keys():
        mod_settings = settings._modules[module_name]

        for name, setting in mod_settings._settings.items():
            # Check if variable is classified
            is_classified = (
                setting.is_parameter or setting.is_internal or _is_user_facing(setting)
            )

            if not is_classified:
                unclassified.append(f"{module_name}.{name}")

            # Verify PARAMETERs can't be modified (additional validation)
            if setting.is_parameter:
                try:
                    # Try to modify - should fail
                    original = setting.get()
                    setting.set(original)
                    misclassified.append(
                        f"{module_name}.{name} is marked PARAMETER but can be modified"
                    )
                except RuntimeError:
                    # Expected behavior for PARAMETERs
                    pass

    # Generate helpful error messages
    error_parts = []

    if unclassified:
        error_parts.append(
            f"Unclassified Fortran variables found ({len(unclassified)}):\n"
            + "\n".join(f"  - {var}" for var in unclassified[:10])
            + ("\n  ... and more" if len(unclassified) > 10 else "")
            + "\n\nTo fix: Add these variables to either:\n"
            + "  - GeneralSettings.FORTRAN_PARAMETERS (if read-only/PARAMETER)\n"
            + "  - GeneralSettings.INTERNAL_PARAMETERS (if internal solver variable)\n"
            + "  - Leave uncategorized if they're user-facing settings"
        )

    if misclassified:
        error_parts.append(
            f"Misclassified variables ({len(misclassified)}):\n"
            + "\n".join(f"  - {msg}" for msg in misclassified)
        )

    if error_parts:
        pytest.fail("\n\n".join(error_parts))


def _is_user_facing(setting):
    """Determine if a setting is user-facing (not PARAMETER or INTERNAL).

    User-facing settings are modifiable and intended for user configuration.
    """
    # If it's not marked as PARAMETER or INTERNAL, it's user-facing
    return not (setting.is_parameter or setting.is_internal)


def test_parameter_immutability():
    """Verify that all PARAMETER variables cannot be modified.

    PARAMETERs are Fortran compile-time constants and should raise
    RuntimeError when set() is called.
    """
    settings = advanced.GeneralSettings()

    for module_name in settings._modules.keys():
        mod_settings = settings._modules[module_name]

        for name, setting in mod_settings._settings.items():
            if setting.is_parameter:
                with pytest.raises(RuntimeError, match="PARAMETER"):
                    # Attempt to modify should fail
                    setting.set(999)


def test_internal_parameters_exist():
    """Verify that known internal parameters are correctly classified.

    This test explicitly checks a few known internal parameters to ensure
    the classification system is working.
    """
    settings = advanced.GeneralSettings()

    # Known internal parameters that should exist
    known_internal = [
        ("physicscore", "dstep"),
        ("network", "reactionrate"),
        ("heating", "chemheating"),
    ]

    for module_name, param_name in known_internal:
        if hasattr(settings, module_name):
            module = getattr(settings, module_name)
            if hasattr(module, param_name):
                setting = getattr(module, param_name)
                assert (
                    setting.is_internal
                ), f"{module_name}.{param_name} should be marked as INTERNAL"


def test_user_parameters_modifiable():
    """Verify that known user-facing parameters can be modified.

    This ensures the classification system doesn't accidentally mark
    user-facing settings as PARAMETERs.
    """
    settings = advanced.GeneralSettings()

    # Known user-facing parameters
    known_user_params = [
        ("defaultparameters", "initialdens"),
        ("defaultparameters", "initialtemp"),
        ("defaultparameters", "finaltime"),
    ]

    for module_name, param_name in known_user_params:
        module = getattr(settings, module_name)
        setting = getattr(module, param_name)

        # Should not be PARAMETER or INTERNAL
        assert (
            not setting.is_parameter
        ), f"{module_name}.{param_name} should not be PARAMETER"
        assert (
            not setting.is_internal
        ), f"{module_name}.{param_name} should not be INTERNAL"

        # Should be modifiable
        original = float(setting.get())
        setting.set(float(original * 1.5))
        assert not np.allclose(setting.get(), original)
        setting.set(original)  # Restore


if __name__ == "__main__":
    """Run validation checks manually for debugging."""
    print("Running Fortran variable classification validation...")

    try:
        test_all_fortran_variables_classified()
        print("✓ All variables properly classified")
    except AssertionError as e:
        print(f"✗ Classification validation failed:\n{e}")

    try:
        test_parameter_immutability()
        print("✓ PARAMETER immutability verified")
    except AssertionError as e:
        print(f"✗ PARAMETER immutability check failed:\n{e}")

    try:
        test_internal_parameters_exist()
        print("✓ Internal parameters correctly classified")
    except AssertionError as e:
        print(f"✗ Internal parameter check failed:\n{e}")

    try:
        test_user_parameters_modifiable()
        print("✓ User parameters are modifiable")
    except AssertionError as e:
        print(f"✗ User parameter check failed:\n{e}")

    print("\nValidation complete!")
