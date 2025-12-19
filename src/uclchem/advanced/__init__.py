"""
UCLCHEM Advanced Module

Runtime configuration and introspection for UCLCHEM's Fortran modules.

This package provides class-based interfaces for:
- HeatingSettings: Control heating and cooling mechanisms
- NetworkState: Access and modify the chemical network state durin runtime
- GeneralSettings: Access and modify all UCLCHEM settings

**Thread Safety Warning:**
All classes in this module modify global Fortran module state and are **NOT thread-safe**.
Do not use these classes with multiprocessing, multithreading, or concurrent model runs.
Settings should only be modified during initialization, before running models.
"""

# Import from package modules
from .advanced_heating import HeatingSettings
from .advanced_network import NetworkState, RuntimeReaction, RuntimeSpecies
from .advanced_settings import GeneralSettings, ModuleSettings, Setting
from .runtime_network import RuntimeNetwork

__all__ = [
    "HeatingSettings",
    "NetworkState",
    "GeneralSettings",
    "RuntimeNetwork",
]
