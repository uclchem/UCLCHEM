# TODO v4.0: Remove this module and all its usages.
"""
Compatibility layer for old Network and LoadedNetwork APIs.

This module provides backward compatibility wrappers for the old Network
and LoadedNetwork classes. It's kept separate and NOT imported by default
to allow for a clean breaking change in v4.0.

When ready to deprecate, this can be imported in __init__.py to provide
warnings and migration paths.
"""

import warnings
from pathlib import Path
from typing import Union

from .network import Network as NewNetwork
from .network import build_network
from .reaction import Reaction
from .species import Species


def Network(species: list[Species] = None, reactions: list[Reaction] = None, **kwargs):
    """Backward compatible Network constructor.

    This function provides compatibility with the old Network class constructor.
    It redirects to the appropriate new factory method based on the arguments.

    Deprecated:
        Use Network.build() or build_network() instead for new code.

    Args:
        species: List of Species objects
        reactions: List of Reaction objects
        **kwargs: Build options (gas_phase_extrapolation, etc.)

    Returns:
        Network: Network instance created via build()

    Examples:
        >>> # Old style (deprecated)
        >>> network = Network(species, reactions, gas_phase_extrapolation=True)

        >>> # New style (recommended)
        >>> network = Network.build(species, reactions, gas_phase_extrapolation=True)
        >>> # or
        >>> network = build_network(species, reactions, gas_phase_extrapolation=True)
    """
    if species is None or reactions is None:
        raise ValueError(
            "Network() requires species and reactions. "
            "For loading from CSV, use Network.from_csv() or load_network_from_csv() instead."
        )

    warnings.warn(
        "Calling Network(species, reactions, ...) is deprecated. "
        "Use Network.build(species, reactions, ...) or build_network(...) instead. "
        "This will be removed in v4.0.",
        DeprecationWarning,
        stacklevel=2,
    )

    return build_network(species, reactions, **kwargs)


class LoadedNetwork:
    """Backward compatible LoadedNetwork class.

    This class provides compatibility with the old LoadedNetwork API.
    It redirects to appropriate Network factory methods.

    Deprecated:
        Use Network.from_csv(), Network.from_lists(), or the module-level
        factory functions instead.
    """

    def __new__(
        cls,
        *,
        species: list[Species] = None,
        reactions: list[Reaction] = None,
        species_filepath: Union[str, Path] = None,
        reactions_filepath: Union[str, Path] = None,
    ):
        """Create a network using old LoadedNetwork API.

        Args:
            species: List of Species objects (use with reactions)
            reactions: List of Reaction objects (use with species)
            species_filepath: Path to species CSV (use with reactions_filepath)
            reactions_filepath: Path to reactions CSV (use with species_filepath)

        Returns:
            Network: Network instance created via appropriate factory method

        Examples:
            >>> # Old style (deprecated)
            >>> network = LoadedNetwork()
            >>> network = LoadedNetwork(species_filepath='s.csv', reactions_filepath='r.csv')
            >>> network = LoadedNetwork(species=sp_list, reactions=rx_list)

            >>> # New style (recommended)
            >>> network = Network.from_csv()
            >>> network = Network.from_csv('s.csv', 'r.csv')
            >>> network = Network.from_lists(sp_list, rx_list)
            >>> # or use module-level functions
            >>> network = load_network_from_csv()
            >>> network = create_network(sp_list, rx_list)
        """
        # Check for invalid combinations
        has_objects = species is not None or reactions is not None
        has_filepaths = species_filepath is not None or reactions_filepath is not None

        if has_objects and has_filepaths:
            raise ValueError(
                "Cannot provide both species/reactions objects and file paths. "
                "Use either (species=..., reactions=...) OR "
                "(species_filepath=..., reactions_filepath=...)"
            )

        # If objects are provided, ensure both are provided
        if has_objects:
            if species is None or reactions is None:
                raise ValueError(
                    "Both species and reactions must be provided together."
                )

            warnings.warn(
                "LoadedNetwork(species=..., reactions=...) is deprecated. "
                "Use Network.from_lists(species, reactions) or "
                "create_network(species, reactions) instead. "
                "This will be removed in v4.0.",
                DeprecationWarning,
                stacklevel=2,
            )

            return NewNetwork.from_lists(species, reactions)
        else:
            # Loading from files
            warnings.warn(
                "LoadedNetwork(species_filepath=..., reactions_filepath=...) is deprecated. "
                "Use Network.from_csv(species_path, reactions_path) or "
                "load_network_from_csv(species_path, reactions_path) instead. "
                "This will be removed in v4.0.",
                DeprecationWarning,
                stacklevel=2,
            )

            return NewNetwork.from_csv(species_filepath, reactions_filepath)


class NetworkState:
    """Backward compatible NetworkState class.

    This class provides compatibility with the old NetworkState API from
    uclchem.advanced. It redirects to Network.from_fortran().

    Deprecated:
        Use Network.from_fortran() or load_network_from_fortran() instead.
    """

    def __new__(cls):
        """Create a network with Fortran interface using old NetworkState API.

        Returns:
            Network: Network instance with Fortran interface

        Examples:
            >>> # Old style (deprecated)
            >>> from uclchem.advanced import NetworkState
            >>> network = NetworkState()
            >>> network._network.alpha[0] = 999.0

            >>> # New style (recommended)
            >>> from uclchem.makerates import Network
            >>> network = Network.from_fortran()
            >>> network.modify_fortran_alpha(0, 999.0)
            >>> # or for direct access
            >>> network.fortran.raw.alpha[0] = 999.0
            >>> # or use module-level function
            >>> network = load_network_from_fortran()
        """
        warnings.warn(
            "NetworkState is deprecated. "
            "Use Network.from_fortran() or load_network_from_fortran() instead. "
            "This will be removed in v4.0.",
            DeprecationWarning,
            stacklevel=2,
        )

        return NewNetwork.from_fortran()


# Compatibility exports for when this module is used
__all__ = [
    "Network",
    "LoadedNetwork",
    "NetworkState",
]
