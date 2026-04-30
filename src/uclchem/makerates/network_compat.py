# TODO v4.0: Remove this module and all its usages.
"""Compatibility layer for old Network and LoadedNetwork APIs.

This module provides backward compatibility wrappers for the old Network
and LoadedNetwork classes. It's kept separate and NOT imported by default
to allow for a clean breaking change in v4.0.

When ready to deprecate, this can be imported in __init__.py to provide
warnings and migration paths.

"""

import warnings
from pathlib import Path
from typing import Any

from uclchem.makerates.network import Network as NewNetwork
from uclchem.makerates.network import build_network
from uclchem.makerates.reaction import Reaction
from uclchem.makerates.species import Species


def Network(  # noqa: N802
    species: list[Species] | None = None,
    reactions: list[Reaction] | None = None,
    **kwargs: Any,
) -> NewNetwork:
    """Backward compatible Network constructor.

    This function provides compatibility with the old Network class constructor.
    It redirects to the appropriate new factory method based on the arguments.

    Deprecated:
        Use Network.build() or build_network() instead for new code.

    Parameters
    ----------
    species : list[Species] | None
        List of Species objects. Default = None.
    reactions : list[Reaction] | None
        List of Reaction objects. Default = None.
    **kwargs : Any
        Build options (gas_phase_extrapolation, etc.)

    Returns
    -------
    NewNetwork
        Network instance created via build_network()

    Raises
    ------
    ValueError
        If `species` or `reactions` is None.

    Examples
    --------
    >>> # Build with validation
    >>> from uclchem.makerates.io_functions import read_species_file, read_reaction_file
    >>> from uclchem.utils import UCLCHEM_ROOT_DIR
    >>>
    >>> species_list, user_defined_bulk = read_species_file(
    ...     UCLCHEM_ROOT_DIR / "../../Makerates/data/default/default_species.csv"
    ... )
    >>> reactions_list, dropped_reactions = read_reaction_file(
    ...     UCLCHEM_ROOT_DIR / "../../Makerates/data/default/default_grain_network.csv",
    ...     species_list,
    ...     "UCL",
    ... )

    >>> # Old style (deprecated)
    >>> network = Network(species_list, reactions_list, gas_phase_extrapolation=True)

    >>> # New style (recommended)
    >>> from uclchem.makerates.network import Network, build_network
    >>> network = Network.build(species_list, reactions_list, gas_phase_extrapolation=True)
    >>> # or
    >>> network = build_network(species_list, reactions_list, gas_phase_extrapolation=True)

    """
    if species is None or reactions is None:
        msg = (
            "Network() requires species and reactions. "
            "For loading from CSV, use Network.from_csv() or load_network_from_csv() instead."
        )
        raise ValueError(msg)

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

    def __new__(  # type: ignore[misc]
        cls,
        *,
        species: list[Species] | None = None,
        reactions: list[Reaction] | None = None,
        species_filepath: str | Path | None = None,
        reactions_filepath: str | Path | None = None,
    ) -> NewNetwork:
        """Create a network using old LoadedNetwork API.

        Parameters
        ----------
        cls : Any
            Not used
        species : list[Species] | None
            List of Species objects (use with reactions).
            Default = None.
        reactions : list[Reaction] | None
            List of Reaction objects (use with species).
            Default = None.
        species_filepath : str | Path | None
            Path to species CSV (use with reactions_filepath).
            Default = None.
        reactions_filepath : str | Path | None
            Path to reactions CSV (use with species_filepath).
            Default = None.

        Returns
        -------
        NewNetwork
            Network instance created via appropriate factory method

        Raises
        ------
        ValueError
            If both `species` and `reactions` and file paths are specified,
            or if only species or reactions is provided.

        Examples
        --------
        >>> # Build with validation
        >>> from uclchem.makerates.io_functions import read_species_file, read_reaction_file
        >>> from uclchem.utils import UCLCHEM_ROOT_DIR
        >>>
        >>> species_path = UCLCHEM_ROOT_DIR / "species.csv"
        >>> reactions_path = UCLCHEM_ROOT_DIR / "reactions.csv"

        >>> species_list, user_defined_bulk = read_species_file(species_path)
        >>> reactions_list, dropped_reactions = read_reaction_file(
        ...     reactions_path, species_list, "UCL"
        ... )

        >>> # Old style (deprecated)
        >>> network = LoadedNetwork()
        >>> network = LoadedNetwork(
        ...     species_filepath=species_path,
        ...     reactions_filepath=reactions_path
        ... )
        >>> network = LoadedNetwork(species=species_list, reactions=reactions_list)

        >>> # New style (recommended)
        >>> from uclchem.makerates.network import Network
        >>> network = Network.from_csv()
        >>> network = Network.from_csv(species_path, reactions_path)
        >>> network = Network.from_lists(species_list, reactions_list)
        >>> # or use module-level functions
        >>> from uclchem.makerates.network import load_network_from_csv, create_network
        >>> network = load_network_from_csv()
        >>> network = create_network(species_list, reactions_list)

        """
        # Check for invalid combinations
        has_objects = species is not None or reactions is not None
        has_filepaths = species_filepath is not None or reactions_filepath is not None

        if has_objects and has_filepaths:
            msg = (
                "Cannot provide both species/reactions objects and file paths. "
                "Use either (species=..., reactions=...) OR "
                "(species_filepath=..., reactions_filepath=...)"
            )
            raise ValueError(msg)

        # If objects are provided, ensure both are provided
        if has_objects:
            if species is None or reactions is None:
                msg = "Both species and reactions must be provided together."
                raise ValueError(msg)

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
    uclchem.advanced. It redirects to ``RuntimeNetwork()`.

    Deprecated:
        Use ``RuntimeNetwork()`` instead.

    """

    def __new__(cls):
        """Create a network with Fortran interface using old NetworkState API.

        Parameters
        ----------
        cls : Any
            Not used

        Returns
        -------
        Network
            Network instance with Fortran interface

        Examples
        --------
        >>> # Old style (deprecated)
        >>> from uclchem.advanced import NetworkState
        >>> network = NetworkState()
        >>> network._network.alpha[0] = 999.0

        >>> # New style (recommended)
        >>> from uclchem.advanced.runtime_network import RuntimeNetwork
        >>>
        >>> network = RuntimeNetwork()
        >>> network.modify_reaction_parameters(0, alpha=999.0)
        >>> # or for direct access
        >>> network._fortran.alpha[0] = 999.0
        >>>

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
    "LoadedNetwork",
    "Network",
    "NetworkState",
]
