"""Unified Network implementation for UCLCHEM.

This module provides the Network class and factory functions for creating networks
in different contexts:
- ``load_network_from_csv()``: Load compiled networks for analysis
- ``build_network()``: Build new networks with full validation
- ``create_network()``: Direct instantiation from objects

The Network class implements a complete interface for accessing and modifying
chemical reaction networks, suitable for build-time and analysis-time use.

For runtime parameter modification during model execution, use RuntimeNetwork
from uclchem.advanced.runtime_network instead.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from copy import deepcopy
from typing import TYPE_CHECKING, Any

import pandas as pd

from uclchem.makerates.reaction import Reaction, reaction_types
from uclchem.makerates.species import Species
from uclchem.utils import UCLCHEM_ROOT_DIR, check_expected_type

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from pathlib import Path

logger = logging.getLogger(__name__)

# ============================================================================
# Abstract Base Classes
# ============================================================================


class NetworkABC(ABC):
    """Base abstract class defining the read-only network interface.

    Defines operations common to ALL network types: reading data, querying,
    and modifying existing parameters. Does NOT include add/remove operations
    since some implementations (like RuntimeNetwork) have fixed structure.

    This interface is implemented by:
    - Network: Full interactive network (via MutableNetworkABC)
    - RuntimeNetwork: Fortran-backed runtime network (read + parameter modification only)
    """

    # Core Properties
    @property
    @abstractmethod
    def species(self) -> dict[str, Species]:
        """Get the species collection."""
        pass

    @property
    @abstractmethod
    def reactions(self) -> dict[int, Reaction]:
        """Get the reactions collection."""
        pass

    # Species Read Interface
    @abstractmethod
    def get_species_list(self) -> list[Species]:
        """Get all species in the network."""
        pass

    @abstractmethod
    def get_species_dict(self) -> dict[str, Species]:
        """Get the species dictionary."""
        pass

    @abstractmethod
    def get_specie(self, specie_name: str) -> Species:
        """Get a species by name (copy).

        Args:
            specie_name (str): species name

        Returns:
            Species: copy of Species instance.

        """
        pass

    # Reaction Read Interface
    @abstractmethod
    def get_reaction_list(self) -> list[Reaction]:
        """Get all reactions in the network."""
        pass

    @abstractmethod
    def get_reaction_dict(self) -> dict[int, Reaction]:
        """Get the reactions dictionary."""
        pass

    @abstractmethod
    def get_reaction(self, reaction_idx: int) -> Reaction:
        """Get a reaction by index (copy).

        Args:
            reaction_idx (int): Index of the reaction

        Returns:
            Reaction: copy of reaction with index reaction_idx.

        """
        pass

    # Query Methods
    @abstractmethod
    def get_reactions_by_types(self, reaction_type: str | list[str]) -> list[Reaction]:
        """Get all reactions of specific type(s).

        Args:
            reaction_type (str | list[str]): Single type or list of types to filter by

        Returns:
            list[Reaction]: List of reactions matching the type(s)
        """
        pass

    @abstractmethod
    def find_similar_reactions(self, reaction: Reaction) -> dict[int, Reaction]:
        """Find reactions with same reactants and products.

        Args:
            reaction (Reaction): Reaction to find similar reactions for

        Returns:
            similar (dict[int, Reaction]): Dictionary of {index: Reaction}
                for matching reactions
        """
        pass

    @abstractmethod
    def get_reaction_index(self, reaction: Reaction) -> int:
        """Get the index of a reaction in the network.

        Args:
            reaction (Reaction): Reaction to find

        Returns:
            int: Index of the reaction

        """
        pass

    # Parameter Modification Methods (modify existing, don't add/remove)
    @abstractmethod
    def change_binding_energy(self, specie: str, new_binding_energy: float) -> None:
        """Change binding energy of a species.

        Args:
            specie (str): string representation of species
            new_binding_energy (float): new binding energy in K

        """
        pass

    @abstractmethod
    def change_reaction_barrier(self, reaction: Reaction, barrier: float) -> None:
        """Change activation barrier of a reaction.

        Args:
            reaction (Reaction): Reaction to change.
            barrier (float): New reaction barrier in K

        """
        pass

    def __repr__(self) -> str:
        """String representation of the network.

        Returns:
            str: string representation of network.

        """
        n_species = len(self.get_species_list())
        n_reactions = len(self.get_reaction_list())
        return (
            f"{self.__class__.__name__}(species={n_species},\n reactions={n_reactions})"
        )


class MutableNetworkABC(NetworkABC):
    """Extended interface for networks that support full CRUD operations.

    Adds add/remove/set operations for species and reactions on top of the
    base NetworkABC interface. Used by the Network class for build-time and
    analysis-time operations.

    NOT implemented by RuntimeNetwork since the Fortran-backed network is read
    directly from the fortran files, only allowing for the editing of existing reaction
    and species parameters.
    """

    # Species Modification Interface
    @abstractmethod
    def add_species(self, species: Species | Sequence[Species | list]) -> None:
        """Add species to network.

        Args:
            species (Species | Sequence[Species | list]): Species object, list of Species,
                or CSV-style entries

        """
        pass

    @abstractmethod
    def remove_species(self, species: str | Species) -> None:
        """Remove a species from network.

        Args:
            species (str | Species): name of species to be deleted, or the Species instance

        Raises:
            ValueError: If no species ``species`` is in the Network.

        """
        pass

    @abstractmethod
    def set_specie(self, species_name: str, species: Species) -> None:
        """Set/update a species.

        Args:
            species_name (str): Name of species
            species (Species): Species to replace the old Species with

        """
        pass

    @abstractmethod
    def set_species_dict(self, new_species_dict: dict[str, Species]) -> None:
        """Replace entire species dictionary.

        Args:
            new_species_dict (dict[str, Species]): dictionary with keys the
                names of the species and values the Species instances.

        """
        pass

    @abstractmethod
    def sort_species(self) -> None:
        """Sort species by type and mass."""
        pass

    # Reaction Modification Interface
    @abstractmethod
    def add_reactions(self, reactions: Reaction | Sequence[list | Reaction]) -> None:
        """Add reactions to network.

        Args:
            reactions (Reaction | Sequence[list | Reaction]): Reaction object, list of Reactions,
                or CSV-style entries

        Raises:
            ValueError: If there is an error when converting the CSV-style entries to
                Reaction instances
            TypeError: If input was not a Reaction, list of Reaction instances, or
                CSV-style entries.
        """
        pass

    @abstractmethod
    def remove_reaction(self, reaction: Reaction) -> None:
        """Remove a reaction from network.

        Args:
            reaction (Reaction): Reaction to remove from the network.

        Raises:
            ValueError: If no matching reaction ``reaction`` is found in the network.
            RuntimeError: If multiple matching reactions are found in the network.

        """
        pass

    @abstractmethod
    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        """Remove a reaction by its index.

        Args:
            reaction_idx (int): Index of reaction to remove.

        """
        pass

    @abstractmethod
    def set_reaction(self, reaction_idx: int, reaction: Reaction) -> None:
        """Set/update a reaction at specific index.

        Args:
            reaction_idx (int): index of reaction to set.
            reaction (Reaction): new reaction instance.

        Raises:
            RuntimeError: If the number of reactions changes.

        """
        pass

    @abstractmethod
    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        """Replace the entire reaction dictionary.

        Args:
            new_dict (dict[int, Reaction]): Dictionary with keys indices and values
                Reaction instances.

        """
        pass

    @abstractmethod
    def sort_reactions(self) -> None:
        """Sort reactions by type and reactants."""
        pass


# ============================================================================
# Base Network Implementation
# ============================================================================


class BaseNetwork(NetworkABC):
    """Base implementation providing common network operations.

    Implements all read and query operations that are common between
    Network and RuntimeNetwork. Both classes store data in _species_dict
    and _reactions_dict, so this base class can implement all the shared logic.

    Subclasses only need to:
    1. Initialize _species_dict and _reactions_dict
    2. Implement modification methods (change_binding_energy, change_reaction_barrier)
    3. Optionally implement add/remove operations (MutableNetworkABC)
    """

    # Subclasses must define these
    _species_dict: dict[str, Species]
    _reactions_dict: dict[int, Reaction]

    # ========================================================================
    # Properties (NetworkABC Implementation)
    # ========================================================================

    @property
    def species(self) -> dict[str, Species]:
        """Get species dictionary.

        Returns:
            dict[str, Species]: species dictionary
        """
        return self._species_dict

    @property
    def reactions(self) -> dict[int, Reaction]:
        """Get reactions dictionary.

        Returns:
            dict[int, Reaction]: Reaction dictionary

        """
        return self._reactions_dict

    # ========================================================================
    # Species Read Interface (NetworkABC Implementation)
    # ========================================================================

    def get_species_list(self) -> list[Species]:
        """Get all species as a list.

        Returns:
            list[Species]: list of all species in the Network.

        """
        return list(self._species_dict.values())

    def get_species_dict(self) -> dict[str, Species]:
        """Get species dictionary (copy).

        Returns:
            dict[str, Species]: copy of species dictionary, with keys of the names
                of the species, and values their Species instances.

        """
        return deepcopy(self._species_dict)

    def get_specie(self, specie_name: str) -> Species:
        """Get a species by name (copy).

        Args:
            specie_name (str): species name

        Returns:
            Species: copy of Species instance.

        """
        return deepcopy(self._species_dict[specie_name])

    # ========================================================================
    # Reaction Read Interface (NetworkABC Implementation)
    # ========================================================================

    def get_reaction_list(self) -> list[Reaction]:
        """Get all reactions as a list.

        Returns:
            list[Reaction]: list of all reactions in the Network.

        """
        return list(self._reactions_dict.values())

    def get_reaction_dict(self) -> dict[int, Reaction]:
        """Get reactions dictionary (copy).

        Returns:
            dict[int, Reaction]: copy of reaction dictionary, with keys of the
                indices and values of the reactions.

        """
        return deepcopy(self._reactions_dict)

    def get_reaction(self, reaction_idx: int) -> Reaction:
        """Get a reaction by index (copy).

        Args:
            reaction_idx (int): Index of the reaction

        Returns:
            Reaction: copy of reaction with index reaction_idx.

        """
        return deepcopy(self._reactions_dict[reaction_idx])

    # ========================================================================
    # Query Methods (NetworkABC Implementation)
    # ========================================================================

    def get_reactions_by_types(self, reaction_type: str | list[str]) -> list[Reaction]:
        """Get all reactions of specific type(s).

        Args:
            reaction_type (str | list[str]): Single type or list of types to filter by

        Returns:
            list[Reaction]: List of reactions matching the type(s)

        """
        if isinstance(reaction_type, str):
            reaction_type = [reaction_type]

        return [
            reaction
            for reaction in self._reactions_dict.values()
            if reaction.get_reaction_type() in reaction_type
        ]

    def find_similar_reactions(self, reaction: Reaction) -> dict[int, Reaction]:
        """Find reactions with same reactants and products.

        Args:
            reaction (Reaction): Reaction to find similar reactions for

        Returns:
            similar (dict[int, Reaction]): Dictionary of {index: Reaction}
                for matching reactions

        """
        check_expected_type(reaction, Reaction)

        similar = {}

        target_reactants = set(reaction.get_reactants()) - {"NAN"}
        target_products = set(reaction.get_products()) - {"NAN"}

        for idx, other_reaction in self._reactions_dict.items():
            other_reactants = set(other_reaction.get_reactants()) - {"NAN"}
            other_products = set(other_reaction.get_products()) - {"NAN"}

            if other_reactants == target_reactants and other_products == target_products:
                similar[idx] = other_reaction

        return similar

    def get_reaction_index(self, reaction: Reaction) -> int:
        """Get the index of a reaction in the network.

        Args:
            reaction (Reaction): Reaction to find

        Returns:
            int: Index of the reaction

        Raises:
            ValueError: If reaction not found or multiple matches exist

        """
        similar = self.find_similar_reactions(reaction)

        if len(similar) == 0:
            msg = f"Reaction {reaction} not found in network"
            raise ValueError(msg)
        elif len(similar) > 1:
            msg = (
                f"Multiple reactions match {reaction}. "
                f"Found indices: {list(similar.keys())}"
            )
            raise ValueError(msg)

        return list(similar.keys())[0]


# ============================================================================
# Network Class
# ============================================================================


class Network(BaseNetwork, MutableNetworkABC):
    """Universal network representation for build and analysis contexts.

    A single Network class that serves all use cases:
    - Build time: Full validation and automatic reaction generation
    - Analysis time: Fast loading of compiled networks from CSV

    The Network class can be created via:
    - Direct instantiation: Network(species_dict, reaction_dict)
    - Factory methods: from_csv(), from_lists(), build()
    - Factory functions: load_network_from_csv(), build_network(), etc.

    All creation methods produce a Network instance that implements the
    full NetworkABC interface, ensuring consistent access patterns.

    For runtime parameter modification during model execution, use RuntimeNetwork
    from uclchem.advanced.runtime_network instead.

    Examples:
        >>> # Load for analysis
        >>> network = Network.from_csv()

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
        >>> network = Network.build(species_list, reactions_list, gas_phase_extrapolation=True)

    Attributes:
        _species_dict (dict[str, Species]): Internal species storage {name: Species}
        _reactions_dict (dict[int, Reaction]): Internal reaction storage {index: Reaction}

    """

    def __init__(
        self, species_dict: dict[str, Species], reaction_dict: dict[int, Reaction]
    ):
        """Initialize network with species and reactions.

        This is the low-level constructor. Most users should prefer factory methods:
        - Network.from_csv() for analysis
        - Network.from_lists() for direct instantiation
        - Network.build() for full build with validation

        Or use the module-level factory functions for clearer documentation.

        For runtime parameter modification, use RuntimeNetwork from
        uclchem.advanced.runtime_network instead.

        Args:
            species_dict (dict[str, Species]): Species dictionary {name: Species}
            reaction_dict (dict[int, Reaction]): Reaction dictionary {index: Reaction}

        """
        self._species_dict = species_dict
        self._reactions_dict = reaction_dict

    # ========================================================================
    # Factory Methods (Class Methods)
    # ========================================================================

    @classmethod
    def from_csv(
        cls,
        species_path: str | Path | None = None,
        reactions_path: str | Path | None = None,
    ) -> Network:
        """Load network from CSV files.

        Loads a pre-compiled network from CSV files without any validation
        or automatic generation. This is the primary method for loading
        networks for analysis purposes.

        Args:
            species_path (str | Path | None): Path to species CSV (None = use default installation).
                Default = None.
            reactions_path (str | Path | None): Path to reactions CSV (None = use default installation).
                Default = None.

        Returns:
            Network: Loaded network instance

        Examples:
            >>> # Load default compiled network
            >>> network = Network.from_csv()

            >>> # Load old/custom network for analysis
            >>> network = Network.from_csv('old/species.csv', 'old/reactions.csv') # doctest: +SKIP

        """
        # Use defaults if not provided
        if species_path is None:
            species_path = UCLCHEM_ROOT_DIR / "species.csv"
        if reactions_path is None:
            reactions_path = UCLCHEM_ROOT_DIR / "reactions.csv"

        logger.debug(f"Loading network from {species_path} and {reactions_path}")

        # Load CSVs
        species_data = pd.read_csv(species_path)
        reactions_data = pd.read_csv(reactions_path)

        # Parse into objects
        species_list = [Species(list(spec)) for idx, spec in species_data.iterrows()]
        reactions_list = [Reaction(list(reac)) for idx, reac in reactions_data.iterrows()]

        # Create dictionaries
        species_dict = {s.get_name(): s for s in species_list}
        reaction_dict = dict(enumerate(reactions_list))

        return cls(species_dict, reaction_dict)

    @classmethod
    def from_lists(
        cls,
        species: list[Species],
        reactions: list[Reaction],
    ) -> Network:
        """Create network directly from lists.

        Direct instantiation from species and reaction lists without any
        validation or automatic generation. Useful for programmatic network
        construction or as a base for NetworkBuilder.

        Args:
            species (list[Species]): List of Species objects
            reactions (list[Reaction]): List of Reaction objects

        Returns:
            Network: Network instance

        Example:
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
            >>> network = Network.from_lists(species_list, reactions_list)

        """
        species_dict = {s.get_name(): s for s in species}
        reaction_dict = dict(enumerate(reactions))
        return cls(species_dict, reaction_dict)

    @classmethod
    def build(
        cls,
        species: list[Species],
        reactions: list[Reaction],
        **build_options: Any,
    ) -> Network:
        """Build network with full validation and automatic generation.

        This is the primary method for building new networks with full validation,
        automatic reaction generation (freeze-out, desorption, bulk), branching
        ratio checks, and all build-time operations. Delegates to NetworkBuilder.

        Args:
            species (list[Species]): List of Species objects
            reactions (list[Reaction]): List of Reaction objects
            **build_options (Any): Options passed to NetworkBuilder:
                - user_defined_bulk: List of user-defined bulk species
                - gas_phase_extrapolation: bool (default False)
                - add_crp_photo_to_grain: bool (default False)
                - derive_reaction_exothermicity: list[str] or None
                - database_reaction_exothermicity: list[Union[str, Path]] or None

        Returns:
            Network: Fully built and validated network

        Examples:
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
            >>> network = Network.build(
            ...     species=species_list,
            ...     reactions=reactions_list,
            ...     gas_phase_extrapolation=True,
            ...     add_crp_photo_to_grain=True
            ... )

        """
        from uclchem.makerates.network_builder import NetworkBuilder  # noqa: PLC0415

        builder = NetworkBuilder(species, reactions, **build_options)
        return builder.build()

    # ========================================================================
    # Properties
    # ========================================================================

    @property
    def species(self) -> dict[str, Species]:
        """Get species dictionary.

        Returns:
            dict[str, Species]: species dictionary

        """
        return self._species_dict

    # Note: Read operations (get_species_list, get_reaction_list, etc.)
    # are inherited from BaseNetwork

    # ========================================================================
    # Species Mutation Interface (MutableNetworkABC Implementation)
    # ========================================================================

    def set_specie(self, species_name: str, species: Species) -> None:
        """Set/update a species.

        Args:
            species_name (str): Name of species
            species (Species): Species to replace the old Species with

        """
        self._species_dict[species_name] = species

    def set_species_dict(self, new_species_dict: dict[str, Species]) -> None:
        """Replace entire species dictionary.

        Args:
            new_species_dict (dict[str, Species]): dictionary with keys the
                names of the species and values the Species instances.

        """
        self._species_dict = new_species_dict

    def add_species(self, species: Species | Sequence[Species | list]) -> None:
        """Add species to network.

        Args:
            species (Species | Sequence[Species | list]): Species object, list of Species,
                or CSV-style entries

        Raises:
            ValueError: If there is an error when converting the CSV-style entries to
                Species instances
            TypeError: If input was not a Species, list of Species instances, or
                CSV-style entries.

        """
        species_list: list[Species]
        # Convert to list of Species objects
        if isinstance(species, list):
            if len(species) == 0:
                logger.warning("Tried to add empty species list, ignoring.")
                return
            elif isinstance(species[0], Species):
                species_list = species  # type: ignore
            elif isinstance(species[0], list):
                try:
                    species_list = [Species(spec) for spec in species]  # type: ignore
                except ValueError as error:
                    msg = "Failed to convert CSV entries to Species objects"
                    raise ValueError(msg) from error
        elif isinstance(species, Species):
            species_list = [species]
        else:
            msg = "Input must be Species object, list of Species, or CSV entries"
            raise TypeError(msg)

        # Add to dictionary
        for specie in species_list:
            # Filter out reaction types
            if specie.get_name() in reaction_types:
                logger.info(f"Ignoring reaction type {specie.get_name()} in species list")
                continue

            # Warn on duplicates
            if specie.get_name() in self._species_dict:
                logger.warning(
                    f"Species {specie.get_name()} already exists, keeping old definition"
                )
                continue

            # Filter out empty species
            if specie.get_name() in {"", "NAN"}:
                continue

            self._species_dict[specie.get_name()] = specie

    def remove_species(self, species: str | Species) -> None:
        """Remove a species from network.

        Args:
            species (str | Species): name of species to be deleted, or the Species instance

        Raises:
            ValueError: If no species ``species`` is in the Network.

        Examples:
            >>> network = Network.from_csv()
            >>>
            >>> network.remove_species("CO2")
            >>> # CO2 is now no longer in our network
            >>> # If we try to remove it again, we will get an error
            >>> network.remove_species("CO2")
            Traceback (most recent call last):
            ...
            ValueError: Species CO2 not found in network.

        """
        if isinstance(species, Species):
            species = species.get_name()
        check_expected_type(species, str)

        if species not in self._species_dict:
            msg = f"Species {species} not found in network."
            raise ValueError(msg)

        del self._species_dict[species]

    def sort_species(self) -> None:
        """Sort species by type and mass, with electron last."""
        species_dict = self.get_species_dict()

        self.set_species_dict(
            dict(
                sorted(
                    species_dict.items(),
                    key=lambda kv: (
                        kv[1].is_ice_species(),
                        kv[1].is_bulk_species(),
                        kv[1].get_mass(),
                    ),
                )
            )
        )

        # Move electron to end if present
        if "E-" in self._species_dict:
            electron = self.get_specie("E-")
            self.remove_species("E-")
            self.add_species(electron)

    # ========================================================================
    # Reaction Mutation Interface (MutableNetworkABC Implementation)
    # ========================================================================

    def set_reaction(self, reaction_idx: int, reaction: Reaction) -> None:
        """Set/update a reaction at specific index.

        Args:
            reaction_idx (int): index of reaction to set.
            reaction (Reaction): new reaction instance.

        Raises:
            RuntimeError: If the number of reactions changes.
        """
        check_expected_type(reaction, Reaction)

        old_length = len(self._reactions_dict)
        self._reactions_dict[reaction_idx] = reaction
        if old_length != len(self._reactions_dict):
            msg = "Setting the reaction caused a change in the number of reactions"
            raise RuntimeError(msg)

    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        """Replace the entire reaction dictionary.

        Args:
            new_dict (dict[int, Reaction]): Dictionary with keys indices and values
                Reaction instances.

        """
        check_expected_type(new_dict, dict)
        for reaction in new_dict.values():
            check_expected_type(reaction, Reaction)

        self._reactions_dict = new_dict

    def add_reactions(self, reactions: Reaction | Sequence[list | Reaction]) -> None:
        """Add reactions to network.

        Args:
            reactions (Reaction | Sequence[list | Reaction]): Reaction object, list of Reactions,
                or CSV-style entries

        Raises:
            ValueError: If there is an error when converting the CSV-style entries to
                Reaction instances
            TypeError: If input was not a Reaction, list of Reaction instances, or
                CSV-style entries.
        """
        # Convert to list of Reaction objects
        reactions_list: list[Reaction]
        if isinstance(reactions, list):
            if len(reactions) == 0:
                logger.warning("Tried to add empty reactions list, ignoring.")
                return
            elif isinstance(reactions[0], Reaction):
                reactions_list = reactions  # type: ignore
            elif isinstance(reactions[0], list):
                try:
                    reactions_list = [Reaction(reac) for reac in reactions]
                except ValueError as error:
                    msg = "Failed to convert CSV entries to Reaction objects"
                    raise ValueError(msg) from error
        elif isinstance(reactions, Reaction):
            reactions_list = [reactions]
        else:
            msg = "Input must be Reaction object, list of Reactions, or CSV entries"
            raise TypeError(msg)

        # Add to dictionary
        for reaction in reactions_list:
            if len(self._reactions_dict) == 0:
                new_idx = 0
            else:
                new_idx = max(self._reactions_dict.keys()) + 1

            self._reactions_dict[new_idx] = reaction

    def remove_reaction(self, reaction: Reaction) -> None:
        """Remove a reaction from network.

        Args:
            reaction (Reaction): Reaction to remove from the network.

        Raises:
            ValueError: If no matching reaction ``reaction`` is found in the network.
            RuntimeError: If multiple matching reactions are found in the network.

        """
        similar_reactions = list(self.find_similar_reactions(reaction).items())

        if len(similar_reactions) == 1:
            reaction_idx, _ = similar_reactions[0]
            del self._reactions_dict[reaction_idx]
        elif len(similar_reactions) == 0:
            msg = f"Reaction {reaction} not found in network"
            raise ValueError(msg)
        else:
            msg = (
                f"Found {len(similar_reactions)} reactions matching {reaction}. "
                "Use remove_reaction_by_index for piecewise reactions."
            )
            raise RuntimeError(msg)

    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        """Remove a reaction by its index.

        Args:
            reaction_idx (int): Index of reaction to remove.

        """
        if reaction_idx in self._reactions_dict:
            del self._reactions_dict[reaction_idx]
        else:
            logger.warning(f"Reaction index {reaction_idx} not found in network")

    def get_reactions_by_types(self, reaction_type: str | list[str]) -> list[Reaction]:
        """Get the union of all reactions of a certain type.

        Args:
            reaction_type (str | list[str]): The reaction type to filter on

        Returns:
            list[Reaction]: A list of reactions of the specified type

        """
        if isinstance(reaction_type, str):
            reaction_type = [reaction_type]
        return [
            r
            for r in self.get_reaction_list()
            if (r.get_reaction_type() in reaction_type)
        ]

    def sort_reactions(self) -> None:
        """Sort reactions by type and first reactant.

        Raises:
            RuntimeError: If the sorting of species causes the number of species
                in the network to change.
        """
        reaction_dict = self.get_reaction_dict()

        self.set_reaction_dict(
            dict(
                sorted(
                    reaction_dict.items(),
                    key=lambda kv: (
                        kv[1].get_reaction_type(),
                        kv[1].get_reactants()[0],
                    ),
                )
            )
        )
        if len(reaction_dict) != len(self.get_reaction_dict()):
            msg = "Sorting the species caused a difference in the number of species"
            raise RuntimeError(msg)

    # Note: Query methods (find_similar_reactions, get_reaction_index, etc.)
    # are inherited from BaseNetwork

    # ========================================================================
    # Parameter Modification Methods (NetworkABC Implementation)
    # ========================================================================

    def change_binding_energy(self, specie: str, new_binding_energy: float) -> None:
        """Change binding energy of a species.

        Handles special case of @H2O which affects other bulk species.

        Args:
            specie (str): string representation of species
            new_binding_energy (float): new binding energy in K

        Raises:
            ValueError: If `specie` is not in the network.
        """
        all_species = self.get_species_list()
        all_species_names = [s.get_name() for s in all_species]

        if specie not in all_species_names:
            msg = f"Species {specie} not found in network"
            raise ValueError(msg)

        # Special handling for @H2O (affects all bulk species)
        if specie == "@H2O":
            old_h2o_be = self._species_dict["@H2O"].get_binding_energy()
            for species_obj in all_species:
                if (
                    "@" in species_obj.get_name()
                    and species_obj.get_binding_energy() == old_h2o_be
                ):
                    species_obj.set_binding_energy(new_binding_energy)
        else:
            # Warn if changing bulk species that was H2O-limited
            if "@" in specie and "@H2O" in self._species_dict:
                h2o_be = self._species_dict["@H2O"].get_binding_energy()
                if self._species_dict[specie].get_binding_energy() == h2o_be:
                    logger.warning(
                        f"Changing binding energy of bulk species {specie} "
                        "that was previously @H2O binding energy limited"
                    )

            self._species_dict[specie].set_binding_energy(new_binding_energy)

    def change_reaction_barrier(self, reaction: Reaction, barrier: float) -> None:
        """Change activation barrier of a reaction.

        Looks up reaction in Network by its reactants and products.
        If Fortran interface is available, also updates Fortran.

        Args:
            reaction (Reaction): Reaction to change.
            barrier (float): New reaction barrier in K

        Raises:
            RuntimeError: If multiple matching reactions are found in the network.
        """
        similar_reactions = list(self.find_similar_reactions(reaction).items())

        if len(similar_reactions) == 1:
            reaction_idx, _ = similar_reactions[0]
            self._reactions_dict[reaction_idx].set_gamma(barrier)

        elif len(similar_reactions) == 0:
            logger.warning(f"Reaction {reaction} not found in network")
        else:
            msg = (
                f"Found {len(similar_reactions)} reactions matching {reaction}. "
                "Cannot uniquely identify which barrier to change."
            )
            raise RuntimeError(msg)

    def set_important_reaction_indices(self, indices: Mapping[str, int | None]) -> None:
        """Set the indices of important reactions that are treated differently in UCLCHEM.

        For example, H2 photodissociation has its own treatment in UCLCHEM.

        Not to be called by the user, but by ``uclchem.makerates.network_builder.NetworkBuilder```

        Args:
            indices (Mapping[str, int | None]): mapping from string representing the
                reaction to the index of that reaction in the network.

        Raises:
            ValueError: If any of the values in ``indices`` are None.

        """
        if any(index is None for index in indices.values()):
            msg = "Important reaction had index None"
            raise ValueError(msg)
        self.important_reactions = indices

    def set_important_species_indices(self, indices: Mapping[str, int | None]) -> None:
        """Set the indices of important species, like H2.

        Not to be called by the user, but by ``uclchem.makerates.network_builder.NetworkBuilder```

        Args:
            indices (Mapping[str, int | None]): mapping from name of the species to the
                index of that species in the network.

        Raises:
            ValueError: If any of the values in ``indices`` are None.

        """
        if any(index is None for index in indices.values()):
            msg = "Important species had index None"
            raise ValueError(msg)
        self.important_species = indices


# ============================================================================
# Factory Functions (Module-Level)
# ============================================================================


def load_network_from_csv(
    species_path: str | Path | None = None,
    reactions_path: str | Path | None = None,
) -> Network:
    """Load a network from CSV files for analysis.

    This is a module-level factory function that provides clear documentation
    and intent. It calls Network.from_csv() internally.

    Use this when analyzing pre-compiled networks, comparing network versions,
    or loading old networks for analysis.

    Args:
        species_path (str | Path | None): Path to species CSV (None = use default installation).
            Default = None.
        reactions_path (str | Path | None): Path to reactions CSV (None = use default installation).
            Default = None.

    Returns:
        Network: Loaded network instance

    Examples:
        >>> # Load default cmpiled network
        >>> network = load_network_from_csv()

        >>> # Load old version for comparison
        >>> old_network = load_network_from_csv(
        ...     'archive/v3.0/species.csv',
        ...     'archive/v3.0/reactions.csv'
        ... ) # doctest: +SKIP
        >>> print(f"Species added: {len(network.get_species_list()) - len(old_network.get_species_list())}") # doctest: +SKIP

    """  # noqa: W505
    return Network.from_csv(species_path, reactions_path)


def build_network(
    species: list[Species],
    reactions: list[Reaction],
    user_defined_bulk: list | None = None,
    gas_phase_extrapolation: bool = False,
    add_crp_photo_to_grain: bool = False,
    derive_reaction_exothermicity: list[str] | None = None,
    database_reaction_exothermicity: list[str | Path] | None = None,
) -> Network:
    """Build a new network with full validation and automatic generation.

    This is a module-level factory function that provides clear documentation
    and intent. It calls Network.build() internally.

    Use this when creating new networks at build time. This performs:
    - Input validation and duplicate checking
    - Automatic freeze-out reaction generation
    - Automatic bulk species and reaction generation
    - Automatic desorption reaction generation
    - Branching ratio validation and correction
    - Temperature range collision detection
    - Optional gas-phase extrapolation
    - Optional reaction exothermicity calculation

    Args:
        species (list[Species]): List of Species objects
        reactions (list[Reaction]): List of Reaction objects
        user_defined_bulk (list | None): User-specified bulk species. Default = None.
        gas_phase_extrapolation (bool): Extrapolate gas-phase reactions temperatures.
            Default = False.
        add_crp_photo_to_grain (bool): Add CRP/PHOTON reactions to grain surface.
            Default = False.
        derive_reaction_exothermicity (list[str] | None): Reaction types to calculate exothermicity for.
            Default = None.
        database_reaction_exothermicity (list[str | Path] | None): Custom exothermicity database files.
            Default = None.

    Returns:
        Network: Fully built and validated network

    Examples:
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

        >>> # Build network with standard options
        >>> network = build_network(
        ...     species=species_list,
        ...     reactions=reactions_list,
        ...     gas_phase_extrapolation=True
        ... )

        >>> # Build with custom exothermicity
        >>> network = build_network(
        ...     species=species_list,
        ...     reactions=reactions_list,
        ...     derive_reaction_exothermicity=['PHOTON', 'CRP'],
        ...     database_reaction_exothermicity=['custom_heating.csv']
        ... ) # doctest: +SKIP

    """
    return Network.build(
        species=species,
        reactions=reactions,
        user_defined_bulk=user_defined_bulk,
        gas_phase_extrapolation=gas_phase_extrapolation,
        add_crp_photo_to_grain=add_crp_photo_to_grain,
        derive_reaction_exothermicity=derive_reaction_exothermicity,
        database_reaction_exothermicity=database_reaction_exothermicity,
    )


def create_network(
    species: list[Species],
    reactions: list[Reaction],
) -> Network:
    """Create a network directly from species and reaction lists.

    This is a module-level factory function that provides clear documentation
    and intent. It calls Network.from_lists() internally.

    Use this when you want direct instantiation without validation or
    automatic generation. Suitable for programmatic network construction
    or when you already have a fully prepared network.

    Args:
        species (list[Species]): List of Species objects
        reactions (list[Reaction]): List of Reaction objects

    Returns:
        Network: Network instance

    Example:
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
        >>>
        >>> network = create_network(species_list, reactions_list)
        >>>
        >>> # Add some additional reactions
        >>> network.add_reactions(additional_reactions) # doctest: +SKIP

    """
    return Network.from_lists(species, reactions)
