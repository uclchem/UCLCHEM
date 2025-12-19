"""
Unified Network implementation for UCLCHEM.

This module provides the Network class and factory functions for creating networks
in different contexts:
- load_network_from_csv(): Load compiled networks for analysis
- build_network(): Build new networks with full validation
- create_network(): Direct instantiation from objects

The Network class implements a complete interface for accessing and modifying
chemical reaction networks, suitable for build-time and analysis-time use.

For runtime parameter modification during model execution, use RuntimeNetwork
from uclchem.advanced.runtime_network instead.
"""

import logging
from abc import ABC, abstractmethod
from copy import deepcopy
from os import path
from pathlib import Path
from typing import Union

import pandas as pd

from uclchem.utils import _ROOT

from .reaction import Reaction, reaction_types
from .species import Species

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
    def species(self):
        """Get the species collection."""
        pass

    @property
    @abstractmethod
    def reactions(self):
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
        """Get a specific species by name."""
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
        """Get a specific reaction by index."""
        pass

    # Query Methods
    @abstractmethod
    def get_reactions_by_types(
        self, reaction_type: Union[str, list[str]]
    ) -> list[Reaction]:
        """Get all reactions of specific type(s)."""
        pass

    @abstractmethod
    def find_similar_reactions(self, reaction: Reaction) -> dict[int, Reaction]:
        """Find reactions with same reactants/products."""
        pass

    @abstractmethod
    def get_reaction_index(self, reaction: Reaction) -> int:
        """Get the unique index of a reaction."""
        pass

    # Parameter Modification Methods (modify existing, don't add/remove)
    @abstractmethod
    def change_binding_energy(self, specie: str, new_binding_energy: float) -> None:
        """Change binding energy of an existing species."""
        pass

    @abstractmethod
    def change_reaction_barrier(self, reaction: Reaction, barrier: float) -> None:
        """Change activation barrier of an existing reaction."""
        pass

    def __repr__(self) -> str:
        """String representation of the network."""
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
    def add_species(self, species: Union[Species, list[Species]]) -> None:
        """Add one or more species to the network."""
        pass

    @abstractmethod
    def remove_species(self, specie_name: str) -> None:
        """Remove a species from the network."""
        pass

    @abstractmethod
    def set_specie(self, species_name: str, species: Species) -> None:
        """Set/update a species in the network."""
        pass

    @abstractmethod
    def set_species_dict(self, new_species_dict: dict[str, Species]) -> None:
        """Replace the entire species dictionary."""
        pass

    @abstractmethod
    def sort_species(self) -> None:
        """Sort species by type and mass."""
        pass

    # Reaction Modification Interface
    @abstractmethod
    def add_reactions(self, reactions: Union[Reaction, list[Reaction]]) -> None:
        """Add one or more reactions to the network."""
        pass

    @abstractmethod
    def remove_reaction(self, reaction: Reaction) -> None:
        """Remove a reaction from the network."""
        pass

    @abstractmethod
    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        """Remove a reaction by its index."""
        pass

    @abstractmethod
    def set_reaction(self, reaction_idx: int, reaction: Reaction) -> None:
        """Set/update a reaction at a specific index."""
        pass

    @abstractmethod
    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        """Replace the entire reaction dictionary."""
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
    def species(self):
        """Get species dictionary."""
        return self._species_dict

    @property
    def reactions(self):
        """Get reactions dictionary."""
        return self._reactions_dict

    # ========================================================================
    # Species Read Interface (NetworkABC Implementation)
    # ========================================================================

    def get_species_list(self) -> list[Species]:
        """Get all species as a list."""
        return list(self._species_dict.values())

    def get_species_dict(self) -> dict[str, Species]:
        """Get species dictionary (copy)."""
        return deepcopy(self._species_dict)

    def get_specie(self, specie_name: str) -> Species:
        """Get a species by name (copy)."""
        return deepcopy(self._species_dict[specie_name])

    # ========================================================================
    # Reaction Read Interface (NetworkABC Implementation)
    # ========================================================================

    def get_reaction_list(self) -> list[Reaction]:
        """Get all reactions as a list."""
        return list(self._reactions_dict.values())

    def get_reaction_dict(self) -> dict[int, Reaction]:
        """Get reactions dictionary (copy)."""
        return deepcopy(self._reactions_dict)

    def get_reaction(self, reaction_idx: int) -> Reaction:
        """Get a reaction by index (copy)."""
        return deepcopy(self._reactions_dict[reaction_idx])

    # ========================================================================
    # Query Methods (NetworkABC Implementation)
    # ========================================================================

    def get_reactions_by_types(
        self, reaction_type: Union[str, list[str]]
    ) -> list[Reaction]:
        """Get all reactions of specific type(s).

        Args:
            reaction_type: Single type or list of types to filter by

        Returns:
            List of reactions matching the type(s)
        """
        if isinstance(reaction_type, str):
            reaction_type = [reaction_type]

        return [
            reaction
            for reaction in self._reactions_dict.values()
            if reaction.get_type() in reaction_type
        ]

    def find_similar_reactions(self, reaction: Reaction) -> dict[int, Reaction]:
        """Find reactions with same reactants and products.

        Args:
            reaction: Reaction to find similar reactions for

        Returns:
            Dictionary of {index: Reaction} for matching reactions
        """
        similar = {}

        target_reactants = set(reaction.get_reactants()) - {"NAN"}
        target_products = set(reaction.get_products()) - {"NAN"}

        for idx, other_reaction in self._reactions_dict.items():
            other_reactants = set(other_reaction.get_reactants()) - {"NAN"}
            other_products = set(other_reaction.get_products()) - {"NAN"}

            if (
                other_reactants == target_reactants
                and other_products == target_products
            ):
                similar[idx] = other_reaction

        return similar

    def get_reaction_index(self, reaction: Reaction) -> int:
        """Get the index of a reaction in the network.

        Args:
            reaction: Reaction to find

        Returns:
            Index of the reaction

        Raises:
            ValueError: If reaction not found or multiple matches exist
        """
        similar = self.find_similar_reactions(reaction)

        if len(similar) == 0:
            raise ValueError(f"Reaction {reaction} not found in network")
        elif len(similar) > 1:
            raise ValueError(
                f"Multiple reactions match {reaction}. "
                f"Found indices: {list(similar.keys())}"
            )

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
        >>> network = Network.build(species, reactions, gas_phase_extrapolation=True)

    Attributes:
        _species_dict: Internal species storage {name: Species}
        _reactions_dict: Internal reaction storage {index: Reaction}
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
            species_dict: Species dictionary {name: Species}
            reaction_dict: Reaction dictionary {index: Reaction}
        """
        self._species_dict = species_dict
        self._reactions_dict = reaction_dict

    # ========================================================================
    # Factory Methods (Class Methods)
    # ========================================================================

    @classmethod
    def from_csv(
        cls,
        species_path: Union[str, Path, None] = None,
        reactions_path: Union[str, Path, None] = None,
    ) -> "Network":
        """Load network from CSV files.

        Loads a pre-compiled network from CSV files without any validation
        or automatic generation. This is the primary method for loading
        networks for analysis purposes.

        Args:
            species_path: Path to species CSV (None = use default installation)
            reactions_path: Path to reactions CSV (None = use default installation)

        Returns:
            Network: Loaded network instance

        Examples:
            >>> # Load default compiled network
            >>> network = Network.from_csv()

            >>> # Load old/custom network for analysis
            >>> network = Network.from_csv('old/species.csv', 'old/reactions.csv')
        """
        # Use defaults if not provided
        if species_path is None:
            species_path = path.join(_ROOT, "species.csv")
        if reactions_path is None:
            reactions_path = path.join(_ROOT, "reactions.csv")

        logging.debug(f"Loading network from {species_path} and {reactions_path}")

        # Load CSVs
        species_data = pd.read_csv(species_path)
        reactions_data = pd.read_csv(reactions_path)

        # Parse into objects
        species_list = [Species(list(spec)) for idx, spec in species_data.iterrows()]
        reactions_list = [
            Reaction(list(reac)) for idx, reac in reactions_data.iterrows()
        ]

        # Create dictionaries
        species_dict = {s.get_name(): s for s in species_list}
        reaction_dict = {k: v for k, v in enumerate(reactions_list)}

        return cls(species_dict, reaction_dict)

    @classmethod
    def from_lists(
        cls,
        species: list[Species],
        reactions: list[Reaction],
    ) -> "Network":
        """Create network directly from lists.

        Direct instantiation from species and reaction lists without any
        validation or automatic generation. Useful for programmatic network
        construction or as a base for NetworkBuilder.

        Args:
            species: List of Species objects
            reactions: List of Reaction objects

        Returns:
            Network: Network instance

        Example:
            >>> network = Network.from_lists(species_list, reactions_list)
        """
        species_dict = {s.get_name(): s for s in species}
        reaction_dict = {k: v for k, v in enumerate(reactions)}
        return cls(species_dict, reaction_dict)

    @classmethod
    def build(
        cls, species: list[Species], reactions: list[Reaction], **build_options
    ) -> "Network":
        """Build network with full validation and automatic generation.

        This is the primary method for building new networks with full validation,
        automatic reaction generation (freeze-out, desorption, bulk), branching
        ratio checks, and all build-time operations. Delegates to NetworkBuilder.

        Args:
            species: List of Species objects
            reactions: List of Reaction objects
            **build_options: Options passed to NetworkBuilder:
                - user_defined_bulk: List of user-defined bulk species
                - gas_phase_extrapolation: bool (default False)
                - add_crp_photo_to_grain: bool (default False)
                - derive_reaction_exothermicity: list[str] or None
                - database_reaction_exothermicity: list[Union[str, Path]] or None

        Returns:
            Network: Fully built and validated network

        Examples:
            >>> network = Network.build(
            ...     species=species_list,
            ...     reactions=reactions_list,
            ...     gas_phase_extrapolation=True,
            ...     add_crp_photo_to_grain=True
            ... )
        """
        from .network_builder import NetworkBuilder

        builder = NetworkBuilder(species, reactions, **build_options)
        return builder.build()

    # ========================================================================
    # Properties
    # ========================================================================

    @property
    def species(self):
        """Get species dictionary."""
        return self._species_dict

    # Note: Read operations (get_species_list, get_reaction_list, etc.)
    # are inherited from BaseNetwork

    # ========================================================================
    # Species Mutation Interface (MutableNetworkABC Implementation)
    # ========================================================================

    def set_specie(self, species_name: str, species: Species) -> None:
        """Set/update a species."""
        self._species_dict[species_name] = species

    def set_species_dict(self, new_species_dict: dict[str, Species]) -> None:
        """Replace entire species dictionary."""
        self._species_dict = new_species_dict

    def add_species(
        self, species: Union[Union[Species, list], list[Union[Species, list]]]
    ) -> None:
        """Add species to network.

        Args:
            species: Species object, list of Species, or CSV-style entries
        """
        # Convert to list of Species objects
        if isinstance(species, list):
            if len(species) == 0:
                logging.warning("Tried to add empty species list, ignoring.")
                return
            elif isinstance(species[0], Species):
                pass  # Already Species objects
            elif isinstance(species[0], list):
                try:
                    species = [Species(spec) for spec in species]
                except ValueError as error:
                    raise ValueError(
                        "Failed to convert CSV entries to Species objects"
                    ) from error
        elif isinstance(species, Species):
            species = [species]
        else:
            raise ValueError(
                "Input must be Species object, list of Species, or CSV entries"
            )

        # Add to dictionary
        for specie in species:
            # Filter out reaction types
            if specie.get_name() in reaction_types:
                logging.info(
                    f"Ignoring reaction type {specie.get_name()} in species list"
                )
                continue

            # Warn on duplicates
            if specie.get_name() in self._species_dict:
                logging.warning(
                    f"Species {specie.get_name()} already exists, keeping old definition"
                )
                continue

            # Filter out empty species
            if specie.get_name() in ["", "NAN"]:
                continue

            self._species_dict[specie.get_name()] = specie

    def remove_species(self, specie_name: str) -> None:
        """Remove a species from network."""
        if specie_name in self._species_dict:
            del self._species_dict[specie_name]
        else:
            logging.warning(f"Species {specie_name} not found in network")

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
        """Set/update a reaction at specific index."""
        old_length = len(self._reactions_dict)
        self._reactions_dict[reaction_idx] = reaction
        assert old_length == len(
            self._reactions_dict
        ), "Setting the reaction caused a change in the number of reactions"

    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        """Replace entire reaction dictionary."""
        self._reactions_dict = new_dict

    def add_reactions(
        self, reactions: Union[Union[Reaction, list], list[Union[Reaction, list]]]
    ) -> None:
        """Add reactions to network.

        Args:
            reactions: Reaction object, list of Reactions, or CSV-style entries
        """
        # Convert to list of Reaction objects
        if isinstance(reactions, list):
            if len(reactions) == 0:
                logging.warning("Tried to add empty reactions list, ignoring.")
                return
            elif isinstance(reactions[0], Reaction):
                pass  # Already Reaction objects
            elif isinstance(reactions[0], list):
                try:
                    reactions = [Reaction(reac) for reac in reactions]
                except ValueError as error:
                    raise ValueError(
                        "Failed to convert CSV entries to Reaction objects"
                    ) from error
        elif isinstance(reactions, Reaction):
            reactions = [reactions]
        else:
            raise ValueError(
                "Input must be Reaction object, list of Reactions, or CSV entries"
            )

        # Add to dictionary
        for reaction in reactions:
            if len(self._reactions_dict) == 0:
                new_idx = 0
            else:
                new_idx = max(self._reactions_dict.keys()) + 1

            self._reactions_dict[new_idx] = reaction

    def remove_reaction(self, reaction: Reaction) -> None:
        """Remove a reaction from network."""
        similar_reactions = list(self.find_similar_reactions(reaction).items())

        if len(similar_reactions) == 1:
            reaction_idx, _ = similar_reactions[0]
            del self._reactions_dict[reaction_idx]
        elif len(similar_reactions) == 0:
            logging.warning(f"Reaction {reaction} not found in network")
        else:
            raise RuntimeError(
                f"Found {len(similar_reactions)} reactions matching {reaction}. "
                "Use remove_reaction_by_index for piecewise reactions."
            )

    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        """Remove a reaction by its index."""
        if reaction_idx in self._reactions_dict:
            del self._reactions_dict[reaction_idx]
        else:
            logging.warning(f"Reaction index {reaction_idx} not found in network")

    def sort_reactions(self) -> None:
        """Sort reactions by type and first reactant."""
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

    # Note: Query methods (find_similar_reactions, get_reaction_index, etc.)
    # are inherited from BaseNetwork

    # ========================================================================
    # Parameter Modification Methods (NetworkABC Implementation)
    # ========================================================================

    def change_binding_energy(self, specie: str, new_binding_energy: float) -> None:
        """Change binding energy of a species.

        Handles special case of @H2O which affects other bulk species.
        """
        all_species = self.get_species_list()
        all_species_names = [s.get_name() for s in all_species]

        if specie not in all_species_names:
            raise ValueError(f"Species {specie} not found in network")

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
                    logging.warning(
                        f"Changing binding energy of bulk species {specie} "
                        "that was previously @H2O binding energy limited"
                    )

            self._species_dict[specie].set_binding_energy(new_binding_energy)

    def change_reaction_barrier(self, reaction: Reaction, barrier: float) -> None:
        """Change activation barrier of a reaction.

        If Fortran interface is available, also updates Fortran.
        """
        similar_reactions = list(self.find_similar_reactions(reaction).items())

        if len(similar_reactions) == 1:
            reaction_idx, _ = similar_reactions[0]
            self._reactions_dict[reaction_idx].set_gamma(barrier)

        elif len(similar_reactions) == 0:
            logging.warning(f"Reaction {reaction} not found in network")
        else:
            raise RuntimeError(
                f"Found {len(similar_reactions)} reactions matching {reaction}. "
                "Cannot uniquely identify which barrier to change."
            )


# ============================================================================
# Factory Functions (Module-Level)
# ============================================================================


def load_network_from_csv(
    species_path: Union[str, Path, None] = None,
    reactions_path: Union[str, Path, None] = None,
) -> Network:
    """Load a network from CSV files for analysis.

    This is a module-level factory function that provides clear documentation
    and intent. It calls Network.from_csv() internally.

    Use this when analyzing pre-compiled networks, comparing network versions,
    or loading old networks for analysis.

    Args:
        species_path: Path to species CSV (None = use default installation)
        reactions_path: Path to reactions CSV (None = use default installation)

    Returns:
        Network: Loaded network instance

    Examples:
        >>> # Load default compiled network
        >>> network = load_network_from_csv()

        >>> # Load old version for comparison
        >>> old_network = load_network_from_csv(
        ...     'archive/v3.0/species.csv',
        ...     'archive/v3.0/reactions.csv'
        ... )
        >>> print(f"Species added: {len(network.get_species_list()) - len(old_network.get_species_list())}")
    """
    return Network.from_csv(species_path, reactions_path)


def build_network(
    species: list[Species],
    reactions: list[Reaction],
    user_defined_bulk: list = None,
    gas_phase_extrapolation: bool = False,
    add_crp_photo_to_grain: bool = False,
    derive_reaction_exothermicity: list[str] = None,
    database_reaction_exothermicity: list[Union[str, Path]] = None,
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
        species: List of Species objects
        reactions: List of Reaction objects
        user_defined_bulk: User-specified bulk species (optional)
        gas_phase_extrapolation: Extrapolate gas-phase reactions to grain
        add_crp_photo_to_grain: Add CRP/PHOTON reactions to grain surface
        derive_reaction_exothermicity: Reaction types to calculate exothermicity for
        database_reaction_exothermicity: Custom exothermicity database files

    Returns:
        Network: Fully built and validated network

    Examples:
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
        ... )
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
        species: List of Species objects
        reactions: List of Reaction objects

    Returns:
        Network: Network instance

    Example:
        >>> network = create_network(species_list, reactions_list)
        >>> network.add_reactions(additional_reactions)
    """
    return Network.from_lists(species, reactions)
