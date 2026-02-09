"""
Runtime network interface for UCLCHEM's compiled chemical network.

This module provides RuntimeNetwork for runtime parameter modification during
model execution. Unlike the Network class in makerates (which supports full CRUD),
RuntimeNetwork provides read access and parameter modification only, as the
compiled Fortran network has fixed structure.

Key capabilities:
- Read all species and reactions
- Modify reaction parameters (alpha, beta, gamma)
- Modify species binding energies
- Disable reactions (set alpha=0)
- Reset to initial state
- Cannot add new species or reactions (Fortran arrays are fixed size)
- Cannot truly remove species or reactions (only disable)

Thread Safety Warning:
    RuntimeNetwork modifies global Fortran module state and is NOT thread-safe.
    Do not use with multiprocessing, multithreading, or concurrent model runs.
"""

import warnings
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from uclchem.utils import UCLCHEM_ROOT_DIR

# Import the base network implementation from makerates
from ..makerates.network import BaseNetwork
from ..makerates.reaction import Reaction, skip_reaction_validation
from ..makerates.species import Species


class RuntimeNetwork(BaseNetwork):
    """Runtime interface to UCLCHEM's compiled Fortran network.

    Provides read access and parameter modification for the compiled chemical
    network during model execution. The network structure (species/reactions)
    is fixed, but parameters can be modified.

    To "remove" a reaction: set its alpha parameter to 0.0 using disable_reaction().
    To reset changes: call reset_to_initial_state().

    Examples:
        >>> # Load runtime network
        >>> network = RuntimeNetwork.from_fortran()
        >>> print(f"Species: {len(network.get_species_list())}")
        >>> print(f"Reactions: {len(network.get_reaction_list())}")

        >>> # Modify parameters
        >>> network.modify_reaction_parameters(0, alpha=1e-10, beta=2.0)
        >>> network.change_binding_energy("H2O", 5773.0)

        >>> # Disable a reaction
        >>> network.disable_reaction(5)

        >>> # Reset when done
        >>> network.reset_to_initial_state()

    Thread Safety Warning:
        Modifies global Fortran module state. NOT thread-safe.
        Do not use with multiprocessing or multithreading.
    """

    def __init__(self):
        """Initialize RuntimeNetwork by loading the compiled Fortran module.

        Automatically imports uclchemwrap.network and validates against
        species.csv and reactions.csv from the installation directory.

        Raises:
            ImportError: If uclchemwrap.network cannot be imported
            ValueError: If CSV data doesn't match Fortran network dimensions
        """
        # Import the compiled Fortran network module
        try:
            from uclchemwrap import network as network_module
        except ImportError:
            raise ImportError(
                "Cannot import Fortran network module. "
                "Ensure UCLCHEM is properly compiled and installed."
            )

        self._fortran = network_module

        # Load CSV files for validation and better indexing
        self._load_csv_files()

        # Load species and reactions from Fortran
        self._species_dict = self._load_species_from_fortran()
        self._reactions_dict = self._load_reactions_from_fortran()

        # Validate that CSV matches Fortran dimensions
        self._validate_dimensions()

        # Cache initial state for reset functionality
        self._cache_initial_state()

    def _load_csv_files(self):
        """Load species and reactions CSV files from installation directory.

        These provide better indexing and validation against the Fortran network.
        """
        species_path = UCLCHEM_ROOT_DIR / "species.csv"
        reactions_path = UCLCHEM_ROOT_DIR / "reactions.csv"

        if not species_path.is_file():
            raise FileNotFoundError(f"Species CSV not found: {species_path}")
        if not reactions_path.is_file():
            raise FileNotFoundError(f"Reactions CSV not found: {reactions_path}")

        self._species_csv = pd.read_csv(species_path)
        self._reactions_csv = pd.read_csv(reactions_path)

    def _validate_dimensions(self):
        """Validate that CSV data matches Fortran network dimensions.

        Raises:
            ValueError: If dimensions don't match
        """
        n_species_csv = len(self._species_csv)
        n_species_fortran = len(self._fortran.specname)

        if n_species_csv != n_species_fortran:
            raise RuntimeError(
                f"Species count mismatch: {n_species_csv} in CSV vs "
                f"{n_species_fortran} in compiled Fortran network. "
                "The installation may be corrupted or out of sync."
            )

        n_reactions_csv = len(self._reactions_csv)
        n_reactions_fortran = len(self._fortran.alpha)

        if n_reactions_csv != n_reactions_fortran:
            raise RuntimeError(
                f"Reaction count mismatch: {n_reactions_csv} in CSV vs "
                f"{n_reactions_fortran} in compiled Fortran network. "
                "The installation may be corrupted or out of sync."
            )

        # Additional validation: check species names match
        for i in range(min(10, n_species_csv)):  # Check first 10 for quick validation
            csv_name = self._species_csv.iloc[i]["NAME"]
            fortran_name = str(np.char.decode(self._fortran.specname[i])).strip()

            if csv_name != fortran_name:
                warnings.warn(
                    f"Species name mismatch at index {i}: '{csv_name}' in CSV vs "
                    f"'{fortran_name}' in Fortran. Network may be out of sync.",
                    RuntimeWarning,
                )
                break

    def _load_species_from_fortran(self) -> dict[str, Species]:
        """Load species from Fortran arrays into Python Species objects.

        Returns:
            Dictionary mapping species names to Species objects
        """
        species_dict = {}

        n_species = len(self._fortran.specname)

        for i in range(n_species):
            # Extract data from Fortran arrays
            name = str(np.char.decode(self._fortran.specname[i])).strip()
            mass = float(self._fortran.mass[i])

            # Optional fields with bounds checking
            binding_energy = (
                float(self._fortran.bindingenergy[i])
                if i < len(self._fortran.bindingenergy)
                else 0.0
            )

            # Formation enthalpy if available
            enthalpy = (
                float(self._fortran.formationenthalpy[i])
                if hasattr(self._fortran, "formationenthalpy")
                and i < len(self._fortran.formationenthalpy)
                else 0.0
            )

            # Create Species object (CSV-style row format)
            # [NAME, MASS, BINDING_ENERGY, SOLID_FRACTION, MONO_FRACTION, VOLCANO_FRACTION, ENTHALPY]
            species_row = [name, mass, binding_energy, 1.0, 1.0, 0.0, enthalpy]
            species = Species(species_row)

            species_dict[name] = species

        return species_dict

    def _load_reactions_from_fortran(self) -> dict[int, Reaction]:
        """Load reactions from Fortran arrays into Python Reaction objects.

        Returns:
            Dictionary mapping reaction indices to Reaction objects
        """
        reactions_dict = {}

        n_reactions = len(self._fortran.alpha)

        # Load reactions without validation - the compiled Fortran network
        # may contain modeling simplifications (e.g., pseudo-hydrogenation)
        # that intentionally don't conserve elements
        with skip_reaction_validation():
            for i in range(n_reactions):
                # Extract reactant indices (1-based in Fortran, 0 for NAN)
                re1 = int(self._fortran.re1[i])
                re2 = int(self._fortran.re2[i])
                re3 = int(self._fortran.re3[i]) if hasattr(self._fortran, "re3") else 0

                # Extract product indices
                p1 = int(self._fortran.p1[i])
                p2 = int(self._fortran.p2[i])
                p3 = int(self._fortran.p3[i])
                p4 = int(self._fortran.p4[i])

                # Convert indices to species names
                reactant1 = self._get_species_name(re1) if re1 > 0 else "NAN"
                reactant2 = self._get_species_name(re2) if re2 > 0 else "NAN"
                reactant3 = self._get_species_name(re3) if re3 > 0 else "NAN"

                product1 = self._get_species_name(p1) if p1 > 0 else "NAN"
                product2 = self._get_species_name(p2) if p2 > 0 else "NAN"
                product3 = self._get_species_name(p3) if p3 > 0 else "NAN"
                product4 = self._get_species_name(p4) if p4 > 0 else "NAN"

                # Extract rate parameters
                alpha = float(self._fortran.alpha[i])
                beta = float(self._fortran.beta[i])
                gamma = float(self._fortran.gama[i])

                # Temperature ranges
                temp_low = (
                    float(self._fortran.mintemps[i])
                    if hasattr(self._fortran, "mintemps")
                    else 0.0
                )
                temp_high = (
                    float(self._fortran.maxtemps[i])
                    if hasattr(self._fortran, "maxtemps")
                    else 1e6
                )

                # Optional fields
                reduced_mass = (
                    float(self._fortran.reducedmasses[i])
                    if hasattr(self._fortran, "reducedmasses")
                    and i < len(self._fortran.reducedmasses)
                    else 0.0
                )
                exothermicity = (
                    float(self._fortran.exothermicities[i])
                    if hasattr(self._fortran, "exothermicities")
                    and i < len(self._fortran.exothermicities)
                    else 0.0
                )

                # Create Reaction object (CSV-style row format)
                # [R1, R2, R3, P1, P2, P3, P4, alpha, beta, gamma, Tmin, Tmax, reduced_mass, extrapolate, exothermicity]
                reaction_row = [
                    reactant1,
                    reactant2,
                    reactant3,
                    product1,
                    product2,
                    product3,
                    product4,
                    alpha,
                    beta,
                    gamma,
                    temp_low,
                    temp_high,
                    reduced_mass,
                    0,  # extrapolate flag
                    exothermicity,
                ]

                reaction = Reaction(reaction_row)
                reactions_dict[i] = reaction

        return reactions_dict

    def _get_species_name(self, index: int) -> str:
        """Get species name from 1-based Fortran index.

        Args:
            index: 1-based species index in Fortran

        Returns:
            Species name
        """
        if index <= 0 or index > len(self._fortran.specname):
            return "NAN"

        array_idx = index - 1  # Convert to 0-based
        return str(np.char.decode(self._fortran.specname[array_idx])).strip()

    def _cache_initial_state(self):
        """Cache the initial state of all modifiable Fortran parameters.

        Allows fast reset without re-reading CSV files. Caches:
        - Reaction parameters: alpha, beta, gama
        - Species parameters: bindingenergy
        """
        self._initial_alpha = np.copy(self._fortran.alpha)
        self._initial_beta = np.copy(self._fortran.beta)
        self._initial_gama = np.copy(self._fortran.gama)
        self._initial_bindingenergy = np.copy(self._fortran.bindingenergy)

    # ========================================================================
    # Properties (NetworkABC Implementation)
    # ========================================================================

    @property
    def species(self):
        """Get species dictionary."""
        return self._species_dict

    # Note: Read operations (get_species_list, get_species_dict, get_specie,
    # get_reaction_list, get_reaction_dict, get_reaction) are inherited from BaseNetwork

    # ========================================================================
    # Species Interface - Unsupported Operations
    # ========================================================================

    def add_species(self, species: Union[Species, List[Species]]) -> None:
        """NOT SUPPORTED: Cannot add species to compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError(
            "Cannot add species to RuntimeNetwork. "
            "The compiled Fortran network has fixed structure. "
            "Use Network class (from makerates) for building new networks."
        )

    def remove_species(self, specie_name: str) -> None:
        """NOT SUPPORTED: Cannot remove species from compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError(
            "Cannot remove species from RuntimeNetwork. "
            "The compiled Fortran network has fixed structure."
        )

    def set_specie(self, species_name: str, species: Species) -> None:
        """NOT SUPPORTED: Cannot replace species in compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError(
            "Cannot set species in RuntimeNetwork. "
            "Use change_binding_energy() to modify species parameters."
        )

    def set_species_dict(self, new_species_dict: dict[str, Species]) -> None:
        """NOT SUPPORTED: Cannot replace species dictionary.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError("Cannot replace species dictionary in RuntimeNetwork.")

    def sort_species(self) -> None:
        """NOT SUPPORTED: Species order is fixed in compiled network.

        Raises:
            NotImplementedError: Always - species order is fixed
        """
        raise NotImplementedError(
            "Cannot sort species in RuntimeNetwork. Species order is fixed."
        )

    # ========================================================================
    # Reaction Interface - Unsupported Operations
    # ========================================================================

    def add_reactions(self, reactions: Union[Reaction, List[Reaction]]) -> None:
        """NOT SUPPORTED: Cannot add reactions to compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError(
            "Cannot add reactions to RuntimeNetwork. "
            "The compiled Fortran network has fixed structure. "
            "Use Network class (from makerates) for building new networks."
        )

    def remove_reaction(self, reaction: Reaction) -> None:
        """NOT SUPPORTED: Cannot remove reactions from compiled Fortran network.

        Use disable_reaction() to set alpha=0 instead.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError(
            "Cannot remove reactions from RuntimeNetwork. "
            "Use disable_reaction() to set alpha=0 to effectively disable a reaction."
        )

    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        """NOT SUPPORTED: Cannot remove reactions from compiled Fortran network.

        Use disable_reaction() to set alpha=0 instead.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError(
            "Cannot remove reactions from RuntimeNetwork. "
            "Use disable_reaction() to set alpha=0 to effectively disable a reaction."
        )

    def set_reaction(self, reaction_idx: int, reaction: Reaction) -> None:
        """NOT SUPPORTED: Cannot replace reactions in compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError(
            "Cannot replace reactions in RuntimeNetwork. "
            "Use modify_reaction_parameters() to modify reaction parameters."
        )

    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        """NOT SUPPORTED: Cannot replace reaction dictionary.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size
        """
        raise NotImplementedError("Cannot replace reaction dictionary in RuntimeNetwork.")

    def sort_reactions(self) -> None:
        """NOT SUPPORTED: Reaction order is fixed in compiled network.

        Raises:
            NotImplementedError: Always - reaction order is fixed
        """
        raise NotImplementedError(
            "Cannot sort reactions in RuntimeNetwork. Reaction order is fixed."
        )

    # Note: Query methods (get_reactions_by_types, find_similar_reactions,
    # get_reaction_index) are inherited from BaseNetwork

    # ========================================================================
    # Parameter Modification Methods (NetworkABC Implementation)
    # ========================================================================

    def change_binding_energy(self, specie: str, new_binding_energy: float) -> None:
        """Change binding energy of a species (modifies Fortran array).

        Args:
            specie: Name of the species
            new_binding_energy: New binding energy in Kelvin

        Raises:
            KeyError: If species not found
        """
        # Find species index (1-based)
        species_names = [
            self._get_species_name(i + 1) for i in range(len(self._fortran.specname))
        ]

        try:
            species_idx = species_names.index(specie)
        except ValueError:
            raise KeyError(f"Species '{specie}' not found in network")

        # Modify Fortran array (0-based)
        self._fortran.bindingenergy[species_idx] = float(new_binding_energy)

        # Update cached species object
        if specie in self._species_dict:
            self._species_dict[specie].set_binding_energy(new_binding_energy)

    def change_reaction_barrier(self, reaction: Reaction, barrier: float) -> None:
        """Change activation barrier of a reaction (modifies Fortran gamma).

        Args:
            reaction: Reaction to modify
            barrier: New activation barrier in Kelvin

        Raises:
            ValueError: If reaction not found or multiple matches
        """
        assert reaction.is_ice_reaction(), "Only ice reactions have modifiable barriers."
        reaction_idx = self.get_reaction_index(reaction)
        self.modify_reaction_parameters(reaction_idx, gamma=barrier)

    # ========================================================================
    # RuntimeNetwork-Specific Methods
    # ========================================================================

    def modify_reaction_parameters(
        self,
        reaction_idx: int,
        alpha: Optional[float] = None,
        beta: Optional[float] = None,
        gamma: Optional[float] = None,
    ) -> None:
        """Modify reaction rate parameters in Fortran arrays.

        Args:
            reaction_idx: Index of reaction to modify (0-based)
            alpha: New alpha value (pre-exponential factor)
            beta: New beta value (temperature exponent)
            gamma: New gamma value (activation energy in K)

        Raises:
            IndexError: If reaction_idx out of range

        Example:
            >>> network.modify_reaction_parameters(0, alpha=1e-10, beta=2.0)
        """
        if reaction_idx < 0 or reaction_idx >= len(self._fortran.alpha):
            raise IndexError(f"Reaction index {reaction_idx} out of range")

        if alpha is not None:
            self._fortran.alpha[reaction_idx] = float(alpha)
            if reaction_idx in self._reactions_dict:
                self._reactions_dict[reaction_idx].set_alpha(float(alpha))

        if beta is not None:
            self._fortran.beta[reaction_idx] = float(beta)
            if reaction_idx in self._reactions_dict:
                self._reactions_dict[reaction_idx].set_beta(float(beta))

        if gamma is not None:
            self._fortran.gama[reaction_idx] = float(gamma)
            if reaction_idx in self._reactions_dict:
                self._reactions_dict[reaction_idx].set_gamma(float(gamma))

    def disable_reaction(self, reaction_idx: int) -> None:
        """Disable a reaction by setting alpha=0.

        This is the only way to "remove" a reaction at runtime since
        the Fortran network has fixed structure. Setting alpha=0 makes
        the reaction have zero rate.

        Args:
            reaction_idx: Index of reaction to disable (0-based)

        Example:
            >>> network.disable_reaction(5)
        """
        self.modify_reaction_parameters(reaction_idx, alpha=0.0)

    def reset_to_initial_state(self) -> None:
        """Reset all Fortran parameters to their initial cached values.

        Restores:
        - All reaction parameters (alpha, beta, gamma)
        - All species binding energies

        Example:
            >>> network.modify_reaction_parameters(0, alpha=999.0)
            >>> network.reset_to_initial_state()  # Restores original alpha
        """
        np.copyto(self._fortran.alpha, self._initial_alpha)
        np.copyto(self._fortran.beta, self._initial_beta)
        np.copyto(self._fortran.gama, self._initial_gama)
        np.copyto(self._fortran.bindingenergy, self._initial_bindingenergy)

        # Reload species and reactions to sync cached objects
        self._species_dict = self._load_species_from_fortran()
        self._reactions_dict = self._load_reactions_from_fortran()

    @property
    def fortran_module(self):
        """Direct access to Fortran module for advanced users.

        Warning: Use with caution. Direct modification bypasses safety checks.
        """
        warnings.warn(
            "Direct access to Fortran module is discouraged, this can break ungracefully. "
            "Use GeneralSettings and instead"
        )
        return self._fortran
