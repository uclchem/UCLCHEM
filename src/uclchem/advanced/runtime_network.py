"""Runtime network interface for UCLCHEM's compiled chemical network.

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

import logging
import warnings
from types import ModuleType

import numpy as np
import pandas as pd

# Import the base network implementation from makerates
from uclchem.makerates.network import BaseNetwork
from uclchem.makerates.reaction import Reaction, skip_reaction_validation
from uclchem.makerates.species import Species
from uclchem.utils import UCLCHEM_ROOT_DIR

logger = logging.getLogger(__name__)


class RuntimeNetwork(BaseNetwork):
    """Runtime interface to UCLCHEM's compiled Fortran network.

    Provides read access and parameter modification for the compiled chemical
    network during model execution. The network structure (species/reactions)
    is fixed, but parameters can be modified.

    To "remove" a reaction: set its alpha parameter to 0.0 using disable_reaction().
    To reset changes: call reset_to_initial_state().

    Examples:
        >>> # Load runtime network
        >>> network = RuntimeNetwork()
        >>> print(f"Species: {len(network.get_species_list())}")
        Species: ...
        >>> print(f"Reactions: {len(network.get_reaction_list())}")
        Reactions: ...

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

    # Fortran uses 9999 as a sentinel value for reaction type keywords
    # (GAR, PHOTON, CRP, CRPHOT, etc.) instead of species indices
    _FORTRAN_KEYWORD_SENTINEL = 9999

    def __init__(self):
        """Initialize RuntimeNetwork by loading the compiled Fortran module.

        Automatically imports uclchemwrap.network and validates against
        species.csv and reactions.csv from the installation directory.

        Raises:
            ImportError: If uclchemwrap.network cannot be imported

        """
        # Import the compiled Fortran network module
        try:
            from uclchemwrap import network as network_module
        except ImportError as e:
            msg = (
                "Cannot import Fortran network module. "
                "Ensure UCLCHEM is properly compiled and installed."
            )
            raise ImportError(msg) from e

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

        Raises:
            FileNotFoundError: If `"UCLCHEM_ROOT_DIR/species.csv"` or
                `"UCLCHEM_ROOT_DIR/reactions.csv"` are not valid files.
        """
        species_path = UCLCHEM_ROOT_DIR / "species.csv"
        reactions_path = UCLCHEM_ROOT_DIR / "reactions.csv"

        if not species_path.is_file():
            msg = f"Species CSV not found: {species_path}"
            raise FileNotFoundError(msg)
        if not reactions_path.is_file():
            msg = f"Reactions CSV not found: {reactions_path}"
            raise FileNotFoundError(msg)

        self._species_csv = pd.read_csv(species_path)
        self._reactions_csv = pd.read_csv(reactions_path)

    def _validate_dimensions(self):
        """Validate that CSV data matches Fortran network dimensions.

        Raises:
            RuntimeError: If dimensions don't match

        """
        logger.debug("Validating dimensions of species in csv and in Fortran")

        n_species_csv = len(self._species_csv)
        n_species_fortran = len(self._fortran.specname)

        if n_species_csv != n_species_fortran:
            msg = (
                f"Species count mismatch: {n_species_csv} in CSV vs "
                f"{n_species_fortran} in compiled Fortran network. "
                "The installation may be corrupted or out of sync."
            )
            raise RuntimeError(msg)
        logger.debug("Ok!")

        logger.debug("Validating dimensions of reactions in csv and in Fortran")
        n_reactions_csv = len(self._reactions_csv)
        n_reactions_fortran = len(self._fortran.alpha)

        if n_reactions_csv != n_reactions_fortran:
            msg = (
                f"Reaction count mismatch: {n_reactions_csv} in CSV vs "
                f"{n_reactions_fortran} in compiled Fortran network. "
                "The installation may be corrupted or out of sync."
            )
            raise RuntimeError(msg)
        logger.debug("Ok!")

        logger.debug(
            "Validating that the species names of the first couple of species match"
        )
        for i in range(min(10, n_species_csv)):  # Check first 10 for quick validation
            csv_name = self._species_csv.iloc[i]["NAME"]
            fortran_name = str(np.char.decode(self._fortran.specname[i])).strip()

            if csv_name != fortran_name:
                warnings.warn(
                    f"Species name mismatch at index {i}: '{csv_name}' in CSV vs "
                    f"'{fortran_name}' in Fortran. Network may be out of sync.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                break
        logger.debug("Ok!")

    def _load_species_from_fortran(self) -> dict[str, Species]:
        """Load species from Fortran arrays into Python Species objects.

        Returns:
            Dictionary mapping species names to Species objects

        """
        species_dict = {}

        n_species = len(self._fortran.specname)
        logger.debug(f"Loading {n_species} species from Fortran")

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
            # [
            #    NAME, MASS, BINDING_ENERGY, SOLID_FRACTION,
            #    MONO_FRACTION, VOLCANO_FRACTION, ENTHALPY,
            # ]
            species_row = [name, mass, binding_energy, 1.0, 1.0, 0.0, enthalpy]
            species = Species(species_row)

            species_dict[name] = species

        return species_dict

    def _load_reactions_from_fortran(self) -> dict[int, Reaction]:
        """Load reactions from Fortran arrays into Python Reaction objects.

        Fortran uses special sentinel values for reactant/product indices:
        - Valid species: 1 to n_species (1-based indexing)
        - Empty slot: 0
        - Reaction keywords (GAR, PHOTON, CRP, etc.): 9999

        Returns:
            Dictionary mapping reaction indices to Reaction objects

        """
        reactions_dict = {}

        n_reactions = len(self._fortran.alpha)
        n_species = len(self._fortran.specname)
        logger.debug(f"Loading {n_reactions} reactions from Fortran")

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
                # Valid species indices: 1 to n_species (1-based)
                # Keyword sentinel (9999): fall back to CSV for GAR, PHOTON, etc.
                reactant1 = self._get_reactant_name(i, re1, n_species, "Reactant 1")
                reactant2 = self._get_reactant_name(i, re2, n_species, "Reactant 2")
                reactant3 = self._get_reactant_name(i, re3, n_species, "Reactant 3")

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
                # ruff: noqa: ERA001
                # [
                #    R1, R2, R3, P1, P2, P3, P4, alpha, beta, gamma,
                #    Tmin, Tmax, reduced_mass, extrapolate, exothermicity,
                # ]
                # ruff: noqa: ERA001
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

    def _get_reactant_name(
        self, reaction_idx: int, fortran_idx: int, n_species: int, csv_column: str
    ) -> str:
        """Get reactant name from Fortran index, handling keyword sentinels.

        Fortran uses different values for reactant indices:
        - 1 to n_species: Valid species index (1-based)
        - 0: Empty reactant slot (returns "NAN")
        - 9999: Reaction type keyword (GAR, PHOTON, CRP, etc.)

        For the keyword sentinel (9999), we fall back to the CSV file to
        get the actual keyword string.

        Args:
            reaction_idx: Index of reaction in reactions list
            fortran_idx: Fortran array value for this reactant
            n_species: Total number of species in network
            csv_column: Column name in reactions CSV ("Reactant 1", etc.)

        Returns:
            Species name or reaction type keyword

        """
        # Check if this is a valid species index (1-based)
        if 0 < fortran_idx <= n_species:
            return self._get_species_name(fortran_idx)

        # For sentinel values (0 or 9999), fall back to CSV
        # Use pd.isna() to properly detect NaN values from pandas
        csv_val = self._reactions_csv.iloc[reaction_idx][csv_column]
        if pd.isna(csv_val) or not isinstance(csv_val, str):
            return "NAN"

        # Return the keyword string (GAR, PHOTON, CRP, etc.)
        return str(csv_val).strip().upper()

    def _cache_initial_state(self):
        """Cache the initial state of all modifiable Fortran parameters.

        Allows fast reset without re-reading CSV files. Caches:
        - Reaction parameters: alpha, beta, gama
        - Species parameters: bindingenergy
        """
        logger.debug("Caching initial state of RuntimeNetwork")
        self._initial_alpha = np.copy(self._fortran.alpha)
        self._initial_beta = np.copy(self._fortran.beta)
        self._initial_gama = np.copy(self._fortran.gama)
        self._initial_bindingenergy = np.copy(self._fortran.bindingenergy)

    # ========================================================================
    # Properties (NetworkABC Implementation)
    # ========================================================================

    @property
    def species(self) -> dict[str, Species]:
        """Get species dictionary."""
        return self._species_dict

    # Note: Read operations (get_species_list, get_species_dict, get_specie,
    # get_reaction_list, get_reaction_dict, get_reaction) are inherited from BaseNetwork

    # ========================================================================
    # Species Interface - Unsupported Operations
    # ========================================================================

    def add_species(self, species: Species | list[Species]) -> None:
        """NOT SUPPORTED: Cannot add species to compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = (
            "Cannot add species to RuntimeNetwork. "
            "The compiled Fortran network has fixed structure. "
            "Use Network class (from makerates) for building new networks."
        )
        raise NotImplementedError(msg)

    def remove_species(self, specie_name: str) -> None:
        """NOT SUPPORTED: Cannot remove species from compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = (
            "Cannot remove species from RuntimeNetwork. "
            "The compiled Fortran network has fixed structure."
        )
        raise NotImplementedError(msg)

    def set_specie(self, species_name: str, species: Species) -> None:
        """NOT SUPPORTED: Cannot replace species in compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = (
            "Cannot set species in RuntimeNetwork. "
            "Use change_binding_energy() to modify species parameters."
        )
        raise NotImplementedError(msg)

    def set_species_dict(self, new_species_dict: dict[str, Species]) -> None:
        """NOT SUPPORTED: Cannot replace species dictionary.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = "Cannot replace species dictionary in RuntimeNetwork."
        raise NotImplementedError(msg)

    def sort_species(self) -> None:
        """NOT SUPPORTED: Species order is fixed in compiled network.

        Raises:
            NotImplementedError: Always - species order is fixed

        """
        msg = "Cannot sort species in RuntimeNetwork. Species order is fixed."
        raise NotImplementedError(msg)

    # ========================================================================
    # Reaction Interface - Unsupported Operations
    # ========================================================================

    def add_reactions(self, reactions: Reaction | list[Reaction]) -> None:
        """NOT SUPPORTED: Cannot add reactions to compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = (
            "Cannot add reactions to RuntimeNetwork. "
            "The compiled Fortran network has fixed structure. "
            "Use Network class (from makerates) for building new networks."
        )
        raise NotImplementedError(msg)

    def remove_reaction(self, reaction: Reaction) -> None:
        """NOT SUPPORTED: Cannot remove reactions from compiled Fortran network.

        Use disable_reaction() to set alpha=0 instead.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = (
            "Cannot remove reactions from RuntimeNetwork. "
            "Use disable_reaction() to set alpha=0 to effectively disable a reaction."
        )
        raise NotImplementedError(msg)

    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        """NOT SUPPORTED: Cannot remove reactions from compiled Fortran network.

        Use disable_reaction() to set alpha=0 instead.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = (
            "Cannot remove reactions from RuntimeNetwork. "
            "Use disable_reaction() to set alpha=0 to effectively disable a reaction."
        )
        raise NotImplementedError(msg)

    def set_reaction(self, reaction_idx: int, reaction: Reaction) -> None:
        """NOT SUPPORTED: Cannot replace reactions in compiled Fortran network.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = (
            "Cannot replace reactions in RuntimeNetwork. "
            "Use modify_reaction_parameters() to modify reaction parameters."
        )
        raise NotImplementedError(msg)

    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        """NOT SUPPORTED: Cannot replace reaction dictionary.

        Raises:
            NotImplementedError: Always - Fortran arrays have fixed size

        """
        msg = "Cannot replace reaction dictionary in RuntimeNetwork."
        raise NotImplementedError(msg)

    def sort_reactions(self) -> None:
        """NOT SUPPORTED: Reaction order is fixed in compiled network.

        Raises:
            NotImplementedError: Always - reaction order is fixed

        """
        msg = "Cannot sort reactions in RuntimeNetwork. Reaction order is fixed."
        raise NotImplementedError(msg)

    # Note: Query methods (get_reactions_by_types, find_similar_reactions,
    # get_reaction_index) are inherited from BaseNetwork

    # ========================================================================
    # Parameter Modification Methods (NetworkABC Implementation)
    # ========================================================================

    def change_binding_energy(self, specie: str, new_binding_energy: int | float) -> None:
        """Change binding energy of a species (modifies Fortran array).

        Args:
            specie (str): Name of the species
            new_binding_energy (int | float): New binding energy in Kelvin

        Raises:
            KeyError: If species not found

        """
        # Find species index (1-based)
        species_names = [
            self._get_species_name(i + 1) for i in range(len(self._fortran.specname))
        ]

        try:
            species_idx = species_names.index(specie)
        except ValueError as e:
            msg = f"Species '{specie}' not found in network"
            raise KeyError(msg) from e

        logger.debug(
            f"Changing the binding energy of species {specie} to {new_binding_energy}"
        )
        # Modify Fortran array (0-based)
        self._fortran.bindingenergy[species_idx] = float(new_binding_energy)

        # Update cached species object
        if specie in self._species_dict:
            self._species_dict[specie].set_binding_energy(new_binding_energy)

    def change_reaction_barrier(self, reaction: Reaction, barrier: float) -> None:
        """Change activation barrier of a reaction (modifies Fortran gamma).

        Args:
            reaction (Reaction): Reaction to modify
            barrier (float): New activation barrier in Kelvin

        Raises:
            RuntimeError: If reaction is not a reaction on the ices.

        """
        if not reaction.is_ice_reaction():
            msg = "Only ice reactions have modifiable barriers."
            raise RuntimeError(msg)
        reaction_idx = self.get_reaction_index(reaction)
        self.modify_reaction_parameters(reaction_idx, gamma=barrier)

    # ========================================================================
    # RuntimeNetwork-Specific Methods
    # ========================================================================

    def modify_reaction_parameters(
        self,
        reaction_idx: int,
        alpha: float | None = None,
        beta: float | None = None,
        gamma: float | None = None,
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
            >>> network = RuntimeNetwork()
            >>> network.modify_reaction_parameters(0, alpha=1e-10, beta=2.0)
            >>> network.reset_to_initial_state()

        """
        if reaction_idx < 0 or reaction_idx >= len(self._fortran.alpha):
            msg = f"Reaction index {reaction_idx} out of range"
            raise IndexError(msg)

        if alpha is not None:
            logger.debug(
                f"Setting alpha of reaction with index {reaction_idx} to {alpha}"
            )
            self._fortran.alpha[reaction_idx] = float(alpha)
            if reaction_idx in self._reactions_dict:
                self._reactions_dict[reaction_idx].set_alpha(float(alpha))

        if beta is not None:
            logger.debug(f"Setting beta of reaction with index {reaction_idx} to {alpha}")
            self._fortran.beta[reaction_idx] = float(beta)
            if reaction_idx in self._reactions_dict:
                self._reactions_dict[reaction_idx].set_beta(float(beta))

        if gamma is not None:
            logger.debug(
                f"Setting gamma of reaction with index {reaction_idx} to {alpha}"
            )
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
            >>> network = RuntimeNetwork()
            >>> network.disable_reaction(5)
            >>> network.reset_to_initial_state()

        """
        logger.debug(f"Disabling reaction with index {reaction_idx}")
        self.modify_reaction_parameters(reaction_idx, alpha=0.0)

    def reset_to_initial_state(self) -> None:
        """Reset all Fortran parameters to their initial cached values.

        Restores:
        - All reaction parameters (alpha, beta, gamma)
        - All species binding energies

        Example:
            >>> network = RuntimeNetwork()
            >>> network.modify_reaction_parameters(0, alpha=999.0)
            >>> network.reset_to_initial_state()  # Restores original alpha

        """
        logger.debug("Resetting RuntimeNetwork to initial state")
        np.copyto(self._fortran.alpha, self._initial_alpha)
        np.copyto(self._fortran.beta, self._initial_beta)
        np.copyto(self._fortran.gama, self._initial_gama)
        np.copyto(self._fortran.bindingenergy, self._initial_bindingenergy)

        # Reload species and reactions to sync cached objects
        self._species_dict = self._load_species_from_fortran()
        self._reactions_dict = self._load_reactions_from_fortran()

    @property
    def fortran_module(self) -> ModuleType:
        """Direct access to Fortran module for advanced users.

        Warning: Use with caution. Direct modification bypasses safety checks.
        """
        warnings.warn(
            "Direct access to Fortran module is discouraged, this can break ungracefully. "
            "Use GeneralSettings instead",
            stacklevel=2,
        )
        return self._fortran
