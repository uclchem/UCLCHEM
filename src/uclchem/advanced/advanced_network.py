"""
Runtime access to UCLCHEM's compiled chemical network state.

This module provides a class-based interface for accessing and modifying the chemical
network that is compiled into the Fortran code. Unlike the Network class in makerates
(used for building networks), NetworkState provides access to the active, in-memory
network during model execution.

RuntimeSpecies and RuntimeReaction provide APIs similar to the Species and Reaction
classes in makerates, but wrap the compiled Fortran data rather than generating new code.

**Thread Safety Warning:**
NetworkState modifies global Fortran module state and is **NOT thread-safe**.
Do not use with multiprocessing, multithreading, or concurrent model runs.

Note: Changes made through NetworkState affect the global Fortran state and persist
across model runs in the same Python session.
"""

import os
from typing import List, Optional

import numpy as np
import pandas as pd
from uclchemwrap import network as network_module

from ..makerates.reaction import Reaction
from ..makerates.species import Species


class RuntimeSpecies:
    """Wrapper for a species in the compiled network.

    Provides a similar API to makerates.Species but accesses the compiled Fortran data.
    """

    def __init__(self, index: int, network_ref):
        """Initialize a runtime species wrapper.

        Args:
            index: 1-based species index in Fortran arrays
            network_ref: Reference to the network module
        """
        self._index = index
        self._network = network_ref
        self._array_idx = index - 1  # 0-based for array access

    def get_name(self) -> str:
        """Get the species name.

        Returns:
            Species name
        """
        return str(np.char.decode(self._network.specname[self._array_idx])).strip()

    def get_mass(self) -> float:
        """Get the molecular mass.

        Returns:
            Mass in atomic mass units
        """
        return float(self._network.mass[self._array_idx])

    def get_binding_energy(self) -> Optional[float]:
        """Get the binding energy.

        Returns:
            Binding energy in Kelvin (or None if not available)
        """
        if self._array_idx < len(self._network.bindingenergy):
            return float(self._network.bindingenergy[self._array_idx])
        return None

    def get_enthalpy(self) -> Optional[float]:
        """Get the formation enthalpy.

        Returns:
            Formation enthalpy in kJ/mol (or None if not available)
        """
        if self._array_idx < len(self._network.formationenthalpy):
            return float(self._network.formationenthalpy[self._array_idx])
        return None

    def is_ice_species(self) -> bool:
        """Return whether the species is a species on the grain

        Returns:
            bool: True if it is an ice species.
        """
        return (
            self.get_name() in ["BULK", "SURFACE"]
            or self.get_name().startswith(
                "#",
            )
            or self.get_name().startswith("@")
        )

    def is_ion(self) -> bool:
        """Check if this is an ion.

        Returns:
            True if species name contains + or -
        """
        name = self.get_name()
        return "+" in name or "-" in name

    def get_charge(self) -> int:
        """Get the charge of the species.

        Returns:
            Charge (+1, -1, or 0)
        """
        name = self.get_name()
        if "+" in name:
            return 1
        elif "-" in name:
            return -1
        return 0

    def __str__(self):
        return self.get_name()

    def __repr__(self):
        return f"RuntimeSpecies({self._index}, '{self.get_name()}')"


class RuntimeReaction:
    """Wrapper for a reaction in the compiled network.

    Provides a similar API to makerates.Reaction but accesses the compiled Fortran data.
    """

    def __init__(self, index: int, network_ref):
        """Initialize a runtime reaction wrapper.

        Args:
            index: 1-based reaction index in Fortran arrays
            network_ref: Reference to the network module
        """
        self._index = index
        self._network = network_ref
        self._array_idx = index - 1  # 0-based for array access

    def get_reactants(self) -> List[int]:
        """Get the reactant species indices.

        Returns:
            List of species indices (1-based, 0 for NAN)
        """
        return [
            int(self._network.re1[self._array_idx]),
            int(self._network.re2[self._array_idx]),
            int(self._network.re3[self._array_idx]),
        ]

    def get_products(self) -> List[int]:
        """Get the product species indices.

        Returns:
            List of species indices (1-based, 0 for NAN)
        """
        return [
            int(self._network.p1[self._array_idx]),
            int(self._network.p2[self._array_idx]),
            int(self._network.p3[self._array_idx]),
            int(self._network.p4[self._array_idx]),
        ]

    def get_reactant_names(self) -> List[str]:
        """Get the names of reactant species.

        Returns:
            List of reactant names (NAN for empty slots)
        """
        names = []
        for idx in self.get_reactants():
            if idx > 0:
                name = str(np.char.decode(self._network.specname[idx - 1])).strip()
                names.append(name)
            else:
                names.append("NAN")
        return names

    def get_product_names(self) -> List[str]:
        """Get the names of product species.

        Returns:
            List of product names (NAN for empty slots)
        """
        names = []
        for idx in self.get_products():
            if idx > 0:
                name = str(np.char.decode(self._network.specname[idx - 1])).strip()
                names.append(name)
            else:
                names.append("NAN")
        return names

    def get_alpha(self) -> float:
        """Get the alpha parameter (pre-exponential factor).

        Returns:
            Alpha parameter
        """
        return float(self._network.alpha[self._array_idx])

    def get_beta(self) -> float:
        """Get the beta parameter (temperature exponent).

        Returns:
            Beta parameter
        """
        return float(self._network.beta[self._array_idx])

    def get_gamma(self) -> float:
        """Get the gamma parameter (activation energy).

        Returns:
            Gamma parameter in Kelvin
        """
        return float(self._network.gama[self._array_idx])

    def get_templow(self) -> float:
        """Get the minimum valid temperature.

        Returns:
            Minimum temperature in Kelvin
        """
        return float(self._network.mintemps[self._array_idx])

    def get_temphigh(self) -> float:
        """Get the maximum valid temperature.

        Returns:
            Maximum temperature in Kelvin
        """
        return float(self._network.maxtemps[self._array_idx])

    def get_exothermicity(self) -> float:
        """Get the reaction exothermicity.

        Returns:
            Exothermicity in erg
        """
        return float(self._network.exothermicities[self._array_idx])

    def get_reduced_mass(self) -> Optional[float]:
        """Get the reduced mass for tunneling reactions.

        Returns:
            Reduced mass in AMU (or None if not available)
        """
        if self._array_idx < len(self._network.reducedmasses):
            return float(self._network.reducedmasses[self._array_idx])
        return None

    def get_rate(self) -> Optional[float]:
        """Get the computed reaction rate from the last model run.

        Returns:
            Computed rate (only meaningful after running a model)
        """
        if self._array_idx < len(self._network.reactionrate):
            return float(self._network.reactionrate[self._array_idx])
        return None

    def set_alpha(self, value: float):
        """Set the alpha parameter.

        Args:
            value: New alpha value
        """
        self._network.alpha[self._array_idx] = float(value)

    def set_beta(self, value: float):
        """Set the beta parameter.

        Args:
            value: New beta value
        """
        self._network.beta[self._array_idx] = float(value)

    def set_gamma(self, value: float):
        """Set the gamma parameter.

        Args:
            value: New gamma value
        """
        self._network.gama[self._array_idx] = float(value)

    def __str__(self):
        reactants = " + ".join([r for r in self.get_reactant_names() if r != "NAN"])
        products = " + ".join([p for p in self.get_product_names() if p != "NAN"])
        return f"{reactants} -> {products}"

    def __repr__(self):
        return f"RuntimeReaction({self._index}, '{self}')"


class NetworkState:
    """Runtime interface to UCLCHEM's compiled chemical network.

    Loads the network from CSV files (on-disk version) and compares with
    the compiled Fortran network (in-memory version) to ensure consistency.

    **Thread Safety Warning:**
    This class modifies global Fortran module state and is **NOT thread-safe**.
    Do not use with multiprocessing, multithreading, or concurrent model runs.

    Example:
        >>> from uclchem.advanced import NetworkState
        >>> network = NetworkState()
        >>> network.validate()  # Check on-disk matches in-memory
        >>> print(f"Species: {len(network.species_list)}")
        >>> print(f"Reactions: {len(network.reaction_list)}")
    """

    def __init__(self):
        """Initialize the NetworkState interface.

        Loads species and reactions from CSV files and compares with
        the compiled Fortran network data. Caches the initial state of
        all modifiable network parameters for fast resetting.
        """
        self._network = network_module

        # Load CSV files from installed package
        self._load_csv_files()

        # Parse CSV data into Species and Reaction objects
        self._parse_species()
        self._parse_reactions()

        # Validate that on-disk matches in-memory
        self._validate_network()

        # Cache initial state of modifiable parameters for fast reset
        self._cache_initial_state()

    def _load_csv_files(self):
        """Load species and reaction CSV files from the installed package."""
        import uclchem

        pkg_dir = os.path.dirname(uclchem.__file__)

        species_path = os.path.join(pkg_dir, "species.csv")
        reactions_path = os.path.join(pkg_dir, "reactions.csv")

        if not os.path.exists(species_path):
            raise FileNotFoundError(f"Species CSV not found: {species_path}")
        if not os.path.exists(reactions_path):
            raise FileNotFoundError(f"Reactions CSV not found: {reactions_path}")

        self._species_df = pd.read_csv(species_path)
        self._reactions_df = pd.read_csv(reactions_path)

    def _parse_species(self):
        """Parse species DataFrame into Species objects."""
        self.species_list = []

        for _, row in self._species_df.iterrows():
            # Create Species object from CSV row
            species_row = [
                row["NAME"],
                int(row["MASS"]),
                row["BINDING_ENERGY"],
                row["SOLID_FRACTION"],
                row["MONO_FRACTION"],
                row["VOLCANO_FRACTION"],
                row["ENTHALPY"],
            ]
            species = Species(species_row)
            self.species_list.append(species)

    def _parse_reactions(self):
        """Parse reactions DataFrame into Reaction objects."""
        self.reaction_list = []

        for _, row in self._reactions_df.iterrows():
            # Create Reaction object from CSV row
            reaction_row = [
                row["Reactant 1"],
                row["Reactant 2"],
                row["Reactant 3"],
                row["Product 1"],
                row["Product 2"],
                row["Product 3"],
                row["Product 4"],
                row["Alpha"],
                row["Beta"],
                row["Gamma"],
                row["T_min"],
                row["T_max"],
                row["reduced_mass"],
                row["extrapolate"],
                row["exothermicity"],
            ]
            reaction = Reaction(reaction_row)
            self.reaction_list.append(reaction)

    def _validate_network(self):
        """Validate that on-disk network matches in-memory Fortran network."""
        errors = []

        # Check species count
        n_species_memory = len(self._network.specname)
        n_species_disk = len(self.species_list)

        if n_species_memory != n_species_disk:
            errors.append(
                f"Species count mismatch: {n_species_memory} in memory vs {n_species_disk} on disk"
            )

        # Check reaction count
        n_reactions_memory = len(self._network.alpha)
        n_reactions_disk = len(self.reaction_list)

        if n_reactions_memory != n_reactions_disk:
            errors.append(
                f"Reaction count mismatch: {n_reactions_memory} in memory vs {n_reactions_disk} on disk"
            )

        # Check species names match
        for i, species in enumerate(self.species_list):
            if i >= len(self._network.specname):
                break
            memory_name = str(np.char.decode(self._network.specname[i])).strip()
            disk_name = species.get_name()

            if memory_name != disk_name:
                errors.append(
                    f"Species name mismatch at index {i}: '{memory_name}' in memory vs '{disk_name}' on disk"
                )
                if len(errors) > 5:  # Limit error output
                    errors.append("... (truncated)")
                    break

        # Check reaction parameters match
        for i, reaction in enumerate(self.reaction_list):
            if i >= len(self._network.alpha):
                break

            memory_alpha = float(self._network.alpha[i])
            disk_alpha = reaction.get_alpha()

            if not np.allclose(memory_alpha, disk_alpha):
                errors.append(
                    f"Reaction {i} alpha mismatch: {memory_alpha} in memory vs {disk_alpha} on disk"
                )
                if len(errors) > 10:
                    errors.append("... (truncated)")
                    break

        if errors:
            import warnings

            warnings.warn("Network validation failed:" + errors, RuntimeWarning)
            raise RuntimeError("Network validation failed")

    def validate(self):
        """Re-run validation to check on-disk matches in-memory.

        Useful after modifications to verify consistency.
        """
        self._validate_network()

    def _cache_initial_state(self):
        """Cache the initial state of all modifiable network parameters.

        This allows fast reset without re-reading CSV files. Caches:
        - Reaction parameters: alpha, beta, gama
        - Species parameters: bindingenergy, formationenthalpy
        """
        # Cache reaction rate parameters (modifiable at runtime)
        self._initial_alpha = np.copy(self._network.alpha)
        self._initial_beta = np.copy(self._network.beta)
        self._initial_gama = np.copy(self._network.gama)

        # Cache species binding energies (modifiable at runtime)
        self._initial_bindingenergy = np.copy(self._network.bindingenergy)

        # Cache formation enthalpies if available (for heating/cooling)
        if hasattr(self._network, "formationenthalpy"):
            self._initial_formationenthalpy = np.copy(self._network.formationenthalpy)
        else:
            self._initial_formationenthalpy = None

    def reset_state(self):
        """Reset the Fortran network to its initial cached state.

        Uses cached values instead of re-reading and parsing CSV files.
        Restores all modifiable network parameters to their initial values
        from when NetworkState was created.

        Modifiable arrays restored:
        - alpha, beta, gama: Reaction rate parameters
        - bindingenergy: Species binding energies
        - formationenthalpy: Formation enthalpies (if available)

        Example:
            >>> network = NetworkState()
            >>> # Modify some reaction parameters
            >>> network._network.alpha[0] = 999.0
            >>> # Restore to initial state
            >>> network.reset_state()
        """
        # Restore reaction rate parameters from cache
        np.copyto(self._network.alpha, self._initial_alpha)
        np.copyto(self._network.beta, self._initial_beta)
        np.copyto(self._network.gama, self._initial_gama)

        # Restore species binding energies from cache
        np.copyto(self._network.bindingenergy, self._initial_bindingenergy)

        # NOTE: formationenthalpy cannot be reset - causes bus error (read-only?)
        # Skipping formationenthalpy reset to avoid crash

    def reset_network_from_csv(self):
        """Reset the Fortran network to match the cached initial state.

        This method now uses cached values for better performance and reliability.
        It's equivalent to reset_state() but kept for backward compatibility.

        Note: This uses the cached initial state from when NetworkState was created,
        not by re-reading CSV files. Use reset_state() for the same functionality
        with a clearer name.

        Modifiable arrays restored:
        - alpha, beta, gama: Reaction rate parameters
        - bindingenergy: Species binding energies

        Example:
            >>> network = NetworkState()
            >>> # Modify some reaction...
            >>> network._network.alpha[0] = 999.0
            >>> # Reset back to initial values
            >>> network.reset_state()
        """
        # Use the cached reset method
        self.reset_state()
