"""UCLCHEM Species."""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Iterator
from typing import Any
from warnings import warn

import pandas as pd

from uclchem.utils import find_number_of_consecutive_digits

element_list = [
    "H",
    "D",
    "HE",
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "CL",
    "LI",
    "NA",
    "MG",
    "SI",
    "PAH",
    "15N",
    "13C",
    "18O",
    "E-",
    "FE",
]
element_mass = [
    1,
    2,
    4,
    12,
    14,
    16,
    19,
    31,
    32,
    35,
    3,
    23,
    24,
    28,
    420,
    15,
    13,
    18,
    0,
    56,
]
symbols = ["#", "@", "*", "+", "-", "(", ")"]


def normalize_species_name(name: str) -> str:
    """Normalize a species name to a canonical form.

    Empty strings are preserved as empty strings. Other falsy values (like None)
    are converted to "NAN". Grain prefixes (#/@) are preserved as-is.
    A chemical isomer prefix — a single alphabetic character followed by a hyphen
    (e.g. 'o-', 'p-', 'a-', 'l-') — is lowercased so that input is case-insensitive.
    The chemical formula part is uppercased. All other names are simply uppercased.

    Args:
        name (str): name to normalize.

    Returns:
        str: Normalized species name

    Examples:
        >>> normalize_species_name("o-H2")
        'o-H2'
        >>> normalize_species_name("O-H2") # case-normalised prefix
        'o-H2'
        >>> normalize_species_name("#o-H2")
        '#o-H2'
        >>> normalize_species_name("C-") # negative ion: len == 2, not a prefix
        'C-'
        >>> normalize_species_name("E-") # electron: same rule
        'E-'
        >>> normalize_species_name("H2O")
        'H2O'
        >>> normalize_species_name("") # empty string
        ''
        >>> normalize_species_name(None) # falsy non-string value
        'NAN'

    """
    # Preserve empty strings; convert other falsy values to "NAN"
    if name == "":
        return ""
    if not name:
        return "NAN"
    grain_prefix = ""
    rest = name
    if rest[0] in ("#", "@"):
        grain_prefix = rest[0]
        rest = rest[1:]
    # A chemical prefix is exactly one alpha char + '-' with more formula after it.
    if len(rest) > 2 and rest[1] == "-" and rest[0].isalpha():
        return grain_prefix + rest[0].lower() + "-" + rest[2:].upper()
    return grain_prefix + rest.upper()


species_header = (
    "NAME",
    "MASS",
    "BINDING_ENERGY",
    "SOLID_FRACTION",
    "MONO_FRACTION",
    "VOLCANO_FRACTION",
    "ENTHALPY",
    "DESORPTION_PREF",
    "DIFFUSION_BARRIER",
    "DIFFUSION_PREF",
    "IX",
    "IY",
    "IZ",
    "SYMMETRY_NUMBER",
)


def is_number(s: Any) -> bool:
    """Try to convert input to a float, if it succeeds, return True.

    Args:
        s (Any): Input object to check

    Returns:
        bool: True if a number, False if not.

    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def sanitize_input_float(row: list[Any], index: int, default: Any = 0.0) -> float:
    """Sanitize the input. If the index is out of bounds of the row or the value
    from the row cannot be turned into a float, use the ``default`` value.
    Otherwise, just gets the value from the row.

    Args:
        row (list[Any]): list of objects
        index (int): index within list to use
        default (Any): default value to use. Default = 0.0.

    Returns:
        float: sanitized value.

    """
    output = default
    if len(row) > index and is_number(row[index]):
        output = float(row[index])
    return output


class Species:
    """Species is a class that holds all the information about an individual species in the
    network. It also has convenience functions to check whether the species is a gas or grain
    species and to help compare between species.
    """

    def __init__(self, input_row: list[str] | pd.Series):
        """Parse a species row.

        Uses the new extended order:

        NAME,MASS,BINDING_ENERGY,SOLID_FRACTION,MONO_FRACTION,VOLCANO_FRACTION,ENTHALPY,
        DESORPTION_PREF,DIFFUSION_BARRIER,DIFFUSION_PREF,Ix,Iy,Iz,SYMMETRYFACTOR

        Falls back to sensible defaults when fields are missing.
        """
        if isinstance(input_row, pd.Series):
            input_row = [input_row[field] for field in species_header]

        self.name = normalize_species_name(str(input_row[0]))
        # Detect chemical isomer prefix (e.g. 'o' from 'o-H2' or '#o-H2').
        _rest = self.name[1:] if self.name and self.name[0] in ("#", "@") else self.name
        self.prefix = (
            _rest[0]
            if (len(_rest) > 2 and _rest[1] == "-" and _rest[0].islower())
            else ""
        )
        self.mass = int(input_row[1])

        # binding energy and refractory handling
        try:
            self.is_refractory = str(input_row[2]).lower() == "inf"
            if self.is_refractory:
                self.binding_energy = 99.9e9
            else:
                self.binding_energy = float(input_row[2])
        except (IndexError, ValueError):
            self.is_refractory = False
            self.binding_energy = 0.0

        # core solid/mono/volcano/enthalpy
        self.solidFraction = sanitize_input_float(input_row, 3, 0.0)
        self.monoFraction = sanitize_input_float(input_row, 4, 0.0)
        self.volcFraction = sanitize_input_float(input_row, 5, 0.0)
        self.enthalpy = sanitize_input_float(input_row, 6, 0.0)

        # extended Tobias fields after enthalpy
        # vdes (desorption attempt frequency)
        self.vdes = sanitize_input_float(input_row, 7, 0.0)

        self.diffusion_barrier = sanitize_input_float(input_row, 8, 0.0)
        # vdiff (diffusion attempt frequency)
        self.vdiff = sanitize_input_float(input_row, 9, 0.0)

        # TST prefactors (optional)
        try:
            self.Ix = float(input_row[10])
        except (IndexError, ValueError):
            self.Ix = -999.0
        try:
            self.Iy = float(input_row[11])
        except (IndexError, ValueError):
            self.Iy = -999.0
        try:
            self.Iz = float(input_row[12])
        except (IndexError, ValueError):
            self.Iz = -999.0
        try:
            self.symmetry_factor = int(input_row[13])
        except (IndexError, ValueError, TypeError):
            self.symmetry_factor = -1

        self.set_n_atoms(0)

        # in first instance, assume species freeze/desorb unchanged
        # this is updated by ``check_freeze_desorbs()`` later.
        if self.is_ice_species():
            # this will make any excited species desorb as their base counterparts
            if "*" in self.get_name():
                self.desorb_products = [self.get_name()[1:-1], "NAN", "NAN", "NAN"]
            else:
                self.desorb_products = [self.get_name()[1:], "NAN", "NAN", "NAN"]
        else:
            self.freeze_products = {}

    def get_name(self) -> str:
        """Get the name of the chemical species.

        Returns:
            str: The name

        """
        return self.name

    def get_mass(self) -> int:
        """Get the molecular mass of the chemical species.

        Returns:
            int: The molecular mass in atomic mass units.

        """
        return self.mass

    def get_charge(self) -> int:
        """Get the charge of the chemical species in e.
        Positive integer indicates positive ion, negative indicates negative ion.
        Assumes species are at most charged +1 or -1.

        Returns:
            int: The charge of the species

        """
        if not self.is_ion():
            return 0

        if "+" in self.get_name():
            return 1
        elif self.get_name().endswith("-"):
            return -1

    def get_solid_fraction(self) -> float:
        """Get the solid fraction of the species.

        Returns:
            float: The solid fraction

        """
        return self.solidFraction

    def get_mono_fraction(self) -> float:
        """Get the monolayer fraction of the species.

        Returns:
            float: The monolayer fraction

        """
        return self.monoFraction

    def get_volcano_fraction(self) -> float:
        """Get the volcano fraction of the species.

        Returns:
            float: The volcano fraction

        """
        return self.volcFraction

    def get_enthalpy(self) -> float:
        """Get the ice enthalpy of the species.

        Returns:
            float: The ice enthalpy

        """
        return self.enthalpy

    def set_name(self, name: str) -> None:
        """Set the name of the chemical species.

        Args:
            name (str): The new name for the species

        """
        self.name = normalize_species_name(name)

    def set_mass(self, mass: int) -> None:
        """Set the molecular mass of the chemical species in atomic mass units.

        Args:
            mass (int): The new molecular mass

        """
        self.mass = int(mass)

    def set_binding_energy(self, binding_energy: float) -> None:
        """Set the binding energy of the species in K.

        Args:
            binding_energy (float): The new binding energy in K

        """
        self.binding_energy = float(binding_energy)

    def set_solid_fraction(self, solid_fraction: float) -> None:
        """Set the solid fraction of the species.

        Args:
            solid_fraction (float): The new solid fraction

        """
        self.solidFraction = float(solid_fraction)

    def set_mono_fraction(self, mono_fraction: float) -> None:
        """Set the monolayer fraction of the species.

        Args:
            mono_fraction (float): The new monolayer fraction

        """
        self.monoFraction = float(mono_fraction)

    def set_volcano_fraction(self, volcano_fraction: float) -> None:
        """Set the volcano fraction of the species.

        Args:
            volcano_fraction (float): The new volcano fraction

        """
        self.volcFraction = float(volcano_fraction)

    def set_enthalpy(self, enthalpy: float) -> None:
        """Set the enthalpy of the species in kcal per mole.

        Args:
            enthalpy (float): The new ice enthalpy

        """
        self.enthalpy = float(enthalpy)

    def set_desorb_products(self, new_desorbs: list[str]) -> None:
        """Set the desorption products for species on the surface or in the bulk.
        It is assumed that there is only one desorption pathway.

        Args:
            new_desorbs (list[str]): The new desorption products

        """
        self.desorb_products = new_desorbs

    def get_standard_desorb_products(self) -> list[str]:
        """Return the 1:1 gas-phase counterpart by stripping the grain prefix.

        For #CH4 returns [CH4, NAN, NAN, NAN]. Always a single product — used
        for gasIceList construction and auto-generated THERM/DESOH2/DESCR/DEUVCR
        reactions when the user has not provided explicit ones.

        Returns:
            list[str]: [base_gas_species, NAN, NAN, NAN]
        """
        return [self.get_name()[1:], "NAN", "NAN", "NAN"]

    def get_desorb_products(self) -> list[str]:
        """Obtain the desorbtion products of ice species.

        Returns:
            list[str]: The desorption products

        Raises:
            AttributeError: If the species has no attribute ``desorb_products``.
                This can occur if the species is a gas-phase species.

        """
        if not hasattr(self, "desorb_products"):
            msg = f"Species {self} has no attribute 'desorb_products'"
            raise AttributeError(msg)
        return self.desorb_products

    def set_freeze_products(self, product_list: list[str], freeze_alpha: float) -> None:
        """Add the freeze products of the species, one species can have
        several freeze products.

        Args:
            product_list (list[str]): The list of freeze out products
            freeze_alpha (float): The freeze out ratio.

        It is called alpha, since it is derived from the alpha column
        in the UCLCHEM reaction format:
        https://github.com/uclchem/UCLCHEM/blob/08d37f8c3063f8ff8a9a7aa16d9eff0ed4f99538/Makerates/src/network.py#L160

        """
        self.freeze_products[",".join(product_list)] = freeze_alpha

    def get_freeze_products(self) -> Iterator[dict[list[str], float]]:
        """Obtain the product to which the species freeze out.

        Yields:
            Iterator[dict[str, float]]: Iterator that returns all of the
                freeze out reactions with ratios

        """
        keys = self.freeze_products.keys()
        values = self.freeze_products.values()
        logging.debug(f"freeze keys: {keys}, products {values}")
        for key, value in zip(keys, values):
            yield key.split(","), value

    def get_freeze_products_list(self) -> list[list[str]]:
        """Get all the freeze products without their ratios.

        Returns:
            list[list[str]]: List of freeze products

        """
        # TODO: Write an unit test for get_freeze_product_behaviour
        return [key.split(",") for key in self.freeze_products]

    def get_freeze_alpha(self, product_list: list[str]) -> float:
        """Obtain the freeze out ratio of a species for a certain reaction.

        Args:
            product_list (list[str]): For a specific reaction, get the freezeout ratio

        Returns:
            float: The freezeout ratio

        """
        return self.freeze_products[",".join(product_list)]

    def is_grain_species(self) -> bool:
        """Return whether the species is a species on the grain.

        Returns:
            bool: True if it is a grain species.

        """
        warn(
            "This method is deprecated in favour of is_ice_species.",
            DeprecationWarning,
            stacklevel=2,
        )
        return (
            self.get_name() in ["BULK", "SURFACE"]
            or self.get_name().startswith(
                "#",
            )
            or self.get_name().startswith("@")
        )

    def is_ice_species(self) -> bool:
        """Return whether the species is a species on the grain.

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

    def is_surface_species(self) -> bool:
        """Check if the species is on the surface.

        Returns:
            bool: True if a surface species

        """
        return self.get_name().startswith("#")

    def is_bulk_species(self) -> bool:
        """Check if the species is in the bulk.

        Returns:
            bool: True if a bulk species

        """
        return self.get_name().startswith("@")

    def is_ion(self) -> bool:
        """Check if the species is ionized, either positively or negatively.

        Returns:
            bool: True if it is ionized

        """
        return self.get_name().endswith("+") or self.get_name().endswith("-")

    def add_default_freeze(self) -> None:
        """Add a default freezeout, which is freezing out to the species itself,
        but with no ionization.
        """
        freeze = "#" + self.get_name()
        if freeze[-1] in ["+", "-"]:
            freeze = freeze[:-1]
        if self.get_name() == "E-":
            freeze = ""
        self.set_freeze_products([freeze, "NAN", "NAN", "NAN"], 1.0)

    def find_constituents(self, quiet: bool = False) -> Counter[str, int]:
        """Loop through the species' name and work out what its consituent
        atoms are. Then calculate mass and alert user if it doesn't match
        input mass.

        Returns:
            counter (Counter): Counter of how many times each element is in the molecule.

        Examples:
            >>> species = Species(['H2'] + [0] * 10)
            >>> constituents = species.find_constituents()
            >>> # Has the right number of H atoms
            >>> constituents['H']
            2
            >>> # And 0 of the other atoms
            >>> constituents['O']
            0

            >>> species = Species(['(CH3)2'] + [0] * 10)
            >>> constituents = species.find_constituents()
            >>> constituents['C'], constituents["H"]
            (2, 6)

            >>> species = Species(['C60'] + [0] * 10)
            >>> constituents = species.find_constituents()
            >>> constituents['C']
            60

        Raises:
            ValueError: If the molecular formula is not valid, for example it has an
                element not in the element list, has no closing bracket, or starts with a digit.
        """
        # Adapted from https://github.com/uclchem/UCLCHEM/blob/main/src/uclchem/makerates/species.py
        name = self.name
        # Strip chemical isomer prefix (e.g. 'o-' from 'o-H2' or '#o-H2') so the
        # element parser only sees the plain formula.
        if self.prefix:
            if name and name[0] in ("#", "@"):
                # keep the grain prefix, remove 'x-' immediately after it
                name = name[0] + name[len(self.prefix) + 2 :]
            else:
                name = name[len(self.prefix) + 1 :]
        if name[0].isdigit():
            msg = f"First character of formula {name} was a digit. Please put repeated parts in a bracket with number after, e.g. (CH3)2"
            raise ValueError(msg)

        char_idx = 0
        atoms = []
        currently_in_bracket = False
        j = None
        # loop over characters in species name to work out what it is made of
        while char_idx < len(name):
            # if character isn't a + or - then check it, otherwise move on
            if name[char_idx] not in symbols:
                if (
                    char_idx + 1 < len(name)
                    and name[char_idx : char_idx + 2] in element_list
                ):
                    # if next two characters are (eg) 'MG' then atom is Mg not M and G
                    j = char_idx + 2
                # if there aren't two characters left just try next one
                elif name[char_idx] in element_list:
                    j = char_idx + 1

                # if we've found a new element check for numbers otherwise print error
                if j is None or j <= char_idx:
                    msg = f"formula {name} contains element(s) not in element list"
                    raise ValueError(msg)

                num_digits = find_number_of_consecutive_digits(name, j)
                if num_digits == 0:
                    nrepeat = 1
                else:
                    nrepeat = int(name[j : j + num_digits])
                for _ in range(nrepeat):
                    if currently_in_bracket:
                        bracket_content.append(name[char_idx:j])  # noqa: F821
                    else:
                        atoms.append(name[char_idx:j])
                char_idx = j + num_digits
            else:
                # if symbol is start of a bracketed part of molecule, keep track
                if name[char_idx] == "(":
                    currently_in_bracket = True
                    bracket_content = []
                # if it's the end then add bracket contents to list
                elif name[char_idx] == ")":
                    if not currently_in_bracket:
                        msg = f"Found closing bracket before opening bracket in formula {name}"
                        raise ValueError(msg)
                    currently_in_bracket = False
                    num_digits = find_number_of_consecutive_digits(name, char_idx + 1)
                    if num_digits == 0:
                        nrepeat = 1
                    else:
                        nrepeat = int(name[char_idx + 1 : char_idx + 1 + num_digits])
                    for _ in range(nrepeat):
                        atoms.extend(bracket_content)
                    char_idx += num_digits
                char_idx += 1
        if currently_in_bracket:
            msg = f"Opening bracket was not closed in formula {name}"
            raise ValueError(msg)
        counter = Counter()
        for element in element_list:
            counter[element] = atoms.count(element)

        mass = 0
        for atom in atoms:
            mass += element_mass[element_list.index(atom)]
        if mass != int(self.get_mass()):
            if not quiet:
                logging.warning(
                    f"Input mass of {self.get_name()} ({self.get_mass()}) does not match calculated mass of constituents, using calculated mass: {int(mass)}"
                )
            self.set_mass(int(mass))
        self.n_atoms = counter.total()

        return counter

    def get_n_atoms(self) -> int:
        """Get the number of atoms in the molecule.

        Returns:
            int: The number of atoms

        """
        return self.n_atoms

    def set_n_atoms(self, new_n_atoms: int) -> None:
        """Set the number of atoms.

        Args:
            new_n_atoms (int): The new number of atoms

        """
        self.n_atoms = new_n_atoms

    def get_binding_energy(self) -> float:
        """Get the binding energy of the species in K.

        Returns:
            float: The binding energy in K

        """
        return self.binding_energy

    def get_vdes(self) -> float:
        """Get the desorption prefactor.

        Returns:
            float: The desorption prefactor in s-1

        """
        return float(self.vdes)

    def get_desorption_pref(self) -> float:
        """Get the desorption prefactor.

        Alias getter matching CSV column name ``desorption_pref``.

        Returns:
            float: The desorption prefactor in s-1

        """
        return self.get_vdes()

    def get_diffusion_barrier(self) -> float:
        """Get the diffusion barrier for the species.

        Returns:
            float: The diffusion barrier in K

        """
        return self.diffusion_barrier

    def get_vdiff(self) -> float:
        """Get the diffusion prefactor.

        Returns:
            float: The diffusion prefactor in s-1

        """
        return float(self.vdiff)

    def get_diffusion_pref(self) -> float:
        """Set the diffusion prefactor.

        Alias getter matching CSV column name ``diffusion_pref``.

        Returns:
            float: The diffusion prefactor in s-1

        """
        return self.get_vdiff()

    def get_Ix(self) -> float:  # noqa: N802
        """Set the moment of inertia along the first principal axis.

        Returns:
            float: moment of inertia in amu/Angstrom^2

        """
        return self.Ix

    def get_Iy(self) -> float:  # noqa: N802
        """Set the moment of inertia along the second principal axis.

        Returns:
            float: moment of inertia in amu/Angstrom^2

        """
        return self.Iy

    def get_Iz(self) -> float:  # noqa: N802
        """Set the moment of inertia along the third principal axis.

        Returns:
            float: moment of inertia in amu/Angstrom^2

        """
        return self.Iz

    def get_symmetry_factor(self) -> int:
        """Get the symmetry factor of the species.

        Returns:
            int: Symmetry factor

        """
        return self.symmetry_factor

    def set_vdes(self, vdes: float) -> None:
        """Set the desorption prefactor.

        Args:
            vdes (float): The desorption prefactor in s-1

        """
        self.vdes = float(vdes)

    def set_desorption_pref(self, vdes: float) -> None:
        """Set the desorption prefactor.

        Alias setter matching CSV column name ``desorption_pref``.

        Args:
            vdes (float): The desorption prefactor in s-1

        """
        self.set_vdes(vdes)

    def set_diffusion_barrier(self, barrier: float) -> None:
        """Set the diffusion barrier for the species.

        Args:
            barrier (float): Diffusion barrier in K

        """
        self.diffusion_barrier = float(barrier)

    def set_vdiff(self, vdiff: float) -> None:
        """Set the diffusion prefactor (vdiff) for the species.

        Args:
            vdiff (float): The diffusion prefactor in s-1

        """
        self.vdiff = float(vdiff)

    def set_diffusion_pref(self, vdiff: float) -> None:
        """Set the diffusion prefactor.

        Alias setter matching CSV column name ``diffusion_pref``.

        Args:
            vdiff (float): The diffusion prefactor in s-1

        """
        self.set_vdiff(vdiff)

    def set_Ix(self, Ix: float) -> None:  # noqa: N802, N803
        """Set the moment of inertia along the first principal axis.

        Args:
            Ix (float): desired moment of inertia (in amu/Angstrom^2)

        """
        self.Ix = float(Ix)

    def set_Iy(self, Iy: float) -> None:  # noqa: N802, N803
        """Set the moment of inertia along the second principal axis.

        Args:
            Iy (float): desired moment of inertia (in amu/Angstrom^2)

        """
        self.Iy = float(Iy)

    def set_Iz(self, Iz: float) -> None:  # noqa: N802, N803
        """Set the moment of inertia along the third principal axis.

        Args:
            Iz (float): desired moment of inertia (in amu/Angstrom^2)

        """
        self.Iz = float(Iz)

    def set_symmetry_factor(self, sym: int | str) -> None:
        """Set the symmetry factor of the species.

        Sets the symmetry factor to -1 if ``sym`` cannot be turned into an
        integer.

        Args:
            sym (int | str): Symmetry factor

        """
        try:
            self.symmetry_factor = int(sym)
        except (ValueError, TypeError):
            self.symmetry_factor = -1

    def __eq__(self, other: str | Species) -> bool:
        """Check for equality based on either a string or another Species instance.

        Args:
            other (str, Species): Another species

        Returns:
            bool: True if two species are identical.

        Raises:
            NotImplementedError: We can only compare between species or strings of species.

        """
        if isinstance(other, Species):
            return self.get_name() == other.get_name()
        elif isinstance(other, str):
            return self.get_name() == other
        else:
            msg = "We can only compare between species or strings of species"
            raise NotImplementedError(msg)

    def __lt__(self, other: Species) -> bool:
        """Compare the mass of the species.

        Args:
            other (Species): Another species instance

        Returns:
            bool: True if less than the other species

        """
        return self.get_mass() < other.get_mass()

    def __gt__(self, other: Species) -> bool:
        """Compare the mass of the species.

        Args:
            other (Species): Another species instance

        Returns:
            bool: True if larger than than the other species

        """
        return self.get_mass() > other.get_mass()

    def __repr__(self) -> str:
        return f"Specie: {self.get_name()}"

    def __str__(self) -> str:
        return self.get_name()

    def calculate_rotational_partition_factor(self) -> float:
        """Calculate 1/sigma*(SQRT(IxIyIz)) for non-linear molecules, and
        1/sigma*(SQRT(IyIz)) for linear molecules.

        Returns -999.0 if molecular inertia data is not available (backward compatibility).
        This signals that TST prefactors cannot be used for this species.

        Returns:
            float: Rotational partition factor scaled by 1e50, or -999.0 if unavailable

        """
        if self.n_atoms == 1:
            # For atoms, this is undefined, just return a value such that
            # it is clearly an atomic species.
            return -1.0
        if self.Ix == -999.0 or self.Iy == -999.0 or self.Iz == -999.0:
            # For species without custom input Ix, Iy and Iz, we cannot do this,
            # Return sentinel value for backward compatibility
            return -999.0
        if self.symmetry_factor <= 0:
            # No symmetry factor provided
            return -999.0

        # Ix, Iy and Iz are in units of amu Angstrom^2,
        # so need to convert to kg m2
        import numpy as np

        amu = 1.66053907e-27  # kg/amu
        scalingFactor = 1e50  # noqa: N806

        if not self.is_linear():
            return (
                (1.0 / self.symmetry_factor)
                * np.sqrt(self.Ix * self.Iy * self.Iz * amu**3 / 1e60)
                * scalingFactor
            )
        else:
            return (
                (1.0 / self.symmetry_factor)
                * np.sqrt(self.Iy * self.Iz * amu**2 / 1e40)
                * scalingFactor
            )

    def is_linear(self) -> bool:
        """Check if molecule is linear based on moment of inertia.

        For linear molecules, Ix = 0 (rotation axis along molecular axis has no inertia).

        Returns:
            bool: True if linear, False otherwise

        """
        if self.n_atoms == 1:
            # Atomic species are not linear (doesn't matter, filtered out anyway)
            return False
        if self.n_atoms == 2:
            # Diatomic molecules are always linear
            return True
        if self.Ix == -999.0 or self.Iy == -999.0 or self.Iz == -999.0:
            # No inertia data available
            return False
        if not self.is_ice_species():
            # Only implement for grain species
            return False
        return self.Ix == 0.0

    def check_symmetry_factor(self) -> None:
        """Check the symmetry factor provided by the user.

        Checks if n_atoms == 2, that if its homoatomic (e.g. H2), that
        sigma == 2, and if it is heteroatomic, (e.g. OH), sigma == 1

        """
        if self.n_atoms == 1:  # Nothing to check
            return
        if self.n_atoms > 2:  # Can not correctly check everything
            return
        constituents = self.find_constituents(quiet=True)
        if (
            len(constituents) == 2
        ):  # Only one constituent, i.e. both atoms are the same element.
            if self.symmetry_factor == 1:
                return
            msg = f"For diatomic molecule consisting of two different atoms (in this case {self.name}), the symmetry factor should be 1, but was given to be {self.symmetry_factor}. Correcting to 1."
            logging.warning(msg)
            self.symmetry_factor = 1
            return
        if self.symmetry_factor == 2:
            return
        msg = f"For diatomic molecule consisting of two of the same atoms (in this case {self.name}), the symmetry factor should be 2, but was given to be {self.symmetry_factor}. Correcting to 2."
        logging.warning(msg)
        self.symmetry_factor = 2
