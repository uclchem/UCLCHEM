"""UCLCHEM Species."""

from __future__ import annotations

import logging
from collections import Counter
from typing import TYPE_CHECKING
from warnings import warn

import numpy as np
import pandas as pd

from uclchem.makerates.utils import (
    find_number_of_consecutive_digits,
    normalize_species_name,
    sanitize_input_float,
)
from uclchem.utils import (
    MISSING_VALUE_FLOAT,
    MISSING_VALUE_INTEGER,
)

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger(__name__)

elementList = [
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
elementMass = [
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
IGNORED_SPECIES_PARSING_SYMBOLS = ["#", "@", "*", "+", "-", "(", ")"]


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


def get_element_counts_per_species(
    species_list: list[Species],
) -> tuple[tuple[str, ...], np.typing.NDArray[np.int32]]:
    """Get the number of times each element occurs in each species.

    Generic element-count 2D array for runtime conservation checking.
    Covers every element that appears at least once across all species.

    Parameters
    ----------
    species_list : list[Species]
        List of species

    Returns
    -------
    unique_elements : tuple[str]
        A sorted list of all the unique elements in the list of species.
    elem_count_2d : np.typing.NDArray[np.int_]
        2D array with shape ``n_spec * n_elem`` of how often each element
        occurs in each species.

    """
    all_constituents: list[Counter[str]] = []
    unique_elements_set: set[str] = set()
    for species in species_list:
        try:
            constituents = species.find_constituents(quiet=True)
            all_constituents.append(constituents)
            unique_elements_set.update(
                element for element, count in constituents.items() if count > 0
            )
        except (ValueError, Exception):
            all_constituents.append(Counter())

    unique_elements = tuple(
        sorted(e for e in unique_elements_set if e.upper() not in {"E", "E-"})
    )

    n_elems = len(unique_elements)
    elem_count_2d = np.zeros((len(species_list), n_elems), dtype=int)
    for si, elem_constituents in enumerate(all_constituents):
        for ei, elem in enumerate(unique_elements):
            elem_count_2d[si, ei] = elem_constituents.get(elem, 0)

    return unique_elements, elem_count_2d


def determine_molecular_mass(constituents: Counter[str]) -> int:
    """Determine the molecular mass for a set of atoms.

    Parameters
    ----------
    constituents : Counter[str]
        Counter of elemental makeup, for example obtained from
        :func:`determine_constituents`.

    Returns
    -------
    mass : int
        Mass of the species.

    """
    mass = 0
    for element, count in constituents.items():
        mass += elementMass[elementList.index(element)] * count
    return mass


def determine_constituents(normalized_species_name: str) -> Counter[str]:
    """Loop through a species' name and work out what its constituent atoms are.

    Parameters
    ----------
    normalized_species_name : str
        Species name that has been normalized by :func:`normalize_species_name`
        if necessary.

    Returns
    -------
    Counter[str]
        Counter of how many times each element is in the molecule.

    Raises
    ------
    ValueError
        If the species name is invalid, meaning contains species not in the
        elements list, has an opening but no closing bracket, etc.

    Examples
    --------
    >>> constituents = determine_constituents('H2')
    >>> # Has the right number of H atoms
    >>> constituents['H']
    2
    >>> # And 0 of the other atoms
    >>> constituents['O']
    0

    >>> constituents = determine_constituents('(CH3)2')
    >>> constituents['C'], constituents["H"]
    (2, 6)

    >>> constituents = determine_constituents('C60')
    >>> constituents['C']
    60

    """
    if normalized_species_name[0].isdigit():
        msg = f"First character of formula {normalized_species_name} was a digit. Please put repeated parts in a bracket with number after, e.g. (CH3)2"
        raise ValueError(msg)

    char_idx = 0
    atoms: list[str] = []
    currently_in_bracket = False
    bracket_content: list[str] = []
    j = None
    # loop over characters in species name to work out what it is made of
    while char_idx < len(normalized_species_name):
        # if character isn't a + or - then check it, otherwise move on
        if normalized_species_name[char_idx] not in IGNORED_SPECIES_PARSING_SYMBOLS:
            if (
                char_idx + 1 < len(normalized_species_name)
                and normalized_species_name[char_idx : char_idx + 2] in elementList
            ):
                # if next two characters are (eg) 'MG' then atom is Mg not M and G
                j = char_idx + 2
            # if there aren't two characters left just try next one
            elif normalized_species_name[char_idx] in elementList:
                j = char_idx + 1

            # if we've found a new element check for numbers otherwise print error
            if j is None or j <= char_idx:
                msg = f"formula {normalized_species_name} contains element(s) not in element list"
                raise ValueError(msg)

            num_digits = find_number_of_consecutive_digits(normalized_species_name, j)
            if num_digits == 0:
                nrepeat = 1
            else:
                nrepeat = int(normalized_species_name[j : j + num_digits])
            for _ in range(nrepeat):
                if currently_in_bracket:
                    bracket_content.append(normalized_species_name[char_idx:j])
                else:
                    atoms.append(normalized_species_name[char_idx:j])
            char_idx = j + num_digits
        else:
            # if symbol is start of a bracketed part of molecule, keep track
            if normalized_species_name[char_idx] == "(":
                currently_in_bracket = True
                bracket_content = []
            # if it's the end then add bracket contents to list
            elif normalized_species_name[char_idx] == ")":
                if not currently_in_bracket:
                    msg = f"Found closing bracket before opening bracket in formula {normalized_species_name}"
                    raise ValueError(msg)
                currently_in_bracket = False
                num_digits = find_number_of_consecutive_digits(
                    normalized_species_name, char_idx + 1
                )
                if num_digits == 0:
                    nrepeat = 1
                else:
                    nrepeat = int(
                        normalized_species_name[char_idx + 1 : char_idx + 1 + num_digits]
                    )
                for _ in range(nrepeat):
                    atoms.extend(bracket_content)
                char_idx += num_digits
            char_idx += 1
    if currently_in_bracket:
        msg = f"Opening bracket was not closed in formula {normalized_species_name}"
        raise ValueError(msg)

    counter: Counter[str] = Counter()
    for element in elementList:
        counter[element] = atoms.count(element)

    return counter


class Species:
    """Species is a class that holds all the information about an individual species in the.

    network. It also has convenience functions to check whether the species is a gas or grain
    species and to help compare between species.

    """

    def __init__(self, input_row: list[str | float] | pd.Series):
        """Parse a species row.

        Uses the new extended order:

        NAME,MASS,BINDING_ENERGY,SOLID_FRACTION,MONO_FRACTION,VOLCANO_FRACTION,ENTHALPY,
        DESORPTION_PREF,DIFFUSION_BARRIER,DIFFUSION_PREF,Ix,Iy,Iz,SYMMETRYFACTOR

        Falls back to sensible defaults when fields are missing.

        Parameters
        ----------
        input_row : list[str | float] | pd.Series
            Row from the species CSV, as a list of values.

        """
        if isinstance(input_row, pd.Series):
            input_row = [input_row[field] for field in species_header]

        self.name = normalize_species_name(str(input_row[0]))
        # Detect chemical isomer prefix (e.g. 'o' from 'o-H2' or '#o-H2').
        rest = self.name[1:] if self.name and self.name[0] in {"#", "@"} else self.name
        self.prefix = (
            rest[0]
            if (len(rest) > 2 and rest[1] == "-" and rest[0].islower())  # ruff: ignore[magic-value-comparison]
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
            self.binding_energy = MISSING_VALUE_FLOAT

        # core solid/mono/volcano/enthalpy
        self.solidFraction = sanitize_input_float(input_row, 3, MISSING_VALUE_FLOAT)
        self.monoFraction = sanitize_input_float(input_row, 4, MISSING_VALUE_FLOAT)
        self.volcFraction = sanitize_input_float(input_row, 5, MISSING_VALUE_FLOAT)
        self.enthalpy = sanitize_input_float(input_row, 6, MISSING_VALUE_FLOAT)

        # extended Tobias fields after enthalpy
        # vdes (desorption attempt frequency)
        self.vdes = sanitize_input_float(input_row, 7, MISSING_VALUE_FLOAT)

        self.diffusion_barrier = sanitize_input_float(input_row, 8, MISSING_VALUE_FLOAT)
        # vdiff (diffusion attempt frequency)
        self.vdiff = sanitize_input_float(input_row, 9, MISSING_VALUE_FLOAT)

        # TST prefactors (optional)
        try:
            self.Ix = float(input_row[10])
        except (IndexError, ValueError):
            self.Ix = MISSING_VALUE_FLOAT
        try:
            self.Iy = float(input_row[11])
        except (IndexError, ValueError):
            self.Iy = MISSING_VALUE_FLOAT
        try:
            self.Iz = float(input_row[12])
        except (IndexError, ValueError):
            self.Iz = MISSING_VALUE_FLOAT
        try:
            self.symmetry_factor = int(input_row[13])
        except (IndexError, ValueError, TypeError):
            self.symmetry_factor = MISSING_VALUE_INTEGER

        self.set_n_atoms(0)

        # ODE term strings, populated by build_ode_string() in io_functions.py
        self.gains: str = ""
        self.losses: str = ""

        # in first instance, assume species freeze/desorb unchanged
        # this is updated by `check_freeze_desorbs()` later.
        if self.is_ice_species():
            # this will make any excited species desorb as their base counterparts
            if "*" in self.get_name():
                self.desorb_products = [self.get_name()[1:-1], "NAN", "NAN", "NAN"]
            else:
                self.desorb_products = [self.get_name()[1:], "NAN", "NAN", "NAN"]
        else:
            self.freeze_products: dict[str, float] = {}

    def get_name(self) -> str:
        """Get the name of the chemical species.

        Returns
        -------
        str
            The name

        """
        return self.name

    def get_mass(self) -> int:
        """Get the molecular mass of the chemical species.

        Returns
        -------
        int
            The molecular mass in atomic mass units.

        """
        return self.mass

    def get_charge(self) -> int:
        """Get the charge of the chemical species in e.

        Positive integer indicates positive ion, negative indicates negative ion.
        Assumes species are at most charged +1 or -1.

        Returns
        -------
        int
            The charge of the species

        """
        if not self.is_ion():
            return 0

        if "+" in self.get_name():
            return 1
        if self.get_name().endswith("-"):
            return -1
        return 0

    def get_solid_fraction(self) -> float:
        """Get the solid fraction of the species.

        Returns
        -------
        float
            The solid fraction

        """
        return self.solidFraction

    def get_mono_fraction(self) -> float:
        """Get the monolayer fraction of the species.

        Returns
        -------
        float
            The monolayer fraction

        """
        return self.monoFraction

    def get_volcano_fraction(self) -> float:
        """Get the volcano fraction of the species.

        Returns
        -------
        float
            The volcano fraction

        """
        return self.volcFraction

    def get_enthalpy(self) -> float:
        """Get the ice enthalpy of the species.

        Returns
        -------
        float
            The ice enthalpy

        """
        return self.enthalpy

    def set_name(self, name: str) -> None:
        """Set the name of the chemical species.

        Parameters
        ----------
        name : str
            The new name for the species

        """
        self.name = normalize_species_name(name)

    def set_mass(self, mass: int) -> None:
        """Set the molecular mass of the chemical species in atomic mass units.

        Parameters
        ----------
        mass : int
            The new molecular mass

        """
        self.mass = int(mass)

    def set_binding_energy(self, binding_energy: float) -> None:
        """Set the binding energy of the species in K.

        Parameters
        ----------
        binding_energy : float
            The new binding energy in K

        """
        self.binding_energy = float(binding_energy)

    def set_solid_fraction(self, solid_fraction: float) -> None:
        """Set the solid fraction of the species.

        Parameters
        ----------
        solid_fraction : float
            The new solid fraction

        """
        self.solidFraction = float(solid_fraction)

    def set_mono_fraction(self, mono_fraction: float) -> None:
        """Set the monolayer fraction of the species.

        Parameters
        ----------
        mono_fraction : float
            The new monolayer fraction

        """
        self.monoFraction = float(mono_fraction)

    def set_volcano_fraction(self, volcano_fraction: float) -> None:
        """Set the volcano fraction of the species.

        Parameters
        ----------
        volcano_fraction : float
            The new volcano fraction

        """
        self.volcFraction = float(volcano_fraction)

    def set_enthalpy(self, enthalpy: float) -> None:
        """Set the enthalpy of the species in kcal per mole.

        Parameters
        ----------
        enthalpy : float
            The new ice enthalpy

        """
        self.enthalpy = float(enthalpy)

    def set_desorb_products(self, new_desorbs: list[str]) -> None:
        """Set the desorption products for species on the surface or in the bulk.

        It is assumed that there is only one desorption pathway.

        Parameters
        ----------
        new_desorbs : list[str]
            The new desorption products

        """
        self.desorb_products = new_desorbs

    def get_standard_desorb_products(self) -> list[str]:
        """Return the 1:1 gas-phase counterpart by stripping the grain prefix.

        For #CH4 returns [CH4, NAN, NAN, NAN]. Always a single product — used
        for gasIceList construction and auto-generated THERM/DESOH2/DESCR/DEUVCR
        reactions when the user has not provided explicit ones.

        Returns
        -------
        list[str]
            [base_gas_species, NAN, NAN, NAN]

        """
        return [self.get_name()[1:], "NAN", "NAN", "NAN"]

    def get_desorb_products(self) -> list[str]:
        """Obtain the desorbtion products of ice species.

        Returns
        -------
        list[str]
            The desorption products

        Raises
        ------
        AttributeError
            If the species has no attribute `desorb_products`.
            This can occur if the species is a gas-phase species.

        """
        if not hasattr(self, "desorb_products"):
            msg = f"Species {self} has no attribute 'desorb_products'"
            raise AttributeError(msg)
        return self.desorb_products

    def set_freeze_products(self, product_list: list[str], freeze_alpha: float) -> None:
        """Add the freeze products of the species, one species can have.

        several freeze products.

        Parameters
        ----------
        product_list : list[str]
            The list of freeze out products
        freeze_alpha : float
            The freeze out ratio.

        """
        self.freeze_products[",".join(product_list)] = freeze_alpha

    def get_freeze_products(self) -> Iterator[tuple[list[str], float]]:
        """Obtain the product to which the species freeze out.

        Yields
        ------
        tuple[list[str], float]
            Iterator that returns all of the
            freeze out reactions with ratios

        """
        keys = self.freeze_products.keys()
        values = self.freeze_products.values()
        logger.debug(f"freeze keys: {keys}, products {values}")
        for key, value in zip(keys, values, strict=False):
            yield key.split(","), value

    def get_freeze_products_list(self) -> list[list[str]]:
        """Get all the freeze products without their ratios.

        Returns
        -------
        list[list[str]]
            List of freeze products

        """
        # TODO: Write an unit test for get_freeze_product_behavior
        return [key.split(",") for key in self.freeze_products]

    def get_freeze_alpha(self, product_list: list[str]) -> float:
        """Obtain the freeze out ratio of a species for a certain reaction.

        Parameters
        ----------
        product_list : list[str]
            For a specific reaction, get the freezeout ratio

        Returns
        -------
        float
            The freezeout ratio

        """
        return self.freeze_products[",".join(product_list)]

    def is_grain_species(self) -> bool:
        """Return whether the species is a species on the grain.

        Returns
        -------
        bool
            True if it is a grain species.

        """
        warn(
            "This method is deprecated in favor of is_ice_species.",
            DeprecationWarning,
            stacklevel=2,
        )
        return (
            self.get_name() in {"BULK", "SURFACE"}
            or self.get_name().startswith(
                "#",
            )
            or self.get_name().startswith("@")
        )

    def is_ice_species(self) -> bool:
        """Return whether the species is a species on the grain.

        Returns
        -------
        bool
            True if it is an ice species.

        """
        return (
            self.get_name() in {"BULK", "SURFACE"}
            or self.get_name().startswith(
                "#",
            )
            or self.get_name().startswith("@")
        )

    def is_surface_species(self) -> bool:
        """Check if the species is on the surface.

        Returns
        -------
        bool
            True if a surface species

        """
        return self.get_name().startswith("#")

    def is_bulk_species(self) -> bool:
        """Check if the species is in the bulk.

        Returns
        -------
        bool
            True if a bulk species

        """
        return self.get_name().startswith("@")

    def is_ion(self) -> bool:
        """Check if the species is ionized, either positively or negatively.

        Returns
        -------
        bool
            True if it is ionized

        """
        return self.get_name().endswith("+") or self.get_name().endswith("-")

    def add_default_freeze(self) -> None:
        """Add a default freezeout, which is freezing out to the species itself,.

        but with no ionization.

        """
        freeze = "#" + self.get_name()
        if freeze[-1] in {"+", "-"}:
            freeze = freeze[:-1]
        if self.get_name() == "E-":
            freeze = ""
        self.set_freeze_products([freeze, "NAN", "NAN", "NAN"], 1.0)

    def find_constituents(self, *, quiet: bool = False) -> Counter[str]:
        """Loop through the species' name and work out what its constituent atoms are.

        Then calculate mass and alert user if it doesn't match the input mass.

        Parameters
        ----------
        quiet : bool
            If ``True``, suppress warnings about unknown constituents. Defaults to ``False``.

        Returns
        -------
        Counter[str]
            Counter of how many times each element is in the molecule.

        Examples
        --------
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

        """
        name = self.name
        # Strip chemical isomer prefix (e.g. 'o-' from 'o-H2' or '#o-H2') so the
        # element parser only sees the plain formula.
        if self.prefix:
            if name and name[0] in {"#", "@"}:
                # keep the grain prefix, remove 'x-' immediately after it
                name = name[0] + name[len(self.prefix) + 2 :]
            else:
                name = name[len(self.prefix) + 1 :]

        counter = determine_constituents(name)
        mass = determine_molecular_mass(counter)
        if mass != int(self.get_mass()):
            if not quiet:
                logger.warning(
                    f"Input mass of {self.get_name()} ({self.get_mass()}) does not match calculated mass of constituents, using calculated mass: {int(mass)}"
                )
            self.set_mass(int(mass))
        self.n_atoms = counter.total()

        return counter

    def get_n_atoms(self) -> int:
        """Get the number of atoms in the molecule.

        Returns
        -------
        int
            The number of atoms

        """
        return self.n_atoms

    def set_n_atoms(self, new_n_atoms: int) -> None:
        """Set the number of atoms.

        Parameters
        ----------
        new_n_atoms : int
            The new number of atoms

        """
        self.n_atoms = new_n_atoms

    def to_UCL_format(self) -> str:
        """Serialize to the extended UCLCHEM species CSV order.

        Order: NAME,MASS,BINDING_ENERGY,SOLID_FRACTION,MONO_FRACTION,VOLCANO_FRACTION,ENTHALPY,
               DESORPTION_PREF,DIFFUSION_BARRIER,DIFFUSION_PREF,Ix,Iy,Iz,SYMMETRYFACTOR

        Returns
        -------
        str
            String with species values shown in format shown above.

        """
        return f"{self.get_name()},{self.get_mass()},{self.get_binding_energy()},{self.get_solid_fraction()},{self.get_mono_fraction()},{self.get_volcano_fraction()},{self.get_enthalpy()},{self.get_vdes()},{self.get_diffusion_barrier()},{self.get_vdiff()},{self.get_Ix()},{self.get_Iy()},{self.get_Iz()},{self.get_symmetry_factor()}"

    def get_binding_energy(self) -> float:
        """Get the binding energy of the species in K.

        Returns
        -------
        float
            The binding energy in K

        """
        return self.binding_energy

    def get_vdes(self) -> float:
        """Get the desorption prefactor.

        Returns
        -------
        float
            The desorption prefactor in s-1

        """
        return float(self.vdes)

    def get_desorption_pref(self) -> float:
        """Get the desorption prefactor.

        Alias getter matching CSV column name `desorption_pref`.

        Returns
        -------
        float
            The desorption prefactor in s-1

        """
        return self.get_vdes()

    def get_diffusion_barrier(self) -> float:
        """Get the diffusion barrier for the species.

        Returns
        -------
        float
            The diffusion barrier in K

        """
        return self.diffusion_barrier

    def get_vdiff(self) -> float:
        """Get the diffusion prefactor.

        Returns
        -------
        float
            The diffusion prefactor in s-1

        """
        return float(self.vdiff)

    def get_diffusion_pref(self) -> float:
        """Set the diffusion prefactor.

        Alias getter matching CSV column name `diffusion_pref`.

        Returns
        -------
        float
            The diffusion prefactor in s-1

        """
        return self.get_vdiff()

    def get_Ix(self) -> float:
        """Set the moment of inertia along the first principal axis.

        Returns
        -------
        float
            moment of inertia in amu/Angstrom^2

        """
        return self.Ix

    def get_Iy(self) -> float:
        """Set the moment of inertia along the second principal axis.

        Returns
        -------
        float
            moment of inertia in amu/Angstrom^2

        """
        return self.Iy

    def get_Iz(self) -> float:
        """Set the moment of inertia along the third principal axis.

        Returns
        -------
        float
            moment of inertia in amu/Angstrom^2

        """
        return self.Iz

    def get_symmetry_factor(self) -> int:
        """Get the symmetry factor of the species.

        Returns
        -------
        int
            Symmetry factor

        """
        return self.symmetry_factor

    def set_vdes(self, vdes: float) -> None:
        """Set the desorption prefactor.

        Parameters
        ----------
        vdes : float
            The desorption prefactor in s-1

        """
        self.vdes = float(vdes)

    def set_desorption_pref(self, vdes: float) -> None:
        """Set the desorption prefactor.

        Alias setter matching CSV column name `desorption_pref`.

        Parameters
        ----------
        vdes : float
            The desorption prefactor in s-1

        """
        self.set_vdes(vdes)

    def set_diffusion_barrier(self, barrier: float) -> None:
        """Set the diffusion barrier for the species.

        Parameters
        ----------
        barrier : float
            Diffusion barrier in K

        """
        self.diffusion_barrier = float(barrier)

    def set_vdiff(self, vdiff: float) -> None:
        """Set the diffusion prefactor (vdiff) for the species.

        Parameters
        ----------
        vdiff : float
            The diffusion prefactor in s-1

        """
        self.vdiff = float(vdiff)

    def set_diffusion_pref(self, vdiff: float) -> None:
        """Set the diffusion prefactor.

        Alias setter matching CSV column name `diffusion_pref`.

        Parameters
        ----------
        vdiff : float
            The diffusion prefactor in s-1

        """
        self.set_vdiff(vdiff)

    def set_Ix(self, Ix: float) -> None:
        """Set the moment of inertia along the first principal axis.

        Parameters
        ----------
        Ix : float
            desired moment of inertia (in amu/Angstrom^2)

        """
        self.Ix = float(Ix)

    def set_Iy(self, Iy: float) -> None:
        """Set the moment of inertia along the second principal axis.

        Parameters
        ----------
        Iy : float
            desired moment of inertia (in amu/Angstrom^2)

        """
        self.Iy = float(Iy)

    def set_Iz(self, Iz: float) -> None:
        """Set the moment of inertia along the third principal axis.

        Parameters
        ----------
        Iz : float
            desired moment of inertia (in amu/Angstrom^2)

        """
        self.Iz = float(Iz)

    def set_symmetry_factor(self, sym: int | str) -> None:
        """Set the symmetry factor of the species.

        Sets the symmetry factor to -1 if `sym` cannot be turned into an
        integer.

        Parameters
        ----------
        sym : int | str
            Symmetry factor

        """
        try:
            self.symmetry_factor = int(sym)
        except (ValueError, TypeError):
            self.symmetry_factor = MISSING_VALUE_INTEGER

    def __eq__(self, other: object) -> bool:
        """Check equality based on species name.

        Parameters
        ----------
        other : object
            Another species or species name string.

        Returns
        -------
        bool
            True if two species are identical.

        """
        if isinstance(other, Species):
            return self.get_name() == other.get_name()
        if isinstance(other, str):
            return self.get_name() == other
        return NotImplemented

    def __hash__(self) -> int:
        """Hash based on species name, consistent with __eq__.

        Returns
        -------
        int
            Hash of the species name.

        """
        return hash(self.get_name())

    def __lt__(self, other: Species) -> bool:
        """Compare the mass of the species.

        Parameters
        ----------
        other : Species
            Another species instance

        Returns
        -------
        bool
            True if less than the other species

        """
        return self.get_mass() < other.get_mass()

    def __gt__(self, other: Species) -> bool:
        """Compare the mass of the species.

        Parameters
        ----------
        other : Species
            Another species instance

        Returns
        -------
        bool
            True if larger than than the other species

        """
        return self.get_mass() > other.get_mass()

    def __repr__(self) -> str:
        """Return a detailed string representation of the species.

        Returns
        -------
        str
            Detailed species string.

        """
        return f"Specie: {self.get_name()}"

    def __str__(self) -> str:
        """Return the species name as a string.

        Returns
        -------
        str
            The species name.

        """
        return self.get_name()

    def calculate_rotational_partition_factor(self) -> float:
        """Calculate 1/sigma*(SQRT(IxIyIz)) for non-linear molecules, and.

        1/sigma*(SQRT(IyIz)) for linear molecules.

        Returns -999.0 if molecular inertia data is not available (backward compatibility).
        This signals that TST prefactors cannot be used for this species.

        Returns
        -------
        float
            Rotational partition factor scaled by 1e50, or -999.0 if unavailable

        """
        if self.n_atoms == 1:
            # For atoms, this is undefined, just return a value such that
            # it is clearly an atomic species.
            return -1.0
        if MISSING_VALUE_FLOAT in {self.Ix, self.Iy, self.Iz}:
            # For species without custom input Ix, Iy and Iz, we cannot do this,
            # Return sentinel value for backward compatibility
            return MISSING_VALUE_FLOAT
        if self.symmetry_factor == MISSING_VALUE_INTEGER:
            # No symmetry factor provided
            return MISSING_VALUE_FLOAT

        # Ix, Iy and Iz are in units of amu Angstrom^2,
        # so need to convert to kg m2
        amu = 1.66053907e-27  # kg/amu
        scaling_factor = 1e50

        if not self.is_linear():
            return (
                (1.0 / self.symmetry_factor)
                * np.sqrt(self.Ix * self.Iy * self.Iz * amu**3 / 1e60)
                * scaling_factor
            )
        return (
            (1.0 / self.symmetry_factor)
            * np.sqrt(self.Iy * self.Iz * amu**2 / 1e40)
            * scaling_factor
        )

    def is_linear(self) -> bool:
        """Check if molecule is linear based on moment of inertia.

        For linear molecules, Ix = 0 (rotation axis along molecular axis has no inertia).

        Returns
        -------
        bool
            True if linear, False otherwise

        """
        if self.n_atoms == 1:
            # Atomic species are not linear (doesn't matter, filtered out anyway)
            return False
        if self.n_atoms == 2:  # ruff: ignore[magic-value-comparison]
            # Diatomic molecules are always linear
            return True
        if MISSING_VALUE_FLOAT in {self.Ix, self.Iy, self.Iz}:
            # No inertia data available
            return False
        if not self.is_ice_species():
            # Only implement for grain species
            return False
        return self.Ix == 0

    def check_symmetry_factor(self) -> None:
        """Check the symmetry factor provided by the user.

        Checks if n_atoms == 2, that if its homoatomic (e.g. H2), that
        sigma == 2, and if it is heteroatomic, (e.g. OH), sigma == 1

        """
        if self.n_atoms == 1:  # Nothing to check
            return
        if self.n_atoms > 2:  # ruff: ignore[magic-value-comparison]  # Can not correctly check everything
            return
        constituents = self.find_constituents(quiet=True)
        if len(constituents) == 2:  # ruff: ignore[magic-value-comparison] # Two constituents, i.e. two different atoms.
            if self.symmetry_factor == 1:
                return
            msg = f"For diatomic molecule consisting of two different atoms (in this case {self.name}), the symmetry factor should be 1, but was given to be {self.symmetry_factor}. Correcting to 1."
            logger.warning(msg)
            self.symmetry_factor = 1
            return
        if self.symmetry_factor == 2:  # ruff: ignore[magic-value-comparison]
            return
        msg = f"For diatomic molecule consisting of two of the same atoms (in this case {self.name}), the symmetry factor should be 2, but was given to be {self.symmetry_factor}. Correcting to 2."
        logger.warning(msg)
        self.symmetry_factor = 2
