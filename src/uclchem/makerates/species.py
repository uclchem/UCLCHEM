import logging
from collections import Counter
from warnings import warn

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
symbols = ["#", "@", "*", "+", "-", "(", ")"]

species_header = (
    "name",
    "mass",
    "binding_energy",
    "solid_fraction",
    "mono_fraction",
    "volcano_fraction",
    "enthalpy",
    "Ix",
    "Iy",
    "Iz",
    "symmetry_number",
)


def is_number(s) -> bool:
    """Try to convert input to a float, if it succeeds, return True.

    Args:
        s: Input element to check for

    Returns:
        bool: True if a number, False if not.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def sanitize_input_float(row, index, default=0.0) -> float:
    output = default
    if len(row) > index and is_number(row[index]):
        output = float(row[index])
    return output


class Species:
    """Species is a class that holds all the information about an individual species in the
    network. It also has convenience functions to check whether the species is a gas or grain
    species and to help compare between species.
    """

    def __init__(self, inputRow):
        """A class representing chemical species, it reads in rows which are formatted as follows:
        NAME,MASS,BINDING ENERGY,SOLID FRACTION,MONO FRACTION,VOLCANO FRACTION,ENTHALPY,Ix,Iy,Iz,SYMMETRY_NUMBER

        The last 4 columns (Ix, Iy, Iz, SYMMETRY_NUMBER) are optional for backward compatibility.
        If not provided, TST prefactors will not be available for that species.

        Args:
            inputRow (list): Row from species CSV file
        """
        self.name = inputRow[0].upper()
        self.mass = int(inputRow[1])
        self.binding_energy = 0.0
        self.solidFraction = 0.0
        self.monoFraction = 0.0
        self.volcFraction = 0.0
        self.enthalpy = 0.0

        if len(inputRow) > 2:
            self.is_refractory = str(inputRow[2]).lower() == "inf"
            if self.is_refractory:
                self.set_binding_energy(99.9e9)
            else:
                self.set_binding_energy(inputRow[2])
        else:
            self.is_refractory = False
            self.set_binding_energy(0.0)

        self.set_solid_fraction(sanitize_input_float(inputRow, 3, 0.0))
        self.set_mono_fraction(sanitize_input_float(inputRow, 4, 0.0))
        self.set_volcano_fraction(sanitize_input_float(inputRow, 5, 0.0))
        self.set_enthalpy(sanitize_input_float(inputRow, 6, 0.0))

        # TST prefactor support - backward compatible with old species files
        # If Ix, Iy, Iz, SYMMETRY_NUMBER columns are missing, use sentinel values
        self.Ix = sanitize_input_float(inputRow, 7, -999.0)
        self.Iy = sanitize_input_float(inputRow, 8, -999.0)
        self.Iz = sanitize_input_float(inputRow, 9, -999.0)
        try:
            if len(inputRow) > 10:
                sym_val = inputRow[10]
                if isinstance(sym_val, (int, float)):
                    self.symmetry_factor = int(sym_val) if sym_val != 0 else -1
                elif isinstance(sym_val, str) and sym_val.strip():
                    self.symmetry_factor = int(sym_val)
                else:
                    self.symmetry_factor = -1
            else:
                self.symmetry_factor = -1
        except (ValueError, IndexError, TypeError):
            self.symmetry_factor = -1

        self.set_n_atoms(0)

        # in first instance, assume species freeze/desorb unchanged
        # this is updated by `check_freeze_desorbs()` later.
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
        """Get the molecular mass of the chemical species

        Returns:
            int: The molecular mass
        """
        return self.mass

    def get_charge(self) -> int:
        """Get the charge of the chemical species in e. Positive integer indicates positive ion,
        negative indicates negative ion. Assumes species are at most charged +1 or -1.

        Returns:
            int: The charge of the species
        """
        if not self.is_ion():
            return 0
        elif "+" in self.get_name():
            return 1
        elif "-" in self.get_name():
            return -1

    def get_binding_energy(self) -> float:
        """Get the binding energy of the species in K

        Returns:
            float: The binding energy in K
        """
        return self.binding_energy

    def get_solid_fraction(self) -> float:
        """Get the solid fraction of the species

        Returns:
            float: The solid fraction
        """
        return self.solidFraction

    def get_mono_fraction(self) -> float:
        """Get the monolayer fraction of the species

        Returns:
            float: The monolayer fraction
        """
        return self.monoFraction

    def get_volcano_fraction(self) -> float:
        """Get the volcano fraction of the species

        Returns:
            float: The volcano fraction
        """
        return self.volcFraction

    def get_enthalpy(self) -> float:
        """Get the ice enthalpy of the species

        Returns:
            float: The ice enthalpy
        """
        return self.enthalpy

    def set_name(self, name: str) -> None:
        """Set the name of the chemical species.

        Args:
            name (str): The new name for the species
        """
        self.name = name.upper()

    def set_mass(self, mass: int) -> None:
        """Set the molecular mass of the chemical species

        Args:
            mass (int): The new molecular mass
        """
        self.mass = int(mass)

    def set_binding_energy(self, binding_energy: float) -> None:
        """Set the binding energy of the species in K

        Args:
            binding_energy (float): The new binding energy in K
        """
        self.binding_energy = float(binding_energy)

    def set_solid_fraction(self, solid_fraction: float) -> None:
        """Set the solid fraction of the species

        Args:
            solid_fraction (float): The new solid fraction
        """
        self.solidFraction = float(solid_fraction)

    def set_mono_fraction(self, mono_fraction: float) -> None:
        """Set the monolayer fraction of the species

        Args:
            mono_fraction (float): The new monolayer fraction
        """
        self.monoFraction = float(mono_fraction)

    def set_volcano_fraction(self, volcano_fraction: float) -> None:
        """Set the volcano fraction of the species

        Args:
            volcano_fraction (float): The new volcano fraction
        """
        self.volcFraction = float(volcano_fraction)

    def set_enthalpy(self, enthalpy: float) -> None:
        """Set the enthalpy of the species in kcal per

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

    def get_desorb_products(self) -> list[str]:
        """Obtain the desorbtion products of ice species

        Returns:
            list[str]: The desorption products
        """
        return self.desorb_products

    def set_freeze_products(self, product_list: list[str], freeze_alpha: float) -> None:
        """Add the freeze products of the species, one species can have several freeze products.

        Args:
            product_list (list[str]): The list of freeze out products
            freeze_alpha (float): The freeze out ratio.

        It is called alpha, since it is derived from the alpha column in the UCLCHEM reaction format:
        https://github.com/uclchem/UCLCHEM/blob/08d37f8c3063f8ff8a9a7aa16d9eff0ed4f99538/Makerates/src/network.py#L160
        """

        self.freeze_products[",".join(product_list)] = freeze_alpha

    def get_freeze_products(self) -> dict[list[str], float]:
        """Obtain the product to which the species freeze out

        Returns:
            dict[str, float]: Reactions and their respective freeze out ratios.

        Yields:
            Iterator[dict[str, float]]: Iterator that returns all of the freeze out reactions with ratios
        """
        keys = self.freeze_products.keys()
        values = self.freeze_products.values()
        logging.debug(f"freeze keys: {keys}, products {values}")
        for key, value in zip(keys, values):
            yield key.split(","), value

    def get_freeze_products_list(self) -> list[list[str]]:
        """Returns all the freeze products without their ratios

        Returns:
            list[list[str]]: List of freeze products
        """
        # TODO: Write an unit test for get_freeze_product_behaviour
        return [key.split(",") for key in self.freeze_products.keys()]

    def get_freeze_alpha(self, product_list: list[str]) -> float:
        """Obtain the freeze out ratio of a species for a certain reaction

        Args:
            product_list (list[str]): For a specific reaction, get the freezeout ratio

        Returns:
            float: The freezeout ratio
        """
        return self.freeze_products[",".join(product_list)]

    def is_grain_species(self) -> bool:
        """Return whether the species is a species on the grain

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

    def is_surface_species(self) -> bool:
        """Checks if the species is on the surface

        Returns:
            bool: True if a surface species
        """
        return self.get_name().startswith("#")

    def is_bulk_species(self) -> bool:
        """Checks if the species is in the bulk

        Returns:
            bool: True if a bulk species
        """
        return self.get_name().startswith("@")

    def is_ion(self) -> bool:
        """Checks if the species is ionized, either postively or negatively.

        Returns:
            bool: True if it is an ionized
        """
        return self.get_name().endswith("+") or self.get_name().endswith("-")

    def add_default_freeze(self) -> None:
        """Adds a defalt freezeout, which is freezing out to the species itself, but with no ionization."""
        freeze = "#" + self.get_name()
        if freeze[-1] in ["+", "-"]:
            freeze = freeze[:-1]
        if self.get_name() == "E-":
            freeze = ""
        self.set_freeze_products([freeze, "NAN", "NAN", "NAN"], 1.0)

    def find_constituents(self, quiet=False):
        """Loop through the species' name and work out what its consituent
        atoms are. Then calculate mass and alert user if it doesn't match
        input mass.
        """
        speciesName = self.get_name()[:]
        i = 0
        atoms = []
        bracket = False
        bracketContent = []
        # loop over characters in species name to work out what it is made of
        while i < len(speciesName):
            # if character isn't a #,+ or - then check it otherwise move on
            if speciesName[i] not in symbols:
                if i + 1 < len(speciesName):
                    # if next two characters are (eg) 'MG' then atom is Mg not M and G
                    if speciesName[i : i + 3] in elementList:
                        j = i + 3
                    elif speciesName[i : i + 2] in elementList:
                        j = i + 2
                    # otherwise work out which element it is
                    elif speciesName[i] in elementList:
                        j = i + 1

                # if there aren't two characters left just try next one
                elif speciesName[i] in elementList:
                    j = i + 1
                # if we've found a new element check for numbers otherwise print error
                if j > i:
                    if bracket:
                        bracketContent.append(speciesName[i:j])
                    else:
                        atoms.append(speciesName[i:j])  # add element to list
                    if j < len(speciesName):
                        if is_number(speciesName[j]):
                            if int(speciesName[j]) > 1:
                                for k in range(1, int(speciesName[j])):
                                    if bracket:
                                        bracketContent.append(speciesName[i:j])
                                    else:
                                        atoms.append(speciesName[i:j])
                                i = j + 1
                            else:
                                i = j
                        else:
                            i = j
                    else:
                        i = j
                else:
                    raise ValueError(
                        f"Contains elements not in element list: {speciesName}"
                    )
                    logging.warning(speciesName[i])
                    logging.warning(
                        "\t{0} contains elements not in element list:".format(
                            speciesName
                        )
                    )
                    logging.warning(elementList)
            else:
                # if symbol is start of a bracketed part of molecule, keep track
                if speciesName[i] == "(":
                    bracket = True
                    bracketContent = []
                    i += 1
                # if it's the end then add bracket contents to list
                elif speciesName[i] == ")":
                    if is_number(speciesName[i + 1]):
                        for k in range(0, int(speciesName[i + 1])):
                            atoms.extend(bracketContent)
                        i += 2
                    else:
                        atoms.extend(bracketContent)
                        i += 1
                # otherwise move on
                else:
                    i += 1

        self.n_atoms = len(atoms)
        mass = 0
        for atom in atoms:
            mass += elementMass[elementList.index(atom)]
        if mass != int(self.get_mass()):
            if not quiet:
                logging.warning(
                    f"Input mass of {self.get_name()} ({self.get_mass()}) does not match calculated mass of constituents, using calculated mass: {int(mass)}"
                )
            self.set_mass(int(mass))
        counter = Counter()
        for element in elementList:
            counter[element] = atoms.count(element)
        return counter

    def get_n_atoms(self) -> int:
        """Obtain the number of atoms in the molecule

        Returns:
            int: The number of atoms
        """
        return self.n_atoms

    def set_n_atoms(self, new_n_atoms: int) -> None:
        """Set the number of atoms

        Args:
            new_n_atoms (int): The new number of atoms
        """
        self.n_atoms = new_n_atoms

    def to_UCL_format(self) -> str:
        return f"{self.get_name()},{self.get_mass()},{self.get_binding_energy()},{self.get_solid_fraction()},{self.get_mono_fraction()},{self.get_volcano_fraction()},{self.get_enthalpy()}"

    def __eq__(self, other):
        """Check for equality based on either a string or another Species instance.

        Args:
            other (str, Species): Another species

        Raises:
            NotImplementedError: We can only compare between species or strings of species.

        Returns:
            bool: True if two species are identical.
        """
        if isinstance(other, Species):
            return self.get_name() == other.get_name()
        elif isinstance(other, str):
            return self.get_name() == other
        else:
            raise NotImplementedError(
                "We can only compare between species or strings of species"
            )

    def __lt__(self, other) -> bool:
        """Compare the mass of the species

        Args:
            other (Species): Another species instance

        Returns:
            bool: True if less than the other species
        """
        return self.get_mass() < other.get_mass()

    def __gt__(self, other) -> bool:
        """Compare the mass of the species

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
        scalingFactor = 1e50

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
