"""UCLCHEM Reaction."""

from __future__ import annotations

import logging
from collections import Counter
from contextlib import contextmanager
from copy import deepcopy
from typing import TYPE_CHECKING

from uclchem.makerates.species import (
    Species,
    determine_constituents,
    determine_molecular_mass,
    elementList,
    elementMass,
    species_header,
    strip_prefix_from_species_name,
)
from uclchem.makerates.utils import normalize_species_name
from uclchem.utils import MISSING_VALUE_FLOAT

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

logger = logging.getLogger(__name__)

# Global flag for validation control
_skip_reaction_validation = False


@contextmanager
def skip_reaction_validation() -> Iterator[None]:
    """Context manager to temporarily disable reaction validation.

    This is useful when loading pre-validated networks where you do not necessarily
    want to check element and charge conservation.

    Yields
    ------
    None
        Control is yielded to the ``with`` block.

    Examples
    --------
    >>> with skip_reaction_validation():
    ...     reaction = Reaction(["#C2N", "#H", "LH", "#CH3CNH", "NAN", "NAN", "NAN"]+ [0] * 10)
    >>> reaction = Reaction(["#C2N", "#H", "LH", "#CH3CNH", "NAN", "NAN", "NAN"] + [0] * 10)
    Traceback (most recent call last):
    ...
    ValueError: Elements not conserved in a reaction.
    The following reaction caused this error: #C2N + #H + LH -> #CH3CNH.
    ...

    """
    global _skip_reaction_validation
    old_value = _skip_reaction_validation
    _skip_reaction_validation = True
    try:
        yield
    finally:
        _skip_reaction_validation = old_value


REACTION_TYPES = frozenset(
    [
        "PHOTON",
        "CRP",
        "CRPHOT",
        "FREEZE",
        "DESORB",
        "THERM",
        "DESOH2",
        "DESCR",
        "DEUVCR",
        "H2FORM",
        "ER",
        "ERDES",
        "LH",
        "LHDES",
        "BULKSWAP",
        "SURFSWAP",
        "IONOPOL1",
        "IONOPOL2",
        "CRS",
        "EXSOLID",
        "EXRELAX",
        "GAR",
        "TWOBODY",
        "ED",
    ]
)


LH_REACTION_TYPES = frozenset({"LH", "LHDES"})

ER_REACTION_TYPES = frozenset({"ER", "ERDES"})

TUNNELING_REACTION_TYPES = LH_REACTION_TYPES | ER_REACTION_TYPES

BULK_REACTION_TYPES = frozenset({"CRP", "CRPHOT", "PHOTON", "LH", "EXSOLID", "EXRELAX"})

reaction_header = [
    "REACTANT 1",
    "REACTANT 2",
    "REACTANT 3",
    "PRODUCT 1",
    "PRODUCT 2",
    "PRODUCT 3",
    "PRODUCT 4",
    "ALPHA",
    "BETA",
    "GAMMA",
    "T_MIN",
    "T_MAX",
    "REDUCED_MASS",
    "EXTRAPOLATE",
    "EXOTHERMICITY",
]


def _infer_reaction_type(reactants: list[str]) -> str:
    if reactants[2] in REACTION_TYPES:
        return reactants[2]
    if reactants[1] in REACTION_TYPES:
        return reactants[1]
    return "TWOBODY"


def check_duplicate_reactions(
    reactions: Sequence[Reaction],
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Check for any duplicate reactions in a list of reactions.

    Parameters
    ----------
    reactions : Sequence[Reaction]
        List of reactions

    Returns
    -------
    overlapping : list[tuple[int, int]]
        List of indices of duplicate reaction pairs.
    overlapping_umist : list[tuple[int, int]]
        List of indices of duplicate reaction pairs that have their source set as ``UMIST``.

    """
    overlapping: list[tuple[int, int]] = []
    overlapping_umist: list[tuple[int, int]] = []

    duplicates = [False] * len(reactions)
    for i, reaction1 in enumerate(reactions):
        for j, reaction2 in enumerate(reactions):
            if j <= i or duplicates[j]:
                continue

            if reaction1 != reaction2:
                continue

            if (reaction1.get_templow() >= reaction2.get_temphigh()) or (
                reaction1.get_temphigh() <= reaction2.get_templow()
            ):
                continue

            if reaction1.get_source() == reaction2.get_source() == "UMIST":
                overlapping_umist.append((i, j))
                continue

            overlapping.append((i, j))
            duplicates[i], duplicates[j] = True, True
    return overlapping, overlapping_umist


class Reaction:
    """Representation of reactions."""

    def __init__(
        self, input_row: list[str | float] | Reaction, reaction_source: str | None = None
    ):
        """Initialize a Reaction object.

        Parameters
        ----------
        input_row : list[str | float] | Reaction
            Either a Reaction object to copy,
            or a list/array with reaction data
        reaction_source : str | None
            Optional source identifier for the reaction.
            Default = None.

        Raises
        ------
        ValueError
            If the length of `input_row` is not long enough.

        """
        if isinstance(input_row, Reaction):
            self.set_reactants(input_row.get_reactants())
            self.set_products(input_row.get_products())
            self.set_alpha(input_row.get_alpha())
            self.set_beta(input_row.get_beta())
            self.set_gamma(input_row.get_gamma())
            self.set_templow(input_row.get_templow())
            self.set_temphigh(input_row.get_temphigh())
            self.set_reduced_mass(input_row.get_reduced_mass())
            self.set_extrapolation(enabled=input_row.get_extrapolation())
            self.set_exothermicity(input_row.get_exothermicity())
        else:
            try:
                self.set_reactants(
                    [
                        normalize_species_name(str(input_row[0])),
                        normalize_species_name(str(input_row[1])),
                        normalize_species_name(str(input_row[2])),
                    ]
                )
                self.set_products(
                    [
                        normalize_species_name(str(input_row[3])),
                        normalize_species_name(str(input_row[4])),
                        normalize_species_name(str(input_row[5])),
                        normalize_species_name(str(input_row[6])),
                    ]
                )

                self.set_alpha(float(input_row[7]))
                self.set_beta(float(input_row[8]))
                self.set_gamma(float(input_row[9]))
                self.set_templow(float(input_row[10]))
                self.set_temphigh(float(input_row[11]))
                if len(input_row) > 12:  # ruff: ignore[magic-value-comparison]
                    self.set_reduced_mass(float(input_row[12]))
                else:
                    self.set_reduced_mass(MISSING_VALUE_FLOAT)
                self.set_extrapolation(
                    enabled=bool(input_row[13]) if len(input_row) > 13 else False  # ruff: ignore[magic-value-comparison]
                )
                self.set_exothermicity(
                    float(input_row[14])
                    if (len(input_row) > 14 and input_row[14])  # ruff: ignore[magic-value-comparison]
                    else 0.0
                )

            except IndexError as e:
                msg = "Input for Reaction should be a list of length 12 with optional 13th entry for reduced mass and 14th for extrapolation flag."
                raise ValueError(msg) from e

        if not _skip_reaction_validation:
            self.check_element_conservation()
            self.check_charge_conservation()
            self.check_reaction_type_is_possible()

        self.duplicate = False
        self.source = reaction_source  # The source of the reaction, e.g. UMIST, KIDA or user defined
        self._reaction_type = _infer_reaction_type(self.get_reactants())

        # body_count is the number of factors of density to include in ODE
        # we drop a factor of density from both the LHS and RHS of ODES
        # So reactions with 1 body have no factors of density
        # which we manage by counting from -1
        self.body_count = -1
        for reactant in self.get_reactants():
            if (reactant not in REACTION_TYPES) and reactant != "NAN":
                self.body_count += 1
            if reactant in {"DESOH2", "FREEZE"}:
                self.body_count += 1
            if reactant in LH_REACTION_TYPES:
                self.body_count -= 1

        if (self.get_reaction_type() == "FREEZE") and (
            self.get_reactants()[0][-1] == "+"
        ):
            self.beta = 1

        if self.get_reaction_type() in TUNNELING_REACTION_TYPES and (
            self._reduced_mass in {MISSING_VALUE_FLOAT, 0.0}
        ):
            # If the reaction is tunneling based, and no reduced mass was supplied,
            # try to predict it.
            self.predict_reduced_mass()

    # Simple getters and setters for parsing the inputrow or changing parameters
    def get_reactants(self) -> list[str]:
        """Get the four reactants present in the reaction,.

        padded with NAN for nonexistent entries.

        Returns
        -------
        list[str]
            The four reactants names

        """
        return self._reactants[:]

    def get_pure_reactants(self) -> list[str]:
        """Get only the pure species, no reaction types and NAN entries.

        Returns
        -------
        list[str]
            The list of reacting species.

        """
        return [r for r in self._reactants[:] if r not in {*REACTION_TYPES, "NAN"}]

    def get_sorted_reactants(self) -> list[str]:
        """Get the four reactants present in the reaction,.

        sorted for fast comparisons.

        Returns
        -------
        list[str]
            The four sorted reactant names

        """
        return self._sorted_reactants

    def set_reactants(self, reactants: list[str]) -> None:
        """Set the four reactants present in the reaction,.

        padded with NAN for nonexistent entries.

        Parameters
        ----------
        reactants : list[str]
            The four reactants names

        """
        self._reactants = reactants
        self._reaction_type = _infer_reaction_type(reactants)
        # Store a sorted version for comparisons
        self._sorted_reactants = sorted(self._reactants)

    def get_products(self) -> list[str]:
        """Get the four products present in the reaction,.

        padded with NAN for nonexistent entries.

        Returns
        -------
        list[str]
            The four products names

        """
        return self._products[:]

    def get_pure_products(self) -> list[str]:
        """Get only the pure species that are products,.

        no reaction types and NAN entries.

        Returns
        -------
        list[str]
            The list of produced species.

        """
        return [r for r in self._products[:] if r not in {*REACTION_TYPES, "NAN"}]

    def get_sorted_products(self) -> list[str]:
        """Get the four products present in the reaction,.

        sorted for fast comparisons.

        Returns
        -------
        list[str]
            The four sorted products names

        """
        return self._sorted_products

    def set_products(self, products: list[str]) -> None:
        """Set the four products present in the reaction, padded with NAN for nonexistent entries.

        Parameters
        ----------
        products : list[str]
            The four products names

        """
        self._products = products
        # Store a sorted version for comparisons
        self._sorted_products = sorted(self._products)

    def get_alpha(self) -> float:
        """Get the alpha parameter from the Kooij-Arrhenius equation.

        Returns
        -------
        float
            the alpha parameter of the reaction

        """
        return self._alpha

    def set_alpha(self, alpha: float) -> None:
        """Set the alpha parameter from the Kooij-Arrhenius equation.

        Parameters
        ----------
        alpha : float
            the alpha parameter of the reaction

        """
        self._alpha = alpha

    def get_beta(self) -> float:
        """Get the beta parameter from the Kooij-Arrhenius equation.

        Returns
        -------
        float
            the beta parameter of the reaction

        """
        return self._beta

    def set_beta(self, beta: float) -> None:
        """Set the beta parameter from the Kooij-Arrhenius equation.

        Parameters
        ----------
        beta : float
            the beta parameter of the reaction

        """
        self._beta = beta

    def set_gamma(self, gamma: float) -> None:
        """Set the gamma parameter from the Kooij-Arrhenius equation.

        Parameters
        ----------
        gamma : float
            the gamma parameter of the reaction

        """
        self._gamma = gamma

    def get_gamma(self) -> float:
        """Get the gamma parameter from the Kooij-Arrhenius equation.

        Returns
        -------
        float
            the gamma parameter of the reaction

        """
        return self._gamma

    def set_templow(self, templow: float) -> None:
        """Set the lower temperature boundary of the reaction in Kelvin.

        Parameters
        ----------
        templow : float
            the lower temperature boundary

        """
        self._templow = templow

    def get_templow(self) -> float:
        """Get the lower temperature boundary of the reaction in Kelvin.

        Returns
        -------
        float
            the lower temperature boundary

        """
        return self._templow

    def set_temphigh(self, temphigh: float) -> None:
        """Set the higher temperature boundary of the reaction in Kelvin.

        Parameters
        ----------
        temphigh : float
            the higher temperature boundary

        """
        self._temphigh = temphigh

    def get_temphigh(self) -> float:
        """Get the higher temperature boundary of the reaction in Kelvin.

        Returns
        -------
        float
            the higher temperature boundary

        """
        return self._temphigh

    def set_exothermicity(self, rate: float) -> None:
        """Set the cooling/heating for the reaction in erg s^-1.

        Parameters
        ----------
        rate : float
            the reaction enthalpy change

        """
        self._exothermicity = rate

    def get_exothermicity(self) -> float:
        """Get the cooling/heating for the reaction in erg s^-1.

        Returns
        -------
        float
            the reaction enthalpy change

        """
        return self._exothermicity

    def predict_reduced_mass(self) -> None:
        """Predict the reduced mass of the tunneling particle in this reaction.

        This is used in the calculation of the tunneling rates.

        Raises
        ------
        RuntimeError
            If the reaction only has one reactant.

        Examples
        --------
        >>> reaction = Reaction(["#CH3OH", "#H", "LH", "#CH3O", "#H2", "NAN", "NAN"] + [0] * 10)
        >>> # Setting a custom reduced mass
        >>> reaction.set_reduced_mass(20.0)
        >>>
        >>> # The custom reduced mass that we set.
        >>> reaction.get_reduced_mass()
        20.0
        >>> # Predicting the reduced mass of the reaction
        >>> reaction.predict_reduced_mass()
        >>> reaction.get_reduced_mass()
        1.0

        >>> # It is called upon Reaction instantiation
        >>> reaction = Reaction(["#CH3OH", "#OH", "LH", "#CH3O", "#H2O", "NAN", "NAN"] + [0] * 10)
        >>> reaction.get_reduced_mass() # mass of atomic hydrogen
        1.0

        """
        reactants = self.get_pure_reactants()
        if len(reactants) == 1:
            msg = f"Tried to predict reduced mass of reaction '{self}' with only one reactant."
            raise RuntimeError(msg)

        reac_constituents = []
        reac_masses = []
        # Get all reactant species and their elemental buildup
        for reac in reactants:
            reac = strip_prefix_from_species_name(reac)
            atoms = determine_constituents(reac)
            reac_constituents.append(atoms)
            reac_masses.append(determine_molecular_mass(atoms))

        products = self.get_pure_products()
        # Get all product species and their elemental buildup
        prod_constituents = []
        prod_masses = []
        for prod in products:
            prod = strip_prefix_from_species_name(prod)
            atoms = determine_constituents(prod)
            prod_constituents.append(atoms)
            prod_masses.append(determine_molecular_mass(atoms))

        # Get mass and number of reactants and products
        naive_reduced_mass = (
            reac_masses[0] * reac_masses[1] / (reac_masses[0] + reac_masses[1])
        )
        n_reacs = len(reac_constituents)
        n_prods = len(prod_constituents)
        if n_reacs == n_prods:
            for reac_constituent in reac_constituents:
                # For each reactant, find which product is
                # closest (most similar in buildup) to it.
                min_total = int(1e10)
                min_diff = None
                for prod_constituent in prod_constituents:
                    diff = deepcopy(reac_constituent)
                    diff.subtract(prod_constituent)
                    total_change = 0
                    for element in elementList:
                        total_change += abs(diff[element])
                    if total_change < min_total:
                        min_total = total_change
                        min_diff = diff
                if min_diff is None:
                    continue
                changing_species = Counter({k: c for k, c in min_diff.items() if c != 0})

                items = changing_species.items()
                if len(items) == 1:
                    # Exchange reaction
                    tuple_items = next(iter(items))
                    if abs(tuple_items[1]) == 1:
                        # One element is switched
                        element_index = elementList.index(tuple_items[0])
                        # Set reduced mass to mass of switched element
                        reduced_mass: float = float(elementMass[element_index])
                        self.set_reduced_mass(reduced_mass)
                        logger.debug(
                            f"Predicted reduced mass of '{self}' to be {self._reduced_mass} (would have been {naive_reduced_mass})"
                        )
                        return
        elif n_reacs == 2 and n_prods == 1:  # ruff: ignore[magic-value-comparison]
            # Addition reaction
            if reactants[0].strip("#@") == reactants[1].strip("#@"):
                # If the two species are the same (e.g. #H+#H-> #H2), set reduced mass to m/2
                mass = reac_masses[0]
                reduced_mass = float(mass) / 2.0
                self.set_reduced_mass(reduced_mass)
                logger.debug(
                    f"Predicted reduced mass of '{self}' to be {self._reduced_mass} (would have been {naive_reduced_mass})"
                )
                return
            if any(species == Counter({"H": 1}) for species in reac_constituents):
                # If one of the species is #H, set reduced mass to 1
                self.set_reduced_mass(1.0)
                logger.debug(
                    f"Predicted reduced mass of '{self}' to be {self._reduced_mass} (would have been {naive_reduced_mass})"
                )
                return
        elif n_reacs == 1 and n_prods == 2:  # ruff: ignore[magic-value-comparison]
            # Splitting reaction. Not in network
            # (also not LH or ER type, so would never get here)
            pass
        msg = f"Could not predict reduced mass of '{self}'. "
        msg += f"Instead, using regular definition with masses of two reactants (mu={naive_reduced_mass:.3})."
        if self._gamma == 0:
            msg += " (Reaction is barrierless anyway)"
        logger.warning(msg)
        self.set_reduced_mass(naive_reduced_mass)

    def set_reduced_mass(self, reduced_mass: float) -> None:
        """Set the reduced mass to be used to calculate tunneling rate in.

        atomic mass units.

        Parameters
        ----------
        reduced_mass : float
            reduced mass of moving atoms

        """
        self._reduced_mass = reduced_mass

    def get_reduced_mass(self) -> float:
        """Get the reduced mass to be used to calculate tunneling rate in.

        atomic mass units.

        Returns
        -------
        float
            reduced mass of moving atoms

        """
        return self._reduced_mass

    # C

    def get_reaction_type(self) -> str:
        """Get the type of a reaction from the reactants.

        First check the third reactant for a reaction type, then the second. If there are none
        in there, it will be regarded as a two body reaction.

        Returns
        -------
        str
            reaction type

        """
        return self._reaction_type

    def get_source(self) -> str | None:
        """Get the source of the reaction.

        Returns
        -------
        str | None
            The source of the reaction

        """
        return self.source

    def set_source(self, source: str) -> None:
        """Set the source of the reaction.

        Parameters
        ----------
        source : str
            The source of the reaction

        """
        self.source = source

    def set_extrapolation(self, *, enabled: bool = True) -> None:
        """Set whether extrapolation is applied for this reaction.

        Parameters
        ----------
        enabled : bool
            whether extrapolation is applied. Default=True.

        """
        logger.info(f"Setting for {self} extrapolation to {enabled}")
        self.extrapolate = enabled

    def get_extrapolation(self) -> bool:
        """Get whether extrapolation is applied for this reaction.

        Returns
        -------
        bool
            whether extrapolation is applied.

        """
        return self.extrapolate

    def convert_surf_to_bulk(self) -> None:
        """Convert the surface species to bulk species in place for this reaction."""
        self.set_reactants([reac.replace("#", "@") for reac in self.get_reactants()])
        self.set_products([prod.replace("#", "@") for prod in self.get_products()])

    def check_element_conservation(self) -> None:
        """Check the conservation of elements.

        Raises
        ------
        ValueError
            If the elements are not conserved by the reaction.

        """
        if self.get_reaction_type() in {"FREEZE", "DESORB"}:
            return

        counter_reactants: Counter[str] = Counter()
        for reac in self._reactants:
            if reac in {*REACTION_TYPES, "NAN", "E-"}:
                continue
            atoms_counter_specie = determine_constituents(
                strip_prefix_from_species_name(reac)
            )
            counter_reactants += atoms_counter_specie

        counter_products: Counter[str] = Counter()
        for prod in self._products:
            if prod in {*REACTION_TYPES, "NAN", "E-"}:
                continue
            atoms_counter_specie = determine_constituents(
                strip_prefix_from_species_name(prod)
            )
            counter_products += atoms_counter_specie

        if counter_products != counter_reactants:
            msg = "Elements not conserved in a reaction.\n"
            msg += f"The following reaction caused this error: {self}.\n"
            msg += f"Reactants: {counter_reactants}. Products: {counter_products}"
            raise ValueError(msg)

    def check_charge_conservation(self) -> None:
        """Check that the charge is conserved by this reaction.

        Grain reactions don't need to conserve charge, because grains can
        absorb/release electrons, so they are ignored.

        Raises
        ------
        ValueError
            If charge is not conserved by the reaction.

        """
        # Grain reactions don't need to conserve charge (grains can absorb/release electrons)
        if self.is_ice_reaction(
            include_reactants=True, include_products=True, strict=False
        ):
            return
        charge_reactants = 0
        for reac in self._reactants:
            if reac == "NAN":
                continue
            specie = Species([reac] + [0] * len(species_header))
            charge_reactants += specie.get_charge()
        charge_products = 0
        for prod in self._products:
            if prod == "NAN":
                continue
            specie = Species([prod] + [0] * len(species_header))
            charge_products += specie.get_charge()

        if charge_products != charge_reactants:
            # GAR and FREEZE reactions do not conserve charge
            if self.get_reaction_type() in {"GAR", "FREEZE"}:
                # GAR reactions rely on grain electrons
                # FREEZE reactions can have electron freeze-out (E- -> nothing with alpha=0)
                return
            msg = "Charges not conserved in a reaction.\n"
            msg += f"The following reaction caused this error: {self}.\n"
            msg += f"Reactants: {charge_reactants}. Products: {charge_products}"
            raise ValueError(msg)

    def check_reaction_type_is_possible(self) -> None:
        """Check that the reaction type is valid for the reactants and products.

        Check that the combination of the number of species (reactant and product)
        and their phases, matches with the given reaction type.

        Raises
        ------
        ValueError
            - If a TWOBODY reaction has a species on the ice
            - If an LH reaction has a species in the gas-phase
            - If an ER reaction does not have one reactant on the ice and one in the gas-phase
            - If a THERM reaction does not have all reactants on the ice
                and all products in the gas-phase
            - If a FREEZE reaction does not have all reactants in the gas-phase
                and all products on the ice
            - If an LH, LHDES, ER, ERDES or TWOBODY reaction does not have two reagents
            - If a THERM, FREEZE, DESCR, DESOH2, DEUVCR, SURFSWAP or BULKSWAP reaction
                does not have only one reagent.

        """
        reaction_type = self.get_reaction_type()
        if reaction_type == "TWOBODY":
            if self.is_ice_reaction(
                include_reactants=True, include_products=True, strict=False
            ):
                msg = f"TWOBODY reactions must happen in the gas-phase, but reaction '{self}' has reactants or products on the ice."
                raise ValueError(msg)
        elif reaction_type == "LH":
            if not self.is_ice_reaction(
                include_reactants=True, include_products=True, strict=True
            ):
                msg = f"LH reactions must happen fully on the ice, but '{self}' has reactants or products in the gas."
                raise ValueError(msg)
        elif reaction_type == "ER":
            if not self.is_ice_reaction(
                include_reactants=True, include_products=False, strict=False
            ) or not self.is_gas_reaction(
                include_reactants=True, include_products=False, strict=False
            ):
                msg = f"ER reactions must have one reactant on ice, and one in the gas, but '{self}' did not."
                raise ValueError(msg)
            if not self.is_ice_reaction(
                include_reactants=False, include_products=True, strict=True
            ):
                msg = f"ER reactions must have all products on ice, but '{self}' did not."
                raise ValueError(msg)
        elif reaction_type == "THERM":
            if not self.is_ice_reaction(
                include_reactants=True, include_products=False, strict=True
            ) or not self.is_gas_reaction(
                include_reactants=False, include_products=True, strict=True
            ):
                msg = f"THERM reactions must have all reactants on ice and products in gas, but '{self}' did not."
                raise ValueError(msg)
        elif reaction_type == "FREEZE" and self.get_pure_reactants() != [
            "E-"
        ]:  # Skip electron freeze-out because it is absorbed by grain
            if not self.is_gas_reaction(
                include_reactants=True, include_products=False, strict=True
            ) or not self.is_ice_reaction(
                include_reactants=False, include_products=True, strict=True
            ):
                msg = f"FREEZE reactions must have all reactants in gas and products on ice, but '{self}' did not."
                raise ValueError(msg)

        if reaction_type in {"LH", "LHDES", "ER", "ERDES", "TWOBODY"}:
            if len(self.get_pure_reactants()) != 2:  # ruff: ignore[magic-value-comparison]
                msg = f"Reactions with type '{reaction_type}' should have two reactants, but reaction '{self}' had {len(self.get_pure_reactants())}"
                raise ValueError(msg)
        elif reaction_type in {
            "THERM",
            "FREEZE",
            "DESCR",
            "DESOH2",
            "DEUVCR",
            "SURFSWAP",
            "BULKSWAP",
        }:
            if len(self.get_pure_reactants()) != 1:
                msg = f"Reactions with type '{reaction_type}' should have only one reactant, but reaction '{self}' had {len(self.get_pure_reactants())}"
                raise ValueError(msg)

    def convert_gas_to_surf(self) -> None:
        """Convert the gas-phase species to surface species in place for this reaction.

        If any ions are produced, the ion is assumed to become neutral because it is
        on the surface. If any electrons are produced, they are assumed to be
        absorbed by the grain.

        """
        do_not_convert = [*REACTION_TYPES, "E-", "NAN"]
        self.set_reactants(
            [
                "#" + reac if reac not in do_not_convert else reac
                for reac in self.get_reactants()
            ]
        )
        self.set_products(
            [
                "#" + prod.replace("+", "")
                if prod not in do_not_convert
                else prod.replace("E-", "NAN")
                for prod in self.get_products()
            ]
        )

    def __eq__(self, other: object) -> bool:
        """Check equality based on sorted reactants and products.

        Parameters
        ----------
        other : object
            Another reaction.

        Returns
        -------
        bool
            equality

        """
        if not isinstance(other, Reaction):
            return NotImplemented
        if self._reaction_type != other._reaction_type:
            return False
        return (
            self.get_sorted_reactants() == other.get_sorted_reactants()
            and self.get_sorted_products() == other.get_sorted_products()
        )

    def check_temperature_collision(self, other: Reaction) -> bool:
        """Check if two reactions have overlapping temperature ranges.

        Returning True means there is a collision.

        Parameters
        ----------
        other : Reaction
            Another reaction

        Returns
        -------
        bool
            Whether there is a collision (True), or not (False)

        Raises
        ------
        NotImplementedError
            If `other` is not a `Reaction` instance.
            Currently we can only compare against instantiated Reaction objects.

        """
        if not isinstance(other, Reaction):
            msg = "Equality is not implemented for anything but comparing to other reactions."
            raise NotImplementedError(msg)
        if (other.get_templow() > self.get_templow()) and (
            other.get_templow() < self.get_temphigh()
        ):
            return True
        return (other.get_temphigh() > self.get_templow()) and (
            other.get_temphigh() < self.get_temphigh()
        )

    def changes_surface_count(self) -> bool:
        """Check whether a grain reaction changes number of particles on the surface.

        2 reactants to 2 products won't but two reactants combining to one will.

        Returns
        -------
        bool
            whether the number of ice molecules changes by this reaction.

        """
        if len([x for x in self.get_reactants() if "#" in x]) != len(
            [x for x in self.get_products() if "#" in x]
        ):
            return True
        return len([x for x in self.get_reactants() if "@" in x]) != len(
            [x for x in self.get_products() if "@" in x]
        )

    def changes_total_mantle(self) -> bool:
        """Check if the total grains on the mantle are changed by the reaction.

        Returns
        -------
        bool
            Whether the total ice abundance is affected by this reaction.

        """
        # If it's not just a movement between ice phases
        if ("BULK" not in self.get_reactants()[1]) and (
            "SWAP" not in self.get_reactants()[1]
        ):
            # if the number of ice species changes
            return self.changes_surface_count()
        return False

    def generate_ode_bit(self, i: int, species_names: list[str]) -> None:
        """Generate the ODE string of this reaction.

        Parameters
        ----------
        i : int
            index of reaction in network in python format (counting from 0)
        species_names : list[str]
            List of species names so we can find index
            of reactants in species list

        """
        self.ode_bit = _generate_reaction_ode_bit(
            i,
            species_names,
            self.body_count,
            self.get_reactants(),
            self.get_reaction_type(),
        )

    def to_UCL_format(self) -> str:
        """Convert a reaction to UCLCHEM reaction file format.

        Returns
        -------
        str
            string representing species

        """
        reactants = self.get_reactants()
        joined_reactants = ",".join(
            [reactant if reactant != "NAN" else "" for reactant in reactants]
        )

        products = self.get_products()
        joined_products = ",".join(
            [product if product != "NAN" else "" for product in products]
        )
        reactants_products = joined_reactants + "," + joined_products
        alpha, beta, gamma = (
            self.get_alpha(),
            self.get_beta(),
            self.get_gamma(),
        )
        str_alpha, str_beta, str_gamma = (
            str(alpha).replace("e", "E"),
            str(beta).replace("e", "E"),
            str(gamma).replace("e", "E"),
        )
        if alpha == 0:
            str_alpha = "0"
        if beta == 0:
            str_beta = "0"
        if gamma == 0:
            str_gamma = "0"
        reaction_parameters = f"{str_alpha},{str_beta},{str_gamma}"
        formatted_reaction = reactants_products + "," + reaction_parameters + ",,,,,"
        formatted_reaction += str(int(self.get_extrapolation()))
        return formatted_reaction

    def _is_reaction_wrap(
        self, *, include_reactants: bool = True, include_products: bool = True
    ) -> list[str]:
        if not (include_reactants or include_products):
            msg = "Either include reactants or products"
            raise AssertionError(msg)
        species_to_check = []
        if include_reactants:
            species_to_check += self.get_pure_reactants()
        if include_products:
            species_to_check += self.get_pure_products()
        return species_to_check

    def is_gas_reaction(
        self,
        *,
        include_reactants: bool = True,
        include_products: bool = True,
        strict: bool = True,
    ) -> bool:
        """Check whether it is a gas reaction.

        By default it is strict, meaning that all reactants must be in the gas-phase.
        If strict=False; any reaction in the gas-phase returns true.

        Parameters
        ----------
        include_reactants : bool
            Include the reactants. Defaults to True.
        include_products : bool
            Include the products. Defaults to True.
        strict : bool
            Choose between all (true) or any (false) must be gas phase.
            Defaults to True.

        Returns
        -------
        bool
            Is it a gas phase reaction?

        """
        checklist = [
            not s.startswith(("#", "@"))
            for s in self._is_reaction_wrap(
                include_reactants=include_reactants, include_products=include_products
            )
        ]
        return all(checklist) if strict else any(checklist)

    def is_ice_reaction(
        self,
        *,
        include_reactants: bool = True,
        include_products: bool = True,
        strict: bool = True,
    ) -> bool:
        """Check whether it is an ice (surface OR bulk) reaction.

        By default it is strict (strict=True); all species must be in the ice phase
        If strict=False; any species in ice phase returns True

        Parameters
        ----------
        include_reactants : bool
            Include the reactants. Defaults to True.
        include_products : bool
            Include the products. Defaults to True.
        strict : bool
            Choose between all (true) or any (false) must be ice phase.
            Defaults to True.

        Returns
        -------
        bool
            Is it an ice phase reaction?

        """
        checklist = [
            s.startswith(("#", "@"))
            for s in self._is_reaction_wrap(
                include_reactants=include_reactants, include_products=include_products
            )
        ]
        return all(checklist) if strict else any(checklist)

    def is_surface_reaction(
        self,
        *,
        include_reactants: bool = True,
        include_products: bool = True,
        strict: bool = False,
    ) -> bool:
        """Check whether it is a surface reaction. Defaults to non-strict since many.

        important surface reactions can lead to desorption in some way.

        By default it is NOT strict (strict=False); any species on the surface returns true
        If strict=True; all species must be on the ice phase

        Parameters
        ----------
        include_reactants : bool
            Include the reactants. Defaults to True.
        include_products : bool
            Include the products. Defaults to True.
        strict : bool
            Choose between all (true) or any (false) must be on the surface.
            Defaults to False.

        Returns
        -------
        bool
            Is it a surface reaction?

        """
        checklist = [
            s.startswith("#")
            for s in self._is_reaction_wrap(
                include_reactants=include_reactants, include_products=include_products
            )
        ]
        return all(checklist) if strict else any(checklist)

    def is_bulk_reaction(
        self,
        *,
        include_reactants: bool = True,
        include_products: bool = True,
        strict: bool = False,
    ) -> bool:
        """Check whether it is a bulk reaction. Defaults to non-strict since many.

        important bulk reactions interact with the surface.

        By default it is NOT strict (strict=False); any species in the bulk returns true
        If strict=True; all species must be on the ice phase

        Parameters
        ----------
        include_reactants : bool
            Include the reactants. Defaults to True.
        include_products : bool
            Include the products. Defaults to True.
        strict : bool
            Choose between all (true) or any (false) must in the bulk.
            Defaults to False.

        Returns
        -------
        bool
            Is it a bulk reaction?

        """
        checklist = [
            s.startswith("@")
            for s in self._is_reaction_wrap(
                include_reactants=include_reactants, include_products=include_products
            )
        ]
        return all(checklist) if strict else any(checklist)

    def has_unknown_species(self, species: Iterable[str]) -> bool:
        """Determine whether this reaction involves any species not in `species`.

        Parameters
        ----------
        species : Iterable[str]
            Iterable of available species

        Returns
        -------
        bool
            True if there is a reactant or product that is not in `species`,
            False otherwise

        """
        for spec in self.get_pure_reactants() + self.get_pure_products():
            if spec not in species:
                return True
        return False

    def __str__(self) -> str:
        """Return a human-readable string of the reaction.

        Returns
        -------
        str
            Human-readable reaction string.

        """
        return (
            " + ".join(filter(lambda r: r != "NAN", self.get_reactants()))
            + " -> "
            + " + ".join(filter(lambda p: p != "NAN", self.get_products()))
        )

    def __repr__(self) -> str:
        """Return a detailed string representation of the reaction.

        Returns
        -------
        str
            Detailed reaction string including reaction type.

        """
        return (
            self.get_reaction_type()
            + " reaction: "
            + " + ".join(
                filter(
                    lambda r: (r != "NAN") and (r not in REACTION_TYPES),
                    self.get_reactants(),
                )
            )
            + " -> "
            + " + ".join(filter(lambda p: p != "NAN", self.get_products()))
        )

    def __hash__(self) -> int:
        """Return hash based on reaction parameters and species.

        Returns
        -------
        int
            Hash value for this reaction.

        """
        return hash(
            f"{self.get_alpha(), self.get_beta(), self.get_gamma(), self.get_reactants(), self.get_products(), self.get_templow(), self.get_temphigh()}"
        )


def _get_original_partner(reaction: Reaction | CoupledReaction) -> Reaction:
    while isinstance(reaction, CoupledReaction):
        # If the current loop reaction is also coupled, get its partner.
        partner = reaction.get_partner()
        if partner is None:
            msg = f"CoupledReaction '{reaction}' has no partner"
            raise RuntimeError(msg)
        reaction = partner
    return reaction


class CoupledReaction(Reaction):
    """Representation of reactions that are coupled to another Reaction instance.

    This means that if a reaction has a parameter changed by, for example,
    `network.change_energy_barrier()`, every CoupledReaction that has that instance
    as its partner also has its binding energy changed to that value.

    """

    def __init__(self, input_row: list[str | float] | Reaction):
        """Initialize the CoupledReaction.

        Parameters
        ----------
        input_row : list[str | float] | Reaction
            Either a Reaction object to copy, or a list with reaction data.

        """
        super().__init__(input_row)
        self.partner: Reaction | None = None

    def set_partner(self, partner: Reaction | CoupledReaction) -> None:
        """Set the partner.

        Parameters
        ----------
        partner : Reaction | CoupledReaction
            partner of this reaction. If a CoupledReaction, will walk down the tree of partners
            until it finds an uncoupled reaction (i.e. just a :class:`Reaction`).

        Raises
        ------
        TypeError
            If `partner` is not an instance of :class:`Reaction`.

        """
        if not isinstance(partner, Reaction):
            msg = f"partner should be of type Reaction, but got type {type(partner)}"
            raise TypeError(msg)
        self.partner = _get_original_partner(partner)

    def get_partner(self) -> Reaction | None:
        """Get the partner.

        Returns
        -------
        Reaction | None
            partner of this reaction.

        """
        return self.partner


def _generate_reaction_ode_bit(
    i: int,
    species_names: list[str],
    body_count: int,
    reactants: list[str],
    reaction_type: str | None = None,
) -> str:
    """Every reaction contributes a fixed rate of change to whatever species it.

    affects. We create the string of fortran code describing that change here.

    Parameters
    ----------
    i : int
        index of reaction in network in python format (counting from 0)
    species_names : list[str]
        List of species names so we can find index
        of reactants in species list
    body_count : int
        Number of bodies in the reaction.
        Used to determine how many factors of density to include
    reactants : list[str]
        The reactants of the reaction
    reaction_type : str | None
        reaction type. Only important if
        it is GAR or not. Defaults to ``None``.

    Returns
    -------
    str
        String fragment of the ODE term for this reaction.

    """
    ode_bit = f"+RATE_CONSTANTS({i + 1})"
    # every body after the first requires a factor of density
    for _ in range(body_count):
        ode_bit += "*D"

    # GAR needs an additional factor of density:
    if reaction_type == "GAR":
        ode_bit += "*D"

    # then bring in factors of abundances
    for species in reactants:
        if species in species_names:
            ode_bit += f"*Y({species_names.index(species) + 1})"
        elif species == "BULKSWAP":
            ode_bit += "*ratioSurfaceToBulk"
        elif species == "SURFSWAP":
            ode_bit += "*totalSwap/safeMantle"
        elif species in {"DEUVCR", "DESCR", "DESOH2", "ER", "ERDES"}:
            ode_bit += "/safeMantle"
            if species == "DESOH2":
                ode_bit += f"*Y({species_names.index('H') + 1})"
        elif species == "ED":
            ode_bit += f"*Y({species_names.index('#H2') + 1})"

        if "H2FORM" in reactants:
            # only 1 factor of H abundance in Cazaux & Tielens 2004 H2 formation
            # so stop looping after first iteration
            break

    if "LH" in reactants[2] and "@" in reactants[0]:
        ode_bit += "*bulkLayersReciprocal"
    return ode_bit
