import logging
from copy import deepcopy
from typing import Union

reaction_types = [
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
    "TWOBODY"
]

tunneling_reaction_types = [
    "LH",
    "LHDES",
    "ER",
    "ERDES",
]

from collections import Counter
from copy import deepcopy

from uclchem.makerates.species import Species, elementList, elementMass, species_header


class Reaction:
    def __init__(self, inputRow, reaction_source=None):
        if isinstance(inputRow, Reaction):
            self.set_reactants(inputRow.get_reactants())
            self.set_products(inputRow.get_products())
            self.check_element_conservation()
            self.check_charge_conservation()
            self.set_alpha(inputRow.get_alpha())
            self.set_beta(inputRow.get_beta())
            self.set_gamma(inputRow.get_gamma())
            self.set_templow(inputRow.get_templow())
            self.set_temphigh(inputRow.get_temphigh())
            self.set_reduced_mass(inputRow.get_reduced_mass())
            self.set_extrapolation(inputRow.get_extrapolation())
            self.set_exothermicity(inputRow.get_exothermicity())
        else:
            try:
                self.set_reactants(
                    [
                        str(inputRow[0]).upper(),
                        str(inputRow[1]).upper(),
                        self.NANCheck(str(inputRow[2])).upper(),
                    ]
                )
                self.set_products(
                    [
                        self.NANCheck(str(inputRow[3])).upper(),
                        self.NANCheck(str(inputRow[4])).upper(),
                        self.NANCheck(str(inputRow[5])).upper(),
                        self.NANCheck(str(inputRow[6])).upper(),
                    ]
                )
                self.check_element_conservation()
                self.check_charge_conservation()
                
                self.set_alpha(float(inputRow[7]))
                self.set_beta(float(inputRow[8]))
                self.set_gamma(float(inputRow[9]))
                self.set_templow(float(inputRow[10]))
                self.set_temphigh(float(inputRow[11]))
                if len(inputRow) > 12:
                    self.set_reduced_mass(float(inputRow[12]))
                self.set_extrapolation(bool(inputRow[13]) if len(inputRow) > 13 else False)
                self.set_exothermicity(float(inputRow[14]) if len(inputRow) > 14 else 0.0)
                
            except IndexError as error:
                raise ValueError(
                    "Input for Reaction should be a list of length 12 with optional 13th entry for reduced mass and 14th for extrapolation flag."
                ) from error
        self.duplicate = False
        self.source = reaction_source  # The source of the reaction, e.g. UMIST, KIDA or user defined

        # body_count is the number of factors of density to include in ODE
        # we drop a factor of density from both the LHS and RHS of ODES
        # So reactions with 1 body have no factors of density which we manage by counting from -1
        self.body_count = -1
        for reactant in self.get_reactants():
            if (reactant not in reaction_types) and reactant != "NAN":
                self.body_count += 1
            if reactant in ["DESOH2", "FREEZE"]:
                self.body_count += 1
            if reactant in ["LH", "LHDES"]:
                self.body_count -= 1

        if (self.get_reaction_type() == "FREEZE") and (
            self.get_reactants()[0][-1] == "+"
        ):
            self.beta = 1

        if (
            self.get_reaction_type() in tunneling_reaction_types
            and self._reduced_mass == 0.0
        ):
            # If the reaction is tunneling based, and no reduced mass was supplied, try to predict it.
            self.predict_reduced_mass()

    # Simple getters and setters for parsing the inputrow or changing parameters
    def get_reactants(self) -> list[str]:
        """Get the four reactants present in the reaction, padded with NAN for nonexistent

        Returns:
            list[str]: The four reactants names
        """
        return self._reactants[:]

    def get_pure_reactants(self) -> list[str]:
        """Get only the pure species, no reaction types and NAN entries

        Returns:
            list[str]: The list of reacting species.
        """
        return [
            r
            for r in self._reactants[:]
            if r
            not in reaction_types
            + [
                "NAN",
            ]
        ]

    def get_sorted_reactants(self) -> list[str]:
        """Get the four reactants present in the reaction, sorted for fast comparisons

        Args:
            reactants (list[str]): The four sorted reactant names
        """
        return self._sorted_reactants

    def set_reactants(self, reactants: list[str]) -> None:
        """Set the four reactants present in the reaction, padded with NAN for nonexistent

        Args:
            reactants (list[str]): The four reactants names
        """
        self._reactants = reactants
        # Store a sorted version for comparisons
        self._sorted_reactants = sorted(self._reactants)
        
    def get_products(self) -> list[str]:
        """Get the four products present in the reaction, padded with NAN for nonexistent

        Args:
            reactants (list[str]): The four products names
        """
        return self._products[:]

    def get_pure_products(self) -> list[str]:
        """Get only the pure species that are products, no reaction types and NAN entries

        Returns:
            list[str]: The list of produced species.
        """
        return [
            r
            for r in self._products[:]
            if r
            not in reaction_types
            + [
                "NAN",
            ]
        ]

    def get_sorted_products(self) -> list[str]:
        """Get the four products present in the reaction, sorted for fast comparisons

        Args:
            products (list[str]): The four sorted products names
        """
        return self._sorted_products

    def set_products(self, products: list[str]) -> None:
        """Set the four products present in the reaction, padded with NAN for nonexistent

        Args:
            products (list[str]): The four products names
        """
        self._products = products
        # Store a sorted version for comparisons
        self._sorted_products = sorted(self._products)

    def get_alpha(self) -> float:
        """Get the alpha parameter from the Kooij-Arrhenius equation

        Returns:
            float: the alpha parameter of the reaction
        """
        return self._alpha

    def set_alpha(self, alpha: float) -> None:
        """Set the alpha parameter from the Kooij-Arrhenius equation

        Args:
            alpha (float): the alpha parameter of the reaction
        """
        self._alpha = alpha

    def get_beta(self) -> float:
        """Get the beta parameter from the Kooij-Arrhenius equation

        Returns:
            float: the beta parameter of the reaction
        """
        return self._beta

    def set_beta(self, beta: float) -> None:
        """Set the beta parameter from the Kooij-Arrhenius equation

        Args:
            beta (float): the beta parameter of the reaction
        """
        self._beta = beta

    def set_gamma(self, gamma: float) -> None:
        """Set the gamma parameter from the Kooij-Arrhenius equation

        Args:
            gamma (float): the gamma parameter of the reaction
        """
        self._gamma = gamma

    def get_gamma(self) -> float:
        """Get the gamma  parameter from the Kooij-Arrhenius equation

        Returns:
            float: the gamma parameter of the reaction
        """
        return self._gamma

    def set_templow(self, templow: float) -> None:
        """Set the lower temperature boundary of the reaction in Kelvin

        Args:
            templow (float): the lower temperature boundary
        """
        self._templow = templow

    def get_templow(self) -> float:
        """Get the lower temperature boundary of the reaction in Kelvin

        Returns:
            float: the lower temperature boundary
        """
        return self._templow

    def set_temphigh(self, temphigh: float) -> None:
        """Set the higher temperature boundary of the reaction in Kelvin

        Args:
            templow (float): the higher temperature boundary
        """
        self._temphigh = temphigh

    def get_temphigh(self) -> float:
        """Get the higher temperature boundary of the reaction in Kelvin

        Returns:
            float: the higher temperature boundary
        """
        return self._temphigh
    
    def set_exothermicity(self, rate: float) -> None:
        """Set the cooling/heating for the reaction in erg s^-1

        Args:
            delta_h (float): the reaction enthalpy change
        """
        self._exothermicity = rate
    
    def get_exothermicity(self) -> float:
        """Get the cooling/heating for the reaction in erg s^-1

        Returns:
            float: the reaction enthalpy change
        """
        return self._exothermicity       

    def predict_reduced_mass(self) -> None:
        """Predict the reduced mass of the tunneling particle in this reaction.
        This is used in the calculation of the tunneling rates.
        """
        reac_constituents = []
        reac_species = []
        # Get all reactant species and their elemental buildup
        for reac in self._reactants:
            if reac in reaction_types:
                continue
            specie = Species([reac] + [0] * len(species_header))
            atoms = specie.find_constituents(quiet=True)
            reac_species.append(specie)
            reac_constituents.append(atoms)

        prod_constituents = []
        prod_species = []
        # Get all product species and their elemental buildup
        for prod in self._products:
            if prod in "NAN":
                continue
            specie = Species([prod] + [0] * len(species_header))
            atoms = specie.find_constituents(quiet=True)
            prod_species.append(specie)
            prod_constituents.append(atoms)

        # Get mass and number of reactants and products
        m_reacs = [reac_specie.get_mass() for reac_specie in reac_species]
        naive_reduced_mass = m_reacs[0] * m_reacs[1] / (m_reacs[0] + m_reacs[1])
        n_reacs = len(reac_constituents)
        n_prods = len(prod_constituents)
        if n_reacs == n_prods:
            for i, reac_constituent in enumerate(reac_constituents):
                # For each reactant, find which product is closest (most similar in buildup) to it.
                min_total = int(1e10)
                min_copy = None
                for j, prod_constituent in enumerate(prod_constituents):
                    diff = deepcopy(reac_constituent)
                    diff.subtract(prod_constituent)
                    total_change = 0
                    for element in elementList:
                        total_change += abs(diff[element])
                    if total_change < min_total:
                        min_total = total_change
                        min_index = j
                        min_diff = diff
                changing_species = Counter(
                    {k: c for k, c in min_diff.items() if c != 0}
                )

                items = changing_species.items()
                if len(items) == 1:
                    # Exchange reaction
                    tuple_items = tuple(items)[0]
                    if abs(tuple_items[1]) == 1:
                        # One element is switched
                        element_index = elementList.index(tuple_items[0])
                        # Set reduced mass to mass of switched element
                        reduced_mass = elementMass[element_index]
                        self.set_reduced_mass(float(reduced_mass))
                        logging.debug(
                            f"Predicted reduced mass of '{self}' to be {self._reduced_mass} (would have been {naive_reduced_mass})"
                        )
                        return
        elif n_reacs == 2 and n_prods == 1:
            # Addition reaction
            if reac_species[0].get_name().strip("#@") == reac_species[1].get_name().strip("#@"):
                # If the two species are the same (e.g. #H+#H-> #H2), set reduced mass to m/2
                mass = reac_species[0].get_mass()
                # mass = elementMass[elementList.index(reac_species[0].get_name().strip("#@"))]
                reduced_mass = float(mass) / 2.0
                self.set_reduced_mass(reduced_mass)
                logging.debug(
                    f"Predicted reduced mass of '{self}' to be {self._reduced_mass} (would have been {naive_reduced_mass})"
                )
                return
            elif any(species == Counter({"H": 1}) for species in reac_constituents):
                # If one of the species is #H, set reduced mass to 1
                self.set_reduced_mass(1.0)
                logging.debug(
                    f"Predicted reduced mass of '{self}' to be {self._reduced_mass} (would have been {naive_reduced_mass})"
                )
                return
            else:
                pass
        elif n_reacs == 1 and n_prods == 2:
            # Splitting reaction. Not in network (also not LH or ER type, so would never get here)
            pass
        msg = f"Could not predict reduced mass of '{self}' cleverly.\n"
        msg += f"Instead, using regular definition with masses of two reactants (mu={naive_reduced_mass:.3})."
        if self._gamma == 0.0:
            msg += " (Reaction is barrierless anyway)"
        logging.warning(msg)
        self.set_reduced_mass(naive_reduced_mass)

    def set_reduced_mass(self, reduced_mass: float) -> None:
        """Set the reduced mass to be used to calculate tunneling rate in AMU

        Args:
            reduced_mass (float): reduced mass of moving atoms
        """
        self._reduced_mass = reduced_mass

    def get_reduced_mass(self) -> float:
        """Get the reduced mass to be used to calculate tunneling rate in AMU

        Returns:
            float: reduced mass of moving atoms
        """
        return self._reduced_mass

    ## C

    def NANCheck(self, a):
        """Convert any Falsy statement to a NAN string

        Args:
            a: thing to check for falsiness

        Returns:
            bool: input a if truthy, otherwise NAN
        """
        return a if a else "NAN"

    def get_reaction_type(self) -> str:
        """Get the type of a reaction from the reactants
        First check the third reactant for a reaction type, then the second. If there are none
        in there, it will be regarded as a two body reaction.

        Returns:
            str:
        """
        if self.get_reactants()[2] in reaction_types:
            return self.get_reactants()[2]
        elif self.get_reactants()[1] in reaction_types:
            return self.get_reactants()[1]
        else:
            return "TWOBODY"

    def get_source(self) -> str:
        """Get the source of the reaction

        Returns:
            str: The source of the reaction
        """
        return self.source

    def set_source(self, source: str) -> None:
        """Set the source of the reaction

        Args:
            source (str): The source of the reaction
        """
        self.source = source
       
    def set_extrapolation(self, flag: bool) -> None:
        logging.info(f"Setting for {self} extrapolation to {flag}")  
        assert isinstance(flag, bool)
        self.extrapolate = flag
        
    def get_extrapolation(self) -> bool:
        return self.extrapolate

    def convert_surf_to_bulk(self) -> None:
        """Convert the surface species to bulk species in place for this reaction."""
        self.set_reactants([reac.replace("#", "@") for reac in self.get_reactants()])
        self.set_products([prod.replace("#", "@") for prod in self.get_products()])

    def check_element_conservation(self) -> None:
        if self.get_reaction_type() in ["FREEZE", "DESORB"]:
            return

        counter_reactants = Counter()
        for reac in self._reactants:
            if reac in reaction_types:
                continue
            if reac in ["NAN", "E-"]:
                continue
            specie = Species([reac] + [0] * len(species_header))
            atoms_counter_specie = specie.find_constituents(quiet=True)
            counter_reactants += atoms_counter_specie

        counter_products = Counter()
        for prod in self._products:
            if prod in reaction_types:
                continue
            if prod in ["NAN", "E-"]:
                continue
            specie = Species([prod] + [0] * len(species_header))
            atoms_counter_specie = specie.find_constituents(quiet=True)
            counter_products += atoms_counter_specie

        if counter_products != counter_reactants:
            msg = "Elements not conserved in a reaction.\n"
            msg += f"The following reaction caused this error: {self}.\n"
            msg += f"Reactants: {counter_reactants}. Products: {counter_products}"
            raise ValueError(msg)
        
    def check_charge_conservation(self) -> None:
        if self.get_reaction_type() in [
            "FREEZE",
            "DESORB",
            "DESOH2",
            "DESCR",
            "DEUVCR",
            "THERM",
        ]:
            return
        charge_reactants = 0
        for reac in self._reactants:
            if reac in ["NAN"]:
                continue
            specie = Species([reac] + [0] * len(species_header))
            charge_reactants += specie.get_charge()
        charge_products = 0
        for prod in self._products:
            if prod in ["NAN"]:
                continue
            specie = Species([prod] + [0] * len(species_header))
            charge_products += specie.get_charge()

        if charge_products != charge_reactants:
            if self.get_reaction_type() == "GAR":
                # GAR reactions do not conserve charge since it relies on grain electrons.
                return
            msg = "Charges not conserved in a reaction.\n"
            msg += f"The following reaction caused this error: {self}.\n"
            msg += f"Reactants: {charge_reactants}. Products: {charge_products}"
            raise ValueError(msg)


    def convert_gas_to_surf(self) -> None:
        """Convert the gas-phase species to surface species in place for this reaction.
        If any ions are produced, the ion is assumed to become neutral because it is on the surface.
        If any electrons are produced, they are assumed to be absorbed by the grain."""
        do_not_convert = reaction_types + ["E-", "NAN"]
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

    def __eq__(self, other) -> bool:
        """Check for equality against another reaction based on the products and reactants.
        Note that it does not check for the temperature ranges that the reactions might have!
        The Reaction.check_temperature_collision can be used for this purpose.

        Args:
            other: Another reaction set.

        Returns:
            bool: equality
        """
        if not isinstance(other, Reaction) and not isinstance(other, CoupledReaction):
            raise NotImplementedError(
                "Equality is not implemented for anything but comparing to other reactions."
            )
        if self.get_sorted_reactants() == other.get_sorted_reactants():
            if self.get_sorted_products() == other.get_sorted_products():
                return True
        return False

    def check_temperature_collision(self, other) -> bool:
        """Check if two reactions have overlapping temperature ranges, returning True means there is a collision.

        Args:
            other: Another reaction

        Raises:
            NotImplementedError: Currently we can only compare against instantiated Reaction objects.

        Returns:
            bool: Whether there is a collision (True), or not (False)
        """
        if not isinstance(other, Reaction) and not isinstance(other, CoupledReaction):
            raise NotImplementedError(
                "Equality is not implemented for anything but comparing to other reactions."
            )
        if (other.get_templow() > self.get_templow()) and (
            other.get_templow() < self.get_temphigh()
        ):
            return True
        if (other.get_temphigh() > self.get_templow()) and (
            other.get_temphigh() < self.get_temphigh()
        ):
            return True
        return False

    def changes_surface_count(self):
        """
        This checks whether a grain reaction changes number of particles on the surface
        2 reactants to 2 products won't but two reactants combining to one will.
        """
        if len([x for x in self.get_reactants() if "#" in x]) != len(
            [x for x in self.get_products() if "#" in x]
        ):
            return True
        if len([x for x in self.get_reactants() if "@" in x]) != len(
            [x for x in self.get_products() if "@" in x]
        ):
            return True
        return False

    def changes_total_mantle(self):
        """Check if the total grains on the mantle are changed by the reaction."""
        # If it's not just a movement between ice phases
        if ("BULK" not in self.get_reactants()[1]) and (
            "SWAP" not in self.get_reactants()[1]
        ):
            # if the number of ice species changes
            if self.changes_surface_count():
                return True
            else:
                return False
        else:
            return False

    def generate_ode_bit(self, i: int, species_names: list):
        self.ode_bit = _generate_reaction_ode_bit(i, species_names, self.body_count, self.get_reactants(), self.get_reaction_type())

    def to_UCL_format(self):
        """Convert a reaction to UCLCHEM reaction file format"""
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

    def _is_reaction_wrap(self, include_reactants=True, include_products=True):
        assert include_reactants or include_products, (
            "Either include reactants or products"
        )
        species_to_check = []
        if include_reactants:
            species_to_check += self.get_pure_reactants()
        if include_products:
            species_to_check += self.get_pure_products()
        return species_to_check

    def is_gas_reaction(
        self, include_reactants=True, include_products=True, strict=True
    ) -> bool:
        """Check whether it is a gas reaction, by default it is strict - all
        reactions must be in the gas-phase - if strict=False; any reaction in
        the gas-phase returns true.

        Args:
            include_reactants (bool, optional): Include the reactants. Defaults to True.
            include_products (bool, optional): Include the products. Defaults to True.
            strict (bool, optional): Choose between all (true) or any (false) must be gas phase . Defaults to True.

        Returns:
            bool: Is it a gas phase reaction?
        """
        checklist = [
            not (s.startswith("#") or s.startswith("@"))
            for s in self._is_reaction_wrap(include_reactants, include_products)
        ]
        return all(checklist) if strict else any(checklist)

    def is_ice_reaction(
        self, include_reactants=True, include_products=True, strict=True
    ) -> bool:
        """Check whether it is an ice (surface OR bulk) reaction
        
        By default it is strict (strict=True); all species must be in the ice phase
        If strict=False; any species in ice phase returns True

        Args:
            include_reactants (bool, optional): Include the reactants. Defaults to True.
            include_products (bool, optional): Include the products. Defaults to True.
            strict (bool, optional): Choose between all (true) or any (false) must be ice phase . Defaults to True.

        Returns:
            bool: Is it an ice phase reaction?
        """
        checklist = [
            (s.startswith("#") or s.startswith("@"))
            for s in self._is_reaction_wrap(include_reactants, include_products)
        ]
        return all(checklist) if strict else any(checklist)

    def is_surface_reaction(
        self, include_reactants=True, include_products=True, strict=False
    ) -> bool:
        """Check whether it is a surface reaction, defaults to non-strict since many
        important surface reactions can lead to desorption in some way.
        
        By default it is NOT strict (strict=False); any species on the surface returns true
        If strict=True; all species must be on the ice phase

        Args:
            include_reactants (bool, optional): Include the reactants. Defaults to True.
            include_products (bool, optional): Include the products. Defaults to True.
            strict (bool, optional): Choose between all (true) or any (false) must be on the surface . Defaults to False.

        Returns:
            bool: Is it a surface reaction?
        """
        checklist = [
            s.startswith("#")
            for s in self._is_reaction_wrap(include_reactants, include_products)
        ]
        return all(checklist) if strict else any(checklist)

    def is_bulk_reaction(
        self, include_reactants=True, include_products=True, strict=False
    ) -> bool:
        """Check whether it is a bulk reaction, defaults to non-strict since many
        important bulk reactions interact with the surface.
        
        By default it is NOT strict (strict=False); any species in the bulk returns true
        If strict=True; all species must be on the ice phase

        Args:
            include_reactants (bool, optional): Include the reactants. Defaults to True.
            include_products (bool, optional): Include the products. Defaults to True.
            strict (bool, optional): Choose between all (true) or any (false) must in the bulk . Defaults to False.

        Returns:
            bool: Is it a bulk reaction?
        """
        checklist = [
            s.startswith("@")
            for s in self._is_reaction_wrap(include_reactants, include_products)
        ]
        return all(checklist) if strict else any(checklist)

    def __str__(self):
        return (
            " + ".join(filter(lambda r: r != "NAN", self.get_reactants()))
            + " -> "
            + " + ".join(filter(lambda p: p != "NAN", self.get_products()))
        )

    def __repr__(self):
        return (
            self.get_reaction_type()
            + " reaction: "
            + " + ".join(
                filter(
                    lambda r: (r != "NAN") and (r not in reaction_types),
                    self.get_reactants(),
                )
            )
            + " -> "
            + " + ".join(filter(lambda p: p != "NAN", self.get_products()))
        )
        
    def __hash__(self):
        return hash(f"{self.get_alpha(), self.get_beta(), self.get_gamma(), self.get_reactants(), self.get_products(), self.get_templow(), self.get_temphigh()}")


class CoupledReaction(Reaction):
    def __init__(self, input):
        super().__init__(input)
        self.partner = None

    def set_partner(self, partner: Reaction):
        self.partner = partner

    def get_partner(self):
        return self.partner

def _generate_reaction_ode_bit(i: int, species_names: list, body_count: int, reactants: list[str], reaction_type: str=None):
        """Every reaction contributes a fixed rate of change to whatever species it
        affects. We create the string of fortran code describing that change here.

        Args:
            i (int): index of reaction in network in python format (counting from 0)
            species_names (list): List of species names so we can find index of reactants in species list
            body_count (bool): Number of bodies in the reaction, used to determine how many factors of density to include
            reactants (list[str]): The reactants of the reaction
        """
        ode_bit = f"+RATE({i + 1})"
        # every body after the first requires a factor of density
        for body in range(body_count):
            ode_bit = ode_bit + "*D"
        
        # GAR needs an additional factor of density: 
        if reaction_type == "GAR":
            ode_bit = ode_bit + "*D"

        # then bring in factors of abundances
        for species in reactants:
            if species in species_names:
                ode_bit += f"*Y({species_names.index(species) + 1})"

            elif species == "BULKSWAP":
                ode_bit += "*bulkLayersReciprocal"
            elif species == "SURFSWAP":
                ode_bit += "*totalSwap/safeMantle"
            elif species in ["DEUVCR", "DESCR", "DESOH2", "ER", "ERDES"]:
                ode_bit = ode_bit + "/safeMantle"
                if species == "DESOH2":
                    ode_bit = ode_bit + f"*Y({species_names.index('H') + 1})"

            if "H2FORM" in reactants:
                # only 1 factor of H abundance in Cazaux & Tielens 2004 H2 formation so stop looping after first iteration
                break

        if "LH" in reactants[2]:
            if "@" in reactants[0]:
                ode_bit += "*bulkLayersReciprocal"
        return ode_bit