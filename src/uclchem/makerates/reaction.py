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
]


class Reaction:
    def __init__(self, inputRow, reaction_source=None):
        if isinstance(inputRow, Reaction):
            self.set_reactants(inputRow.get_reactants())
            self.set_products(inputRow.get_products())
            self.set_alpha(inputRow.get_alpha())
            self.set_beta(inputRow.get_beta())
            self.set_gamma(inputRow.get_gamma())
            self.set_templow(inputRow.get_templow())
            self.set_temphigh(inputRow.get_temphigh())
            self.set_extrapolation(inputRow.get_extrapolation())
        else:
            try:
                self.set_reactants(
                    [
                        inputRow[0].upper(),
                        inputRow[1].upper(),
                        self.NANCheck(inputRow[2]).upper(),
                    ]
                )
                self.set_products(
                    [
                        self.NANCheck(inputRow[3].upper()),
                        self.NANCheck(inputRow[4]).upper(),
                        self.NANCheck(inputRow[5]).upper(),
                        self.NANCheck(inputRow[6]).upper(),
                    ]
                )
                self.set_alpha(float(inputRow[7]))
                self.set_beta(float(inputRow[8]))
                self.set_gamma(float(inputRow[9]))
                self.set_templow(float(inputRow[10]))
                self.set_temphigh(float(inputRow[11]))
                self.set_extrapolation(bool(inputRow[12]) if len(inputRow) > 12 else False)
            except IndexError as error:
                raise ValueError(
                    "Input for Reaction should be a list of length 12"
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
        assert isinstance(flag, bool)
        self.extrapolate = flag
        
    def get_extrapolation(self) -> bool:
        return self.extrapolate

    def convert_surf_to_bulk(self) -> None:
        """Convert the surface species to bulk species in place for this reaction."""
        self.set_reactants([reac.replace("#", "@") for reac in self.get_reactants()])
        self.set_products([prod.replace("#", "@") for prod in self.get_products()])

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
        self.ode_bit = _generate_reaction_ode_bit(i, species_names, self.body_count, self.get_reactants())

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

def _generate_reaction_ode_bit(i: int, species_names: list, body_count: int, reactants: list[str]):
        """Every reaction contributes a fixed rate of change to whatever species it
        affects. We create the string of fortran code describing that change here.

        Args:
            i (int): index of reaction in network in python format (counting from 0)
            species_names (list): List of species names so we can find index of reactants in species list
            three_phase (bool): Bool indicating whether this is three phase network
        """
        ode_bit = f"+RATE({i + 1})"
        # every body after the first requires a factor of density
        for body in range(body_count):
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