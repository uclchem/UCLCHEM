reaction_types = [
    "PHOTON",
    "CRP",
    "CRPHOT",
    "FREEZE",
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
]


class Reaction:
    def __init__(self, inputRow):
        self.reactants = [
            inputRow[0].upper(),
            inputRow[1].upper(),
            self.NANCheck(inputRow[2]).upper(),
        ]
        self.products = [
            inputRow[3].upper(),
            self.NANCheck(inputRow[4]).upper(),
            self.NANCheck(inputRow[5]).upper(),
            self.NANCheck(inputRow[6]).upper(),
        ]
        self.alpha = float(inputRow[7])
        self.beta = float(inputRow[8])
        self.gamma = float(inputRow[9])
        self.templow = float(inputRow[10])
        self.temphigh = float(inputRow[11])
        self.reac_type = self.get_reaction_type()
        self.duplicate = False

        # body_count is the number of factors of density to include in ODE
        # we drop a factor of density from both the LHS and RHS of ODES
        # So reactions with 1 body have no factors of density which we manage by counting from -1
        self.body_count = -1
        for reactant in self.reactants:
            if (reactant not in reaction_types) and reactant != "NAN":
                self.body_count += 1
            if reactant in ["DESOH2", "FREEZE"]:
                self.body_count += 1
            if reactant in ["LH", "LHDES"]:
                self.body_count -= 1

        if (self.get_reaction_type() == "FREEZE") and (self.reactants[0][-1] == "+"):
            self.beta = 1

    def NANCheck(self, a):
        aa = a if a else "NAN"
        return aa

    def get_reaction_type(self):
        if self.reactants[2] in reaction_types:
            return self.reactants[2]
        else:
            if self.reactants[1] in reaction_types:
                return self.reactants[1]
            else:
                return "TWOBODY"

    def convert_to_bulk(self):
        for i in range(len(self.reactants)):
            self.reactants[i] = self.reactants[i].replace("#", "@")
        for i in range(len(self.products)):
            self.products[i] = self.products[i].replace("#", "@")

    def __eq__(self, other):
        if set(self.reactants) == set(other.reactants):
            if set(self.products) == set(other.products):
                return True
        return False

    def changes_surface_count(self):
        """
        This checks whether a grain reaction changes number of particles on the surface
        2 reactants to 2 products won't but two reactants combining to one will.
        """
        if len([x for x in self.reactants if "#" in x]) != len(
            [x for x in self.products if "#" in x]
        ):
            return True
        if len([x for x in self.reactants if "@" in x]) != len(
            [x for x in self.products if "@" in x]
        ):
            return True
        return False

    def changes_total_mantle(self):
        # If it's not just a movement between ice phases
        if ("BULK" not in self.reactants[1]) and ("SWAP" not in self.reactants[1]):
            # if the number of ice species changes
            if self.changes_surface_count():
                return True
            else:
                return False
        else:
            return False

    def generate_ode_bit(self, i, species_names, three_phase):
        """Every reaction contributes a fixed rate of change to whatever species it
        affects. We create the string of fortran code describing that change here.

        Args:
            i (int): index of reaction in network in python format (counting from 0)
            species_names (list): List of species names so we can find index of reactants in species list
            three_phase (bool): Bool indicating whether this is three phase network
        """
        ode_bit = f"+RATE({i+1})"
        # every body after the first requires a factor of density
        for body in range(self.body_count):
            ode_bit = ode_bit + f"*D"

        # then bring in factors of abundances
        for species in self.reactants:
            if species in species_names:
                ode_bit += f"*Y({species_names.index(species)+1})"
            elif species == "BULKSWAP":
                ode_bit += "*bulkLayersReciprocal"
            elif species == "SURFSWAP":
                ode_bit += "*totalSwap/safeMantle"
            elif species in ["DEUVCR", "DESCR", "DESOH2", "ER", "ERDES"]:
                ode_bit = ode_bit + f"/safeMantle"
                if species == "DESOH2":
                    ode_bit = ode_bit + f"*Y({species_names.index('H')+1})"
            elif (species in ["THERM"]) and not (three_phase):
                ode_bit += f"*D/safeMantle"
            if "H2FORM" in self.reactants:
                # only 1 factor of H abundance in Cazaux & Tielens 2004 H2 formation so stop looping after first iteration
                break

        if "LH" in self.reactants[2]:
            if "@" in self.reactants[0]:
                ode_bit += "*bulkLayersReciprocal"
        self.ode_bit = ode_bit

    def __str__(self):
        return " + ".join(self.reactants) + " -> " + " + ".join(self.products)
