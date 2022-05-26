"""
This python file contains all functions for de-duplicating species and reaction lists,
checking for common errors, and automatic addition of reactions such as freeze out,
desorption and bulk reactions for three phase models.
"""
from .species import Species, elementList
from .reaction import Reaction
from copy import deepcopy
from numpy import unique
from numpy import any as np_any


class Network:
    def __init__(self, species, reactions, three_phase=False):
        """
        Simple class to store network information such as indices of important reactions.
        Also logical home of functions meant to make network sensible.
        """

        self.species_list = species
        self.remove_duplicate_species()
        self.excited_species =  self.check_for_excited_species()
        self.three_phase = three_phase
        if self.three_phase:
            self.add_bulk_species()
        self.species_list.sort()
        self.species_list.append(Species(["E-", 0, 0, 0, 0, 0, 0]))
        self.species_list[-1].n_atoms = 1
        self.reaction_list = reactions
        self.add_grain_reactions()

    def remove_duplicate_species(self):
        """Alerts user if the same species appears twice in species_list
        then de-duplicates list

        """
        for species in self.species_list:
            if self.species_list.count(species) > 1:
                print(f"\t {species.name} appears twice in input species list")

        self.species_list = list(unique(self.species_list))

    def check_and_filter_species(self):
        """Check every speces in network appears in at least one reaction.
        Remove any that do not and alert user.
        """
        species_names = [species.name for species in self.species_list]
        # check for species not involved in any reactions
        lostSpecies = []
        for species in self.species_list:

            # keep species that appear in a reaction
            reac_keeps = False
            for reaction in self.reaction_list:
                if species.name in reaction.reactants or species.name in reaction.products:
                    reac_keeps = True
                    break

            # remove the species if it didn't make it into either keep list
            if not (reac_keeps):
                lostSpecies.append(species.name)
        for species in lostSpecies:
            self.species_list.remove(species)

        # then alert user to changes
        if len(lostSpecies) > 0:
            print("\tSpecies in input list that do not appear in final list:")
            print("\t", lostSpecies)
            print("\n")
        else:
            print("\tAll input species in final network")
        for species in self.species_list:
            species.find_constituents()

        # add in pseudo-species to track mantle
        mantle_specs = []
        new_spec = [999] * 7
        new_spec[0] = "BULK"
        mantle_specs.append(Species(new_spec))
        new_spec[0] = "SURFACE"
        mantle_specs.append(Species(new_spec))
        self.species_list = self.species_list + mantle_specs

    def add_bulk_species(self):
        """For three phase models, MakeRates will produce the version of the species in the bulk
        so that the user doesn't have to endlessly relist the same species
        """
        speciesNames = [species.name for species in self.species_list]
        new_species = []
        try:
            h2o_binding_energy = speciesNames.index("#H2O")
            h2o_binding_energy = self.species_list[h2o_binding_energy].binding_energy
        except:
            error = "You are trying to create a three phase model but #H2O is not in your network"
            error += "\nThis is likely an error so Makerates will not complete"
            error += "\nTry adding #H2O or switching to three_phase=False in Makerates.py"
            raise RuntimeError(error)
        for species in self.species_list:
            if species.is_surface_species():
                if not species.name.replace("#", "@") in speciesNames:
                    new_spec = deepcopy(species)
                    new_spec.name = new_spec.name.replace("#", "@")
                    new_spec.binding_energy = h2o_binding_energy
                    new_species.append(new_spec)
        self.species_list = self.species_list + new_species

    def add_grain_reactions(self):
        # check which species are changed on freeze or desorb
        self.check_freeze_and_desorbs()

        # Need additional grain reactions including non-thermal desorption and chemically induced desorption
        self.add_freeze_reactions()
        self.add_desorb_reactions()
        self.add_chemdes_reactions()
        if self.excited_species:
            reaction_list = self.add_excited_surface_reactions()
        if self.three_phase:
            reaction_list = self.add_bulk_reactions()
        self.check_and_filter_species()

        self.reaction_list = sorted(self.reaction_list, key=lambda x: (x.reac_type, x.reactants[0]))

    def check_freeze_and_desorbs(self):
        """`add_freeze_reactions()` and `add_desorb_reactions()` automatically generate
        all desorption and freeze out reactions. However, user may want to change a species on freeze out
        eg C+ becomes #C rather than #C+. This function checks for that and updates species so they'll
        freeze or desorb correctly when reactions are generated.

        Args:
            species_list (list): list of species objects including all species in network
            reaction_list (list): list of reaction objects including all reactions in network

        Returns:
            list: species and reaction lists with user specified freeze and desorb reactions removed (but species updated)
        """
        species_names = [species.name for species in self.species_list]
        desorbs = [x for x in self.reaction_list if x.reac_type == "DESORB"]
        for desorb in desorbs:
            species_index = species_names.index(desorb.reactants[0])
            self.species_list[species_index].desorb_products = desorb.products
            #also modify bulk species desorb_products
            if "@"+desorb.reactants[0][1:] in species_names:
                bulk_index = species_names.index("@"+desorb.reactants[0][1:])
                self.species_list[bulk_index].desorb_products = desorb.products

        #for all listed freeze out reactions, add them to correct species
        freezes = [x for x in self.reaction_list if x.reac_type == "FREEZE"]
        for freeze in freezes:
            species_index = species_names.index(freeze.reactants[0])
            if self.species_list[species_index].freeze_products == None:
                self.species_list[species_index].freeze_products = {}
            self.species_list[species_index].freeze_products[
                ",".join(freeze.products)
            ] = freeze.alpha
            
        #then add default freeze out for species without a listed freeze out
        for species in self.species_list:
            if not species.is_grain_species():
                if species.freeze_products is None:
                    species.add_default_freeze()

        self.reaction_list = [
            reaction for reaction in self.reaction_list if reaction not in freezes + desorbs
        ]

    def add_freeze_reactions(self):
        """Save the user effort by automatically generating freeze out reactions"""

        for species in self.species_list:
            if not species.is_grain_species():
                for products, alpha in species.freeze_products.items():
                    newReaction = Reaction(
                        [species.name, "FREEZE", "NAN"]
                        + products.split(",")
                        + [alpha, 0, species.binding_energy, 0.0, 10000.0]
                    )
                    self.reaction_list.append(newReaction)

    def add_desorb_reactions(self):
        """Save the user effort by automatically generating desorption reactions"""
        desorb_reacs = ["DESOH2", "DESCR", "DEUVCR", "THERM"]

        for species in self.species_list:
            if species.is_surface_species():
                for reacType in desorb_reacs:
                    newReaction = Reaction(
                        [species.name, reacType, "NAN"]
                        + species.desorb_products
                        + [1, 0, species.binding_energy, 0.0, 10000.0]
                    )
                    self.reaction_list.append(newReaction)
            if species.is_bulk_species() and not species.is_refractory:
                newReaction = Reaction(
                    [species.name, "THERM", "NAN"]
                    + species.desorb_products
                    + [1, 0, species.binding_energy, 0.0, 10000.0]
                )
                self.reaction_list.append(newReaction)

    def add_chemdes_reactions(self):
        """We have the user list all Langmuir-Hinshelwood and Eley-Rideal
        reactions once. Then we duplicate so that the reaction branches
        with products on grain and products desorbing.
        """
        new_reacs = []
        for reaction in self.reaction_list:
            if reaction.reac_type in ["LH", "ER"]:
                new_reac = deepcopy(reaction)
                new_reac.reac_type = new_reac.reac_type + "DES"
                new_reac.reactants[2] = new_reac.reactants[2] + "DES"
                for i, product in enumerate(new_reac.products):
                    if ("#" in product) or ("@" in product):
                        new_reac.products[i] = new_reac.products[i][1:]
                    else:
                        if product != "NAN":
                            print(
                                "All Langmuir-Hinshelwood and Eley-Rideal reactions should be input with products on grains only."
                            )
                            print(
                                "The fraction of products that enter the gas is dealt with by Makerates and UCLCHEM."
                            )
                            print("the following reaction caused this warning")
                            print("\t", reaction)
                new_reacs.append(new_reac)

        self.reaction_list = self.reaction_list + new_reacs

    def check_for_excited_species(self):
        check = False
        for species in self.species_list:
            if '*' in species.name:
                check=True 
        return check

    def add_excited_surface_reactions(self):
        """ All excited species will relax to the ground state if they do not react
        the vibrational frequency of the species is used as a pseudo approximation of the rate coefficient
        We assume all grain reactions have an excited variant. For example:
        #A, #B LH #C will have the variants:
        #A*, #B EXSOLID #C  and  #A, #B* EXSOLID #C
        If only one of the reactants in the base reaction has an excited counterpart then
        only one excited version of that reaction is created.
        """
        excited_species = [x for x in self.species_list if "*" in x.name]
        lh_reactions = [x for x in self.reaction_list if "LH" in x.reactants]
        lh_reactions = lh_reactions + [x for x in self.reaction_list if "LHDES" in x.reactants]

        new_reactions = []
        
        # add relaxation of excited species
        for spec in excited_species:
            relax_reac = [spec.name,'EXRELAX','NAN',spec.name[:-1],'NAN','NAN','NAN',1.0,0.0,0.0,0.0,10000]
            new_react = Reaction(relax_reac)
            new_reactions.append(new_react)


        for reaction in lh_reactions:
            # if both #A and #B have excited counterparts
            if reaction.reactants[0]+'*' in [excited_species[i].name for i in range(len(excited_species))] and reaction.reactants[1]+'*' in [excited_species[i].name for i in range(len(excited_species))]:
                new_reac_A_list = [reaction.reactants[0] + '*', reaction.reactants[1], 'EXSOLID']
                new_reac_A_list = new_reac_A_list + reaction.products + [reaction.alpha, 0, 0, 0, 10000]
                new_reac_B_list = [reaction.reactants[0],reaction.reactants[1]+'*','EXSOLID']
                new_reac_B_list = new_reac_B_list + reaction.products + [reaction.alpha, 0, 0, 0, 10000]

                new_reac_A = Reaction(new_reac_A_list)
                new_reac_B = Reaction(new_reac_B_list)
                
                # stops duplicate reactions e.g. #H* + #H and #H + #H*
                if new_reac_A != new_reac_B:
                    new_reactions.append(new_reac_A)
                    new_reactions.append(new_reac_B)
                else:
                    new_reactions.append(new_reac_A)
                    
            # if only #A has an excited counterpart
            elif reaction.reactants[0]+'*' in [excited_species[i].name for i in range(len(excited_species))]:
                new_reac_A_list = [reaction.reactants[0] + '*', reaction.reactants[1], 'EXSOLID']
                new_reac_A_list = new_reac_A_list + reaction.products + [reaction.alpha, 0, 0, 0, 10000]
                new_reac_A = Reaction(new_reac_A_list)
                new_reactions.append(new_reac_A)

            # if only #B has an excited counterpart
            elif reaction.reactants[1]+'*' in [excited_species[i].name for i in range(len(excited_species))]:
                new_reac_B_list = [reaction.reactants[0], reaction.reactants[1]+'*', 'EXSOLID']
                new_reac_B_list = new_reac_B_list + reaction.products + [reaction.alpha, 0, 0, 0, 10000]
                new_reac_B = Reaction(new_reac_B_list)
                new_reactions.append(new_reac_B)
    
        self.reaction_list = self.reaction_list + new_reactions


    def add_bulk_reactions(self):
        """We assume any reaction that happens on the surface of grains can also happen
        in the bulk (just more slowly due to binding energy). The user therefore only
        lists surface reactions in their input reaction file and we duplicate here.
        """
        lh_reactions = [x for x in self.reaction_list if "LH" in x.reactants]
        lh_reactions = lh_reactions + [x for x in self.reaction_list if "LHDES" in x.reactants]
        ex_reactions = [x for x in self.reaction_list if "CRS" in x.reactants]
        ex_reactions = ex_reactions + [x for x in self.reaction_list if "EXSOLID" in x.reactants]
        ex_reactions = ex_reactions + [x for x in self.reaction_list if "EXRELAX" in x.reactants]
        surface_reactions = lh_reactions + ex_reactions

        new_reactions = []
        for reaction in surface_reactions:
            new_reac = deepcopy(reaction)
            new_reac.convert_to_bulk()
            new_reactions.append(new_reac)
        new_reactions = [reac for reac in new_reactions if reac not in self.reaction_list]

        bulk_species = [x for x in self.species_list if "@" in x.name]
        for species in bulk_species:
            # add individual swapping
            if not species.is_refractory:
                new_reac_list = [species.name, "BULKSWAP", "NAN", species.name.replace("@", "#")]
                new_reac_list = new_reac_list + ["NAN", "NAN", "NAN", 1, 0, 0, 0, 10000]
                new_reac = Reaction(new_reac_list)
                new_reactions.append(new_reac)

            # and the reverse
            new_reac_list[0] = species.name.replace("@", "#")
            new_reac_list[1] = "SURFSWAP"
            new_reac_list[3] = species.name
            new_reac = Reaction(new_reac_list)
            new_reactions.append(new_reac)

        self.reaction_list = self.reaction_list + new_reactions

    # check reactions to alert user of potential issues including repeat reactions
    # and multiple freeze out routes
    def check_network(self):
        """Run through the list of reactions and check for obvious errors such
        as duplicate reactions, multiple freeze out routes (to warn, not necessarily
        an error), etc.
        """
        self.freeze_checks()
        self.duplicate_checks()
        self.index_important_reactions()
        self.index_important_species()

    def freeze_checks(self):
        """Check that every species freezes out and alert the user if a
        species freezes out via mutiple routes. This isn't necessarily an
        error so best just print.
        """
        print("\tSpecies with multiple freeze outs, check alphas:")
        for spec in self.species_list:
            freezes = 0
            for reaction in self.reaction_list:
                if spec.name in reaction.reactants and "FREEZE" in reaction.reactants:
                    freezes += 1
            if freezes > 1:
                print(f"\t{spec.name} freezes out through {freezes} routes")
            elif freezes < 1 and not spec.is_grain_species():
                print(f"\t{spec.name} does not freeze out")

    def duplicate_checks(self):
        """
        Check reaction network to make sure no reaction appears twice unless
        they have different temperature ranges.
        """
        print("\n\tPossible duplicate reactions for manual removal:")
        duplicates = False
        for i, reaction1 in enumerate(self.reaction_list):
            if not reaction1.duplicate:
                for j, reaction2 in enumerate(self.reaction_list):
                    if i != j:
                        if reaction1 == reaction2:
                            if not (
                                (reaction1.templow >= reaction2.temphigh)
                                or (reaction1.temphigh <= reaction2.templow)
                            ):
                                print(f"\tReactions {i+1} and {j+1} are possible duplicates")
                                print(reaction1)
                                print(reaction2)
                                duplicates = True
                                # adjust temperatures so temperature ranges are adjacent
                                if reaction1.temphigh > reaction2.temphigh:
                                    if reaction1.templow < reaction2.temphigh:
                                        print(
                                            f"\tReactions {i+1} and {j+1} have non-adjacent temperature ranges"
                                        )
                                reaction1.duplicate = True
                                reaction2.duplicate = True
        if not duplicates:
            print("\tNone")

    def index_important_reactions(self):
        """We have a whole bunch of important reactions and we want to store
        their indices. We find them all here.
        """

        # Any None values in dictionary will raise an error
        # therefore these reactions are mandatory and makerates will not complete if the user doesn't supply them.
        self.important_reactions = {
            "nR_H2Form_CT": None,
            "nR_H2Form_ERDes": None,
            "nR_H2Form_ER": None,
            "nR_H2Form_LH": None,
            "nR_H2Form_LHDes": None,
            "nR_HFreeze": None,
            "nR_EFreeze": None,
            "nR_H2_hv": None,
        }

        # this looks complex but each if statement just uniquely identifies a special reaction
        # if found, it is added to the dictionary with its fortran index as the value
        for i, reaction in enumerate(self.reaction_list):
            if ("CO" in reaction.reactants) and ("PHOTON" in reaction.reactants):
                if "O" in reaction.products and "C" in reaction.products:
                    self.important_reactions["nR_CO_hv"] = i + 1
            if ("C" in reaction.reactants) and ("PHOTON" in reaction.reactants):
                self.important_reactions["nR_C_hv"] = i + 1
            if "H2FORM" in reaction.reactants:
                self.important_reactions["nR_H2Form_CT"] = i + 1
            if ("H" in reaction.reactants) and ("#H" in reaction.reactants):
                if "H2" in reaction.products:
                    self.important_reactions["nR_H2Form_ERDes"] = i + 1
                elif "#H2" in reaction.products:
                    self.important_reactions["nR_H2Form_ER"] = i + 1
            if (reaction.reactants.count("#H") == 2) and ("LH" in reaction.reactants):
                self.important_reactions["nR_H2Form_LH"] = i + 1
            if (reaction.reactants.count("#H") == 2) and ("LHDES" in reaction.reactants):
                self.important_reactions["nR_H2Form_LHDes"] = i + 1
            if ("H" in reaction.reactants) and ("FREEZE" in reaction.reactants):
                self.important_reactions["nR_HFreeze"] = i + 1
            if ("H2" in reaction.reactants) and ("FREEZE" in reaction.reactants):
                self.important_reactions["nR_H2Freeze"] = i + 1
            if ("E-" in reaction.reactants) and ("FREEZE" in reaction.reactants):
                self.important_reactions["nR_EFreeze"] = i + 1
            if ("H2" in reaction.reactants) and ("PHOTON" in reaction.reactants):
                self.important_reactions["nR_H2_hv"] = i + 1
            if (
                ("H2" in reaction.reactants)
                and ("CRP" in reaction.reactants)
                and (reaction.products.count("H") == 2)
            ):
                self.important_reactions["nR_H2_crp"] = i + 1
        if np_any([value is None for value in self.important_reactions.values()]):
            missing_reac_error = "Input reaction file is missing mandatory reactions"
            missing_reac_error += (
                "\nH and E- freeze out as well as H2 formation and photodissociation"
            )
            missing_reac_error += " must all be included in user reaction list. Check default_grain_network.csv for example"
            raise RuntimeError(missing_reac_error)

    def index_important_species(self):
        self.species_indices = {}
        names = [species.name for species in self.species_list]
        for element in [
            "C+",
            "H+",
            "H2",
            "SI+",
            "S+",
            "CL+",
            "CO",
            "HE+",
            "#H",
            "#H2",
            "#N",
            "#O",
            "#OH",
            "SURFACE",
            "BULK",
        ] + elementList:
            try:
                species_index = names.index(element) + 1
            except:
                print(f"\t{element} not in network, adding dummy index")
                species_index = len(self.species_list) + 1
            name = "n" + element.lower().replace("+", "x").replace("e-", "elec").replace("#", "g")
            self.species_indices[name] = species_index
