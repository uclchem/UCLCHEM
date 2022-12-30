"""
This python file contains all functions for de-duplicating species and reaction lists,
checking for common errors, and automatic addition of reactions such as freeze out,
desorption and bulk reactions for three phase models.
"""
from .species import Species, elementList
from .reaction import Reaction, reaction_types
import logging
from copy import deepcopy
from numpy import unique
from numpy import any as np_any
from typing import Union


class Network:
    def __init__(self, species, reactions, three_phase=False, user_defined_bulk=[]):
        """
        Simple class to store network information such as indices of important reactions.
        Also logical home of functions meant to make network sensible.
        """

        assert len(set([s.name for s in species])) == len(
            species
        ), "Cannot have duplicate species in the species list."
        self._species_dict = {s.name: s for s in species}
        # self.remove_duplicate_species()
        self.excited_species = self.check_for_excited_species()
        self.user_defined_bulk = user_defined_bulk
        self.three_phase = three_phase
        # We used to add bulk here?
        # if self.three_phase:
        #     self.add_bulk_species()
        # self.species_list.sort()
        electron_specie = Species(["E-", 0, 0.0, 0, 0, 0, 0])
        electron_specie.n_atoms = 1
        self.add_species(electron_specie)
        self._reactions_dict = {k: v for k, v in enumerate(reactions)}

        #### Add reactions and species   ####
        # check which species are changed on freeze or desorb
        self.check_freeze_and_desorbs()

        # Need additional grain reactions including non-thermal desorption and chemically induced desorption
        self.add_freeze_reactions()
        if self.three_phase:
            self.add_bulk_species()
            self.add_bulk_reactions()
        self.add_desorb_reactions()
        self.add_chemdes_reactions()
        if self.excited_species:
            self.add_excited_surface_reactions()

        # Sort the reactions before returning them, this is important for convergence of the ODE
        self.sort_reactions()
        self.sort_species()

        self.check_and_filter_species()

    @property
    def species(self):
        return self._species_dict

    @property
    def reactions(self):
        return self._reactions_dict.values()

    # Make sure that the reaction_list or species_list are never explicitely accessed.

    @property
    def reaction_list(self):
        logging.warning(
            "Reaction list should not be accessed directly, but using the get and set methods"
        )
        return None

    @reaction_list.setter
    def reaction_list(self, value):
        raise Exception(
            "Do not set reaction lists explicitely, use the add and remove interfaces for reactions"
        )

    @property
    def species_list(self):
        logging.warning(
            "Species list should not be accessed directly, but using the get and set methods"
        )
        return None

    @species_list.setter
    def species_list(self, value):
        raise Exception(
            "Do not set reaction lists explicitely, use the add and remove interfaces for reactions"
        )

    def add_reactions(
        self, reactions: Union[Union[Reaction, str], list[Union[Reaction, str]]]
    ):
        if isinstance(reactions, list):
            if isinstance(reactions[0], Reaction):
                # if it is a list of reactions, no action is needed.
                pass
            elif isinstance(reactions[0], list):
                try:
                    reactions = [Reaction(reac) for reac in reactions]
                except ValueError as error:
                    raise ValueError(
                        "Failed to convert the list of csv entries to a reaction"
                    ) from error
        elif isinstance(reactions, Reaction):
            reactions = [reactions]
        else:
            raise ValueError(
                "Input must either be (a list of) Reaction class or csv entries"
            )
        current_reaction_list = self.get_reaction_list()
        for reaction in reactions:
            if reaction in current_reaction_list:
                # See if we have a collision with the any reactions with identical reactants and
                # products, but different temperature ranges to avoid double definitions.
                similar_reactions = self.find_similar_reactions(reaction)
                for similar_reaction in similar_reactions:
                    if reaction.check_temperature_collision(similar_reaction):
                        raise RuntimeError(
                            f"There already is a {reaction} present that has overlapping temperature ranges. Not adding the reaction."
                        )

            # quick check to make sure all species in the reaction are in the species list.
            species_to_add = filter(
                lambda spec: (spec not in self.get_species_list())
                and (spec not in ["NAN", "", "E-"] + reaction_types),
                reaction.get_reactants() + reaction.get_products(),
            )
            # if any(
            #     specie not in self.get_species_list()
            #     for specie in reaction.get_reactants() + reaction.get_products()
            # ):
            for specie in species_to_add:
                logging.debug(f"Trying to add specie {specie}")
                # TODO: get more sensible mass
                self.add_species(Species([specie, 1000, 0.0, 0.0, 0.0, 0.0, 0.0]))
            # Index and add the new reaction.
            new_idx = list(self._reactions_dict.keys())[-1] + 1
            self._reactions_dict[new_idx] = reaction

    def find_similar_reactions(self, reaction: Reaction) -> list[int]:
        return list(
            filter(
                lambda x: x != None,
                [k if v == reaction else None for k, v in self._reactions_dict.items()],
            )
        )

    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        del self._reactions_dict[reaction_idx]

    def remove_reaction(self, reaction: Reaction) -> None:
        # In this reaction we use equality as defined for reactions, this does as of now not include
        # checking temperature ranges, so more than one key could be returned.
        reaction_index = self.find_similar_reactions(reaction)
        if len(reaction_index) == 1:
            logging.debug(
                f"Trying to remove index: {reaction_index[0]}: {self._reactions_dict[reaction_index[0]]} "
            )
            del self._reactions_dict[reaction_index[0]]
        elif len(reaction_index) == 0:
            logging.warning(
                f"The reaction {reaction} is not present in the reaction set, so cannot remove it"
            )
        elif len(reaction_index) > 1:
            raise (
                ValueError(
                    "found more than one indices for the reaction {reaction}, remove by index instead of by reaction."
                )
            )

    def get_reaction(self, reaction_idx: int) -> Reaction:
        return deepcopy(self._reactions_dict[reaction_idx])

    def set_reaction(self, reaction_idx: int, reaction: Reaction) -> None:
        self._reactions_dict[reaction_idx] = Reaction

    def get_reaction_dict(self) -> dict[int, Reaction]:
        return self._reactions_dict

    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        self._reactions_dict = new_dict

    def get_reaction_list(self) -> list[Reaction]:
        return list(self._reactions_dict.values())

    def sort_reactions(self) -> None:
        """Sort the reactions by the reaction type first, reactants second."""
        reaction_dict = self.get_reaction_dict()
        logging.debug(
            f"Before sorting reactions {[(k, v) for i, (k, v) in enumerate(self.get_reaction_dict().items()) if i < 5]}"
        )
        self.set_reaction_dict(
            dict(
                sorted(
                    reaction_dict.items(),
                    key=lambda kv: (
                        kv[1].get_reaction_type(),
                        kv[1].get_reactants()[0],
                    ),
                )
            )
        )
        logging.debug(
            f"After sorting reactions {[(k,v ) for i, (k, v) in enumerate(self.get_reaction_dict().items()) if i < 5]}"
        )
        assert len(reaction_dict) == len(
            self.get_reaction_dict()
        ), "Sorting the species caused a difference in the number of species"

    def add_species(
        self, species: Union[Union[Species, str], list[Union[Species, str]]]
    ):
        if isinstance(species, list):
            if isinstance(species[0], Species):
                # if it is a list of Species, no action is needed.
                pass
            elif isinstance(species[0], list):
                try:
                    species = [Species(spec) for spec in species]
                except ValueError as error:
                    raise ValueError(
                        "Failed to convert the list of csv entries to a species"
                    ) from error
        elif isinstance(species, Species):
            species = [species]
        else:
            raise ValueError(
                "Input must either be (a list of) Species class or csv entries"
            )
        for specie in species:
            if specie in self.get_species_list():
                logging.warning(
                    "Trying to add a duplicate specie, ignoring the new definition and keeping the old one"
                )
                continue
            if specie in reaction_types:
                logging.info(
                    f"Trying to add a reaction type {specie}, reactions are not supposed to be in the species file! "
                )
                continue

            # Filter out reactants that react into nothing, i.e. electrons.
            if specie.name:
                # Filter out ice species with zero binding energy as they tend to give problems
                if specie.name[0] in ["@", "#"] and specie.binding_energy == 0.0:
                    raise ValueError(
                        f"Trying to add an ice specie {specie.name} with zero binding energy, this is not possible. Make sure this specie was added manually to the species file."
                    )
                self._species_dict[specie.name] = specie
            else:
                logging.warning(
                    f"You try to add a falsy specie called '{specie.name}', this cannot be done."
                )

    def remove_species(self, specie_name):
        del self._species_dict[specie_name]

    def get_species_list(self):
        return list(self._species_dict.values())

    def get_species_dict(self):
        return self._species_dict

    def get_specie(self, specie_name: str):
        return deepcopy(self._species_dict[specie_name])

    def set_specie(self, specie_name: str, specie: Species):
        self._species_dict[specie_name] = specie

    def set_species_dict(self, new_species_dict: dict[str, Species]):
        self._species_dict = new_species_dict

    def sort_species(self):
        species_dict = self.get_species_dict()
        logging.debug(
            f"Before sorting species {[(k,v ) for i, (k, v) in enumerate(species_dict.items()) if i < 5]}"
        )
        self.set_species_dict(
            dict(
                sorted(
                    species_dict.items(),
                    # key=lambda kv: (kv[1].get_mass(),),
                    key=lambda kv: custom_lookup[kv[1].name],
                )
            )
        )
        logging.debug(
            f"After sorting species {[(k,v ) for i, (k, v) in enumerate(self.get_species_dict().items()) if i < 5]}"
        )
        assert len(species_dict) == len(
            self.get_species_dict()
        ), "Sorting the species caused a difference in the number of species"
        electron = self.get_specie("E-")
        self.remove_species("E-")
        self.add_species(electron)

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

    def check_and_filter_species(self):
        """Check every speces in network appears in at least one reaction.
        Remove any that do not and alert user.
        """
        species_names = [species.name for species in self.get_species_list()]
        # check for species not involved in any reactions
        lostSpecies = []
        for species in self.get_species_list():

            # keep species that appear in a reaction
            reac_keeps = False
            for reaction in self.get_reaction_list():
                if (
                    species.name in reaction.get_reactants()
                    or species.name in reaction.get_products()
                ):
                    reac_keeps = True
                    break

            # remove the species if it didn't make it into either keep list
            if not (reac_keeps):
                lostSpecies.append(species.name)
        for species in lostSpecies:
            logging.warning(
                f"Trying to remove {species} as it is not present in the reactions"
            )
            self.remove_species(species)

        # then alert user to changes
        if len(lostSpecies) > 0:
            logging.warning(
                "\tSpecies in input list that do not appear in final list:\t"
                + str(lostSpecies)
            )
        else:
            logging.info("\tAll input species in final network")
        logging.debug(f"The network consists of species: {self.get_species_list()}")
        for species in self.get_species_list():
            species.find_constituents()

        # add in pseudo-species to track mantle
        mantle_specs = []
        new_spec = [999] * 7
        new_spec[0] = "BULK"
        mantle_specs.append(Species(new_spec))
        new_spec[0] = "SURFACE"
        mantle_specs.append(Species(new_spec))
        self.add_species(mantle_specs)

    def add_bulk_species(self):
        """For three phase models, MakeRates will produce the version of the species in the bulk
        so that the user doesn't have to endlessly relist the same species
        """
        logging.debug("Adding bulk species")
        speciesNames = [species.name for species in self.get_species_list()]
        userSpecies = [manualSpec.name for manualSpec in self.user_defined_bulk]
        new_species = []
        try:
            h2o_binding_energy = speciesNames.index("#H2O")
            h2o_binding_energy = self.get_species_list()[
                h2o_binding_energy
            ].binding_energy
        except:
            error = "You are trying to create a three phase model but #H2O is not in your network"
            error += "\nThis is likely an error so Makerates will not complete"
            error += (
                "\nTry adding #H2O or switching to three_phase=False in Makerates.py"
            )
            raise RuntimeError(error)
        for species in self.get_species_list():
            if species.is_surface_species():
                if not species.name.replace("#", "@") in speciesNames:
                    new_spec = deepcopy(species)
                    new_spec.name = new_spec.name.replace("#", "@")
                    if new_spec.name in userSpecies:
                        definedBinding = [
                            userSpec.binding_energy
                            for userSpec in self.user_defined_bulk
                            if userSpec.name == new_spec.name
                        ]
                        new_spec.binding_energy = definedBinding[0]
                    else:
                        new_spec.binding_energy = h2o_binding_energy
                    new_species.append(new_spec)
                    logging.debug(
                        f"Adding a bulk partner species for {species}, new {new_spec}"
                    )
        self.add_species(new_species)

    def check_freeze_and_desorbs(self) -> None:
        """`add_freeze_reactions()` and `add_desorb_reactions()` automatically generate
        all desorption and freeze out reactions. However, user may want to change a species on freeze out
        eg C+ becomes #C rather than #C+. This function checks for that and updates species so they'll
        freeze or desorb correctly when reactions are generated.

        LEGACY DOCUMENTATION ?????
        Args:
            species_list (list): list of species objects including all species in network
            reaction_list (list): list of reaction objects including all reactions in network

        Returns:
            list: species and reaction lists with user specified freeze and desorb reactions removed (but species updated)
        """
        desorbs = [
            x for x in self.get_reaction_list() if x.get_reaction_type() == "DESORB"
        ]
        for desorb in desorbs:
            specie = self.get_specie(desorb.get_reactants()[0])
            specie.set_desorb_products(desorb.get_products())
            self.set_specie(desorb.get_reactants()[0], specie)
            # also modify bulk species desorb_products
            bulk_name = "@" + desorb.get_reactants()[0][1:]
            if bulk_name in self.species:
                specie = self.get_specie(bulk_name)
                specie.set_desorb_products(desorb.get_products())
                self.set_specie(bulk_name, specie)

        # for all listed freeze out reactions, add them to correct species
        freezes = [
            x for x in self.get_reaction_list() if x.get_reaction_type() == "FREEZE"
        ]
        for freeze in freezes:
            logging.debug(freeze)
            specie = self.get_specie(freeze.get_reactants()[0])
            specie.set_freeze_products(freeze.get_products(), freeze.get_alpha())
            self.set_specie(freeze.get_reactants()[0], specie)

        # then add default freeze out for species without a listed freeze out
        for species_name, specie in self.get_species_dict().items():
            if (not specie.is_grain_species()) and (
                not specie.get_freeze_products_list()
            ):
                logging.warning(
                    f"Adding a default freezeout for {specie} to the specie"
                )
                specie.add_default_freeze()
                self.set_specie(species_name, specie)

        # WHY: Here we filter all the freeze and desorb reactions for some reason?
        [self.remove_reaction(reaction) for reaction in desorbs + freezes]

    def add_freeze_reactions(self):
        """Save the user effort by automatically generating freeze out reactions"""
        logging.debug("Adding the freeze out reactions!")
        new_reactions = []
        new_species = []
        for species in self.get_species_list():
            logging.debug(f"Checking if {species} needs to have its freezeout added")
            if not species.is_grain_species():
                for products, alpha in species.get_freeze_products():
                    new_reactions.append(
                        Reaction(
                            [species.name, "FREEZE", "NAN"]
                            + products
                            + [alpha, 0.0, species.binding_energy, 0.0, 10000.0]
                        )
                    )
                    # Check if the product is in the species list
                    if products[0] not in self.get_species_list():
                        logging.info(f"Trying to add new specie {products}")
                        new_species.append(
                            Species(
                                [
                                    products[0],
                                    species.mass,
                                    species.binding_energy,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                ]
                            )
                        )
        self.add_reactions(new_reactions)
        self.add_species(new_species)

    def add_desorb_reactions(self):
        """Save the user effort by automatically generating desorption reactions"""
        desorb_reacs = ["DESOH2", "DESCR", "DEUVCR", "THERM"]
        logging.debug("Adding desorbtion reactions!")
        new_reactions = []
        for species in self.get_species_list():
            if species.is_surface_species():
                for reacType in desorb_reacs:
                    new_reactions.append(
                        Reaction(
                            [species.name, reacType, "NAN"]
                            + species.get_desorb_products()
                            + [1, 0, species.binding_energy, 0.0, 10000.0]
                        )
                    )
            if species.is_bulk_species() and not species.is_refractory:
                new_reactions.append(
                    Reaction(
                        [species.name, "THERM", "NAN"]
                        + species.get_desorb_products()
                        + [1, 0, species.binding_energy, 0.0, 10000.0]
                    )
                )
        self.add_reactions(new_reactions)

    def add_chemdes_reactions(self):
        """We have the user list all Langmuir-Hinshelwood and Eley-Rideal
        reactions once. Then we duplicate so that the reaction branches
        with products on grain and products desorbing.
        """
        logging.debug("Adding desascociation reactions for LH and ER mechanisms")
        new_reactions = []
        for reaction in self.get_reaction_list():
            if reaction.get_reaction_type() in ["LH", "ER"]:
                new_reaction = deepcopy(reaction)
                # Convert to disassociation reaction
                new_reactants = new_reaction.get_reactants()
                new_reactants[2] += "DES"
                new_reaction.set_reactants(new_reactants)

                # Replace the species on the grain/bulk with species in gas
                new_products = new_reaction.get_products()
                # Check there are no products incompatible with LH/ER reactions by counting them
                if new_products.count("NAN") + sum(
                    [
                        prod.startswith("#") or prod.startswith("@")
                        for prod in new_products
                    ]
                ) != len(new_products):
                    logging.warning(
                        "All Langmuir-Hinshelwood and Eley-Rideal reactions should be input with products on grains only.\n"
                        + "The fraction of products that enter the gas is dealt with by Makerates and UCLCHEM.\n"
                        + "the following reaction caused this warning\t"
                        + str(reaction)
                    )
                # Replace all grain or bulk products with gas phase counterparts
                new_products = [
                    product.replace("#", "").replace("@", "")
                    for product in new_products
                ]
                new_reaction.set_products(new_products)
                logging.debug(
                    f"Adding deassociation reaction for {reaction}, new reaction {new_reaction}"
                )
                new_reactions.append(new_reaction)
        self.add_reactions(new_reactions)

    def check_for_excited_species(self):
        check = False
        for species in self.get_species_list():
            if "*" in species.name:
                check = True
        return check

    def add_excited_surface_reactions(self) -> None:
        """All excited species will relax to the ground state if they do not react
        the vibrational frequency of the species is used as a pseudo approximation of the rate coefficient
        We assume all grain reactions have an excited variant. For example:
        #A, #B LH #C will have the variants:
        #A*, #B EXSOLID #C  and  #A, #B* EXSOLID #C
        If only one of the reactants in the base reaction has an excited counterpart then
        only one excited version of that reaction is created.
        """
        logging.debug("Adding excited surface reactions")
        excited_species = [x for x in self.get_species_list() if "*" in x.name]
        lh_reactions = [
            x for x in self.get_reaction_list() if "LH" in x.get_reactants()
        ]
        lh_reactions += [
            x for x in self.get_reaction_list() if "LHDES" in x.get_reactants()
        ]

        new_reactions = []

        # add relaxation of excited species
        for spec in excited_species:
            relax_reac = [
                spec.name,
                "EXRELAX",
                "NAN",
                spec.name[:-1],
                "NAN",
                "NAN",
                "NAN",
                1.0,
                0.0,
                0.0,
                0.0,
                10000,
            ]
            new_react = Reaction(relax_reac)
            new_reactions.append(new_react)

        for reaction in lh_reactions:
            # if both #A and #B have excited counterparts
            if reaction.get_reactants()[0] + "*" in [
                specie.name for specie in excited_species
            ] and reaction.get_reactants()[1] + "*" in [
                specie.name for specie in excited_species
            ]:
                new_reac_A_list = [
                    reaction.get_reactants()[0] + "*",
                    reaction.get_reactants()[1],
                    "EXSOLID",
                ]
                new_reac_A_list = (
                    new_reac_A_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000]
                )
                new_reac_B_list = [
                    reaction.get_reactants()[0],
                    reaction.get_reactants()[1] + "*",
                    "EXSOLID",
                ]
                new_reac_B_list = (
                    new_reac_B_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000]
                )

                new_reac_A = Reaction(new_reac_A_list)
                new_reac_B = Reaction(new_reac_B_list)

                # stops duplicate reactions e.g. #H* + #H and #H + #H*
                if new_reac_A != new_reac_B:
                    new_reactions.append(new_reac_A)
                    new_reactions.append(new_reac_B)
                else:
                    new_reactions.append(new_reac_A)

            # if only #A has an excited counterpart
            elif reaction.get_reactants()[0] + "*" in [
                specie.name for specie in excited_species
            ]:
                new_reac_A_list = [
                    reaction.get_reactants()[0] + "*",
                    reaction.get_reactants()[1],
                    "EXSOLID",
                ]
                new_reac_A_list = (
                    new_reac_A_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000]
                )
                new_reac_A = Reaction(new_reac_A_list)
                new_reactions.append(new_reac_A)

            # if only #B has an excited counterpart
            elif reaction.get_reactants()[1] + "*" in [
                specie.name for specie in excited_species
            ]:
                new_reac_B_list = [
                    reaction.get_reactants()[0],
                    reaction.get_reactants()[1] + "*",
                    "EXSOLID",
                ]
                new_reac_B_list = (
                    new_reac_B_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000]
                )
                new_reac_B = Reaction(new_reac_B_list)
                new_reactions.append(new_reac_B)
        self.add_reactions(new_reactions)

    def add_bulk_reactions(self):
        """We assume any reaction that happens on the surface of grains can also happen
        in the bulk (just more slowly due to binding energy). The user therefore only
        lists surface reactions in their input reaction file and we duplicate here.
        """
        logging.debug("Adding bulk reactions")
        current_reaction_list = self.get_reaction_list()
        lh_reactions = [x for x in current_reaction_list if "LH" in x.get_reactants()]
        lh_reactions = lh_reactions + [
            x for x in current_reaction_list if "LHDES" in x.get_reactants()
        ]
        ex_reactions = [x for x in current_reaction_list if "CRS" in x.get_reactants()]
        ex_reactions = ex_reactions + [
            x for x in current_reaction_list if "EXSOLID" in x.get_reactants()
        ]
        ex_reactions = ex_reactions + [
            x for x in current_reaction_list if "EXRELAX" in x.get_reactants()
        ]
        surface_reactions = lh_reactions + ex_reactions

        new_reactions = []
        for reaction in surface_reactions:
            new_reac = deepcopy(reaction)
            new_reac.convert_to_bulk()
            new_reactions.append(new_reac)
        new_reactions = [
            reac for reac in new_reactions if reac not in current_reaction_list
        ]

        bulk_species = [x for x in self.get_species_list() if "@" in x.name]
        for species in bulk_species:
            # add individual swapping
            if not species.is_refractory:
                new_reac_list = [
                    species.name,
                    "BULKSWAP",
                    "NAN",
                    species.name.replace("@", "#"),
                ]
                new_reac_list = new_reac_list + ["NAN", "NAN", "NAN", 1, 0, 0, 0, 10000]
                new_reactions.append(Reaction(new_reac_list))

            # and the reverse
            new_reac_list[0] = species.name.replace("@", "#")
            new_reac_list[1] = "SURFSWAP"
            new_reac_list[3] = species.name
            new_reactions.append(Reaction(new_reac_list))
        logging.debug(
            f"The following bulk reactions are added to the reactions: {new_reactions}"
        )
        self.add_reactions(new_reactions)

    def freeze_checks(self):
        """Check that every species freezes out and alert the user if a
        species freezes out via mutiple routes. This isn't necessarily an
        error so best just print.
        """
        logging.info(
            "\tCheck that species have surface counterparts or if they have multiple freeze outs/check alphas:\n"
        )
        for spec in self.get_species_list():
            if not spec.is_grain_species() and spec.name[-1] not in ["+", "-"]:
                exist_check = 0
                for checkSpeck in self.get_species_list():
                    if checkSpeck.name == "#" + spec.name:
                        exist_check += 1
                if exist_check == 0:
                    logging.warning(
                        f"{spec.name} does not have a surface counterpart in given default species file."
                        + "\n\tThis sets the binding energy to zero, it might cause species conservation errors."
                    )
            freezes = 0
            freezeout_reactions = []
            for reaction in self.get_reaction_list():
                if (spec.name in reaction.get_reactants()) and (
                    "FREEZE" in reaction.get_reactants()
                ):
                    freezes += 1
                    freezeout_reactions.append(reaction)
            if freezes == 1:
                logging.info(
                    f"\t{spec.name} freezes out through {freezeout_reactions[0]}"
                )
            if freezes > 1:
                logging.info(f"\t{spec.name} freezes out through {freezes} routes")
            elif freezes < 1 and not spec.is_grain_species():
                logging.info(f"\t{spec.name} does not freeze out")

    def duplicate_checks(self):
        """
        Check reaction network to make sure no reaction appears twice unless
        they have different temperature ranges.
        """
        logging.info("\tPossible duplicate reactions for manual removal:")
        duplicates = False
        for i, reaction1 in enumerate(self.get_reaction_list()):
            if not reaction1.duplicate:
                for j, reaction2 in enumerate(self.get_reaction_list()):
                    # Save half the checks by only doing half the comparisons
                    if j > i:
                        if reaction1 == reaction2:
                            if not (
                                (reaction1.get_templow() >= reaction2.get_temphigh())
                                or (reaction1.get_temphigh() <= reaction2.get_templow())
                            ):
                                logging.warning(
                                    f"\tReactions {i+1} and {j+1} are possible duplicates\n\t\t"
                                    + str(reaction1)
                                    + f"with temperature range [{reaction1.get_templow()}, {reaction1.get_temphigh()}]"
                                    + "\n\t\t"
                                    + str(reaction2)
                                    + f"with temperature range [{reaction2.get_templow()}, {reaction2.get_temphigh()}]"
                                )
                                duplicates = True
                                # adjust temperatures so temperature ranges are adjacent
                                if reaction1.get_temphigh() > reaction2.get_temphigh():
                                    if (
                                        reaction1.get_templow()
                                        < reaction2.get_temphigh()
                                    ):
                                        logging.warning(
                                            f"\tReactions {reaction1} and {reaction2} have non-adjacent temperature ranges"
                                        )
                                reaction1.duplicate = True
                                reaction2.duplicate = True
        if not duplicates:
            logging.info("\tNone")

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
        for i, reaction in enumerate(self.get_reaction_list()):
            # CO + PHOTON -> O + C
            if ("CO" in reaction.get_reactants()) and (
                "PHOTON" in reaction.get_reactants()
            ):
                if "O" in reaction.get_products() and "C" in reaction.get_products():
                    self.important_reactions["nR_CO_hv"] = i + 1
            # C + PHOTON -> ???
            if ("C" in reaction.get_reactants()) and (
                "PHOTON" in reaction.get_reactants()
            ):
                self.important_reactions["nR_C_hv"] = i + 1
            # H2FORM -> ???
            if "H2FORM" in reaction.get_reactants():
                self.important_reactions["nR_H2Form_CT"] = i + 1
            # H + #H -> H2 AND H + #H -> #H2
            if ("H" in reaction.get_reactants()) and ("#H" in reaction.get_reactants()):
                if "H2" in reaction.get_products():
                    self.important_reactions["nR_H2Form_ERDes"] = i + 1
                elif "#H2" in reaction.get_products():
                    self.important_reactions["nR_H2Form_ER"] = i + 1
            # #H + #H + LH -> ???
            if (reaction.get_reactants().count("#H") == 2) and (
                "LH" in reaction.get_reactants()
            ):
                self.important_reactions["nR_H2Form_LH"] = i + 1
            # #H + #H + LHDES -> ???
            if (reaction.get_reactants().count("#H") == 2) and (
                "LHDES" in reaction.get_reactants()
            ):
                self.important_reactions["nR_H2Form_LHDes"] = i + 1
            # H + FREEZE -> ???
            if ("H" in reaction.get_reactants()) and (
                "FREEZE" in reaction.get_reactants()
            ):
                self.important_reactions["nR_HFreeze"] = i + 1
            # H2 + FREEZE -> ???
            if ("H2" in reaction.get_reactants()) and (
                "FREEZE" in reaction.get_reactants()
            ):
                self.important_reactions["nR_H2Freeze"] = i + 1
            # E- + FREEZE -> ???
            if ("E-" in reaction.get_reactants()) and (
                "FREEZE" in reaction.get_reactants()
            ):
                self.important_reactions["nR_EFreeze"] = i + 1
            # H2 + PHOTON -? ...
            if ("H2" in reaction.get_reactants()) and (
                "PHOTON" in reaction.get_reactants()
            ):
                self.important_reactions["nR_H2_hv"] = i + 1
            # H2 + CRP -> H + H
            if (
                ("H2" in reaction.get_reactants())
                and ("CRP" in reaction.get_reactants())
                and (reaction.get_products().count("H") == 2)
            ):
                self.important_reactions["nR_H2_crp"] = i + 1
        if np_any([value is None for value in self.important_reactions.values()]):
            logging.debug(self.important_reactions)
            missing_reac_error = "Input reaction file is missing mandatory reactions"
            missing_reac_error += (
                "\nH and E- freeze out as well as H2 formation and photodissociation"
            )
            missing_reac_error += " must all be included in user reaction list. Check default_grain_network.csv for example"
            raise RuntimeError(missing_reac_error)

    def index_important_species(self):
        self.species_indices = {}
        names = [species.name for species in self.get_species_list()]
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
                logging.info(f"\t{element} not in network, adding dummy index")
                species_index = len(self.get_species_list()) + 1
            name = "n" + element.lower().replace("+", "x").replace(
                "e-", "elec"
            ).replace("#", "g")
            self.species_indices[name] = species_index

    def __repr__(self):
        return (
            "Reaction network with \nSpecies:\n"
            + ", ".join(map(str, self.get_species_list()))
            + "\nReactions:\n"
            + "\n".join(map(repr, self.get_reaction_list()))
        )


custom_lookup = {
    "H": 0,
    "#H": 1,
    "H+": 2,
    "@H": 3,
    "H2": 4,
    "#H2": 5,
    "H2+": 6,
    "@H2": 7,
    "H3+": 8,
    "HE": 9,
    "#HE": 10,
    "HE+": 11,
    "@HE": 12,
    "HEH+": 13,
    "C+": 14,
    "#C": 15,
    "C": 16,
    "@C": 17,
    "CH": 18,
    "CH+": 19,
    "#CH": 20,
    "@CH": 21,
    "CH2": 22,
    "CH2+": 23,
    "#CH2": 24,
    "N": 25,
    "N+": 26,
    "#N": 27,
    "@CH2": 28,
    "@N": 29,
    "#NH": 30,
    "NH+": 31,
    "NH": 32,
    "#CH3": 33,
    "CH3+": 34,
    "CH3": 35,
    "@NH": 36,
    "@CH3": 37,
    "#CH4": 38,
    "CH4": 39,
    "CH4+": 40,
    "NH2": 41,
    "NH2+": 42,
    "#NH2": 43,
    "O": 44,
    "#O": 45,
    "O+": 46,
    "@CH4": 47,
    "@NH2": 48,
    "@O": 49,
    "#OH": 50,
    "OH+": 51,
    "OH": 52,
    "CH5+": 53,
    "NH3": 54,
    "#NH3": 55,
    "NH3+": 56,
    "@OH": 57,
    "@NH3": 58,
    "#H2O": 59,
    "H2O": 60,
    "H2O+": 61,
    "NH4+": 62,
    "@H2O": 63,
    "H3O+": 64,
    "#MG": 65,
    "#C2": 66,
    "MG+": 67,
    "C2+": 68,
    "C2": 69,
    "MG": 70,
    "@MG": 71,
    "@C2": 72,
    "C2H": 73,
    "C2H+": 74,
    "#C2H": 75,
    "@C2H": 76,
    "CN+": 77,
    "#C2H2": 78,
    "C2H2": 79,
    "C2H2+": 80,
    "CN": 81,
    "#CN": 82,
    "@C2H2": 83,
    "@CN": 84,
    "#C2H3": 85,
    "HNC": 86,
    "HCN+": 87,
    "#HCN": 88,
    "C2H3": 89,
    "#HNC": 90,
    "HCN": 91,
    "@C2H3": 92,
    "@HCN": 93,
    "@HNC": 94,
    "SI+": 95,
    "#H2CN": 96,
    "SI": 97,
    "C2H4": 98,
    "#C2H4": 99,
    "N2+": 100,
    "N2": 101,
    "HCNH+": 102,
    "H2CN": 103,
    "CO+": 104,
    "CO": 105,
    "#N2": 106,
    "#CO": 107,
    "#SI": 108,
    "@H2CN": 109,
    "@C2H4": 110,
    "@N2": 111,
    "@CO": 112,
    "@SI": 113,
    "#SIH": 114,
    "#C2H5": 115,
    "C2H5": 116,
    "SIH": 117,
    "SIH+": 118,
    "HOC+": 119,
    "#HCO": 120,
    "HCO+": 121,
    "HCO": 122,
    "N2H+": 123,
    "@SIH": 124,
    "@C2H5": 125,
    "@HCO": 126,
    "#NO": 127,
    "#H2CO": 128,
    "H2CO": 129,
    "H2CO+": 130,
    "NO": 131,
    "NO+": 132,
    "SIH2": 133,
    "SIH2+": 134,
    "#SIH2": 135,
    "@NO": 136,
    "@H2CO": 137,
    "@SIH2": 138,
    "CH2OH": 139,
    "#CH2OH": 140,
    "#HNO": 141,
    "#SIH3": 142,
    "SIH3": 143,
    "HNO+": 144,
    "HNO": 145,
    "H3CO": 146,
    "H3CO+": 147,
    "SIH3+": 148,
    "#H3CO": 149,
    "@CH2OH": 150,
    "@HNO": 151,
    "@SIH3": 152,
    "@H3CO": 153,
    "SIH4": 154,
    "S+": 155,
    "#S": 156,
    "SIH4+": 157,
    "O2+": 158,
    "#CH3OH": 159,
    "#SIH4": 160,
    "#O2": 161,
    "S": 162,
    "H2NO+": 163,
    "O2": 164,
    "CH3OH": 165,
    "@S": 166,
    "@CH3OH": 167,
    "@SIH4": 168,
    "@O2": 169,
    "HS": 170,
    "#HS": 171,
    "CH3OH2+": 172,
    "SIH5+": 173,
    "O2H+": 174,
    "HS+": 175,
    "O2H": 176,
    "#O2H": 177,
    "@HS": 178,
    "@O2H": 179,
    "#H2S": 180,
    "H2S": 181,
    "H2S+": 182,
    "@H2S": 183,
    "CL": 184,
    "CL+": 185,
    "#CL": 186,
    "H3S+": 187,
    "@CL": 188,
    "HCL+": 189,
    "HCL": 190,
    "C3+": 191,
    "#HCL": 192,
    "@HCL": 193,
    "H2CL+": 194,
    "C2N": 195,
    "#C2N": 196,
    "C2N+": 197,
    "C3H2": 198,
    "#C3H2": 199,
    "@C2N": 200,
    "@C3H2": 201,
    "C2NH+": 202,
    "#SIC": 203,
    "CH3CCH": 204,
    "#CH3CCH": 205,
    "SIC": 206,
    "SIC+": 207,
    "@SIC": 208,
    "@CH3CCH": 209,
    "#CH3CN": 210,
    "C3H5+": 211,
    "CH3CN": 212,
    "@CH3CN": 213,
    "CH2CO": 214,
    "CH3CNH": 215,
    "CH3CNH+": 216,
    "OCN": 217,
    "#OCN": 218,
    "#CH2CO": 219,
    "#CH3CNH": 220,
    "@OCN": 221,
    "@CH2CO": 222,
    "@CH3CNH": 223,
    "#HNCO": 224,
    "HNCO": 225,
    "@HNCO": 226,
    "#CH3CHO": 227,
    "CH3CHO": 228,
    "SIO+": 229,
    "CS": 230,
    "CS+": 231,
    "CO2": 232,
    "#CS": 233,
    "#SIO": 234,
    "SIO": 235,
    "#CO2": 236,
    "@CH3CHO": 237,
    "@CS": 238,
    "@SIO": 239,
    "@CO2": 240,
    "#NH2CHO": 241,
    "#HCOO": 242,
    "HCOO": 243,
    "NH2CHO": 244,
    "SIOH+": 245,
    "HCO2+": 246,
    "HCS+": 247,
    "#HSIO": 248,
    "HCS": 249,
    "#HCS": 250,
    "@NH2CHO": 251,
    "@HCOO": 252,
    "@HSIO": 253,
    "@HCS": 254,
    "#H2SIO": 255,
    "#H2CS": 256,
    "#NS": 257,
    "#HCOOH": 258,
    "HCOOH": 259,
    "#NO2": 260,
    "NS": 261,
    "NO2": 262,
    "H2CS": 263,
    "H2SIO": 264,
    "H2CS+": 265,
    "NS+": 266,
    "@H2SIO": 267,
    "@H2CS": 268,
    "@NS": 269,
    "@HCOOH": 270,
    "@NO2": 271,
    "H3CS+": 272,
    "HNS+": 273,
    "SO": 274,
    "SO+": 275,
    "#SO": 276,
    "@SO": 277,
    "#C4H": 278,
    "C4H": 279,
    "HSO+": 280,
    "@C4H": 281,
    "C3N": 282,
    "#C3N": 283,
    "@C3N": 284,
    "#HC3N": 285,
    "HC3N": 286,
    "@HC3N": 287,
    "#SIC2": 288,
    "SIC2+": 289,
    "C2N2+": 290,
    "SIC2": 291,
    "NCCN": 292,
    "#NCCN": 293,
    "@SIC2": 294,
    "@NCCN": 295,
    "SIS+": 296,
    "#OCS": 297,
    "#SIS": 298,
    "OCS": 299,
    "OCS+": 300,
    "SIS": 301,
    "@OCS": 302,
    "@SIS": 303,
    "HOCS+": 304,
    "HSIS+": 305,
    "C4N": 306,
    "C4N+": 307,
    "#C4N": 308,
    "@C4N": 309,
    "SO2": 310,
    "S2+": 311,
    "SO2+": 312,
    "S2": 313,
    "SIC3": 314,
    "SIC3+": 315,
    "#SIC3": 316,
    "#S2": 317,
    "#SO2": 318,
    "@SIC3": 319,
    "@S2": 320,
    "@SO2": 321,
    "HS2": 322,
    "HS2+": 323,
    "HSO2+": 324,
    "#HS2": 325,
    "@HS2": 326,
    "H2S2+": 327,
    "H2S2": 328,
    "#H2S2": 329,
    "@H2S2": 330,
    "E-": 331,
    "BULK": 332,
    "SURFACE": 333,
    "HSIO": 334,
}
