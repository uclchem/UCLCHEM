"""
This python file contains all functions for de-duplicating species and reaction lists,
checking for common errors, and automatic addition of reactions such as freeze out,
desorption and bulk reactions for three phase models.
"""

import logging
import sys
from copy import deepcopy
from pathlib import Path
from typing import Union

from numpy import any as np_any

from .reaction import CoupledReaction, Reaction, reaction_types
from .species import Species, elementList


class Network:
    """The network class stores all the information about reaction network."""

    def __init__(
        self,
        species: list[Species],
        reactions: list[Reaction],
        user_defined_bulk: list = [],
        gas_phase_extrapolation: bool = False,
        add_crp_photo_to_grain: bool = False,
    ):
        """A class to store network information such as indices of important reactions.

        The class fully utilizes getters and setters, which can be used to add/remove
        reactions and the species involved. Important is that you do not directly edit
        the internal dictionaries that store the species and reactions, unless you
        know what you are doing. The network by default checks for duplicates in species
        and identical reactions that overlap in temperature ranges, potentially causing
        problems.

        Args:
            species (list[Species]): A list of chemical species that are added to the network
            reactions (list[Reaction]): A list of chemical reactions that are added to the network
            user_defined_bulk (list, optional): List of user defined bulk. Defaults to [].
            add_crp_photo_to_grain (bool, optional): Whether to add CRP, CRPHOT and PHOTON reactions from gas-phase into solid phase too.
        """
        assert len(set([s.name for s in species])) == len(species), (
            "Cannot have duplicate species in the species list."
        )
        self.set_species_dict({s.name: s for s in species})
        self.excited_species = self.check_for_excited_species()
        self.user_defined_bulk = user_defined_bulk
        self.add_crp_photo_to_grain = add_crp_photo_to_grain
        electron_specie = Species(["E-", 0, 0.0, 0, 0, 0, 0])
        electron_specie.n_atoms = 1
        self.add_species(electron_specie)
        self.set_reaction_dict({k: v for k, v in enumerate(reactions)})

        #### Add reactions and species   ####
        # check which species are changed on freeze or desorb
        self.check_freeze_and_desorbs()

        # Need additional grain reactions including non-thermal desorption and chemically induced desorption
        self.add_freeze_reactions()
        if self.add_crp_photo_to_grain:
            self.add_CRP_and_PHOTO_reactions_to_grain()
        self.add_bulk_species()
        self.add_bulk_reactions()
        self.add_desorb_reactions()
        self.add_chemdes_reactions()
        if self.excited_species:
            self.add_excited_surface_reactions()

        # Ensure that the branching ratios are correct, if not, edit the network to enforce it.
        self.branching_ratios_checks()

        # Extrapolate Gas phase reactions if needed:
        if gas_phase_extrapolation:
            self.add_gas_phase_extrapolation()

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
        return self.get_reaction_list()

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
        return self.get_species_list()

    @species_list.setter
    def species_list(self, value):
        raise Exception(
            "Do not set species lists explicitely, use the add and remove interfaces for reactions"
        )

    def add_reactions(
        self, reactions: Union[Union[Reaction, str], list[Union[Reaction, str]]]
    ):
        """Add a reaction, list of inputs to the Reaction class or list of reactions to the network.

        Args:
            reactions (Union[Union[Reaction, str], list[Union[Reaction, str]]]): Reaction or list or reactions

        """
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
        logging.debug(
            f"n_reactions {len(current_reaction_list)} before adding new reactions to the internal dict"
        )
        for reaction in reactions:
            if reaction in current_reaction_list:
                # See if we have a collision with the any reactions with identical reactants and
                # products, but different temperature ranges to avoid double definitions.
                similar_reactions = self.find_similar_reactions(reaction)
                for similar_reaction_key in similar_reactions:
                    if reaction.check_temperature_collision(
                        similar_reactions[similar_reaction_key]
                    ):
                        raise RuntimeError(
                            f"There already is a {reaction} present that has overlapping temperature ranges. Check the reaction set."
                        )

            # quick check to make sure all species in the reaction are in the species list.
            species_not_present = [
                spec
                for spec in reaction.get_reactants() + reaction.get_products()
                if (spec not in self.get_species_list())
                and (spec not in ["NAN", "", "E-"] + reaction_types)
            ]
            for specie in species_not_present:
                logging.debug(f"Trying to add specie {specie}")
                # TODO: get more sensible mass
                self.add_species(Species([specie, -1, 0.0, 0.0, 0.0, 0.0, 0.0]))
            # Index and add the new reaction.
            new_idx = max(list(self._reactions_dict.keys())) + 1
            if new_idx in self._reactions_dict.keys():
                raise ValueError(
                    f"Makerates is trying to add a reaction with index {new_idx}, but something is already there, please report this to the developers."
                )
            self._reactions_dict[new_idx] = reaction
        logging.debug(
            f"n_reactions {len(current_reaction_list)} after adding them to the internal dict"
        )

    def find_similar_reactions(self, reaction: Reaction) -> dict[int, Reaction]:
        """Reactions are similar if the reaction has the same reactants and products,
        find all reactions that are similar, returning their index and the reaction itself.

        Args:
            reaction (Reaction): Reaction with possible identical (but for temperature range) reactions in the network

        Returns:
            dict[int, Reaction]: A dict with the identical reactions.
        """
        return {k: v for k, v in self._reactions_dict.items() if v == reaction}

    def get_reaction_index(self, reaction: Reaction) -> int:
        """Get the index of a reaction in the internal _reactions_dict.

        Args:
            reaction (Reaction): The reaction to find the index of

        Returns:
            int: The index of the reaction in the internal _reactions_dict
        """
        similar_reactions = self.find_similar_reactions(reaction)
        if len(similar_reactions) == 1:
            return list(similar_reactions.keys())[0]
        elif len(similar_reactions) == 0:
            raise ValueError(f"The reaction {reaction} is not present in the network.")
        else:
            raise RuntimeError(
                f"Found more than one index for the reaction {reaction}, cannot uniquely identify it; use find_similar_reactions instead."
            )

    def remove_reaction_by_index(self, reaction_idx: int) -> None:
        """Remove a reaction by its index in the internal _reactions_dict, this is the only way
        to remove reactions that are defined piecewise across temperature ranges.

        Args:
            reaction_idx (int): Index of the reaction to remove
        """
        del self._reactions_dict[reaction_idx]

    def remove_reaction(self, reaction: Reaction) -> None:
        """Remove the reaction by giving the object itself, this only works if the reaction is
        not piecewise defined across the temperature ranges.

        Args:
            reaction (Reaction): The reaction you wish to delete.
        """
        # In this reaction we use equality as defined for reactions, this does as of now not include
        # checking temperature ranges, so more than one key could be returned.
        #
        # The find_similar_reaction returns a dict of (index[int], reaction[Reaction]),
        # We make it a list of tuples as it is easier to index and manipulate for this case.
        reaction_idx_dict_as_tuples = list(
            self.find_similar_reactions(reaction).items()
        )
        if len(reaction_idx_dict_as_tuples) == 1:
            reac_idx, reac_value = reaction_idx_dict_as_tuples[0]
            logging.debug(f"Trying to remove index: {reac_idx}: {reac_value} ")
            # Remove the reaction with the index from the reaction set:
            del self._reactions_dict[reac_idx]

            if not isinstance(reaction, CoupledReaction):
                to_pop = []
                for k, v in self._reactions_dict.items():
                    if not isinstance(v, CoupledReaction):
                        continue
                    if v.get_partner() == reaction:
                        logging.debug(f"Coupled reaction {v} will also be removed")
                        to_pop.append(v)
                if to_pop:
                    [self.remove_reaction(coupled_reac) for coupled_reac in to_pop]

        elif len(reaction_idx_dict_as_tuples) == 0:
            logging.warning(
                f"The reaction {reaction} is not present in the reaction set, so cannot remove it"
            )
        elif len(reaction_idx_dict_as_tuples) > 1:
            raise (
                RuntimeError(
                    "found more than one indices for the reaction {reaction}, remove by index instead of by reaction."
                )
            )

    def change_reaction_barrier(self, reaction: Reaction, barrier: float) -> None:
        reaction_idx_dict_as_tuples = list(
            self.find_similar_reactions(reaction).items()
        )
        if len(reaction_idx_dict_as_tuples) == 1:
            reac_idx, reac_value = reaction_idx_dict_as_tuples[0]
            logging.debug(
                f"Trying to change barrier of reaction index: {reac_idx}: {reac_value} "
            )
            # Change the barrier of the reaction with the index from the reaction set:
            self._reactions_dict[reac_idx].set_gamma(barrier)

            if not isinstance(reaction, CoupledReaction):
                to_change = []
                for k, v in self._reactions_dict.items():
                    if not isinstance(v, CoupledReaction):
                        continue
                    if v.get_partner() == reaction:
                        logging.debug(f"Coupled reaction {v} will also be removed")
                        to_change.append(v)
                if to_change:
                    [
                        self.change_reaction_barrier(coupled_reac, barrier)
                        for coupled_reac in to_change
                    ]

        elif len(reaction_idx_dict_as_tuples) == 0:
            logging.warning(
                f"The reaction {reaction} is not present in the reaction set, so cannot change its barrier it"
            )
        elif len(reaction_idx_dict_as_tuples) > 1:
            raise (
                RuntimeError("found more than one indices for the reaction {reaction}.")
            )

    def change_binding_energy(self, specie: str, new_binding_energy: float) -> None:
        all_species = self.get_species_list()
        all_species_names = [specie.get_name() for specie in all_species]
        if specie not in all_species_names:
            error = f"Specie {specie} was not found in the network while attempting to change its binding energy."
            raise ValueError(error)
        old_bulk_h2o_binding_energy = all_species[all_species_names.index("@H2O")]
        old_bulk_h2o_binding_energy = old_bulk_h2o_binding_energy.binding_energy
        if specie == "@H2O":
            # If specie is bulk H2O, we need to change binding energies of all other bulk species,
            # as the diffusion is limited by diffusion of bulk H2O. (Ghesquiere 2015)
            for specie_in_network in all_species:
                if (
                    "@" in specie_in_network.get_name()
                    and specie_in_network.binding_energy == old_bulk_h2o_binding_energy
                ):
                    # If the specie had a different bulk binding energy, do not change it, as it was user specified.
                    specie_in_network.binding_energy = new_binding_energy
            return
        if "@" in specie:
            if (
                all_species[all_species_names.index(specie)].binding_energy
                == old_bulk_h2o_binding_energy
            ):
                # If the bulk species has the same binding energy as bulk H2O,
                # but we are trying to change it directly, give user a warning.
                print(
                    f"WARNING: ATTEMPTING TO CHANGE BINDING ENERGY OF BULK SPECIE {specie} THAT WAS PREVIOUSLY @H2O BINDING ENERGY LIMITED"
                )
        for specie_in_network in all_species:
            if specie_in_network.get_name() == specie:
                specie_in_network.binding_energy = new_binding_energy
                return

    def get_reaction(self, reaction_idx: int) -> Reaction:
        """Obtain a reaction from the reaction set given an index of the internal _reactions_dict.

        Args:
            reaction_idx (int): The reaction index

        Returns:
            Reaction: the desired reaction
        """
        return deepcopy(self._reactions_dict[reaction_idx])

    def set_reaction(self, reaction_idx: int, reaction: Reaction) -> None:
        """This setter explicitely sets the reaction for a certain index.

        Args:
            reaction_idx (int): The index to be written to
            reaction (Reaction): The reaction to be added to the index.
        """
        old_length = len(self._reactions_dict)
        self._reactions_dict[reaction_idx] = reaction
        assert old_length == len(self._reactions_dict), (
            "Setting the reaction caused a change in the number of reactions, use add_reaction and remove_reaction for add and remove operations."
        )

    def get_reaction_dict(self) -> dict[int, Reaction]:
        """Returns the whole internal reaction dictionary.

        Returns:
            dict[int, Reaction]: A copy of the internal reactions dictionary.
        """
        return deepcopy(self._reactions_dict)

    def set_reaction_dict(self, new_dict: dict[int, Reaction]) -> None:
        """Override the reactions dictionary with a new dictionar.

        Args:
            new_dict (dict[int, Reaction]): The new reactions_dictionary.
        """
        self._reactions_dict = new_dict

    def get_reaction_list(self) -> list[Reaction]:
        """Obtain all the reactions in the Network.

        Returns:
            list[Reaction]: A list with all the reaction objects
        """
        return list(self._reactions_dict.values())

    def get_reactions_by_types(
        self, reaction_type: Union[str, list[str]]
    ) -> list[Reaction]:
        """Get the union of all reactions of a certain type.

        Args:
            reaction_type (str): The reaction type to filter on

        Returns:
            list[Reaction]: A list of reactions of the specified type
        """
        if isinstance(reaction_type, str):
            reaction_type = [reaction_type]
        return [
            r
            for r in self.get_reaction_list()
            if (r.get_reaction_type() in reaction_type)
        ]

    def sort_reactions(self) -> None:
        """Sort the reaction dictionary by reaction type first and by the first reactant second."""
        # """Sort the reactions by the reaction type first, reactants second."""
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
            f"After sorting reactions {[(k, v) for i, (k, v) in enumerate(self.get_reaction_dict().items()) if i < 5]}"
        )
        assert len(reaction_dict) == len(self.get_reaction_dict()), (
            "Sorting the species caused a difference in the number of species"
        )

    def add_species(
        self, species: Union[Union[Species, str], list[Union[Species, str]]]
    ):
        """Add species to the network, given a (list of) species. If it is a list of strings,
        it tries to instantiate a species class with it. It also checks for duplicate entries  and
        filters out attempts to add reaction types to the species.

        Args:
            species (Union[Union[Species, str], list[Union[Species, str]]]): A (list of) species or strings.

        Raises:
            ValueError: If we cannot parse the (list of) reactions
            ValueError: If an ice specie with binding energy of zero is added.
        """
        if isinstance(species, list):
            if len(species) == 0:
                logging.warning(
                    "Tried to add a list of species, but it was empty, ignoring."
                )
            elif isinstance(species[0], Species):
                # if it is a list of Species, no action is needed.
                pass
            elif isinstance(species[0], list):
                try:
                    species = [Species(spec) for spec in species]
                except ValueError as error:
                    raise ValueError(
                        "Failed to convert the list of csv style entries to a species"
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
                    specie.binding_energy = 5600.0  # Default to water binding energy
                    logging.warning(
                        f"Setting binding energy of ice specie {specie.name} to default of 5600K as it was zero."
                    )
                    # raise ValueError(
                    #     f"Trying to add an ice specie {specie.name} with zero binding energy, this is not possible. Make sure this specie was added manually to the species file."
                    # )
                self._species_dict[specie.name] = specie
            else:
                logging.info(
                    f"You try to add a falsy specie called '{specie.name}', this cannot be done and will be ignored."
                )

    def remove_species(self, specie_name: str) -> None:
        """Remove a specie from the network

        Args:
            specie_name (str): Species to remove
        """
        del self._species_dict[specie_name]

    def get_species_list(self) -> list[Species]:
        """Obtain a list with all the species in the network

        Returns:
            list[Species]: A list of all the species in the reaction network
        """
        return list(self._species_dict.values())

    def get_species_dict(self) -> dict[str, Species]:
        """Get the internal dictionary that stores all the species, it consists
        of all species' names as key, with the species object as value.

        Returns:
            dict[str, Species]: A dictionary with the species
        """
        return deepcopy(self._species_dict)

    def get_specie(self, specie_name: str) -> Species:
        """Get the species of the reaction network (from the internal dictionary)

        Args:
            specie_name (str): the name of the species as a string

        Returns:
            Species: The species object
        """
        return deepcopy(self._species_dict[specie_name])

    def set_specie(self, species_name: str, species: Species) -> None:
        """Set the species of the reaction network in the internal dictionary

        Args:
            species_name (str): The name of the species as string
            species (Species): The Species object to set
        """
        self._species_dict[species_name] = species

    def set_species_dict(self, new_species_dict: dict[str, Species]) -> None:
        """Set the internal species dict

        Args:
            new_species_dict (dict[str, Species]): The new dictionary to set
        """
        self._species_dict = new_species_dict

    def sort_species(self) -> None:
        """Sort the species based on their mass in ascending order. We always make sure the Electron is last."""
        species_dict = self.get_species_dict()
        logging.debug(
            f"Before sorting species {[(k, v) for i, (k, v) in enumerate(species_dict.items()) if i < 5]}"
        )

        self.set_species_dict(
            dict(
                sorted(
                    species_dict.items(),
                    # key=lambda kv: (kv[1].get_mass(),),
                    key=lambda kv: (
                        kv[1].is_grain_species(),
                        kv[1].is_bulk_species(),
                        kv[1].get_mass(),
                    ),
                )
            )
        )
        logging.debug(
            f"After sorting species {[(k, v) for i, (k, v) in enumerate(self.get_species_dict().items()) if i < 5]}"
        )
        assert len(species_dict) == len(self.get_species_dict()), (
            "Sorting the species caused a difference in the number of species"
        )
        electron = self.get_specie("E-")
        self.remove_species("E-")
        self.add_species(electron)

    # check reactions to alert user of potential issues including repeat reactions
    # and multiple freeze out routes
    def check_network(self) -> None:
        """Run through the list of reactions and check for obvious errors such
        as duplicate reactions, multiple freeze out routes (to warn, not necessarily
        an error), etc.
        """
        self.freeze_checks()
        self.duplicate_checks()
        self.index_important_reactions()
        self.index_important_species()

    def check_and_filter_species(self) -> None:
        """Check every speces in network appears in at least one reaction.
        Remove any that do not and alert user.
        """
        # species_names = [species.name for species in self.get_species_list()]
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

    def add_bulk_species(self) -> None:
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
        except ValueError:
            error = "You are trying to create a three phase model but #H2O is not in your network"
            error += "\nThis is likely an error so Makerates will not complete. Try adding #H2O"
            raise RuntimeError(error)
        for species in self.get_species_list():
            if species.is_surface_species():
                if species.name.replace("#", "@") not in speciesNames:
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
                logging.info(f"Adding a default freezeout for {specie} to the specie")
                specie.add_default_freeze()
                self.set_specie(species_name, specie)

        # Here we filter all the freeze and desorb reactions in order to avoid duplicates
        [self.remove_reaction(reaction) for reaction in desorbs + freezes]

    def add_freeze_reactions(self) -> None:
        """Save the user effort by automatically generating freeze out reactions"""
        logging.debug("Adding the freeze out reactions!")
        new_reactions = []
        new_species = []
        for species in self.get_species_list():
            logging.debug(f"Checking if {species} needs to have its freezeout added")
            if not species.is_grain_species():
                for products, alpha in species.get_freeze_products():
                    if species.name == "E-":
                        # Set electron freeze out to zero:
                        alpha = 0.0
                    new_reactions.append(
                        Reaction(
                            [species.name, "FREEZE", "NAN"]
                            + products
                            + [alpha, 0.0, species.binding_energy, 0.0, 10000.0, 0.0]
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
        if new_reactions:
            self.add_reactions(new_reactions)
        if new_species:
            self.add_species(new_species)

    def add_desorb_reactions(self) -> None:
        """Save the user effort by automatically generating desorption reactions"""
        desorb_reacs = ["DESOH2", "DESCR", "DEUVCR", "THERM"]
        logging.debug("Adding desorption reactions!")
        new_reactions = []
        for species in self.get_species_list():
            if species.is_surface_species():
                for reacType in desorb_reacs:
                    new_reactions.append(
                        Reaction(
                            [species.name, reacType, "NAN"]
                            + species.get_desorb_products()
                            + [1, 0, species.binding_energy, 0.0, 10000.0, 0.0]
                        )
                    )
            if species.is_bulk_species() and not species.is_refractory:
                new_reactions.append(
                    Reaction(
                        [species.name, "THERM", "NAN"]
                        + species.get_desorb_products()
                        + [1, 0, species.binding_energy, 0.0, 10000.0, 0.0]
                    )
                )
        self.add_reactions(new_reactions)

    def add_chemdes_reactions(self) -> None:
        """We have the user list all Langmuir-Hinshelwood and Eley-Rideal
        reactions once. Then we duplicate so that the reaction branches
        with products on grain and products desorbing.
        """
        logging.debug("Adding desorption reactions for LH and ER mechanisms")
        new_reactions = []
        existing_desorption_reactions = [
            x
            for x in self.get_reaction_list()
            if x.get_reaction_type() in ["LHDES", "ERDES"]
        ]
        for reaction in self.get_reaction_list():
            if reaction.get_reaction_type() in ["LH", "ER"]:
                # If either the LH or ER reaction already has a desorption reaction, skip it.
                if any(
                    [
                        (
                            existing_reaction.get_reaction_type() + "DES"
                            == reaction.get_reaction_type()
                        )
                        and (
                            existing_reaction.get_pure_reactants()
                            == reaction.get_pure_reactants()
                        )
                        for existing_reaction in existing_desorption_reactions
                    ]
                ):
                    logging.warning(
                        f"We were trying to add an automatic desorb reaction for {reaction}, but it already exists in the network, so skipping it."
                    )
                    continue
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
                    f"Adding desorption reaction for {reaction}, new reaction {new_reaction}"
                )

                new_reaction = CoupledReaction(new_reaction)
                while isinstance(reaction, CoupledReaction):
                    # If the current loop reaction is also coupled, get its partner.
                    reaction = reaction.get_partner()
                new_reaction.set_partner(reaction)

                new_reactions.append(new_reaction)
        self.add_reactions(new_reactions)

    def check_for_excited_species(self) -> bool:
        """Check if there are any exicted species in the network, true if there are any."""
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
                0.0,
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
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
                )
                new_reac_B_list = [
                    reaction.get_reactants()[0],
                    reaction.get_reactants()[1] + "*",
                    "EXSOLID",
                ]
                new_reac_B_list = (
                    new_reac_B_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
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
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
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
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
                )
                new_reac_B = Reaction(new_reac_B_list)
                new_reactions.append(new_reac_B)
        self.add_reactions(new_reactions)

    def add_bulk_reactions(self) -> None:
        """We assume any reaction that happens on the surface of grains can also happen
        in the bulk (just more slowly due to binding energy). The user therefore only
        lists surface reactions in their input reaction file and we duplicate here.
        """
        logging.debug("Adding bulk reactions")
        surface_reactions = self.get_reactions_on_grain()
        bulk_reaction_types = ["CRP", "CRPHOT", "PHOTON", "LH", "EXSOLID", "EXRELAX"]
        surface_reactions_can_be_bulk = [
            reaction
            for reaction in surface_reactions
            if reaction.get_reaction_type() in bulk_reaction_types
        ]
        current_reaction_list = self.get_reaction_list()

        new_reactions = []
        for reaction in surface_reactions_can_be_bulk:
            new_reac = deepcopy(reaction)
            new_reac.convert_surf_to_bulk()
            new_reac = CoupledReaction(new_reac)
            while isinstance(reaction, CoupledReaction):
                # If the current loop reaction is also coupled, get its partner.
                reaction = reaction.get_partner()
            new_reac.set_partner(reaction)
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
                new_reac_list = new_reac_list + [
                    "NAN",
                    "NAN",
                    "NAN",
                    1,
                    0,
                    0,
                    0,
                    10000,
                    0.0,
                ]
                new_reactions.append(Reaction(new_reac_list))

            # and the reverse, going from surface to bulk
            if not species in [
                "@H2"
            ]:  # If species is H2, do not allow it to go from surface to bulk
                new_reac_list[0] = species.name.replace("@", "#")
                new_reac_list[1] = "SURFSWAP"
                new_reac_list[3] = species.name
                new_reactions.append(Reaction(new_reac_list))
        logging.debug(
            f"The following bulk reactions are added to the reactions: {new_reactions}"
        )
        self.add_reactions(new_reactions)

    def freeze_checks(self) -> None:
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

    def get_reactions_on_grain(self) -> list[Reaction]:
        reactions_on_grain = []
        for reaction in self.get_reaction_list():
            reactants = reaction.get_reactants()
            if any("#" in reactant or "@" in reactant for reactant in reactants):
                reactions_on_grain.append(reaction)
        return reactions_on_grain

    def add_CRP_and_PHOTO_reactions_to_grain(self) -> None:
        """Add all the gas-phase reactions with CRP, CRPHOT or PHOTON to the grain surface too"""
        logging.info("Adding gas-phase reactions with CRP, CRPHOT or PHOTON to grain")
        reactions_on_grain = self.get_reactions_on_grain()
        reactions_on_grain_filtered = [
            reaction
            for reaction in reactions_on_grain
            if reaction.get_reaction_type() in ["CRP", "CRPHOT", "PHOTON"]
        ]
        new_reactions = []
        for reaction in self.get_reaction_list():
            if not reaction.get_reaction_type() in ["CRP", "CRPHOT", "PHOTON"]:
                continue
            if reaction in reactions_on_grain_filtered:
                continue
            reactants = reaction.get_reactants()
            products = reaction.get_products()
            if any(
                "@" + reactants[0] in reaction.get_reactants()
                or "#" + reactants[0] in reaction.get_reactants()
                for reaction in reactions_on_grain_filtered
            ):
                # There is already another version in the network, keep that
                continue

            if any("+" in reactant for reactant in reactants):
                logging.debug(
                    f"Reaction {reaction} had reactant ions, which is not possible on grain. Skipping"
                )
                # Cannot have ions in grain
                continue
            # We have now filtered to have only gas-phase reactions that have type CRP, CRPHOT or PHOTON
            new_reaction = deepcopy(reaction)
            new_reaction.convert_gas_to_surf()
            reactants, products = (
                new_reaction.get_reactants(),
                new_reaction.get_products(),
            )
            if reactants[0] == products[0]:
                # This means the reaction was simply an ionization reaction. We can skip this reaction.
                logging.debug(
                    f"Reaction {reaction} is ionization reaction, skip on grain surface"
                )
                continue
            if not all(
                species in self.get_species_list()
                or species in ["NAN", "CRP", "CRPHOT", "PHOTON"]
                for species in reactants + products
            ):
                # Reaction contains species that were not defined on grain surface.
                # Do not add this reaction to the network.
                logging.debug(
                    f"Reaction {reaction} contains species that were not set on grain. Skipping"
                )
                continue
            new_reactions.append(new_reaction)
        logging.debug(f"Adding new reactions to grain")
        self.add_reactions(new_reactions)
        logging.info(f"Added {len(new_reactions)} reactions to grain")

    def duplicate_checks(self) -> None:
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
                                if (
                                    reaction1.get_source() == reaction2.get_source()
                                    and reaction1.get_source() == "UMIST"
                                ):
                                    logging.info(
                                        f"Detected overlapping UMIST reactions {reaction1} wit indices {i + 1} {j + 1}, this is done in UMIST to provide better rates. "
                                    )
                                else:
                                    logging.warning(
                                        f"\tReactions with indices {i + 1} and {j + 1} are possible duplicates\n\t\t"
                                        + str(reaction1)
                                        + f" with temperature range [{reaction1.get_templow()}, {reaction1.get_temphigh()}] and source {reaction1.get_source()}"
                                        + "\n\t\t"
                                        + str(reaction2)
                                        + f" with temperature range [{reaction2.get_templow()}, {reaction2.get_temphigh()}] and source {reaction2.get_source()}"
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

    def index_important_reactions(self) -> None:
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
            reacs = reaction.get_reactants()
            prods = reaction.get_products()
            reaction_filters = {
                "nR_CO_hv": lambda reacs, prods: ("CO" in reacs)
                and ("PHOTON" in reacs)
                and ("O" in prods)
                and ("C" in prods),
                "nR_C_hv": lambda reacs, prods: ("C" in reacs) and ("PHOTON" in reacs),
                "nR_H2Form_CT": lambda reacs, prods: "H2FORM" in reacs,
                "nR_H2Form_ERDes": lambda reacs, prods: ("H" in reacs)
                and ("#H" in reacs)
                and ("H2" in prods),
                "nR_H2Form_ER": lambda reacs, prods: ("H" in reacs)
                and ("#H" in reacs)
                and ("#H2" in prods),
                "nR_H2Form_LH": lambda reacs, prods: (reacs.count("#H") == 2)
                and ("LH" in reacs),
                "nR_H2Form_LHDes": lambda reacs, prods: (reacs.count("#H") == 2)
                and ("LHDES" in reacs),
                "nR_HFreeze": lambda reacs, prods: ("H" in reacs)
                and ("FREEZE" in reacs),
                "nR_H2Freeze": lambda reacs, prods: ("H2" in reacs)
                and ("FREEZE" in reacs),
                "nR_EFreeze": lambda reacs, prods: ("E-" in reacs)
                and ("FREEZE" in reacs),
                "nR_H2_hv": lambda reacs, prods: ("H2" in reacs)
                and ("PHOTON" in reacs),
                "nR_H2_crp": lambda reacs, prods: ("H2" in reacs)
                and ("CRP" in reacs)
                and (prods.count("H") == 2),
            }

            for key, lambda_filter in reaction_filters.items():
                if lambda_filter(reacs, prods):
                    if (
                        key in self.important_reactions
                        and self.important_reactions[key] is not None
                    ):
                        raise RuntimeError(
                            f"When trying to index the important reactions, we found a disastrous reaction {reaction} is a duplicate of {self.important_reactions[key]}, there can only be one reaction that matches {key}"
                        )
                    self.important_reactions[key] = i + 1

        if np_any([value is None for value in self.important_reactions.values()]):
            logging.debug(self.important_reactions)
            missing_reac_error = "Input reaction file is missing mandatory reactions"
            missing_reac_error += (
                "\nH and E- freeze out as well as H2 formation and photodissociation"
            )
            missing_reac_error += " must all be included in user reaction list. Check default_grain_network.csv for example"
            raise RuntimeError(missing_reac_error)

    def index_important_species(self) -> None:
        """Obtain the indices for all the important reactions."""
        self.species_indices = {}
        names = [species.name for species in self.get_species_list()]
        for element in [
            "C+",
            "H+",
            "H2",
            "HE",
            "HE+",
            "N",
            "N+",
            "O",
            "O+",
            "C",
            "C+",
            "SI+",
            "S+",
            "H2O",
            "CH3OH",
            "CL",
            "CL+",
            "CO",
            "MG",
            "MG+",
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
            except ValueError:
                # TODO: The dummy value is currently SURFACE/BULK; We could handle this better somehow
                logging.info(f"\t{element} not in network, adding dummy index")
                species_index = len(self.get_species_list()) + 1
            name = "n" + element.lower().replace("+", "x").replace(
                "e-", "elec"
            ).replace("#", "g")
            self.species_indices[name] = species_index

    def branching_ratios_checks(self) -> None:
        """Check that the branching ratios for the ice reactions sum to 1.0. If they do not, correct them.
        This needs to be done for LH and LHDES separately since we already added the desorption to the network.
        """
        branching_reactions = {}
        for i, reaction in enumerate(self.get_reaction_list()):
            if reaction.get_reaction_type() in ["LH", "LHDES"]:
                reactant_string = ",".join(reaction.get_reactants())
                if reactant_string in branching_reactions:
                    branching_reactions[reactant_string] += reaction.get_alpha()
                else:
                    branching_reactions[reactant_string] = reaction.get_alpha()
        if not all(branching_reactions.values()) == 1.0:
            logging.warning(
                "Some of the branching ratios do not sum to 1.0, correcting those that do not"
            )
            for i, reaction in enumerate(self.get_reaction_list()):
                if reaction.get_reaction_type() in ["LH", "LHDES"]:
                    reactant_string = ",".join(reaction.get_reactants())
                    # Check if we need to correct the branching ratio (smaller than 0.98 is allowed)
                    if (
                        reactant_string in branching_reactions
                        and branching_reactions[reactant_string] != 1.0
                    ):
                        new_reaction = deepcopy(reaction)
                        if branching_reactions[reactant_string] != 0.0:
                            if branching_reactions[reactant_string] < 0.99:
                                logging.warning(
                                    f"You have reaction {reaction} with a branching ratio {branching_reactions[reactant_string]} we are assuming you set this lower on purpose."
                                )
                                continue
                            new_alpha = (
                                reaction.get_alpha()
                                / branching_reactions[reactant_string]
                            )
                            logging.warning(
                                f"Grain reaction {reaction} has a branching ratio of {reaction.get_alpha()}, dividing it by {branching_reactions[reactant_string]} resulting in BR of {new_alpha}"
                            )
                            # TODO: apply to all partners of the reaction
                            reaction_index = self.get_reaction_index(reaction)
                            reaction.set_alpha(new_alpha)
                            self.set_reaction(
                                reaction_idx=reaction_index, reaction=reaction
                            )
                        else:
                            if isinstance(reaction, CoupledReaction) and (
                                not reaction in self.get_reaction_list()
                            ):
                                logging.info(
                                    f"Tried to remove a coupled reaction {reaction}, but it was already removed by one of its partners."
                                )
                            else:
                                logging.warning(
                                    f"Grain reaction {reaction} has a branching ratio of 0.0, removing the reaction altogether"
                                )
                                self.remove_reaction(reaction)

    def add_gas_phase_extrapolation(self):
        for reaction in self.reactions:
            if reaction.get_reaction_type() in ["TWOBODY", "PHOTON", "CRP", "CRPHOT"]:
                similar_reactions = self.find_similar_reactions(reaction)
                # Only enable extrapolation if we have one or overlapping reactions
                # UMIST uses overlapping reactions to get more correct reaction rates.
                if all(
                    [
                        reaction.check_temperature_collision(
                            similar_reactions[similar_reaction_key]
                        )
                        for similar_reaction_key in similar_reactions
                    ]
                ):
                    reaction.set_extrapolation = True

    def __repr__(self) -> str:
        return (
            "Reaction network with \nSpecies:\n"
            + ", ".join(map(str, self.get_species_list()))
            + "\nReactions:\n"
            + "\n".join(map(repr, self.get_reaction_list()))
        )


class LoadedNetwork(Network):
    """Network version that skips all steps and just loads two lists. This is another
    here be dragons version, use this with exceeding caution as no checks are performed for you.

    Args:
        Network (_type_): _description_
    """

    def __init__(self, species: list[Species], reactions: list[Reaction]) -> None:
        """A loader of networks without any checks.

        Here be dragons.

        Args:
            species (list[Species]): A list of species objects
            reactions (list[Reaction]): A list of reaction objects.
        """
        self.set_species_dict({s.name: s for s in species})
        self.set_reaction_dict({k: v for k, v in enumerate(reactions)})
