"""
NetworkBuilder - Handles complex network construction logic.

This module extracts the build-time complexity from the Network class,
providing a clean separation between:
- Network: Data container with unified interface
- NetworkBuilder: Build-time validation and automatic reaction generation
"""

import logging
from copy import deepcopy
from pathlib import Path
from typing import Union

from numpy import any as np_any

from uclchem.makerates.heating import convert_to_erg

from .reaction import CoupledReaction, Reaction, reaction_types
from .species import Species, elementList


class NetworkBuilder:
    """Builder for constructing complex chemical networks.

    Handles all build-time operations:
    - Input validation
    - Automatic freeze-out reactions
    - Automatic bulk species and reactions
    - Automatic desorption reactions
    - Branching ratio validation and correction
    - Temperature range collision detection
    - Gas-phase reaction extrapolation
    - Reaction exothermicity calculation

    This class separates the complex build logic from the Network data container,
    making the code more maintainable and testable.

    Example:
        >>> builder = NetworkBuilder(
        ...     species=species_list,
        ...     reactions=reactions_list,
        ...     gas_phase_extrapolation=True,
        ...     add_crp_photo_to_grain=True
        ... )
        >>> network = builder.build()
    """

    def __init__(
        self,
        species: list[Species],
        reactions: list[Reaction],
        user_defined_bulk: list = None,
        gas_phase_extrapolation: bool = False,
        add_crp_photo_to_grain: bool = False,
        derive_reaction_exothermicity: list[str] = None,
        database_reaction_exothermicity: list[Union[str, Path]] = None,
    ):
        """Initialize the network builder.

        Args:
            species: List of chemical species
            reactions: List of chemical reactions
            user_defined_bulk: User-specified bulk species (optional)
            gas_phase_extrapolation: Extrapolate gas-phase to grain (default: False)
            add_crp_photo_to_grain: Add CRP/PHOTON to grain (default: False)
            derive_reaction_exothermicity: Reaction types to calculate exothermicity for
            database_reaction_exothermicity: Custom exothermicity database files

        Raises:
            AssertionError: If duplicate species are provided
        """
        # Validate inputs
        assert len({s.get_name() for s in species}) == len(
            species
        ), "Cannot have duplicate species in the species list."

        # Store inputs
        self.input_species = species
        self.input_reactions = reactions

        # Store options
        self.user_defined_bulk = user_defined_bulk or []
        self.gas_phase_extrapolation = gas_phase_extrapolation
        self.add_crp_photo_to_grain = add_crp_photo_to_grain
        self.derive_reaction_exothermicity = derive_reaction_exothermicity
        self.database_reaction_exothermicity = database_reaction_exothermicity

        # Will be set during build
        self.network = None
        self.excited_species = False
        self.enthalpies_present = False

    def build(self):
        """Build the network with all validations and automatic additions.

        This method orchestrates all build steps in the correct order:
        1. Create initial network from inputs
        2. Add electron species
        3. Check and handle freeze/desorb species
        4. Add automatic reactions (freeze, bulk, desorb, chemdes)
        5. Validate branching ratios
        6. Apply optional features (extrapolation, exothermicity)
        7. Sort and filter final network

        Returns:
            Network: Fully built and validated network
        """
        # Import here to avoid circular dependency
        from .network import Network

        # Create initial network from inputs
        logging.info(
            "Building network from %d species and %d reactions",
            len(self.input_species),
            len(self.input_reactions),
        )

        self.network = Network.from_lists(self.input_species, self.input_reactions)

        # Store options on network for methods that need them
        self.network.user_defined_bulk = self.user_defined_bulk
        self.network.add_crp_photo_to_grain = self.add_crp_photo_to_grain
        self.network.derive_reaction_exothermicity = self.derive_reaction_exothermicity
        self.network.database_reaction_exothermicity = (
            self.database_reaction_exothermicity
        )
        self.network.enthalpies_present = False

        # Check for excited species
        self.excited_species = self._check_for_excited_species()
        self.network.excited_species = self.excited_species

        # Add electron if not present
        electron_specie = Species(["E-", 0, 0.0, 0, 0, 0, 0])
        electron_specie.set_n_atoms(1)
        self.network.add_species(electron_specie)

        logging.info("Starting automatic reaction and species generation")

        # Check which species change on freeze or desorb
        self._check_freeze_and_desorbs()

        # Add automatic grain reactions
        self._add_freeze_reactions()

        if self.add_crp_photo_to_grain:
            self._add_CRP_and_PHOTO_reactions_to_grain()

        self._add_bulk_species()
        self._add_bulk_reactions()
        self._add_desorb_reactions()
        self._add_chemdes_reactions()

        if self.excited_species:
            self._add_excited_surface_reactions()

        # Validate and correct branching ratios
        self._branching_ratios_checks()

        # Apply optional features
        if self.gas_phase_extrapolation:
            logging.info("Applying gas-phase extrapolation")
            self._add_gas_phase_extrapolation()

        if self.derive_reaction_exothermicity:
            logging.info(
                "Calculating reaction exothermicities for types: %s",
                self.derive_reaction_exothermicity,
            )
            self._add_reaction_enthalpies(self.derive_reaction_exothermicity)

        if self.database_reaction_exothermicity:
            logging.info(
                "Applying custom exothermicity files: %s",
                self.database_reaction_exothermicity,
            )
            self._apply_custom_exothermicities(self.database_reaction_exothermicity)

        # Final sorting and filtering
        logging.info("Sorting and filtering final network")
        self.network.sort_reactions()
        self.network.sort_species()
        self._check_and_filter_species()

        # Run validation checks and indexing
        self._check_network()

        logging.info(
            "Network building complete: %d species, %d reactions",
            len(self.network.get_species_list()),
            len(self.network.get_reaction_list()),
        )

        return self.network

    # ========================================================================
    # Private Build Methods
    # ========================================================================

    def _check_for_excited_species(self) -> bool:
        """Check if there are any excited species in the network.

        Returns:
            bool: True if any species name contains '*'
        """
        return any(
            "*" in species.get_name() for species in self.network.get_species_list()
        )

    def _check_and_filter_species(self) -> None:
        """Check every species in network appears in at least one reaction.
        Remove any that do not and alert user.
        """
        # check for species not involved in any reactions
        lostSpecies = []
        for species in self.network.get_species_list():
            # keep species that appear in a reaction
            reac_keeps = False
            for reaction in self.network.get_reaction_list():
                if (
                    species.get_name() in reaction.get_reactants()
                    or species.get_name() in reaction.get_products()
                ):
                    reac_keeps = True
                    break

            # remove the species if it didn't make it into either keep list
            if not (reac_keeps):
                lostSpecies.append(species.get_name())
        for species in lostSpecies:
            logging.warning(
                f"Trying to remove {species} as it is not present in the reactions"
            )
            self.network.remove_species(species)

        # then alert user to changes
        if len(lostSpecies) > 0:
            logging.warning(
                "\tSpecies in input list that do not appear in final list:\t"
                + str(lostSpecies)
            )
        else:
            logging.info("\tAll input species in final network")
        logging.debug(
            f"The network consists of species: {self.network.get_species_list()}"
        )
        for species in self.network.get_species_list():
            species.find_constituents()

        # add in pseudo-species to track mantle
        mantle_specs = []
        new_spec = [999] * 7
        new_spec[0] = "BULK"
        mantle_specs.append(Species(new_spec))
        new_spec[0] = "SURFACE"
        mantle_specs.append(Species(new_spec))
        self.network.add_species(mantle_specs)

    def _check_freeze_and_desorbs(self) -> None:
        """`_add_freeze_reactions()` and `_add_desorb_reactions()` automatically generate
        all desorption and freeze out reactions. However, user may want to change a species on freeze out
        eg C+ becomes #C rather than #C+. This function checks for that and updates species so they'll
        freeze or desorb correctly when reactions are generated.
        """
        desorbs = [
            x
            for x in self.network.get_reaction_list()
            if x.get_reaction_type() == "DESORB"
        ]
        for desorb in desorbs:
            specie = self.network.get_specie(desorb.get_reactants()[0])
            specie.set_desorb_products(desorb.get_products())
            self.network.set_specie(desorb.get_reactants()[0], specie)
            # also modify bulk species desorb_products
            bulk_name = "@" + desorb.get_reactants()[0][1:]
            if bulk_name in self.network._species_dict:
                specie = self.network.get_specie(bulk_name)
                specie.set_desorb_products(desorb.get_products())
                self.network.set_specie(bulk_name, specie)

        # for all listed freeze out reactions, add them to correct species
        freezes = [
            x
            for x in self.network.get_reaction_list()
            if x.get_reaction_type() == "FREEZE"
        ]
        for freeze in freezes:
            logging.debug(freeze)
            specie = self.network.get_specie(freeze.get_reactants()[0])
            specie.set_freeze_products(freeze.get_products(), freeze.get_alpha())
            self.network.set_specie(freeze.get_reactants()[0], specie)

        # then add default freeze out for species without a listed freeze out
        for species_name, specie in self.network.get_species_dict().items():
            if (not specie.is_ice_species()) and (not specie.get_freeze_products_list()):
                logging.info(f"Adding a default freezeout for {specie} to the specie")
                specie.add_default_freeze()
                self.network.set_specie(species_name, specie)

        # Here we filter all the freeze and desorb reactions in order to avoid duplicates
        [self.network.remove_reaction(reaction) for reaction in desorbs + freezes]

    def _add_freeze_reactions(self) -> None:
        """Save the user effort by automatically generating freeze out reactions"""
        logging.debug("Adding the freeze out reactions!")
        new_reactions = []
        new_species = []
        for species in self.network.get_species_list():
            logging.debug(f"Checking if {species} needs to have its freezeout added")
            if not species.is_ice_species():
                for products, alpha in species.get_freeze_products():
                    if species.get_name() == "E-":
                        # Set electron freeze out to zero:
                        alpha = 0.0
                    new_reactions.append(
                        Reaction(
                            [species.get_name(), "FREEZE", "NAN"]
                            + products
                            + [
                                alpha,
                                0.0,
                                species.get_binding_energy(),
                                0.0,
                                10000.0,
                                0.0,
                            ]
                        )
                    )
                    # Check if the product is in the species list
                    if products[0] not in self.network.get_species_list():
                        logging.info(f"Trying to add new specie {products}")
                        new_species.append(
                            Species(
                                [
                                    products[0],
                                    species.get_mass(),
                                    species.get_binding_energy(),
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                ]
                            )
                        )
        if new_reactions:
            self.network.add_reactions(new_reactions)
        if new_species:
            self.network.add_species(new_species)

    def _add_bulk_species(self) -> None:
        """For three phase models, MakeRates will produce the version of the species in the bulk
        so that the user doesn't have to endlessly relist the same species
        """
        logging.debug("Adding bulk species")
        speciesNames = [species.get_name() for species in self.network.get_species_list()]
        userSpecies = [manualSpec.get_name() for manualSpec in self.user_defined_bulk]
        new_species = []
        try:
            h2o_binding_energy = speciesNames.index("#H2O")
            h2o_binding_energy = self.network.get_species_list()[
                h2o_binding_energy
            ].get_binding_energy()
        except ValueError:
            error = "You are trying to create a three phase model but #H2O is not in your network"
            error += "\nThis is likely an error so Makerates will not complete. Try adding #H2O"
            raise RuntimeError(error)
        for species in self.network.get_species_list():
            if species.is_surface_species():
                if species.get_name().replace("#", "@") not in speciesNames:
                    new_spec = deepcopy(species)
                    new_spec.set_name(new_spec.get_name().replace("#", "@"))
                    if new_spec.get_name() in userSpecies:
                        definedBinding = [
                            userSpec.get_binding_energy()
                            for userSpec in self.user_defined_bulk
                            if userSpec.get_name() == new_spec.get_name()
                        ]
                        new_spec.set_binding_energy(definedBinding[0])
                    else:
                        new_spec.set_binding_energy(h2o_binding_energy)
                    new_species.append(new_spec)
                    logging.debug(
                        f"Adding a bulk partner species for {species}, new {new_spec}"
                    )
        self.network.add_species(new_species)

    def _add_bulk_reactions(self) -> None:
        """We assume any reaction that happens on the surface of grains can also happen
        in the bulk (just more slowly due to binding energy). The user therefore only
        lists surface reactions in their input reaction file and we duplicate here.
        """
        logging.debug("Adding bulk reactions")
        surface_reactions = self._get_reactions_on_grain()
        bulk_reaction_types = ["CRP", "CRPHOT", "PHOTON", "LH", "EXSOLID", "EXRELAX"]
        surface_reactions_can_be_bulk = [
            reaction
            for reaction in surface_reactions
            if reaction.get_reaction_type() in bulk_reaction_types
        ]
        current_reaction_list = self.network.get_reaction_list()

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

        bulk_species = [x for x in self.network.get_species_list() if "@" in x.get_name()]
        for species in bulk_species:
            # add individual swapping
            if not species.is_refractory:
                new_reac_list = [
                    species.get_name(),
                    "BULKSWAP",
                    "NAN",
                    species.get_name().replace("@", "#"),
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
            if species not in [
                "@H2"
            ]:  # If species is H2, do not allow it to go from surface to bulk
                new_reac_list[0] = species.get_name().replace("@", "#")
                new_reac_list[1] = "SURFSWAP"
                new_reac_list[3] = species.get_name()
                new_reactions.append(Reaction(new_reac_list))
        logging.debug(
            f"The following bulk reactions are added to the reactions: {new_reactions}"
        )
        self.network.add_reactions(new_reactions)

    def _add_desorb_reactions(self) -> None:
        """Save the user effort by automatically generating desorption reactions"""
        desorb_reacs = ["DESOH2", "DESCR", "DEUVCR", "THERM"]
        logging.debug("Adding desorption reactions!")
        new_reactions = []
        for species in self.network.get_species_list():
            if species.is_surface_species():
                for reacType in desorb_reacs:
                    new_reactions.append(
                        Reaction(
                            [species.get_name(), reacType, "NAN"]
                            + species.get_desorb_products()
                            + [1, 0, species.get_binding_energy(), 0.0, 10000.0, 0.0]
                        )
                    )
            if species.is_bulk_species() and not species.is_refractory:
                new_reactions.append(
                    Reaction(
                        [species.get_name(), "THERM", "NAN"]
                        + species.get_desorb_products()
                        + [1, 0, species.get_binding_energy(), 0.0, 10000.0, 0.0]
                    )
                )
        self.network.add_reactions(new_reactions)

    def _add_chemdes_reactions(self) -> None:
        """We have the user list all Langmuir-Hinshelwood and Eley-Rideal
        reactions once. Then we duplicate so that the reaction branches
        with products on grain and products desorbing.
        """
        logging.debug("Adding chemical desorption reactions for LH and ER mechanisms")
        new_reactions = []
        species_list = self.network.get_species_list()
        species_names = [species.name for species in species_list]
        for reaction in self.network.get_reaction_list():
            reactants = reaction.get_reactants()
            if reactants[0][0] == "@" or reactants[1][0] == "@":
                continue
            if reaction.get_reaction_type() in ["LH", "ER"]:
                nProducts = sum(prod != "NAN" for prod in reaction.get_products())

                reactionPartner = deepcopy(reaction)
                while isinstance(reactionPartner, CoupledReaction):
                    # If the current loop reaction is also coupled, get its partner.
                    reactionPartner = reactionPartner.get_partner()

                for i in range(nProducts):
                    # For each of the products, make a new reaction where it is desorbed
                    new_reaction = deepcopy(reaction)

                    # Convert to disassociation reaction
                    new_reactants = new_reaction.get_reactants()
                    new_reactants[2] += "DES"
                    new_reaction.set_reactants(new_reactants)

                    # Replace the species on the grain/bulk with species in gas
                    new_products = new_reaction.get_products()
                    # Check there are no products incompatible with LH/ER reactions by counting them
                    if new_products.count("NAN") + sum(
                        prod.startswith("#") or prod.startswith("@")
                        for prod in new_products
                    ) != len(new_products):
                        logging.warning(
                            "All Langmuir-Hinshelwood and Eley-Rideal reactions should be input with products on grains only.\n"
                            + "The fraction of products that enter the gas is dealt with by Makerates and UCLCHEM.\n"
                            + "the following reaction caused this warning\t"
                            + str(reaction)
                        )

                    # Replace all grain or bulk products with gas phase counterparts
                    desorbProducts = species_list[
                        species_names.index(new_products[i])
                    ].get_desorb_products()

                    if desorbProducts[1] == "NAN":
                        new_products[i] = desorbProducts[0]
                    elif desorbProducts[2] == "NAN":
                        if i < 2 and new_products[i + 1] != "NAN":
                            # Move i+1th product over to i+2th product
                            new_products[i + 2] = new_products[i + 1]
                        new_products[i + 1] = desorbProducts[1]
                        new_products[i] = desorbProducts[0]
                    else:
                        raise NotImplementedError()

                    new_reaction.set_products(new_products)
                    logging.debug(
                        f"Adding chemical desorption reaction for {reaction}, new reaction {new_reaction}"
                    )

                    new_reaction = CoupledReaction(new_reaction)
                    new_reaction.set_partner(reactionPartner)

                    new_reactions.append(new_reaction)
        self.network.add_reactions(new_reactions)

    def _add_excited_surface_reactions(self) -> None:
        """All excited species will relax to the ground state if they do not react
        the vibrational frequency of the species is used as a pseudo approximation of the rate coefficient
        We assume all grain reactions have an excited variant. For example:
        #A, #B LH #C will have the variants:
        #A*, #B EXSOLID #C  and  #A, #B* EXSOLID #C
        If only one of the reactants in the base reaction has an excited counterpart then
        only one excited version of that reaction is created.
        """
        logging.debug("Adding excited surface reactions")
        excited_species = [
            x for x in self.network.get_species_list() if "*" in x.get_name()
        ]
        lh_reactions = [
            x for x in self.network.get_reaction_list() if "LH" in x.get_reactants()
        ]
        lh_reactions += [
            x for x in self.network.get_reaction_list() if "LHDES" in x.get_reactants()
        ]

        new_reactions = []

        # add relaxation of excited species
        for spec in excited_species:
            relax_reac = [
                spec.get_name(),
                "EXRELAX",
                "NAN",
                spec.get_name()[:-1],
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
                specie.get_name() for specie in excited_species
            ] and reaction.get_reactants()[1] + "*" in [
                specie.get_name() for specie in excited_species
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
                specie.get_name() for specie in excited_species
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
                specie.get_name() for specie in excited_species
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
        self.network.add_reactions(new_reactions)

    def _add_CRP_and_PHOTO_reactions_to_grain(self) -> None:
        """Add all the gas-phase reactions with CRP, CRPHOT or PHOTON to the grain surface too"""
        logging.info("Adding gas-phase reactions with CRP, CRPHOT or PHOTON to grain")
        reactions_on_grain = self._get_reactions_on_grain()
        reactions_on_grain_filtered = [
            reaction
            for reaction in reactions_on_grain
            if reaction.get_reaction_type() in ["CRP", "CRPHOT", "PHOTON"]
        ]
        new_reactions = []
        for reaction in self.network.get_reaction_list():
            if reaction.get_reaction_type() not in ["CRP", "CRPHOT", "PHOTON"]:
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
                species in self.network.get_species_list()
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
        logging.debug("Adding new reactions to grain")
        self.network.add_reactions(new_reactions)
        logging.info(f"Added {len(new_reactions)} reactions to grain")

    def _branching_ratios_checks(self) -> None:
        """Check that the branching ratios for the ice reactions sum to 1.0. If they do not, correct them.
        This needs to be done for LH and LHDES separately since we already added the desorption to the network.
        """
        branching_reactions = {}
        for i, reaction in enumerate(self.network.get_reaction_list()):
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
            for i, reaction in enumerate(self.network.get_reaction_list()):
                if reaction.get_reaction_type() in ["LH", "LHDES"]:
                    reactant_string = ",".join(reaction.get_reactants())
                    # Check if we need to correct the branching ratio (smaller than 0.98 is allowed)
                    if (
                        reactant_string in branching_reactions
                        and branching_reactions[reactant_string] != 1.0
                    ):
                        # new_reaction = deepcopy(reaction)  # Unused variable
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
                            reaction_index = self.network.get_reaction_index(reaction)
                            reaction.set_alpha(new_alpha)
                            self.network.set_reaction(
                                reaction_idx=reaction_index, reaction=reaction
                            )
                        else:
                            if isinstance(reaction, CoupledReaction) and (
                                reaction not in self.network.get_reaction_list()
                            ):
                                logging.info(
                                    f"Tried to remove a coupled reaction {reaction}, but it was already removed by one of its partners."
                                )
                            else:
                                logging.warning(
                                    f"Grain reaction {reaction} has a branching ratio of 0.0, removing the reaction altogether"
                                )
                                self.network.remove_reaction(reaction)

    def _add_gas_phase_extrapolation(self):
        """Enable extrapolation for gas-phase reactions that have unique or overlapping temperature ranges."""
        for reaction in self.network._reactions_dict.values():
            if reaction.is_gas_reaction() and (
                reaction.get_reaction_type() in ["TWOBODY", "PHOTON", "CRP", "CRPHOT"]
            ):
                similar_reactions = self.network.find_similar_reactions(reaction)
                # Only enable extrapolation if we have one or overlapping reactions
                # UMIST uses overlapping reactions to get more correct reaction rates.
                if all(
                    (reaction.get_templow() == v.get_templow())
                    and (reaction.get_temphigh() == v.get_temphigh())
                    for k, v in similar_reactions.items()
                ):
                    reaction.set_extrapolation(True)

    def _add_reaction_enthalpies(self, enthalpy_reaction_types):
        """Add reaction enthalpies (exothermicity) to reactions for heating/cooling calculations.

        Args:
            enthalpy_reaction_types: List of reaction types or "ALL" or "GAS"
        """
        exclude_ices = True
        if not isinstance(enthalpy_reaction_types, list):
            enthalpy_reaction_types = [enthalpy_reaction_types]
        if enthalpy_reaction_types[0].upper() == "ALL":
            exclude_ices = False
            enthalpy_reaction_types = list(reaction_types)
        elif enthalpy_reaction_types[0].upper() == "GAS":
            enthalpy_reaction_types = list(reaction_types)
        for reaction in self.network.get_reaction_list():
            logging.debug(f"Checking if we need to add enthalpy to {reaction}")
            if reaction.get_reaction_type() in enthalpy_reaction_types:
                if exclude_ices and reaction.is_ice_reaction(strict=(not exclude_ices)):
                    logging.debug("Skipping ice reaction")
                    continue
                if "E-" in (reaction.get_pure_products() + reaction.get_pure_reactants()):
                    logging.debug(
                        "Reaction involving electrons, skipping enthalpy due to poor estimates"
                    )
                    continue
                delta_h = self._compute_exothermicity(reaction)
                # TODO: add a heating efficiency factor in here:
                reaction.set_exothermicity(convert_to_erg(-delta_h, "kcal/mol"))
                logging.debug(
                    f"Setting reaction enthalpy of {reaction} to {delta_h} kcal/mol"
                )

    def _apply_custom_exothermicities(self, database_reaction_exothermicity: list):
        """Apply custom exothermicity values from CSV files to the network reactions.

        Args:
            database_reaction_exothermicity (list): List of paths to custom exothermicity CSV files.
        """
        from .heating import set_custom_exothermicities

        for csv_path in database_reaction_exothermicity:
            logging.info(f"Applying custom exothermicities from {csv_path}")
            set_custom_exothermicities(
                reactions=self.network.get_reaction_list(),
                csv_path=csv_path,
                overwrite=True,
            )

    def _compute_exothermicity(self, reaction: Reaction) -> float:
        """Compute the reaction enthalpy in eV for a given reaction based on the
        species enthalpies.

        Args:
            reaction (Reaction): The reaction to compute the enthalpy for.
        Returns:
            float: The reaction enthalpy in kcal/mol.
        """
        reactants = reaction.get_pure_reactants()
        products = reaction.get_pure_products()
        return sum(self.network._species_dict[p].get_enthalpy() for p in products) - sum(
            self.network._species_dict[r].get_enthalpy() for r in reactants
        )

    def _get_reactions_on_grain(self) -> list[Reaction]:
        """Get all reactions that occur on grain surfaces (# prefix) or in bulk (@prefix)."""
        reactions_on_grain = []
        for reaction in self.network.get_reaction_list():
            reactants = reaction.get_reactants()
            if any("#" in reactant or "@" in reactant for reactant in reactants):
                reactions_on_grain.append(reaction)
        return reactions_on_grain

    # ========================================================================
    # Validation and Indexing Methods
    # ========================================================================

    def _check_network(self) -> None:
        """Run all validation checks and create important indices."""
        self._freeze_checks()
        self._duplicate_checks()
        self._index_important_reactions()
        self._index_important_species()

    def _freeze_checks(self) -> None:
        """Check that every species freezes out and alert the user if a
        species freezes out via mutiple routes. This isn't necessarily an
        error so best just print.
        """
        logging.info(
            "\tCheck that species have surface counterparts or if they have multiple freeze outs/check alphas:\n"
        )
        for spec in self.network.get_species_list():
            if not spec.is_ice_species() and spec.get_name()[-1] not in ["+", "-"]:
                exist_check = 0
                for checkSpeck in self.network.get_species_list():
                    if checkSpeck.get_name() == "#" + spec.get_name():
                        exist_check += 1
                if exist_check == 0:
                    logging.warning(
                        f"{spec.get_name()} does not have a surface counterpart in given default species file."
                        + "\n\tThis sets the binding energy to zero, it might cause species conservation errors."
                    )
            freezes = 0
            freezeout_reactions = []
            for reaction in self.network.get_reaction_list():
                if (spec.get_name() in reaction.get_reactants()) and (
                    "FREEZE" in reaction.get_reactants()
                ):
                    freezes += 1
                    freezeout_reactions.append(reaction)
            if freezes == 1:
                logging.info(
                    f"\t{spec.get_name()} freezes out through {freezeout_reactions[0]}"
                )
            if freezes > 1:
                logging.info(f"\t{spec.get_name()} freezes out through {freezes} routes")
            elif freezes < 1 and not spec.is_ice_species():
                logging.info(f"\t{spec.get_name()} does not freeze out")

    def _duplicate_checks(self) -> None:
        """
        Check reaction network to make sure no reaction appears twice unless
        they have different temperature ranges.
        """
        logging.info("\tPossible duplicate reactions for manual removal:")
        duplicates = False
        for i, reaction1 in enumerate(self.network.get_reaction_list()):
            if not reaction1.duplicate:
                for j, reaction2 in enumerate(self.network.get_reaction_list()):
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
                                    if reaction1.get_templow() < reaction2.get_temphigh():
                                        logging.warning(
                                            f"\tReactions {reaction1} and {reaction2} have non-adjacent temperature ranges"
                                        )
                                reaction1.duplicate = True
                                reaction2.duplicate = True
        if not duplicates:
            logging.info("\tNone")

    def _index_important_reactions(self) -> None:
        """We have a whole bunch of important reactions and we want to store
        their indices. We find them all here.
        """

        # Any None values in dictionary will raise an error
        # therefore these reactions are mandatory and makerates will not complete if the user doesn't supply them.
        self.network.important_reactions = {
            "nR_H2Form_CT": None,
            "nR_H2Form_ERDes": None,
            "nR_H2Form_ER": None,
            "nR_H2Form_LH": None,
            "nR_H2Form_LHDes": None,
            "nR_HFreeze": None,
            "nR_EFreeze": None,
            "nR_H2_hv": None,
            "nR_H2_ED": None,
            "nR_H_ED": None,
        }

        # this looks complex but each if statement just uniquely identifies a special reaction
        # if found, it is added to the dictionary with its fortran index as the value
        for i, reaction in enumerate(self.network.get_reaction_list()):
            # CO + PHOTON -> O + C
            reacs = reaction.get_reactants()
            prods = reaction.get_products()
            reaction_filters = {
                "nR_CO_hv": lambda reacs, prods: (
                    ("CO" in reacs)
                    and ("PHOTON" in reacs)
                    and ("O" in prods)
                    and ("C" in prods)
                ),
                "nR_C_hv": lambda reacs, prods: ("C" in reacs) and ("PHOTON" in reacs),
                "nR_H2Form_CT": lambda reacs, prods: "H2FORM" in reacs,
                "nR_H2Form_ERDes": lambda reacs, prods: (
                    ("H" in reacs) and ("#H" in reacs) and ("H2" in prods)
                ),
                "nR_H2Form_ER": lambda reacs, prods: (
                    ("H" in reacs) and ("#H" in reacs) and ("#H2" in prods)
                ),
                "nR_H2Form_LH": lambda reacs, prods: (
                    (reacs.count("#H") == 2) and ("LH" in reacs)
                ),
                "nR_H2Form_LHDes": lambda reacs, prods: (
                    (reacs.count("#H") == 2) and ("LHDES" in reacs)
                ),
                "nR_HFreeze": lambda reacs, prods: ("H" in reacs) and ("FREEZE" in reacs),
                "nR_H2Freeze": lambda reacs, prods: (
                    ("H2" in reacs) and ("FREEZE" in reacs)
                ),
                "nR_EFreeze": lambda reacs, prods: (
                    ("E-" in reacs) and ("FREEZE" in reacs)
                ),
                "nR_H2_hv": lambda reacs, prods: ("H2" in reacs) and ("PHOTON" in reacs),
                "nR_H2_crp": lambda reacs, prods: (
                    ("H2" in reacs) and ("CRP" in reacs) and (prods.count("H") == 2)
                ),
                "nR_H2_ED": lambda reacs, prods: (
                    ("#H2" in reacs) and ("ED" in reacs) and ("H2" in prods)
                ),
                "nR_H_ED": lambda reacs, prods: (
                    ("#H" in reacs) and ("ED" in reacs) and ("H" in prods)
                ),
            }

            for key, lambda_filter in reaction_filters.items():
                if lambda_filter(reacs, prods):
                    if (
                        key in self.network.important_reactions
                        and self.network.important_reactions[key] is not None
                    ):
                        raise RuntimeError(
                            f"When trying to index the important reactions, we found a disastrous reaction {reaction} is a duplicate of {self.network.important_reactions[key]}, there can only be one reaction that matches {key}"
                        )
                    self.network.important_reactions[key] = i + 1

        if np_any([value is None for value in self.network.important_reactions.values()]):
            logging.debug(self.network.important_reactions)
            missing_reac_error = "Input reaction file is missing mandatory reactions"
            missing_reac_error += (
                "\nH and E- freeze out as well as H2 formation and photodissociation"
            )
            missing_reac_error += " must all be included in user reaction list. Check default_grain_network.csv for example"
            raise RuntimeError(missing_reac_error)

    def _index_important_species(self) -> None:
        """Obtain the indices for all the important reactions."""
        self.network.species_indices = {}
        names = [species.get_name() for species in self.network.get_species_list()]
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
                species_index = len(self.network.get_species_list()) + 1
            name = "n" + element.lower().replace("+", "x").replace("e-", "elec").replace(
                "#", "g"
            )
            self.network.species_indices[name] = species_index
