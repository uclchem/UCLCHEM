"""NetworkBuilder - Handles complex network construction logic.

This module extracts the build-time complexity from the Network class,
providing a clean separation between:
- Network: Data container with unified interface
- NetworkBuilder: Build-time validation and automatic reaction generation

"""

import logging
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING

from numpy import any as np_any

from uclchem.makerates.heating import convert_to_erg, set_custom_exothermicities

from .reaction import (
    BULK_REACTION_TYPES,
    REACTION_TYPES,
    TUNNELING_REACTION_TYPES,
    CoupledReaction,
    Reaction,
    get_duplicate_reactions,
)
from .species import Species, elementList

if TYPE_CHECKING:
    from .network import Network

logger = logging.getLogger(__name__)


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

    Examples
    --------
    >>> from uclchem.makerates.io_functions import read_species_file, read_reaction_file
    >>> from uclchem.utils import UCLCHEM_ROOT_DIR
    >>>
    >>> species_list, user_defined_bulk = read_species_file(
    ...     UCLCHEM_ROOT_DIR / "../../Makerates/data/default/default_species.csv"
    ... )
    >>> reactions_list, dropped_reactions = read_reaction_file(
    ...     UCLCHEM_ROOT_DIR / "../../Makerates/data/default/default_grain_network.csv",
    ...     species_list,
    ...     "UCL",
    ... )
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
        user_defined_bulk: list | None = None,
        derive_reaction_exothermicity: str | set[str] | None = None,
        database_reaction_exothermicity: list[str | Path] | None = None,
        *,
        gas_phase_extrapolation: bool = False,
        add_crp_photo_to_grain: bool = False,
    ):
        """Initialize the network builder.

        Parameters
        ----------
        species : list[Species]
            List of chemical species
        reactions : list[Reaction]
            List of chemical reactions
        user_defined_bulk : list | None
            User-specified bulk species (optional) (Default value = None)
        derive_reaction_exothermicity : str | set[str] | None
            Reaction types to calculate exothermicity for (Default value = None)
        database_reaction_exothermicity : list[str | Path] | None
            Custom exothermicity database files (Default value = None)
        gas_phase_extrapolation : bool
            Extrapolate gas-phase temperature (default: False)
        add_crp_photo_to_grain : bool
            Add CRP/PHOTON to grain (default: False)

        Raises
        ------
        ValueError
            If duplicate species are provided.

        """
        # Validate inputs
        if len({s.get_name() for s in species}) != len(species):
            msg = "Cannot have duplicate species in the species list."
            raise ValueError(msg)

        # Store inputs as immutable
        self.input_species = tuple(species)
        self.input_reactions = tuple(reactions)

        # Store options
        self.user_defined_bulk = user_defined_bulk or []
        self.gas_phase_extrapolation = gas_phase_extrapolation
        self.add_crp_photo_to_grain = add_crp_photo_to_grain
        self.derive_reaction_exothermicity = derive_reaction_exothermicity
        self.database_reaction_exothermicity = database_reaction_exothermicity

        # Will be set during build
        self._network: Network | None = None
        self.excited_species = False
        self.enthalpies_present = False

    @property
    def network(self) -> "Network":
        """The network under construction.

        Returns
        -------
        Network
            The network being built.

        Raises
        ------
        RuntimeError
            If accessed before :meth:`build` has created the network.

        """
        if self._network is None:
            msg = "Network has not been built yet; call build() first."
            raise RuntimeError(msg)
        return self._network

    @network.setter
    def network(self, value: "Network | None") -> None:
        self._network = value

    def build(self) -> "Network":
        """Build the network with all validations and automatic additions.

        This method orchestrates all build steps in the correct order:
        1. Create initial network from inputs
        2. Add electron species
        3. Check and handle freeze/desorb species
        4. Add automatic reactions (freeze, bulk, desorb, chemdes)
        5. Validate branching ratios
        6. Apply optional features (extrapolation, exothermicity)
        7. Sort and filter final network

        Returns
        -------
        Network
            Fully built and validated network

        """
        # Import here to avoid circular dependency
        from .network import Network  # ruff: ignore[import-outside-top-level] circular

        # Create initial network from inputs
        logger.info(
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

        logger.info("Starting automatic reaction and species generation")

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
            logger.info("Applying gas-phase extrapolation")
            self._add_gas_phase_extrapolation()

        if self.derive_reaction_exothermicity:
            logger.info(
                "Calculating reaction exothermicities for types: %s",
                self.derive_reaction_exothermicity,
            )
            self._add_reaction_enthalpies(self.derive_reaction_exothermicity)

        if self.database_reaction_exothermicity:
            logger.info(
                "Applying custom exothermicity files: %s",
                self.database_reaction_exothermicity,
            )
            self._apply_custom_exothermicities(self.database_reaction_exothermicity)

        # Final sorting and filtering
        logger.info("Sorting and filtering final network")
        self.network.sort_reactions()
        self.network.sort_species()
        self._check_and_filter_species()

        # Run validation checks and indexing
        self._check_network()

        logger.info(
            "Network building complete: %d species, %d reactions",
            len(self.network.species),
            len(self.network.reactions),
        )

        return self.network

    # ========================================================================
    # Private Build Methods
    # ========================================================================

    def _check_for_excited_species(self) -> bool:
        """Check if there are any excited species in the network.

        Returns
        -------
        bool
            True if any species name contains '*'

        """
        return any(
            "*" in species.get_name() for species in self.network.get_species_list()
        )

    def _check_and_filter_species(self) -> None:
        """Check every species in network appears in at least one reaction.

        Remove any that do not and alert user.

        """
        # check for species not involved in any reactions
        lost_species = []
        for species in self.network.species:
            # keep species that appear in a reaction
            reac_keeps = False
            for reaction in self.network.reactions.values():
                if (
                    species in reaction.get_reactants()
                    or species in reaction.get_products()
                ):
                    reac_keeps = True
                    break

            # remove the species if it didn't make it into either keep list
            if not reac_keeps:
                lost_species.append(species)
        for lost_name in lost_species:
            logger.warning(
                f"Trying to remove {lost_name} as it is not present in the reactions"
            )
            self.network.remove_species(lost_name)

        # then alert user to changes
        if len(lost_species) > 0:
            logger.warning(
                "\tSpecies in input list that do not appear in final list:\t"
                + str(lost_species)
            )
        else:
            logger.info("\tAll input species in final network")
        logger.debug(
            f"The network consists of species: {self.network.get_species_list()}"
        )

        for spec in self.network.species.values():
            spec.find_constituents()

        # add in pseudo-species to track mantle
        mantle_specs = []
        new_spec: list[str | float] = [999] * 7
        new_spec[0] = "BULK"
        mantle_specs.append(Species(new_spec))
        new_spec[0] = "SURFACE"
        mantle_specs.append(Species(new_spec))
        self.network.add_species(mantle_specs)

    def _check_freeze_and_desorbs(self) -> None:
        """Correct all freeze and desorb reactions for custom input reactions.

        `_add_freeze_reactions()` and `_add_desorb_reactions()` automatically generate
        all desorption and freeze out reactions. However, user may want to change a
        species on freeze out eg C+ becomes #C rather than #C+.
        This function checks for that and updates species so they'll
        freeze or desorb correctly when reactions are generated.

        """
        desorbs = [
            x
            for x in self.network.reactions.values()
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
            for x in self.network.reactions.values()
            if x.get_reaction_type() == "FREEZE"
        ]
        for freeze in freezes:
            logger.debug(freeze)
            specie = self.network.get_specie(freeze.get_reactants()[0])
            specie.set_freeze_products(freeze.get_products(), freeze.get_alpha())
            self.network.set_specie(freeze.get_reactants()[0], specie)

        # then add default freeze out for species without a listed freeze out
        for species_name, specie in self.network.species.items():
            if (not specie.is_ice_species()) and (not specie.get_freeze_products_list()):
                logger.info(f"Adding a default freezeout for {specie} to the specie")
                specie.add_default_freeze()
                self.network.set_specie(species_name, specie)

        # Here we filter all the freeze and desorb reactions in order to avoid duplicates.
        # DESORB reactions whose reactant starts with '#' are DESORB shorthand reactions
        # (expand to THERM/DESOH2/DESCR/DEUVCR in _add_desorb_reactions). Keep those
        # in the list so _add_desorb_reactions() can expand and then remove them.
        desorbs_to_remove = [
            d for d in desorbs if not d.get_reactants()[0].startswith("#")
        ]
        for reaction in desorbs_to_remove + freezes:
            self.network.remove_reaction(reaction)

    def _add_freeze_reactions(self) -> None:
        """Add freeze-out reactions for all species.

        Takes into account custom defined freeze-out reactions. Otherwise,
        species freeze out their neutral counterparts, i.e. `H+ + FREEZE -> #H`

        """
        logger.debug("Adding the freeze out reactions!")
        new_reactions = []
        new_species = []
        for species in self.network.species.values():
            logger.debug(f"Checking if {species} needs to have its freezeout added")
            if not species.is_ice_species():
                for products, alpha in species.get_freeze_products():
                    if species.get_name() == "E-":
                        # Set electron freeze out to zero:
                        alpha = 0.0
                    new_reactions.append(
                        Reaction(
                            [
                                species.get_name(),
                                "FREEZE",
                                "NAN",
                                *products,
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
                    if products[0] not in self.network.species:
                        logger.info(f"Trying to add new specie {products}")
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
        """Copy all surface species to the bulk.

        Their binding energies and diffusion barriers are set to those of H2O,
        to mimic Ghesquire et al, 2015 (https://doi.org/10.1039/C5CP00558B).

        Raises
        ------
        RuntimeError
            If #H2O is not in the network.

        """
        logger.debug("Adding bulk species")
        species_names = [
            species.get_name() for species in self.network.get_species_list()
        ]
        user_species = [manualSpec.get_name() for manualSpec in self.user_defined_bulk]
        new_species = []
        try:
            h2o_index = species_names.index("#H2O")
            h2o_binding_energy = self.network.get_species_list()[
                h2o_index
            ].get_binding_energy()
        except ValueError:
            error = "You are trying to create a three phase model but #H2O is not in your network"
            error += "\nThis is likely an error so Makerates will not complete. Try adding #H2O"
            raise RuntimeError(error) from None
        for species in self.network.get_species_list():
            if (
                species.is_surface_species()
                and species.get_name().replace("#", "@") not in species_names
            ):
                new_spec = deepcopy(species)
                new_spec.set_name(new_spec.get_name().replace("#", "@"))
                if new_spec.get_name() in user_species:
                    defined_binding = [
                        userSpec.get_binding_energy()
                        for userSpec in self.user_defined_bulk
                        if userSpec.get_name() == new_spec.get_name()
                    ]
                    new_spec.set_binding_energy(defined_binding[0])
                else:
                    new_spec.set_binding_energy(h2o_binding_energy)
                new_species.append(new_spec)
                logger.debug(
                    f"Adding a bulk partner species for {species}, new {new_spec}"
                )
        self.network.add_species(new_species)

    def _add_bulk_reactions(self) -> None:
        """We assume any reaction that happens on the surface of grains can also happen.

        in the bulk (just more slowly due to binding energy). The user therefore only
        lists surface reactions in their input reaction file and we duplicate here.

        Raises
        ------
        RuntimeError
            If a CoupledReaction in the network has no partner set.

        """
        logger.debug("Adding bulk reactions")
        surface_reactions = self._get_reactions_on_grain()
        surface_reactions_can_be_bulk = [
            reaction
            for reaction in surface_reactions
            if reaction.get_reaction_type() in BULK_REACTION_TYPES
        ]

        species_names = self.network.species.keys()
        new_reactions: list[Reaction] = []
        for reaction in surface_reactions_can_be_bulk:
            new_reac = deepcopy(reaction)
            new_reac.convert_surf_to_bulk()
            if new_reac.has_unknown_species(species_names):
                msg = f"New bulk reaction '{new_reac}' has unknown species. Add to bulk, or call 'self._add_bulk_species' first."
                raise RuntimeError(msg)
            new_reac = CoupledReaction(new_reac)
            new_reac.set_partner(reaction)
            new_reactions.append(new_reac)

        current_reactions = self.network.reactions.values()
        new_reactions = [reac for reac in new_reactions if reac not in current_reactions]

        bulk_species = [x for x in self.network.species.values() if x.is_bulk_species()]
        new_reac_list: list[str | float] = []
        for species in bulk_species:
            # add individual swapping
            if not species.is_refractory:
                new_reac_list = [
                    species.get_name(),
                    "BULKSWAP",
                    "NAN",
                    species.get_name().replace("@", "#"),
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
            new_reac_list[0] = species.get_name().replace("@", "#")
            new_reac_list[1] = "SURFSWAP"
            new_reac_list[3] = species.get_name()
            new_reactions.append(Reaction(new_reac_list))
        logger.debug(
            f"The following bulk reactions are added to the reactions: {new_reactions}"
        )
        self.network.add_reactions(new_reactions)

    def _add_desorb_reactions(self) -> None:
        """Save the user effort by automatically generating desorption reactions.

        A DESORB reaction in the user input is a shorthand that expands to all four
        physical desorption mechanisms: THERM, DESOH2, DESCR, DEUVCR. After expansion
        the original DESORB reactions are removed from the network.

        If a surface species already has any explicit THERM/DESOH2/DESCR/DEUVCR
        (thermal / h2formation desorption / crir desorption / crir induced photodesorption)
        reaction defined in the network, auto-generation is skipped for that species.
        The user is then responsible for providing all desorption pathways, and the
        alpha values for each mechanism must sum to 1.0.

        Raises
        ------
        ValueError
            If a species has both a DESORB shorthand and explicit desorption reactions.
        ValueError
            If alpha values for explicit desorption reactions do not sum to 1.0.

        """
        desorb_reacs = ["DESOH2", "DESCR", "DEUVCR", "THERM"]
        logger.debug("Adding desorption reactions")

        # Expand DESORB shorthand into all four physical desorption mechanisms.
        desorb_shorthand = [
            r
            for r in self.network.reactions.values()
            if r.get_reaction_type() == "DESORB" and r.get_reactants()[0].startswith("#")
        ]
        if desorb_shorthand:
            # Check for conflicts: if any explicit THERM/DESOH2/DESCR/DEUVCR reaction
            # already exists for a species that also has a DESORB shorthand, raise an
            # error so the user doesn't accidentally have both.
            shorthand_species = {r.get_reactants()[0] for r in desorb_shorthand}
            for r in self.network.get_reaction_list():
                if (
                    r.get_reaction_type() in desorb_reacs
                    and r.get_reactants()[0] in shorthand_species
                ):
                    msg = (
                        f"{r.get_reactants()[0]} has both a DESORB shorthand reaction and an "
                        f"explicit {r.get_reaction_type()} reaction. Use DESORB (which expands "
                        f"to all four mechanisms) or define each mechanism explicitly — not both."
                    )
                    raise ValueError(msg)
            expanded: list[Reaction] = []
            for r in desorb_shorthand:
                reac, _, third = r.get_reactants()
                prods = r.get_products()
                expanded.extend(
                    Reaction(
                        [
                            reac,
                            mech,
                            third,
                            *prods,
                            r.get_alpha(),
                            r.get_beta(),
                            r.get_gamma(),
                            r.get_templow(),
                            r.get_temphigh(),
                            r.get_reduced_mass(),
                        ]
                    )
                    for mech in desorb_reacs
                )
                self.network.remove_reaction(r)
            self.network.add_reactions(expanded)

        # Species that already have at least one explicit desorption reaction
        existing_desorbs = {
            r.get_reactants()[0]
            for r in self.network.reactions.values()
            if r.get_reaction_type() in desorb_reacs
        }

        # For species with explicit desorption reactions, update their desorb_products
        # to the first product of the first explicit reaction so that gasIceList can
        # resolve them even when the standard gas counterpart is not in the species list.
        species_dict = self.network.species.copy()
        for species_name in existing_desorbs:
            first_rxn = next(
                r
                for r in self.network.get_reaction_list()
                if r.get_reactants()[0] == species_name
                and r.get_reaction_type() in desorb_reacs
            )
            first_product = first_rxn.get_products()[0]
            if species_name in species_dict and first_product not in {"NAN", ""}:
                species_dict[species_name].set_desorb_products(
                    [first_product, "NAN", "NAN", "NAN"]
                )
                # Also update the bulk (@) counterpart if it exists
                bulk_name = "@" + species_name[1:]
                if bulk_name in species_dict:
                    species_dict[bulk_name].set_desorb_products(
                        [first_product, "NAN", "NAN", "NAN"]
                    )

        # Validate that alpha sums to 1.0 per species per mechanism
        for species_name in existing_desorbs:
            for rtype in desorb_reacs:
                rxns = [
                    r
                    for r in self.network.reactions.values()
                    if r.get_reactants()[0] == species_name
                    and r.get_reaction_type() == rtype
                ]
                if rxns:
                    total_alpha = sum(r.get_alpha() for r in rxns)
                    if abs(total_alpha - 1.0) > 1e-6:  # ruff: ignore[magic-value-comparison]
                        msg = (
                            f"{species_name} has {rtype} reactions with alpha summing to "
                            f"{total_alpha:.6f}, expected 1.0. "
                            f"Branching ratios for each desorption mechanism must sum to 1."
                        )
                        raise ValueError(msg)

        new_reactions: list[Reaction] = []
        for species in self.network.species.values():
            if species.is_surface_species():
                if species.get_name() in existing_desorbs:
                    logger.debug(
                        f"Skipping auto-generation of desorption reactions for "
                        f"{species.get_name()} — user-defined explicit reactions found."
                    )
                    continue
                new_reactions.extend(
                    Reaction(
                        [
                            species.get_name(),
                            reacType,
                            "NAN",
                            *species.get_standard_desorb_products(),
                            1,
                            0,
                            species.get_binding_energy(),
                            0.0,
                            10000.0,
                            0.0,
                        ]
                    )
                    for reacType in desorb_reacs
                )
            if species.is_bulk_species() and not species.is_refractory:
                # Bulk species may have custom desorb_products set from FREEZE/DESORB
                # reactions (e.g., @HSIO -> SIOH+ for ion freeze-out). Use the custom
                # pathway if available, otherwise use standard gas desorption (strip @).
                if hasattr(species, "desorb_products"):
                    desorb_prods = species.get_desorb_products()
                else:
                    desorb_prods = species.get_standard_desorb_products()
                new_reactions.append(
                    Reaction(
                        [
                            species.get_name(),
                            "THERM",
                            "NAN",
                            *desorb_prods,
                            1,
                            0,
                            species.get_binding_energy(),
                            0.0,
                            10000.0,
                            0.0,
                        ]
                    )
                )
        self.network.add_reactions(new_reactions)

    def _add_chemdes_reactions(self) -> None:
        """Add chemical desorption reactions from LH and ER reactions.

        We have the user list all Langmuir-Hinshelwood and Eley-Rideal reactions once.
        Then we duplicate so that the reaction branches
        with products on grain and products desorbing.

        Raises
        ------
        NotImplementedError
            ChemDes reaction loading is not yet implemented.

        """
        logger.debug("Adding chemical desorption reactions for LH and ER mechanisms")
        new_reactions = []
        species_list = self.network.get_species_list()
        species_names = [species.name for species in species_list]
        for reaction in self.network.get_reaction_list():
            reactants = reaction.get_reactants()
            if reactants[0][0] == "@" or reactants[1][0] == "@":
                continue
            if reaction.get_reaction_type() in {"LH", "ER"}:
                n_products = sum(prod != "NAN" for prod in reaction.get_products())

                reaction_partner: Reaction = deepcopy(reaction)

                for i in range(n_products):
                    # For each of the products, make a new reaction where it is desorbed
                    new_reaction = deepcopy(reaction)
                    # Each LHDES variant carries 1/n_products of the parent LH rate,
                    # because the reaction energy can desorb any one of the products.
                    new_reaction.set_alpha(reaction.get_alpha() / n_products)

                    # Convert to disassociation reaction
                    new_reactants = new_reaction.get_reactants()
                    new_reactants[2] += "DES"
                    new_reaction.set_reactants(new_reactants)

                    # Replace the species on the grain/bulk with species in gas
                    new_products = new_reaction.get_products()

                    # Replace all grain or bulk products with gas phase counterparts.
                    # CRITICAL: Distinguish between chemical and physical desorption pathways.
                    #
                    # LHDES/ERDES (chemical desorption): Always use standard gas desorption
                    # (strip # or @ prefix). Chemical desorption via reaction exothermicity
                    # produces neutral species only, without ionization.
                    # Example: #HSIO + H (chemical) -> HSIO (gas, neutral)
                    #
                    # THERM/DESOH2/DESCR/DEUVCR (physical desorption): Use custom DESORB
                    # pathway if defined (set via FREEZE/DESORB reactions). This allows
                    # non-standard products (e.g., ionization during freeze-out).
                    # Example: #HSIO -> SIOH+ (via custom FREEZE/DESORB definition)
                    if new_reaction.get_reaction_type() in {"LHDES", "ERDES"}:
                        # Chemical desorption: strip grain/bulk prefix for standard gas product
                        gas_product = new_products[i][1:]  # Remove # or @ prefix
                        desorb_products = [gas_product, "NAN", "NAN", "NAN"]
                    else:
                        # Physical desorption: use custom desorb pathway
                        desorb_products = species_list[
                            species_names.index(new_products[i])
                        ].get_desorb_products()

                    if desorb_products[1] == "NAN":
                        new_products[i] = desorb_products[0]
                    elif desorb_products[2] == "NAN":
                        if i < 2 and new_products[i + 1] != "NAN":  # ruff: ignore[magic-value-comparison]
                            # Move i+1th product over to i+2th product
                            new_products[i + 2] = new_products[i + 1]
                        new_products[i + 1] = desorb_products[1]
                        new_products[i] = desorb_products[0]
                    else:
                        raise NotImplementedError()

                    new_reaction.set_products(new_products)

                    if new_reaction in self.network.get_reaction_list():
                        logger.warning(
                            f"Custom chemical desorption reaction {new_reaction} was added, not adding the default generated one to the network. The custom chemical desorption reaction is not made a CoupledReaction"
                        )
                        break

                    logger.debug(
                        f"Adding chemical desorption reaction for {reaction}, new reaction {new_reaction}"
                    )

                    new_reaction = CoupledReaction(new_reaction)
                    new_reaction.set_partner(reaction_partner)

                    new_reactions.append(new_reaction)
        self.network.add_reactions(new_reactions)

    def _add_excited_surface_reactions(self) -> None:
        """All excited species will relax to the ground state if they do not react.

        the vibrational frequency of the species is used as a pseudo approximation
        of the rate coefficient. We assume all grain reactions have an excited variant.
        For example:
        #A, #B LH #C will have the variants:
        #A*, #B EXSOLID #C  and  #A, #B* EXSOLID #C
        If only one of the reactants in the base reaction has an excited counterpart then
        only one excited version of that reaction is created.

        """
        logger.debug("Adding excited surface reactions")
        excited_species = [
            x for x in self.network.get_species_list() if "*" in x.get_name()
        ]
        lh_reactions = [
            x for x in self.network.get_reaction_list() if "LH" in x.get_reactants()
        ]
        lh_reactions += [
            x for x in self.network.get_reaction_list() if "LHDES" in x.get_reactants()
        ]

        new_reactions: list[Reaction] = []

        # add relaxation of excited species
        for spec in excited_species:
            relax_reac: list[str | float] = [
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
                new_reac_a_list: list[str | float] = [
                    reaction.get_reactants()[0] + "*",
                    reaction.get_reactants()[1],
                    "EXSOLID",
                ]
                new_reac_a_list = (
                    new_reac_a_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
                )
                new_reac_b_list: list[str | float] = [
                    reaction.get_reactants()[0],
                    reaction.get_reactants()[1] + "*",
                    "EXSOLID",
                ]
                new_reac_b_list = (
                    new_reac_b_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
                )

                new_reac_a = Reaction(new_reac_a_list)
                new_reac_b = Reaction(new_reac_b_list)

                # stops duplicate reactions e.g. #H* + #H and #H + #H*
                if new_reac_a != new_reac_b:
                    new_reactions.extend((new_reac_a, new_reac_b))
                else:
                    new_reactions.append(new_reac_a)

            # if only #A has an excited counterpart
            elif reaction.get_reactants()[0] + "*" in [
                specie.get_name() for specie in excited_species
            ]:
                new_reac_a_list = [
                    reaction.get_reactants()[0] + "*",
                    reaction.get_reactants()[1],
                    "EXSOLID",
                ]
                new_reac_a_list = (
                    new_reac_a_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
                )
                new_reac_a = Reaction(new_reac_a_list)
                new_reactions.append(new_reac_a)

            # if only #B has an excited counterpart
            elif reaction.get_reactants()[1] + "*" in [
                specie.get_name() for specie in excited_species
            ]:
                new_reac_b_list = [
                    reaction.get_reactants()[0],
                    reaction.get_reactants()[1] + "*",
                    "EXSOLID",
                ]
                new_reac_b_list = (
                    new_reac_b_list
                    + reaction.get_products()
                    + [reaction.get_alpha(), 0, 0, 0, 10000, 0.0]
                )
                new_reac_b = Reaction(new_reac_b_list)
                new_reactions.append(new_reac_b)
        self.network.add_reactions(new_reactions)

    def _add_CRP_and_PHOTO_reactions_to_grain(self) -> None:
        """Add all the gas-phase reactions with CRP, CRPHOT or PHOTON.

        to the grain surface too.

        """
        logger.info("Adding gas-phase reactions with CRP, CRPHOT or PHOTON to grain")
        reactions_on_grain = self._get_reactions_on_grain()
        reactions_on_grain_filtered = [
            reaction
            for reaction in reactions_on_grain
            if reaction.get_reaction_type() in {"CRP", "CRPHOT", "PHOTON"}
        ]
        new_reactions = []
        for reaction in self.network.get_reaction_list():
            if reaction.get_reaction_type() not in {"CRP", "CRPHOT", "PHOTON"}:
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
                logger.debug(
                    f"There a custom version of reaction '{reaction}' on the grain (with same reactants) already in the network. Skipping adding the copied one."
                )
                continue

            if any("+" in reactant for reactant in reactants):
                logger.debug(
                    f"Reaction '{reaction}' had reactant ions, which is not possible on grain. Skipping"
                )
                continue

            # We have now filtered to have only gas-phase reactions that have type
            # CRP, CRPHOT or PHOTON
            new_reaction = deepcopy(reaction)
            new_reaction.convert_gas_to_surf()
            reactants, products = (
                new_reaction.get_reactants(),
                new_reaction.get_products(),
            )
            if reactants[0] == products[0]:
                # This means the reaction was simply an ionization reaction.
                # We can skip this reaction.
                logger.debug(
                    f"Reaction {reaction} is ionization reaction, skip on grain surface"
                )
                continue
            if not all(
                species in self.network.get_species_list()
                or species in {"NAN", "CRP", "CRPHOT", "PHOTON"}
                for species in reactants + products
            ):
                logger.debug(
                    f"Reaction copied from gas to grain '{reaction}' contains species that were not set on grain. Skipping"
                )
                continue

            new_reaction = CoupledReaction(new_reaction)
            new_reaction.set_partner(reaction)
            new_reactions.append(new_reaction)
        logger.debug("Adding new reactions to grain")
        self.network.add_reactions(new_reactions)
        logger.info(f"Added {len(new_reactions)} reactions to grain")

    def _branching_ratios_checks(self) -> None:
        """Check that the branching ratios for the ice reactions sum to 1.0.

        If they do not, correct them. This needs to be done for LH and LHDES
        separately since we already added the desorption to the network.

        """
        branching_reactions: dict[str, float] = {}
        for reaction in self.network.get_reaction_list():
            if isinstance(reaction, CoupledReaction):
                continue
            if reaction.get_reaction_type() in TUNNELING_REACTION_TYPES:
                reactant_string = ",".join(reaction.get_reactants())
                if reactant_string in branching_reactions:
                    branching_reactions[reactant_string] += reaction.get_alpha()
                else:
                    branching_reactions[reactant_string] = reaction.get_alpha()

        if any(val != 1 for val in branching_reactions.values()):
            logger.warning(
                "Some of the branching ratios do not sum to 1.0, correcting those that do not"
            )
            for reaction in self.network.get_reaction_list():
                if isinstance(reaction, CoupledReaction):
                    continue
                if reaction.get_reaction_type() in TUNNELING_REACTION_TYPES:
                    reactant_string = ",".join(reaction.get_reactants())
                    # Check if we need to correct the branching ratio
                    # (smaller than 0.98 is allowed)
                    if (
                        reactant_string in branching_reactions
                        and branching_reactions[reactant_string] != 1
                    ):
                        if branching_reactions[reactant_string] == 0:
                            logger.warning(
                                f"Grain reaction {reaction} has a branching ratio of 0.0, removing the reaction altogether"
                            )
                            self.network.remove_reaction(reaction)
                            continue

                        if branching_reactions[reactant_string] < 0.99:  # ruff: ignore[magic-value-comparison]
                            logger.warning(
                                f"You have reaction {reaction} with a branching ratio {branching_reactions[reactant_string]} we are assuming you set this lower on purpose."
                            )
                            continue

                        new_alpha = (
                            reaction.get_alpha() / branching_reactions[reactant_string]
                        )
                        logger.warning(
                            f"Grain reaction {reaction} has a branching ratio of {reaction.get_alpha()}, dividing it by {branching_reactions[reactant_string]} resulting in BR of {new_alpha}"
                        )

                        reaction_index = self.network.get_reaction_index(reaction)
                        reaction.set_alpha(new_alpha)
                        self.network.set_reaction(
                            reaction_idx=reaction_index, reaction=reaction
                        )

                        for coupled_reaction in self.network.get_all_partners(reaction):
                            reaction_index = self.network.get_reaction_index(
                                coupled_reaction
                            )
                            coupled_reaction.set_alpha(new_alpha)
                            self.network.set_reaction(
                                reaction_idx=reaction_index, reaction=coupled_reaction
                            )

    def _add_gas_phase_extrapolation(self) -> None:
        """Enable extrapolation for gas-phase reactions that.

        have unique or overlapping temperature ranges.

        """
        for reaction in self.network._reactions_dict.values():
            if reaction.is_gas_reaction() and (
                reaction.get_reaction_type() in {"TWOBODY", "PHOTON", "CRP", "CRPHOT"}
            ):
                similar_reactions = self.network.find_similar_reactions(reaction)
                # Only enable extrapolation if we have one or overlapping reactions
                # UMIST uses overlapping reactions to get more correct reaction rates.
                if all(
                    (reaction.get_templow() == v.get_templow())
                    and (reaction.get_temphigh() == v.get_temphigh())
                    for v in similar_reactions.values()
                ):
                    reaction.set_extrapolation(enabled=True)

    def _add_reaction_enthalpies(self, enthalpy_reaction_types: str | set[str]) -> None:
        """Add reaction enthalpies (exothermicity) to reactions.

        for heating/cooling calculations.

        Parameters
        ----------
        enthalpy_reaction_types : str | set[str]
            Set of reaction types or "ALL" or "GAS"

        """
        exclude_ices = True
        if isinstance(enthalpy_reaction_types, str):
            if enthalpy_reaction_types.upper() == "ALL":
                exclude_ices = False
                enthalpy_reaction_types = set(REACTION_TYPES)
            elif enthalpy_reaction_types.upper() == "GAS":
                enthalpy_reaction_types = set(REACTION_TYPES)

        for reaction in self.network.reactions.values():
            logger.debug(f"Checking if we need to add enthalpy to {reaction}")
            if reaction.get_reaction_type() in enthalpy_reaction_types:
                if exclude_ices and reaction.is_ice_reaction(strict=(not exclude_ices)):
                    logger.debug("Skipping ice reaction")
                    continue
                if "E-" in (reaction.get_reactants() + reaction.get_products()):
                    logger.debug(
                        "Reaction involving electrons, skipping enthalpy due to poor estimates"
                    )
                    continue
                delta_h = self._compute_exothermicity(reaction)
                # TODO: add a heating efficiency factor in here:
                reaction.set_exothermicity(convert_to_erg(-delta_h, "kcal/mol"))
                logger.debug(
                    f"Setting reaction enthalpy of {reaction} to {delta_h} kcal/mol"
                )

    def _apply_custom_exothermicities(
        self, database_reaction_exothermicity: list[str | Path]
    ) -> None:
        """Apply custom exothermicity values from CSV files to the network reactions.

        Parameters
        ----------
        database_reaction_exothermicity : list[str | Path]
            List of paths
            to custom exothermicity CSV files.

        """
        for csv_path in database_reaction_exothermicity:
            logger.info(f"Applying custom exothermicities from {csv_path}")
            set_custom_exothermicities(
                reactions=self.network.get_reaction_list(),
                csv_path=csv_path,
                overwrite=True,
            )

    def _compute_exothermicity(self, reaction: Reaction) -> float:
        """Compute the reaction enthalpy in eV for a given reaction based on the.

        species enthalpies.

        Parameters
        ----------
        reaction : Reaction
            The reaction to compute the enthalpy for.

        Returns
        -------
        float
            The reaction enthalpy in kcal/mol.

        """
        reactants = reaction.get_pure_reactants()
        products = reaction.get_pure_products()
        return sum(self.network._species_dict[p].get_enthalpy() for p in products) - sum(
            self.network._species_dict[r].get_enthalpy() for r in reactants
        )

    def _get_reactions_on_grain(self) -> list[Reaction]:
        """Get all reactions that occur on grain surfaces (# prefix) or in bulk (@prefix).

        Returns
        -------
        reactions_on_grain : list[Reaction]
            All reactions occurring on the grain.

        """
        reactions_on_grain = [
            reaction
            for reaction in self.network.reactions.values()
            if reaction.is_ice_reaction(include_products=False, strict=False)
        ]
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
        """Check that every species freezes out and alert the user if a.

        species freezes out via multiple routes. This isn't necessarily an
        error so best just print.

        """
        logger.info(
            "\tCheck that species have surface counterparts or if they have multiple freeze outs/check alphas:\n"
        )
        for name, spec in self.network.species.items():
            if not spec.is_ice_species() and not spec.is_ion():
                exist_check = 0
                for check_speck in self.network.species:
                    if check_speck == "#" + spec.get_name():
                        exist_check += 1
                        break
                if exist_check == 0:
                    logger.warning(
                        f"{spec.get_name()} does not have a surface counterpart in given default species file."
                        + "\n\tThis sets the binding energy to zero, it might cause species conservation errors."
                    )

            freezeout_reactions = [
                reaction
                for reaction in self.network.reactions.values()
                if (
                    name in reaction.get_reactants()
                    and reaction.get_reaction_type() == "FREEZE"
                )
            ]
            freezes = len(freezeout_reactions)
            if freezes == 1:
                logger.info(
                    f"\t{spec.get_name()} freezes out through {freezeout_reactions[0]}"
                )
            if freezes > 1:
                logger.info(f"\t{spec.get_name()} freezes out through {freezes} routes")
            elif freezes < 1 and not spec.is_ice_species():
                logger.info(f"\t{spec.get_name()} does not freeze out")

    def _duplicate_checks(self) -> None:
        """Check reaction network to make sure no reaction appears twice unless.

        they have different temperature ranges.

        """
        logger.info("\tPossible duplicate reactions for manual removal:")
        duplicates_found = False

        # We could also pass `self.network.reactions.values()` to check_duplcicate_reactions,
        # but that would lead to n_reac*(n_reac-1)/2 comparisons,
        # whereas this splits up the comparison into buckets of reaction types,
        # because reactions with different types are never duplicates.
        # So, fewer comparisons are necessary, so this is much faster.
        for reaction_type in REACTION_TYPES:
            reactions = self.network.get_reactions_by_types(reaction_type)
            duplicates, duplicates_umist = get_duplicate_reactions(reactions)
            for duplicate in duplicates:
                duplicates_found = True
                reaction1 = reactions[duplicate[0]]
                reaction2 = reactions[duplicate[1]]
                logger.warning(
                    "\tFound reactions that are possible duplicates\n\t\t"
                    + str(reaction1)
                    + f" with temperature range [{reaction1.get_templow()}, {reaction1.get_temphigh()}] and source {reaction1.get_source()}"
                    + "\n\t\t"
                    + str(reaction2)
                    + f" with temperature range [{reaction2.get_templow()}, {reaction2.get_temphigh()}] and source {reaction2.get_source()}"
                )
                if (
                    reaction1.get_temphigh() > reaction2.get_temphigh()
                    and reaction1.get_templow() < reaction2.get_temphigh()
                ):
                    logger.warning(
                        f"\tReactions {reaction1} and {reaction2} have non-adjacent temperature ranges"
                    )
            for duplicate in duplicates_umist:
                reaction1 = reactions[duplicate[0]]
                logger.info(
                    f"Detected overlapping UMIST reaction {reaction1}. This is done in UMIST to provide better rates. "
                )
        if not duplicates_found:
            logger.info("\tNone")

    def _index_important_reactions(self) -> None:
        """We have a whole bunch of important reactions and we want to store.

        their indices. We find them all here.

        Raises
        ------
        RuntimeError
            If an important reaction is found twice, or one of the important reactions
            is not found.

        """
        # Any None values in dictionary will raise an error
        # therefore these reactions are mandatory and
        # makerates will not complete if the user doesn't supply them.
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
                "nR_C_hv": lambda reacs, prods: ("C" in reacs) and ("PHOTON" in reacs),  # ruff: ignore[unused-lambda-argument]
                "nR_H2Form_CT": lambda reacs, prods: "H2FORM" in reacs,  # ruff: ignore[unused-lambda-argument]
                "nR_H2Form_ERDes": lambda reacs, prods: (
                    ("H" in reacs) and ("#H" in reacs) and ("H2" in prods)
                ),
                "nR_H2Form_ER": lambda reacs, prods: (
                    ("H" in reacs) and ("#H" in reacs) and ("#H2" in prods)
                ),
                "nR_H2Form_LH": lambda reacs, prods: (  # ruff: ignore[unused-lambda-argument]
                    (reacs.count("#H") == 2) and ("LH" in reacs)  # ruff: ignore[magic-value-comparison]
                ),
                "nR_H2Form_LHDes": lambda reacs, prods: (  # ruff: ignore[unused-lambda-argument]
                    (reacs.count("#H") == 2) and ("LHDES" in reacs)  # ruff: ignore[magic-value-comparison]
                ),
                "nR_HFreeze": lambda reacs, prods: ("H" in reacs) and ("FREEZE" in reacs),  # ruff: ignore[unused-lambda-argument]
                "nR_H2Freeze": lambda reacs, prods: (  # ruff: ignore[unused-lambda-argument]
                    ("H2" in reacs) and ("FREEZE" in reacs)
                ),
                "nR_EFreeze": lambda reacs, prods: (  # ruff: ignore[unused-lambda-argument]
                    ("E-" in reacs) and ("FREEZE" in reacs)
                ),
                "nR_H2_hv": lambda reacs, prods: ("H2" in reacs) and ("PHOTON" in reacs),  # ruff: ignore[unused-lambda-argument]
                "nR_H2_crp": lambda reacs, prods: (
                    ("H2" in reacs) and ("CRP" in reacs) and (prods.count("H") == 2)  # ruff: ignore[magic-value-comparison]
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
                        msg = f"When trying to index the important reactions, we found a disastrous reaction {reaction} is a duplicate of {self.network.important_reactions[key]}, there can only be one reaction that matches {key}"
                        raise RuntimeError(msg)
                    self.network.important_reactions[key] = i + 1

        if np_any([value is None for value in self.network.important_reactions.values()]):
            logger.debug(self.network.important_reactions)
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
            *elementList,
        ]:
            try:
                species_index = names.index(element) + 1
            except ValueError:
                # TODO(Tobias Dijkhuis): The dummy value is currently SURFACE/BULK; https://github.com/uclchem/UCLCHEM/issues/205
                # We could handle this better somehow
                logger.info(f"\t{element} not in network, adding dummy index")
                species_index = len(self.network.get_species_list()) + 1
            name = "n" + element.lower().replace("+", "x").replace("e-", "elec").replace(
                "#", "g"
            )
            self.network.species_indices[name] = species_index
