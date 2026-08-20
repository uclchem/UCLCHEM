import numpy as np
import pytest

from uclchem.makerates.network import Network
from uclchem.makerates.network_builder import NetworkBuilder
from uclchem.makerates.reaction import CoupledReaction, Reaction
from uclchem.makerates.species import Species

_dummy_reaction_data = ["NAN"] * 10
_dummy_species_data = [0] * 10
_rng = np.random.default_rng()

add_bulk_species_data = [
    ([Species(["H", *_dummy_species_data]), Species(["#H2O", *_dummy_species_data])]),
    (
        [
            Species(["H", *_dummy_species_data]),
            Species(["#H", *_dummy_species_data]),
            Species(["#H2O", 0, _rng.random(), *_dummy_species_data]),
        ]
    ),
]


def _initialize_network(species: list[Species], reactions: list[Reaction]) -> Network:
    builder = NetworkBuilder(species, reactions)
    builder.network = Network.from_lists(species, reactions)
    builder._add_bulk_species()
    builder._add_bulk_reactions()
    builder._add_freeze_reactions()
    builder._add_freeze_reactions()
    return builder.network


add_bulk_reactions_species = [
    Species(["H", *_dummy_species_data]),
    Species(["#H", *_dummy_species_data]),
    Species(["#H2", *_dummy_species_data]),
    Species(["#H2O", *_dummy_species_data]),
]


def test_get_all_partners() -> None:
    original_reaction = Reaction(["#H", "#H", "LH", "#H2", *_dummy_reaction_data])
    coupled_reaction_bulk = CoupledReaction(
        ["@H", "@H", "LH", "@H2", *_dummy_reaction_data]
    )
    coupled_reaction_bulk.set_partner(original_reaction)
    coupled_reaction_lhdes = CoupledReaction(
        ["#H", "#H", "LHDES", "H2", *_dummy_reaction_data]
    )
    coupled_reaction_lhdes.set_partner(original_reaction)
    coupled_reaction_bulk_lhdes = CoupledReaction(
        ["@H", "@H", "LHDES", "H2", *_dummy_reaction_data]
    )
    coupled_reaction_bulk_lhdes.set_partner(coupled_reaction_lhdes)

    coupled_reactions: list[CoupledReaction] = [
        coupled_reaction_bulk,
        coupled_reaction_lhdes,
        coupled_reaction_bulk_lhdes,
    ]
    reactions = [original_reaction, *coupled_reactions]

    network = _initialize_network(add_bulk_reactions_species, reactions)
    reactions_in_network = list(network.reactions.values())
    for reaction in reactions_in_network:
        if reaction == original_reaction:
            assert network.get_all_partners(reaction) == coupled_reactions
        elif reaction in coupled_reactions:
            reaction: CoupledReaction = reaction  # type: ignore[no-redef, ty:invalid-assignment] # Making mypy happy
            assert isinstance(reaction, CoupledReaction)
            assert reaction.partner == original_reaction


add_bulk_reactions_data = [
    (
        [
            Reaction(["H", "H", "NAN", "H2", *_dummy_reaction_data]),
            Reaction(
                [
                    "#H",
                    "#H",
                    "LH",
                    "#H2",
                    "NAN",
                    "NAN",
                    "NAN",
                    _rng.random(),
                    _rng.random(),
                    _rng.random(),
                    *_dummy_reaction_data,
                ]
            ),
        ]
    )
]


@pytest.mark.parametrize("reactions", add_bulk_reactions_data)
def test_change_reaction_barrier(reactions):
    network = _initialize_network(add_bulk_reactions_species, reactions)

    reactions_in_network = list(network.reactions.values())
    for reaction in reactions_in_network:
        if not reaction.is_bulk_reaction(strict=True):
            continue
        surface_reactants = [
            reac.replace("@", "#") for reac in reaction.get_pure_reactants()
        ]
        surface_products = [
            reac.replace("@", "#") for reac in reaction.get_pure_products()
        ]
        surface_reaction = [
            reac
            for reac in reactions_in_network
            if reac.get_pure_reactants() == surface_reactants
            and reac.get_pure_products() == surface_products
        ]
        assert len(surface_reaction) == 1
        surface_reaction = surface_reaction[0]

        assert reaction in network.get_all_partners(surface_reaction)

        network.change_reaction_barrier(surface_reaction, _rng.random())
        assert reaction.get_gamma() == surface_reaction.get_gamma()
