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
reactions = [
    Reaction(["H", "H", "NAN", "H2", *_dummy_reaction_data]),
]


def _initialize_network(
    species: list[Species], reactions: list[Reaction]
) -> NetworkBuilder:
    builder = NetworkBuilder(species, reactions)

    # If we call builder.build(), the whole network will automatically be built,
    # but we want to test specific parts individually.
    builder.network = Network.from_lists(builder.input_species, builder.input_reactions)

    return builder


@pytest.mark.parametrize(("species"), add_bulk_species_data)
def test_add_bulk_species(species: list[Species]):
    builder = _initialize_network(species, reactions)
    builder._add_bulk_species()

    n_surface_species = sum(spec.is_surface_species() for spec in species)
    bulk_species = [
        spec for spec in builder.network.species.values() if spec.is_bulk_species()
    ]
    assert len(bulk_species) == n_surface_species

    # Make sure that all bulk species have the same binding energy as water
    water_index = [spec.name for spec in species].index("#H2O")
    water_binding_energy = species[water_index].binding_energy
    for spec in bulk_species:
        assert spec.binding_energy == water_binding_energy


add_bulk_reactions_species = [
    Species(["H", *_dummy_species_data]),
    Species(["#H", *_dummy_species_data]),
    Species(["#H2", *_dummy_species_data]),
    Species(["#H2O", *_dummy_species_data]),
]

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
def test_add_bulk_reactions_without_bulk_species_raises(reactions: list[Reaction]):
    builder = _initialize_network(add_bulk_reactions_species, reactions)

    with pytest.raises(RuntimeError, match="New bulk reaction"):
        builder._add_bulk_reactions()


@pytest.mark.parametrize("reactions", add_bulk_reactions_data)
def test_add_bulk_reactions(reactions: list[Reaction]):
    builder = _initialize_network(add_bulk_reactions_species, reactions)

    builder._add_bulk_species()
    builder._add_bulk_reactions()

    n_surface_reactions = sum(
        reaction.is_surface_reaction(strict=True) for reaction in reactions
    )
    bulk_reactions = [
        reaction
        for reaction in builder.network.reactions.values()
        if reaction.is_bulk_reaction(strict=True)
    ]
    assert len(bulk_reactions) == n_surface_reactions

    reactions_in_network = list(builder.network.reactions.values())
    for reaction in bulk_reactions:
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
        surface_reaction: Reaction = surface_reaction[0]  # type: ignore[no-redef]
        assert isinstance(surface_reaction, Reaction)  # Making mypy happy

        assert isinstance(reaction, CoupledReaction)
        assert reaction.partner == surface_reaction
        assert reaction.get_gamma() == surface_reaction.get_gamma()
        assert reaction.get_alpha() == surface_reaction.get_alpha()
        assert reaction.get_beta() == surface_reaction.get_beta()
