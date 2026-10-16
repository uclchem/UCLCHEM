import numpy as np
import pytest
import uclchemwrap

import uclchem
import uclchem.constants
from uclchem.makerates.reaction import Reaction
from uclchem.makerates.species import Species

_rng = np.random.default_rng()
_network = uclchem.advanced.RuntimeNetwork()

surface_index = uclchemwrap.network.nsurface - 1
bulk_index = uclchemwrap.network.nbulk - 1


assert uclchemwrap.network.specname[surface_index].strip() == b"SURFACE"
assert uclchemwrap.network.specname[bulk_index].strip() == b"BULK"


def test_surfgrowthuncorrected_shrink():
    rate_constants = np.zeros(uclchem.constants.n_reactions)
    abundances = np.zeros(uclchem.constants.n_species + 2)

    desorb_h = Reaction(["#H", "THERM", "NAN", "H", *["NAN"] * 8])
    reaction_idx = _network.get_reaction_index(desorb_h)
    species_list = _network.get_species_list()
    gas_h_index = species_list.index(Species(["H", *[0] * 12]))
    surf_h_index = species_list.index(Species(["#H", *[0] * 12]))

    rate_constants[reaction_idx] = _rng.random()
    abundances[:] = _rng.random(uclchem.constants.n_species + 2)

    surface_index = uclchemwrap.network.nsurface - 1
    bulk_index = uclchemwrap.network.nbulk - 1

    reaction_rate = rate_constants[reaction_idx] * abundances[surf_h_index]

    surface_coverage = uclchemwrap.surfacereactions.bulkgainfrommantlebuildup()
    density = 1e5

    ydot, surfgrowthuncorrected = uclchemwrap.odes.getydot(
        rate_constants,
        abundances,
        surface_coverage,
        density,
    )

    assert ydot[gas_h_index] == reaction_rate
    assert surfgrowthuncorrected == -reaction_rate

    # surfgrowthuncorrected is less than 0, so bulk should shrink to compensate
    for species_idx, species in enumerate(species_list):
        if species.is_bulk_species():
            assert ydot[species_idx] < 0

    assert ydot[surface_index] + ydot[bulk_index] == pytest.approx(-reaction_rate)

    # The surface should not shrink at all, because bulk compensates for its loss
    assert ydot[surface_index] == pytest.approx(0)
    assert ydot[bulk_index] == pytest.approx(-reaction_rate)


def test_surfgrowthuncorrected_growth():
    rate_constants = np.zeros(uclchem.constants.n_reactions)
    abundances = np.zeros(uclchem.constants.n_species + 2)

    freeze_h = Reaction(["H", "FREEZE", "NAN", "#H", *["NAN"] * 8])
    reaction_idx = _network.get_reaction_index(freeze_h)
    species_list = _network.get_species_list()
    gas_h_index = species_list.index(Species(["H", *[0] * 12]))
    surf_h_index = species_list.index(Species(["#H", *[0] * 12]))

    rate_constants[reaction_idx] = _rng.random()
    abundances[gas_h_index] = 1

    surface_coverage = uclchemwrap.surfacereactions.bulkgainfrommantlebuildup()
    density = 1e5

    reaction_rate = rate_constants[reaction_idx] * abundances[gas_h_index] * density

    ydot, surfgrowthuncorrected = uclchemwrap.odes.getydot(
        rate_constants,
        abundances,
        surface_coverage,
        density,
    )

    assert ydot[surf_h_index] == reaction_rate
    assert surfgrowthuncorrected == reaction_rate

    for species_idx, species in enumerate(species_list):
        if species.is_bulk_species():
            assert ydot[species_idx] == 0

    assert ydot[surface_index] + ydot[bulk_index] == pytest.approx(reaction_rate)

    assert ydot[surface_index] == pytest.approx(reaction_rate)
    assert ydot[bulk_index] == pytest.approx(0)
