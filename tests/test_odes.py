import numpy as np
import pytest
import uclchemwrap

import uclchem
import uclchem.constants
from uclchem.makerates.io_functions import MIN_BURIABLE_MANTLE_FRACTION
from uclchem.makerates.reaction import Reaction
from uclchem.makerates.species import Species

_rng = np.random.default_rng()
_network = uclchem.advanced.RuntimeNetwork()

surface_index = uclchemwrap.network.nsurface - 1
bulk_index = uclchemwrap.network.nbulk - 1


assert uclchemwrap.network.specname[surface_index].strip() == b"SURFACE"
assert uclchemwrap.network.specname[bulk_index].strip() == b"BULK"


@pytest.mark.parametrize("bulk_scale", [2.0, 0.5], ids=["bulk_thicker", "bulk_thinner"])
def test_surfgrowthuncorrected_shrink(bulk_scale):
    rate_constants = np.zeros(uclchem.constants.n_reactions)
    abundances = np.zeros(uclchem.constants.n_species + 2)

    desorb_h = Reaction(["#H", "THERM", "NAN", "H", *["NAN"] * 8])
    reaction_idx = _network.get_reaction_index(desorb_h)
    species_list = _network.get_species_list()
    gas_h_index = species_list.index(Species(["H", *[0] * 12]))
    surf_h_index = species_list.index(Species(["#H", *[0] * 12]))
    surface_species = [i for i, s in enumerate(species_list) if s.is_surface_species()]
    bulk_species = [i for i, s in enumerate(species_list) if s.is_bulk_species()]

    rng = np.random.default_rng(seed=195)
    rate_constants[reaction_idx] = rng.random()
    abundances[:] = rng.random(uclchem.constants.n_species + 2)
    # Fix the total bulk relative to the total surface, so that both regimes of the
    # Garrod & Pauly (2011) transfer are tested deterministically.
    total_surface = abundances[surface_species].sum()
    abundances[bulk_species] *= (
        bulk_scale * total_surface / abundances[bulk_species].sum()
    )

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
    for species_idx in bulk_species:
        assert ydot[species_idx] < 0

    assert ydot[surface_index] + ydot[bulk_index] == pytest.approx(-reaction_rate)

    # Garrod & Pauly (2011): the fraction of lost surface that is replenished from
    # the bulk is min(1, N_bulk / N_surface). A bulk thicker than the surface fully
    # compensates the loss; a thinner bulk only compensates it partially.
    replenished_fraction = min(1.0, bulk_scale)
    assert ydot[bulk_index] == pytest.approx(-replenished_fraction * reaction_rate)
    assert ydot[surface_index] == pytest.approx(
        -(1.0 - replenished_fraction) * reaction_rate, abs=1e-12
    )


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


@pytest.mark.parametrize("coverage_times_mantle", [0.5, 2.0])
@pytest.mark.parametrize("h2_fraction", [0.3, 0.95])
def test_surface_growth_does_not_bury_h2(coverage_times_mantle, h2_fraction):
    rate_constants = np.zeros(uclchem.constants.n_reactions)
    abundances = np.zeros(uclchem.constants.n_species + 2)

    freeze_h = Reaction(["H", "FREEZE", "NAN", "#H", *["NAN"] * 8])
    reaction_idx = _network.get_reaction_index(freeze_h)
    species_list = _network.get_species_list()
    gas_h_index = species_list.index(Species(["H", *[0] * 12]))
    surf_h2_index = species_list.index(Species(["#H2", *[0] * 12]))
    bulk_h2_index = species_list.index(Species(["@H2", *[0] * 12]))
    surface_species = [i for i, s in enumerate(species_list) if s.is_surface_species()]
    bulk_species = [i for i, s in enumerate(species_list) if s.is_bulk_species()]
    other_surface_species = [i for i in surface_species if i != surf_h2_index]

    surface_coverage = uclchemwrap.surfacereactions.bulkgainfrommantlebuildup()
    total_surface = coverage_times_mantle / surface_coverage

    rng = np.random.default_rng(seed=195)
    rate_constants[reaction_idx] = rng.random()
    abundances[gas_h_index] = 1
    abundances[other_surface_species] = rng.random(len(other_surface_species))
    abundances[other_surface_species] *= (
        (1 - h2_fraction) * total_surface / abundances[other_surface_species].sum()
    )
    abundances[surf_h2_index] = h2_fraction * total_surface
    abundances[bulk_species] = rng.random(len(bulk_species)) * total_surface
    density = 1e5

    growth = rate_constants[reaction_idx] * abundances[gas_h_index] * density

    ydot, surfgrowthuncorrected = uclchemwrap.odes.getydot(
        rate_constants,
        abundances,
        surface_coverage,
        density,
    )
    assert surfgrowthuncorrected == pytest.approx(growth)

    # H2 is never buried
    assert ydot[surf_h2_index] == 0
    assert ydot[bulk_h2_index] == 0

    # The total amount that is buried is the same as when H2 could be buried, so the
    # surface stays capped, unless the surface is almost entirely H2.
    buried_fraction = min(1.0, coverage_times_mantle)
    buriable_fraction = 1 - h2_fraction
    buried_fraction *= buriable_fraction / max(
        buriable_fraction, MIN_BURIABLE_MANTLE_FRACTION
    )
    assert ydot[bulk_index] == pytest.approx(buried_fraction * growth)
    assert ydot[surface_index] == pytest.approx((1 - buried_fraction) * growth)
