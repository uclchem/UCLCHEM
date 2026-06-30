import uclchem.makerates
import pytest


def test_get_isotope_reactions():
    network = uclchem.advanced.RuntimeNetwork()
    reactions = network.get_reaction_list()
    for reaction in reactions:
        if not reaction.contains_species("H"):
            continue
        reaction.get_reactions_with_other_isotopes("H", "D")
