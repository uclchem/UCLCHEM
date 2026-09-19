from contextlib import nullcontext

import pytest

from uclchem.makerates.reaction import CoupledReaction, Reaction, skip_reaction_validation

dummy_data = ["NAN"] * 8

with skip_reaction_validation():
    conserves_charge_data = [
        (Reaction(["H", "H", "NAN", "H2", *dummy_data]), nullcontext()),
        (Reaction(["H", "H", "NAN", "H", *dummy_data]), nullcontext()),
        (
            Reaction(["H+", "H", "NAN", "H2", *dummy_data]),
            pytest.raises(ValueError, match="Charges not conserved"),
        ),
        (Reaction(["H+", "H-", "NAN", "H2", *dummy_data]), nullcontext()),
    ]

    conserves_elements_data = [
        (Reaction(["H", "H", "NAN", "H2", *dummy_data]), nullcontext()),
        (
            Reaction(["H", "H", "NAN", "H", *dummy_data]),
            pytest.raises(ValueError, match="Elements not conserved"),
        ),
        (Reaction(["H+", "H", "NAN", "H2", *dummy_data]), nullcontext()),
        (Reaction(["H+", "H-", "NAN", "H2", *dummy_data]), nullcontext()),
        (Reaction(["CH3OH", "#H", "LH", "CH2OH", "H2", *dummy_data]), nullcontext()),
    ]

    reaction_type_possible_data = [
        (
            Reaction(["H", "H", "TWOBODY", "H2", *dummy_data]),
            nullcontext(),
        ),
        (
            Reaction(["#H", "H", "TWOBODY", "H2", *dummy_data]),
            pytest.raises(
                ValueError,
                match="TWOBODY reactions must happen in the gas-phase",
            ),
        ),
        (
            Reaction(["H", "H", "TWOBODY", "#H2", *dummy_data]),
            pytest.raises(
                ValueError,
                match="TWOBODY reactions must happen in the gas-phase",
            ),
        ),
        (
            Reaction(["H", "H", "TWOBODY", "H", "#H", *dummy_data]),
            pytest.raises(
                ValueError, match="TWOBODY reactions must happen in the gas-phase"
            ),
        ),
        (
            Reaction(["H2", "TWOBODY", "NAN", "H", "H", *dummy_data]),
            pytest.raises(
                ValueError,
                match="Reactions with type 'TWOBODY' should have two reactants",
            ),
        ),
        (
            Reaction(["#H2", "LH", "NAN", "#H", "#H", *dummy_data]),
            pytest.raises(
                ValueError, match="Reactions with type 'LH' should have two reactants"
            ),
        ),
        (
            Reaction(["#H", "THERM", "NAN", "H", *dummy_data]),
            nullcontext(),
        ),
        (
            Reaction(["#H", "THERM", "NAN", "#H", *dummy_data]),
            pytest.raises(
                ValueError, match="THERM reactions must have all reactants on ice"
            ),
        ),
        (
            Reaction(["H", "FREEZE", "NAN", "#H", *dummy_data]),
            nullcontext(),
        ),
        (
            Reaction(["#H", "FREEZE", "NAN", "#H", *dummy_data]),
            pytest.raises(
                ValueError,
                match="FREEZE reactions must have all reactants in gas and products on ice",
            ),
        ),
        (
            Reaction(["H", "FREEZE", "NAN", "H", *dummy_data]),
            pytest.raises(
                ValueError,
                match="FREEZE reactions must have all reactants in gas and products on ice",
            ),
        ),
        (
            Reaction(["#H", "#H", "LH", "#H2", *dummy_data]),
            nullcontext(),
        ),
        (
            Reaction(["H", "#H", "LH", "#H2", *dummy_data]),
            pytest.raises(ValueError, match="LH reactions must happen fully on the ice"),
        ),
        (
            Reaction(["#H", "#H", "LH", "H2", *dummy_data]),
            pytest.raises(ValueError, match="LH reactions must happen fully on the ice"),
        ),
        (
            Reaction(["#H", "#CH3OH", "LH", "#H2", "CH2O", *dummy_data]),
            pytest.raises(ValueError, match="LH reactions must happen fully on the ice"),
        ),
        (
            Reaction(["H", "#CH3OH", "ER", "#H2", "#CH2O", *dummy_data]),
            nullcontext(),
        ),
        (
            Reaction(["#H", "#CH3OH", "ER", "#H2", "#CH2O", *dummy_data]),
            pytest.raises(ValueError, match="ER reactions must have one reactant"),
        ),
        (
            Reaction(["H", "CH3OH", "ER", "#H2", "#CH2O", *dummy_data]),
            pytest.raises(ValueError, match="ER reactions must have one reactant"),
        ),
        (
            Reaction(["H", "H", "FREEZE", "#H2", *dummy_data]),
            pytest.raises(
                ValueError,
                match="Reactions with type 'FREEZE' should have only one reactant",
            ),
        ),
    ]


@pytest.mark.parametrize(("reaction", "expected"), conserves_charge_data)
def test_reaction_conserves_charge(reaction: Reaction, expected):
    with expected:
        reaction.check_charge_conservation()


@pytest.mark.parametrize(("reaction", "expected"), conserves_elements_data)
def test_reaction_conserves_elements(reaction: Reaction, expected):
    with expected:
        reaction.check_element_conservation()


@pytest.mark.parametrize(("reaction", "expected"), reaction_type_possible_data)
def test_all_on_grain(reaction: Reaction, expected):
    with expected:
        reaction.check_reaction_type_is_possible()


def test_get_original_partner() -> None:
    original_reaction = Reaction(["#H", "#H", "LH", "#H2", *dummy_data])
    coupled_reaction_bulk = CoupledReaction(["@H", "@H", "LH", "@H2", *dummy_data])
    coupled_reaction_bulk.set_partner(original_reaction)
    coupled_reaction_lhdes = CoupledReaction(["#H", "#H", "LHDES", "H2", *dummy_data])
    coupled_reaction_lhdes.set_partner(original_reaction)
    coupled_reaction_bulk_lhdes = CoupledReaction(
        ["@H", "@H", "LHDES", "H2", *dummy_data]
    )
    coupled_reaction_bulk_lhdes.set_partner(coupled_reaction_lhdes)

    coupled_reactions: list[CoupledReaction] = [
        coupled_reaction_bulk,
        coupled_reaction_lhdes,
        coupled_reaction_bulk_lhdes,
    ]

    for reaction in coupled_reactions:
        assert not isinstance(reaction.get_partner(), CoupledReaction)
        assert reaction.get_partner() == original_reaction
