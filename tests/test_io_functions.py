import numpy as np
import pytest

from uclchem.makerates.io_functions import array_to_string, check_reaction
from uclchem.makerates.reaction import Reaction
from uclchem.makerates.species import Species, normalize_species_name


def test_array_to_string_1d_int():
    arr = np.array([1, 2, 3, 4])
    result = array_to_string("arr1", arr, type="int", parameter=True)
    assert "INTEGER, PARAMETER :: arr1 (4)=(/1,2,3,4/)" in result.replace("\n", "")


def test_array_to_string_1d_float():
    arr = np.array([1.0, 2.0, 3.0])
    result = array_to_string("arr2", arr, type="float", parameter=True)
    assert (
        "REAL(dp), PARAMETER :: arr2 (3)=(/1.0000e+00,2.0000e+00,3.0000e+00/)"
        in result.replace("\n", "")
    )


def test_array_to_string_2d_int():
    arr = np.array([[1, 2, 3], [4, 5, 6]])
    result = array_to_string("arr3", arr, type="int", parameter=True)
    expected = """INTEGER, PARAMETER :: arr3(2,3) = RESHAPE((/ 1,4,2,5,3,6 /), (/ 2, 3 /)&
    &)
"""
    assert result == expected


def test_array_to_string_2d_ones():
    arr = np.ones((5, 7), dtype=int)
    result = array_to_string("arr_ones", arr, type="int", parameter=True)
    expected = """INTEGER, PARAMETER :: arr_ones(5,7) = RESHAPE((/ 1,1,1,1,1,1,1,1,1,1,1,1&
    &,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /), (/ 5, 7 /))
"""
    assert result == expected


def test_array_to_string_2d_float():
    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    result = array_to_string("arr4", arr, type="float", parameter=True)
    expected = """REAL(dp), PARAMETER :: arr4(2,2) = RESHAPE((/ 1.0000e+00,3.0000e+00&
    &,2.0000e+00,4.0000e+00 /), (/ 2, 2 /))
"""
    assert result == expected


def test_array_to_string_2d_string():
    arr = np.array([["A", "B"], ["C", "D"]])
    result = array_to_string("arr5", arr, type="string", parameter=True)
    expected = """CHARACTER(Len=1), PARAMETER :: arr5(2,2) = RESHAPE((/ "A","C","B","D" /)&
    &, (/ 2, 2 /))
"""
    assert result == expected


# ---------------------------------------------------------------------------
# normalize_species_name
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw, expected", [
    # plain species – uppercased
    ("H2O",       "H2O"),
    ("ch3oh",     "CH3OH"),
    # chemical isomer prefix – prefix lowercased, formula uppercased
    ("o-H2",      "o-H2"),
    ("p-H2",      "p-H2"),
    ("a-CH3OH",   "a-CH3OH"),
    ("l-C3H2",    "l-C3H2"),
    # prefix case-normalisation (user typed uppercase prefix letter)
    ("O-H2",      "o-H2"),
    ("P-H2",      "p-H2"),
    # grain prefix preserved; chemical prefix lowercased
    ("#o-H2",     "#o-H2"),
    ("@a-CH3OH",  "@a-CH3OH"),
    ("#O-H2",     "#o-H2"),
    # negative ions (len==2 after stripping grain prefix -> NOT a chem prefix)
    ("C-",        "C-"),
    ("OH-",       "OH-"),
    ("E-",        "E-"),
    ("e-",        "E-"),
    # positive ions
    ("H+",        "H+"),
    ("C+",        "C+"),
    # reaction-type tokens and sentinels pass through unchanged
    ("NAN",       "NAN"),
    ("PHOTON",    "PHOTON"),
    ("FREEZE",    "FREEZE"),
    ("",          ""),
])
def test_normalize_species_name(raw, expected):
    assert normalize_species_name(raw) == expected


# ---------------------------------------------------------------------------
# Species with chemical prefix
# ---------------------------------------------------------------------------

def _make(name, mass=0, be=0.0):
    """Convenience: build a Species from minimal data."""
    return Species([name, mass, be, 0, 0, 0, 0])


def test_prefix_species_name_and_prefix_field():
    s = _make("o-H2", mass=2)
    assert s.name == "o-H2"
    assert s.prefix == "o"


def test_prefix_species_case_normalised():
    """Uppercase prefix letter in input must be lowercased."""
    s = _make("O-H2", mass=2)
    assert s.name == "o-H2"
    assert s.prefix == "o"


def test_prefix_grain_species_name_and_prefix_field():
    s = _make("#o-H2", mass=2, be=650)
    assert s.name == "#o-H2"
    assert s.prefix == "o"
    assert s.is_surface_species()


def test_no_prefix_species_has_empty_prefix():
    assert _make("H2O").prefix == ""
    assert _make("C-").prefix == ""
    assert _make("E-").prefix == ""


# ---------------------------------------------------------------------------
# find_constituents with chemical prefix
# ---------------------------------------------------------------------------

def test_find_constituents_ortho_h2():
    s = _make("o-H2", mass=2)
    counter = s.find_constituents(quiet=True)
    assert counter["H"] == 2
    assert s.mass == 2


def test_find_constituents_para_h2():
    s = _make("p-H2", mass=2)
    counter = s.find_constituents(quiet=True)
    assert counter["H"] == 2
    assert s.mass == 2


def test_find_constituents_anti_methanol():
    s = _make("a-CH3OH", mass=32)
    counter = s.find_constituents(quiet=True)
    assert counter["C"] == 1
    assert counter["H"] == 4
    assert counter["O"] == 1
    assert s.mass == 32


def test_find_constituents_surface_prefix_species():
    s = _make("#o-H2", mass=2, be=650)
    counter = s.find_constituents(quiet=True)
    assert counter["H"] == 2
    assert s.mass == 2


# ---------------------------------------------------------------------------
# get_charge / is_ion with chemical prefix
# ---------------------------------------------------------------------------

def test_prefix_species_not_ion():
    assert not _make("o-H2").is_ion()
    assert not _make("p-H2").is_ion()
    assert not _make("a-CH3OH").is_ion()


def test_prefix_species_charge_is_zero():
    assert _make("o-H2").get_charge() == 0
    assert _make("p-H2").get_charge() == 0


def test_negative_ion_charge_unchanged():
    assert _make("C-").get_charge() == -1
    assert _make("OH-").get_charge() == -1
    assert _make("E-").get_charge() == -1


def test_positive_ion_charge_unchanged():
    assert _make("H+").get_charge() == 1
    assert _make("C+").get_charge() == 1


# ---------------------------------------------------------------------------
# Reaction normalisation with chemical prefix
# ---------------------------------------------------------------------------

def _reac(r1, r2, r3, p1, p2, p3, p4):
    """Build a minimal Reaction (alpha=1, beta=0, gamma=0, T 0–1e9)."""
    return Reaction([r1, r2, r3, p1, p2, p3, p4, 1.0, 0.0, 0.0, 0.0, 1e9, 0.0])


def test_reaction_normalises_prefix_reactant():
    r = _reac("o-H2", "CRP", "NAN", "p-H2", "NAN", "NAN", "NAN")
    assert r.get_reactants()[0] == "o-H2"
    assert r.get_products()[0] == "p-H2"


def test_reaction_normalises_uppercase_prefix():
    """Uppercase prefix in a reaction file is lowercased to match species list."""
    r = _reac("O-H2", "CRP", "NAN", "P-H2", "NAN", "NAN", "NAN")
    assert r.get_reactants()[0] == "o-H2"
    assert r.get_products()[0] == "p-H2"


def test_reaction_element_conservation_ortho_para():
    """o-H2 -> p-H2 must pass element conservation (both are H2)."""
    # Should not raise
    _reac("o-H2", "CRP", "NAN", "p-H2", "NAN", "NAN", "NAN")


def test_reaction_charge_conservation_prefix_neutral():
    """Prefixed neutral species must not upset charge conservation."""
    # o-H2 (0) + H+ (1) -> H3+ (1) – charge conserved
    _reac("o-H2", "H+", "NAN", "H3+", "NAN", "NAN", "NAN")


# ---------------------------------------------------------------------------
# check_reaction with chemical prefix (case-insensitive keep_list matching)
# ---------------------------------------------------------------------------

def _keep_list(*names):
    from uclchem.makerates.reaction import reaction_types
    return ["", "NAN", "E-", "e-", "PHOTON", "CRP"] + list(reaction_types) + list(names)


def test_check_reaction_accepts_prefix_species():
    row = ["o-H2", "CRP", "NAN", "p-H2", "NAN", "NAN", "NAN", 1.0, 0.0, 0.0, 0.0, 1e9, 0.0]
    keep = _keep_list("o-H2", "p-H2")
    assert check_reaction(row, keep)


def test_check_reaction_accepts_uppercase_prefix_input():
    """Reaction file may use uppercase prefix; check_reaction normalises before lookup."""
    row = ["O-H2", "CRP", "NAN", "P-H2", "NAN", "NAN", "NAN", 1.0, 0.0, 0.0, 0.0, 1e9, 0.0]
    keep = _keep_list("o-H2", "p-H2")   # keep_list uses canonical form
    assert check_reaction(row, keep)


def test_check_reaction_drops_unknown_prefix_species():
    """A prefixed species not in the network should cause the reaction to be dropped."""
    row = ["x-H2", "CRP", "NAN", "H2", "NAN", "NAN", "NAN", 1.0, 0.0, 0.0, 0.0, 1e9, 0.0]
    keep = _keep_list("H2")   # x-H2 is not in the network
    assert not check_reaction(row, keep)
