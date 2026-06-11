"""Run doctests for uclchem modules as part of the unit test suite."""

import doctest

import uclchem.makerates.io_functions
import uclchem.makerates.network
import uclchem.makerates.reaction
import uclchem.makerates.species
import uclchem.model
import uclchem.plot
import uclchem.utils

_MODULES = [
    uclchem.utils,
    uclchem.makerates.io_functions,
    uclchem.makerates.species,
    uclchem.makerates.reaction,
    uclchem.makerates.network,
    uclchem.model,
    uclchem.plot,
]

_FLAGS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE


def _run(module: object) -> None:
    results = doctest.testmod(module, verbose=False, optionflags=_FLAGS)  # type: ignore[arg-type]
    assert results.failed == 0, f"{results.failed} doctest(s) failed in {module}"  # type: ignore[union-attr]


def test_utils_doctests() -> None:
    _run(uclchem.utils)


def test_io_functions_doctests() -> None:
    _run(uclchem.makerates.io_functions)


def test_species_doctests() -> None:
    _run(uclchem.makerates.species)


def test_reaction_doctests() -> None:
    _run(uclchem.makerates.reaction)


def test_network_doctests() -> None:
    _run(uclchem.makerates.network)


def test_model_doctests() -> None:
    _run(uclchem.model)


def test_plot_doctests() -> None:
    _run(uclchem.plot)
