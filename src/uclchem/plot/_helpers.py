"""Private utilities shared across the plot subpackage."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib as mpl

if TYPE_CHECKING:
    import pandas as pd

_TAB20_COLORS: list[str] = [mpl.colors.to_hex(c) for c in mpl.colormaps["tab20"].colors]  # type: ignore[attr-defined]


def _color_for(name: str, registry: dict[str, str]) -> str:
    if name not in registry:
        registry[name] = _TAB20_COLORS[len(registry) % len(_TAB20_COLORS)]
    return registry[name]


def _reactants_from_rxn(rxn: str) -> list[str]:
    """Return the non-special reactant tokens from a UCLCHEM reaction string.

    Parameters
    ----------
    rxn : str
        Reaction string in UCLCHEM format, e.g. ``"H3+ + CO -> HCO+ + H2"``.

    Returns
    -------
    list[str]
        Species tokens on the left-hand side, excluding radiation/process
        tokens (``PHOTON``, ``CRP``, ``CRPHOT``, ``UV``, ``FREEZE``).

    """
    skip = {"PHOTON", "CRP", "CRPHOT", "UV", "FREEZE"}
    if "->" not in rxn:
        return []
    lhs = rxn.split("->", maxsplit=1)[0]
    return [p.strip() for p in lhs.split("+") if p.strip() and p.strip() not in skip]


def _filter_reactions(
    rates_window: pd.DataFrame,
    threshold: float,
) -> list[str]:
    """Return reaction columns that exceed a relative rate threshold.

    Parameters
    ----------
    rates_window : pd.DataFrame
        Rate values restricted to the filter time window; columns are
        reaction strings.
    threshold : float
        Minimum fraction of the per-timestep maximum rate a reaction must
        reach at least once to be retained.

    Returns
    -------
    list[str]
        Column names (reaction strings) that pass the threshold criterion.

    """
    keep = []
    for rxn in rates_window.columns:
        max_per_step = rates_window.max(axis=1)
        if (rates_window[rxn] > threshold * max_per_step).any():
            keep.append(rxn)
    return keep
