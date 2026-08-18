"""Publication-quality matplotlib style for UCLCHEM.

Applies STIX serif fonts and MNRAS-compatible sizes without requiring a
LaTeX installation.  The style is applied automatically when ``uclchem`` is
imported.  Set the environment variable ``UCLCHEM_NO_STYLE=1`` before
importing to suppress auto-apply.

**Key Functions:**

- :func:`apply` - Apply the UCLCHEM rcParams globally
- :func:`reset` - Restore matplotlib defaults
- :func:`context` - Temporary style context manager
- :func:`format_chemical_formula` - Convert species names to mathtext labels
- :func:`format_reaction_label` - Convert reaction strings to mathtext labels

**Example Usage:**

    >>> import uclchem  # style applied on import
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> _ = ax.set_xlabel(uclchem.style.format_chemical_formula("HCO+"))
    >>> uclchem.style.reset()      # restore matplotlib defaults
    >>> uclchem.style.apply()      # re-apply UCLCHEM style

"""

from __future__ import annotations

import os
import re
from contextlib import contextmanager
from typing import TYPE_CHECKING

import matplotlib as mpl
import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from collections.abc import Generator

UCLCHEM_STYLE: dict[str, object] = {
    # fonts
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 9,
    # axes labels / ticks
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    # ticks: inward, mirrored
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    # lines
    "lines.linewidth": 1.5,
    "axes.linewidth": 0.8,
    # grid
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    # output
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    # layout
    "figure.constrained_layout.use": True,
}

_SPECIAL_TOKENS: dict[str, str] = {
    "PHOTON": r"\gamma",
    "CRP": r"\mathrm{CRP}",
    "CRPHOT": r"\mathrm{CR\gamma}",
    "UV": r"\mathrm{UV}",
    "FREEZE": r"\mathrm{FREEZE}",
}


def apply(style: dict[str, object] = UCLCHEM_STYLE) -> None:
    """Apply *style* to the global matplotlib rcParams.

    Parameters
    ----------
    style : dict[str, object]
        Mapping of rcParam keys to values.  Defaults to
        :data:`UCLCHEM_STYLE`.

    """
    mpl.rcParams.update(style)


def reset() -> None:
    """Restore matplotlib's built-in default rcParams."""
    mpl.rcdefaults()


@contextmanager
def context(
    style: dict[str, object] = UCLCHEM_STYLE,
) -> Generator[None, None, None]:
    """Apply *style* temporarily as a context manager.

    Parameters
    ----------
    style : dict[str, object]
        Mapping of rcParam keys to values.  Defaults to
        :data:`UCLCHEM_STYLE`.

    Yields
    ------
    None
        All ``plt`` calls inside the block use *style*; settings are
        restored on exit.

    Examples
    --------
    >>> import uclchem
    >>> with uclchem.style.context():
    ...     pass  # plots here use UCLCHEM style

    """
    with plt.style.context(style):
        yield


def format_chemical_formula(name: str) -> str:
    r"""Convert a UCLCHEM species name to a mathtext label.

    Subscripts digit runs and superscripts charge suffixes so the result
    renders correctly with matplotlib's built-in mathtext engine (no LaTeX
    required).

    Parameters
    ----------
    name : str
        UCLCHEM species name, e.g. ``"HCO+"``, ``"#H2O"``, ``"@CO2"``,
        ``"E-"``, ``"o-H2"``.

    Returns
    -------
    str
        Mathtext string, e.g. ``r"$\mathrm{HCO^{+}}$"``.

    Examples
    --------
    >>> format_chemical_formula("HCO+")
    '$\\mathrm{HCO^{+}}$'
    >>> format_chemical_formula("H2O")
    '$\\mathrm{H_{2}O}$'
    >>> format_chemical_formula("#H2O")
    '#$\\mathrm{H_{2}O}$'

    >>> # The '$' indicating 'all ice' is escaped
    >>> format_chemical_formula("$H2O")
    '\\$$\\mathrm{H_{2}O}$'

    """
    s = name.strip()

    phase_prefix = ""
    if s.startswith(("#", "@", "$")):
        phase_prefix = s[0]
        if phase_prefix == "$":
            phase_prefix = "\\" + phase_prefix
        s = s[1:]

    iso_prefix = ""
    iso_match = re.match(r"^([op])-", s)
    if iso_match:
        iso_prefix = iso_match.group(0)
        s = s[len(iso_prefix) :]

    charge = ""
    charge_match = re.search(r"([+\-]+)$", s)
    if charge_match:
        charge = f"^{{{charge_match.group(1)}}}"
        s = s[: charge_match.start()]

    s = re.sub(r"(\d+)", r"_{\1}", s)

    inner = f"{iso_prefix}{s}{charge}"
    return f"{phase_prefix}$\\mathrm{{{inner}}}$"


def format_reaction_label(rxn: str, k_mean: float | None = None) -> str:
    r"""Convert a UCLCHEM reaction string to a mathtext legend label.

    Parameters
    ----------
    rxn : str
        Reaction string in UCLCHEM format, e.g.
        ``"H3+ + CO -> HCO+ + H2"``.
    k_mean : float | None
        Mean rate coefficient to append to the label.  If ``None`` no
        rate is shown.  Default: ``None``.

    Returns
    -------
    str
        Mathtext string suitable for use as a matplotlib legend label.

    Notes
    -----
    Whitespace is required around each "+" that indicates a different molecule,
    because otherwise it cannot distinguish between two separate species or
    the charge of the first molecule.

    Examples
    --------
    >>> format_reaction_label("H2 + O -> OH + H")
    '$\\mathrm{H_{2}}$ + $\\mathrm{O}$ $\\rightarrow$ $\\mathrm{OH}$ + $\\mathrm{H}$'
    >>> format_reaction_label("H3+ + CO -> H2 + HCO+")
    '$\\mathrm{H_{3}^{+}}$ + $\\mathrm{CO}$ $\\rightarrow$ $\\mathrm{H_{2}}$ + $\\mathrm{HCO^{+}}$'

    """
    rxn = re.sub(r"__\d+$", "", rxn)

    if "->" not in rxn:
        return rxn

    lhs, rhs = rxn.split("->", 1)

    def _fmt_side(side: str) -> str:
        parts = [p.strip() for p in side.split(" +")]
        formatted = []
        for p in parts:
            if not p:
                continue
            if p in _SPECIAL_TOKENS:
                formatted.append(f"${_SPECIAL_TOKENS[p]}$")
            else:
                formatted.append(format_chemical_formula(p))
        return " + ".join(formatted)

    label = f"{_fmt_side(lhs)} $\\rightarrow$ {_fmt_side(rhs)}"
    if k_mean is not None:
        label += rf" ($k$ = {k_mean:.2e})"
    return label


if not os.environ.get("UCLCHEM_NO_STYLE"):
    apply()
