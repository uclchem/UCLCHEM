"""Atomic plot units — each function draws a single panel onto a provided axis."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import pandas as pd

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

from uclchem.style import format_chemical_formula, format_reaction_label

from ._helpers import _color_for

#: CV threshold above which a rate constant is considered time-varying.
_K_VARY_THRESHOLD = 0.01


def plot_species(
    ax: plt.Axes,
    df: pd.DataFrame,
    species: list[str],
    legend: bool = True,
    plot_kwargs: dict[str, Any] | None = None,
) -> plt.Axes:
    """Plot the abundance of a list of species through time directly onto an axis.

    Parameters
    ----------
    ax : plt.Axes
        An axis object to plot on
    df : pd.DataFrame
        A dataframe created by
        ``uclchem.analysis.read_output_file``, ``uclchem.model.load_model`` or
        ``uclchem.model.Model.get_dataframes``.
    species : list[str]
        A list of species names to be plotted.
        If species name starts with "$" instead of "#" or "@",
        plots the sum of surface and bulk abundances
    legend : bool
        Whether to add a legend to the plot. Default = True.
    plot_kwargs : dict[str, Any] | None
        keyword arguments passed to ``ax.plot``.
        Default = None.

    Returns
    -------
    ax : plt.Axes
        Modified input axis is returned

    Raises
    ------
    KeyError
        if no ``"Time"`` column is present in ``df``.

    """
    if plot_kwargs is None:
        plot_kwargs = {}
    for species_name in species:
        linestyle = "solid"
        if species_name[0] == "$":
            abundances = df[species_name.replace("$", "#")]
            linestyle = "dashed"
            if species_name.replace("$", "@") in df.columns:
                abundances += df[species_name.replace("$", "@")]
        else:
            abundances = df[species_name]
        plot_kwargs["linestyle"] = linestyle
        plot_kwargs["label"] = species_name
        # Support legacy code that use either "age" or "Time" as the time variable
        if "age" in df.columns:
            timecolumn = "age"
        elif "Time" in df.columns:
            timecolumn = "Time"
        else:
            msg = "No time variable in dataframe"
            raise KeyError(msg)
        ax.plot(
            df[timecolumn],
            abundances,
            lw=2,
            **plot_kwargs,
        )
    ax.set(yscale="log")
    if legend:
        ax.legend()
    return ax


def draw_panel_abundances(
    ax: plt.Axes,
    time: pd.Series,
    species: str,
    chem: pd.DataFrame,
    companion: list[str] | None = None,
    *,
    reactant_species: set[str] | None = None,
    color_registry: dict[str, str] | None = None,
) -> plt.Axes:
    """Draw species abundances onto *ax* (Panel A of a deepdive figure).

    The target *species* is drawn in black at full weight.  Each entry in
    *companion* is drawn in a tab20 color; species that appear in
    *reactant_species* get a thicker, more opaque line to signal their
    direct chemical involvement.

    Parameters
    ----------
    ax : plt.Axes
        Axes to draw on.
    time : pd.Series
        Time series (years) for the x-axis.
    species : str
        UCLCHEM name of the primary species.
    chem : pd.DataFrame
        Chemistry (abundance) DataFrame, already filtered to the desired
        time range.
    companion : list[str] | None
        Additional species to overlay.  Pass ``None`` (default) to show
        only *species*.
    reactant_species : set[str] | None
        Species that appear as reactants in the top reactions; these are
        rendered with higher visual weight.  Default: empty set.
    color_registry : dict[str, str] | None
        Shared color map (species name → hex string).  Pass the same dict
        to multiple panel calls to keep colors consistent.  A fresh
        registry is created when ``None`` is passed.  Default: ``None``.

    Returns
    -------
    ax : plt.Axes
        The modified axes.

    """
    if companion is None:
        companion = []
    if reactant_species is None:
        reactant_species = set()
    if color_registry is None:
        color_registry = {}

    sp_color = {sp: _color_for(sp, color_registry) for sp in companion}
    sp_ls = {sp: ["-", "--"][i // 6] for i, sp in enumerate(companion)}

    if species in chem.columns:
        ax.plot(
            time,
            chem[species],
            label=format_chemical_formula(species),
            linewidth=2.0,
            color="black",
            alpha=0.9,
            zorder=10,
        )
    for sp in companion:
        in_rxns = sp in reactant_species
        ax.plot(
            time,
            chem[sp],
            label=format_chemical_formula(sp),
            linewidth=1.5 if in_rxns else 0.8,
            color=sp_color[sp],
            alpha=0.85 if in_rxns else 0.5,
            linestyle=sp_ls[sp],
        )
    ax.set_ylabel("Abundance (w.r.t. H)")
    ax.text(
        0.02,
        0.98,
        "A",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        weight="bold",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=7, loc="lower left", ncol=3, framealpha=0.9)
    ax.tick_params(labelbottom=False)
    return ax


def _draw_reaction_time_series(
    ax: plt.Axes,
    time: pd.Series,
    prod_df: pd.DataFrame,
    dest_df: pd.DataFrame,
    top_prod: list[str],
    top_dest: list[str],
    color_registry: dict[str, str],
    ylabel: str,
    panel_label: str,
) -> plt.Axes:
    """Shared time-series drawing logic for rates and rate-constants panels.

    Draws total formation/destruction envelopes (black), an "Other" aggregate
    for reactions not in *top_prod* / *top_dest* (gray), and individual lines
    for each listed reaction.

    Parameters
    ----------
    ax : plt.Axes
        Axes to draw on.
    time : pd.Series
        Time series (years) for the x-axis.
    prod_df : pd.DataFrame
        Per-reaction production data (rates or rate constants).
    dest_df : pd.DataFrame
        Per-reaction destruction data.
    top_prod : list[str]
        Reaction column names to draw individually for production.
    top_dest : list[str]
        Reaction column names to draw individually for destruction.
    color_registry : dict[str, str]
        Shared color map (reaction string → hex string).
    ylabel : str
        Label for the y-axis.
    panel_label : str
        Single-character panel identifier drawn in the top-left corner.

    Returns
    -------
    ax : plt.Axes
        The modified axes.

    """
    other_prod = [c for c in prod_df.columns if c not in top_prod]
    other_dest = [c for c in dest_df.columns if c not in top_dest]

    ax.plot(
        time,
        prod_df.sum(axis=1),
        lw=1.5,
        color="black",
        alpha=0.45,
        linestyle="-",
        zorder=10,
        label="Total formation",
    )
    ax.plot(
        time,
        dest_df.sum(axis=1),
        lw=1.5,
        color="black",
        alpha=0.45,
        linestyle="--",
        zorder=10,
        label="Total destruction",
    )
    if other_prod:
        ax.plot(
            time,
            prod_df[other_prod].sum(axis=1),
            lw=1.2,
            color="gray",
            alpha=0.6,
            linestyle=":",
            zorder=9,
            label="Other formation",
        )
    if other_dest:
        ax.plot(
            time,
            dest_df[other_dest].sum(axis=1),
            lw=1.2,
            color="gray",
            alpha=0.6,
            linestyle="-.",
            zorder=9,
            label="Other destruction",
        )
    for rxn in top_prod:
        ax.plot(
            time,
            prod_df[rxn],
            lw=1.2,
            color=_color_for(rxn, color_registry),
            alpha=0.85,
            linestyle="-",
            label=format_reaction_label(rxn),
        )
    for rxn in top_dest:
        ax.plot(
            time,
            dest_df[rxn],
            lw=1.2,
            color=_color_for(rxn, color_registry),
            alpha=0.85,
            linestyle="--",
            label=format_reaction_label(rxn),
        )

    leg = ax.legend(
        loc="lower center",
        mode="expand",
        ncol=2,
        fontsize=7,
        framealpha=0.9,
        handlelength=3.5,
    )
    for handle in leg.legend_handles:
        handle.set_linewidth(1.75)  # type: ignore[union-attr]
    ax.set_xlabel("Time / years")
    ax.set_ylabel(ylabel)
    ax.text(
        0.02,
        0.98,
        panel_label,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        weight="bold",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    return ax


def draw_panel_rates(
    ax: plt.Axes,
    time: pd.Series,
    prod_rates: pd.DataFrame,
    dest_rates: pd.DataFrame,
    top_prod: list[str] | None = None,
    top_dest: list[str] | None = None,
    *,
    top_k: int | None = 5,
    color_registry: dict[str, str] | None = None,
) -> plt.Axes:
    """Draw production and destruction rates onto *ax* (Panel B of a deepdive figure).

    Total formation and destruction envelopes are always drawn in black.
    Reactions listed in *top_prod* / *top_dest* are drawn individually in
    tab20 colors; remaining reactions are summed into a gray "Other" line.

    Parameters
    ----------
    ax : plt.Axes
        Axes to draw on.
    time : pd.Series
        Time series (years) for the x-axis.
    prod_rates : pd.DataFrame
        Per-reaction production rates (abundance wrt H s⁻¹), already
        filtered to the desired time range.
    dest_rates : pd.DataFrame
        Per-reaction destruction rates (absolute values), already filtered.
    top_prod : list[str] | None
        Reaction column names to draw individually for production.
        When ``None``, the top *top_k* reactions by mean rate are selected.
        Default: ``None``.
    top_dest : list[str] | None
        Reaction column names to draw individually for destruction.
        When ``None``, the top *top_k* reactions by mean rate are selected.
        Default: ``None``.
    top_k : int | None
        Number of top reactions to show individually when *top_prod* /
        *top_dest* are ``None``.  Pass ``None`` to show all reactions.
        Default: 5.
    color_registry : dict[str, str] | None
        Shared color map (reaction string → hex string).  Pass the same
        dict to multiple panel calls to keep colors consistent.  A fresh
        registry is created when ``None`` is passed.  Default: ``None``.

    Returns
    -------
    ax : plt.Axes
        The modified axes.

    """
    if color_registry is None:
        color_registry = {}
    if top_prod is None:
        cols = prod_rates.columns
        top_prod = (
            list(prod_rates[cols].mean().nlargest(top_k).index)
            if top_k is not None
            else list(cols)
        )
    if top_dest is None:
        cols = dest_rates.columns
        top_dest = (
            list(dest_rates[cols].mean().nlargest(top_k).index)
            if top_k is not None
            else list(cols)
        )

    return _draw_reaction_time_series(
        ax,
        time,
        prod_rates,
        dest_rates,
        top_prod,
        top_dest,
        color_registry,
        ylabel=r"Reaction rate (abundance wrt H s$^{-1}$)",
        panel_label="B",
    )


def draw_panel_rate_constants(
    ax: plt.Axes,
    time: pd.Series,
    prod_k: pd.DataFrame,
    dest_k: pd.DataFrame,
    top_prod: list[str] | None = None,
    top_dest: list[str] | None = None,
    *,
    top_k: int | None = 5,
    bar: bool = False,
    color_registry: dict[str, str] | None = None,
) -> plt.Axes:
    """Draw rate constants onto *ax* (Panel C of a deepdive figure).

    By default draws rate constants as time series so trends are visible.
    Pass ``bar=True`` for a mean bar chart; a warning is emitted for any
    reaction whose rate constant varies significantly over time (CV > 1 %).

    Parameters
    ----------
    ax : plt.Axes
        Axes to draw on.
    time : pd.Series
        Time series (years) for the x-axis.
    prod_k : pd.DataFrame
        Rate-constant DataFrame for production reactions, already filtered
        to the desired time range.
    dest_k : pd.DataFrame
        Rate-constant DataFrame for destruction reactions.
    top_prod : list[str] | None
        Reaction column names to include for production.
        When ``None``, the top *top_k* reactions by mean rate constant are
        selected.  Default: ``None``.
    top_dest : list[str] | None
        Reaction column names to include for destruction.
        When ``None``, the top *top_k* reactions by mean rate constant are
        selected.  Default: ``None``.
    top_k : int | None
        Number of top reactions to show when *top_prod* / *top_dest* are
        ``None``.  Pass ``None`` to show all.  Default: 5.
    bar : bool
        If ``True``, draw a mean bar chart instead of time series.
        Default: ``False``.
    color_registry : dict[str, str] | None
        Shared color map (reaction string → hex string).  A fresh registry
        is created when ``None`` is passed.  Default: ``None``.

    Returns
    -------
    ax : plt.Axes
        The modified axes.

    """
    if color_registry is None:
        color_registry = {}
    if top_prod is None:
        cols = prod_k.columns
        top_prod = (
            list(prod_k[cols].mean().nlargest(top_k).index)
            if top_k is not None
            else list(cols)
        )
    if top_dest is None:
        cols = dest_k.columns
        top_dest = (
            list(dest_k[cols].mean().nlargest(top_k).index)
            if top_k is not None
            else list(cols)
        )

    if not bar:
        return _draw_reaction_time_series(
            ax,
            time,
            prod_k,
            dest_k,
            top_prod,
            top_dest,
            color_registry,
            ylabel=r"Rate constant $k$ (s$^{-1}$)",
            panel_label="C",
        )

    # Bar mode: warn for time-varying rate constants, then plot means.
    import numpy as np  # ruff: ignore[import-outside-top-level]

    top_prod_k = [r for r in top_prod if r in prod_k.columns]
    top_dest_k = [r for r in top_dest if r in dest_k.columns]

    varying = [
        r
        for r in top_prod_k
        if prod_k[r].mean() > 0 and prod_k[r].std() / prod_k[r].mean() > _K_VARY_THRESHOLD
    ] + [
        r
        for r in top_dest_k
        if dest_k[r].mean() > 0 and dest_k[r].std() / dest_k[r].mean() > _K_VARY_THRESHOLD
    ]
    if varying:
        warnings.warn(
            f"{len(varying)} rate constant(s) vary over time (CV > {_K_VARY_THRESHOLD:.0%}); "
            "bar shows time-mean. Use bar=False for a time-resolved view.\n"
            f"Varying reactions: {', '.join(varying)}",
            stacklevel=2,
        )

    prod_k_mean = prod_k[top_prod_k].mean() if top_prod_k else pd.Series(dtype=float)
    dest_k_mean = dest_k[top_dest_k].mean() if top_dest_k else pd.Series(dtype=float)
    n_prod_k = len(top_prod_k)
    n_dest_k = len(top_dest_k)
    x_prod = np.arange(n_prod_k)
    x_dest = np.arange(n_prod_k, n_prod_k + n_dest_k)
    prod_colors = [_color_for(rxn, color_registry) for rxn in top_prod_k]
    dest_colors = [_color_for(rxn, color_registry) for rxn in top_dest_k]
    ax.bar(x_prod, prod_k_mean.values, width=0.75, color=prod_colors, alpha=0.85)  # type: ignore[arg-type]
    ax.bar(x_dest, dest_k_mean.values, width=0.75, color=dest_colors, alpha=0.85)  # type: ignore[arg-type]
    ax.set_xticks([n_prod_k / 2.0 - 0.5, n_prod_k + n_dest_k / 2.0 - 0.5])
    ax.set_xticklabels(["Formation", "Destruction"])
    ax.set_ylabel(r"Mean $k$ (s$^{-1}$)")
    ax.text(
        0.02,
        0.98,
        "C",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        weight="bold",
    )
    ax.set_yscale("log")
    return ax
