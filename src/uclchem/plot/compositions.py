"""Composition functions — each creates a complete figure by assembling panels."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import matplotlib.figure
import matplotlib.pyplot as plt

if TYPE_CHECKING:
    import pandas as pd

from ._helpers import _filter_reactions, _reactants_from_rxn
from .panels import (
    draw_panel_abundances,
    draw_panel_rate_constants,
    draw_panel_rates,
    plot_species,
)

if TYPE_CHECKING:
    from uclchem.makerates.network import Network


def create_abundance_plot(
    df: pd.DataFrame,
    species: list[str],
    figsize: tuple[float, float] = (16, 9),
    plot_file: str | Path | None = None,
    plot_kwargs: dict[str, Any] | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Create a plot of the abundance of a list of species through time.

    Parameters
    ----------
    df : pd.DataFrame
        Pandas dataframe containing the UCLCHEM output, see
        ``uclchem.analysis.read_output_file``, ``uclchem.model.load_model`` or
        ``uclchem.model.Model.get_dataframes``.
    species : list[str]
        list of strings containing species names.
        Using a $ instead of # or @ will plot the sum of surface and bulk abundances.
    figsize : tuple[float, float]
        Size of figure, width by height in inches.
        Defaults to (16, 9).
    plot_file : str | Path | None
        Path to file where figure will be saved.
        If None, figure is not saved. Defaults to None.
    plot_kwargs : dict[str, Any] | None
        keyword arguments passed to ``ax.plot``.
        Default = None.

    Returns
    -------
    fig : plt.Figure
        created Figure object
    ax : plt.Axes
        created axis object

    """
    if plot_kwargs is None:
        plot_kwargs = {}

    fig, ax = plt.subplots(figsize=figsize, tight_layout=True)

    ax = plot_species(ax, df, species, legend=False, **plot_kwargs)
    ax.legend(loc=4, fontsize="small")

    ax.set_xlabel("Time / years")
    ax.set_ylabel("X$_{Species}$")

    ax.set_yscale("log")
    if plot_file is not None:
        fig.savefig(plot_file)
    return fig, ax


def plot_rate_summary(
    production_df: pd.DataFrame,
    destruction_df: pd.DataFrame,
    step: int,
    xlabel: str = "Reaction rate (abundance wrt H / s)",
    top_k_rates: int = 5,
) -> list[plt.Axes]:
    """Create a summary of the top k production and destruction reactions.

    Parameters
    ----------
    production_df : pd.DataFrame
        dataframe with reaction rates of
        formation reactions of species of interest
    destruction_df : pd.DataFrame
        dataframe with reaction rates of
        destruction reactions of species of interest
    step : int
        time index of dataframes to plot.
    xlabel : str
        xlabel. Default: "Reaction rate (abundance wrt H / s)"
    top_k_rates : int
        Plot top k formation and destruction reactions.
        Default: 5

    Returns
    -------
    axs : list[plt.Axes]
        axes of the plot

    """
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(7, top_k_rates))
    production_df.iloc[step].sort_values(ascending=False)[:top_k_rates].plot.barh(
        ax=axs[0], title="Production", logx=True
    )
    destruction_df.iloc[step].sort_values(ascending=False)[:top_k_rates].plot.barh(
        ax=axs[1], title="Destruction"
    )
    axs[1].set_xlabel(xlabel)
    return axs


def plot_rates_deepdive(
    species: str,
    physics_df: pd.DataFrame,
    chemistry_df: pd.DataFrame,
    rate_constants_df: pd.DataFrame,
    network: Network | None = None,
    *,
    filter_threshold: float = 0.01,
    filter_window: tuple[float, float] = (1e4, 1e6),
    filter_freeze: bool = True,
    max_species_show: int = 12,
    figsize: tuple[float, float] = (8, 12),
    output_path: Path | str | None = None,
    fig: plt.Figure | None = None,
    color_registry: dict[str, str] | None = None,
) -> tuple[matplotlib.figure.FigureBase, plt.Axes, plt.Axes, plt.Axes]:
    """Create a three-panel chemical deep-dive figure for *species*.

    **Panel A** (top): Abundances of *species* and the reactant species
    involved in its top production and destruction reactions.

    **Panel B** (bottom): Individual production (solid) and destruction
    (dashed) reaction rates, plus totals.

    **Panel C** (middle): Bar chart of mean rate constants for the top
    reactions, colored to match Panel B.

    Parameters
    ----------
    species : str
        UCLCHEM species name to analyze, e.g. ``"HCO+"``.
    physics_df : pd.DataFrame
        Physics DataFrame from
        :meth:`~uclchem.model.AbstractModel.get_dataframes`.
    chemistry_df : pd.DataFrame
        Chemistry (abundance) DataFrame.
    rate_constants_df : pd.DataFrame
        Rate-constants DataFrame (``with_rate_constants=True``).
    network : Network | None
        Pre-loaded :class:`~uclchem.makerates.network.Network`.  If
        ``None`` the default network is loaded via
        :meth:`~uclchem.makerates.network.Network.from_csv`.
    filter_threshold : float
        Reactions whose rate never exceeds this fraction of the per-step
        maximum within *filter_window* are excluded.  Default: ``0.01``.
    filter_window : tuple[float, float]
        ``(t_min, t_max)`` in years used for reaction filtering and
        ranking.  Default: ``(1e4, 1e6)``.
    filter_freeze : bool
        If ``True`` (default), exclude freeze-out reactions.
    max_species_show : int
        Maximum number of companion species to draw in Panel A.
        Default: ``12``.
    figsize : tuple[float, float]
        Figure width × height in inches.  Ignored when *fig* is provided.
        Default: ``(8, 12)``.
    output_path : Path | str | None
        If provided, save the figure as both ``<output_path>.pdf`` and
        ``<output_path>.png``.  Only meaningful when *fig* is a top-level
        :class:`~matplotlib.figure.Figure`.  Default: ``None``.
    fig : matplotlib.figure.FigureBase | None
        Existing figure or sub-figure to draw into.  Pass a
        :class:`~matplotlib.figure.SubFigure` obtained from
        ``parent.subfigures()`` to embed this plot inside a larger layout.
        If ``None`` (default) a new figure is created.
    color_registry : dict[str, str] | None
        Mutable mapping from species / reaction name to hex color string.
        Pass the same dict to multiple calls to keep colors consistent
        across subfigures.  If ``None`` (default) a fresh registry is
        created internally.

    Returns
    -------
    fig : matplotlib.figure.FigureBase
        The figure (or sub-figure) containing all three panels.
    ax_abundances : plt.Axes
        Panel A — species abundances.
    ax_rates : plt.Axes
        Panel B — production / destruction rates.
    ax_rate_constants : plt.Axes
        Panel C — mean rate-constant bar chart.

    Notes
    -----
    .. todo::
        Return a ``pd.DataFrame`` of per-reaction rates at the final
        timestep so callers can export to CSV / Excel.

    """
    from uclchem import analysis  # noqa: PLC0415
    from uclchem.makerates.network import Network as _Network  # noqa: PLC0415

    if network is None:
        network = _Network.from_csv()

    _reg: dict[str, str] = {} if color_registry is None else color_registry

    # Compute rates
    _, rates = analysis.rate_constants_to_dy_and_rates(
        physics_df, chemistry_df, rate_constants_df, network=network
    )
    prod_rates_full, dest_rates_full = analysis.get_production_and_destruction(
        species, rates
    )
    prod_k_full, dest_k_full = analysis.get_production_and_destruction(
        species, rate_constants_df
    )

    # Restrict to positive time steps
    pos_mask = physics_df["Time"] > 0
    time = physics_df["Time"][pos_mask]

    prod_rates = prod_rates_full[pos_mask]
    dest_rates = dest_rates_full[pos_mask].abs()
    prod_k = prod_k_full[pos_mask]
    dest_k = dest_k_full[pos_mask].abs()

    if filter_freeze:
        prod_rates = prod_rates[
            [c for c in prod_rates.columns if "FREEZE" not in c.upper()]
        ]
        dest_rates = dest_rates[
            [c for c in dest_rates.columns if "FREEZE" not in c.upper()]
        ]
        prod_k = prod_k[[c for c in prod_k.columns if "FREEZE" not in c.upper()]]
        dest_k = dest_k[[c for c in dest_k.columns if "FREEZE" not in c.upper()]]

    # Select top reactions within the filter window
    win_mask = (time >= filter_window[0]) & (time <= filter_window[1])
    top_prod = _filter_reactions(prod_rates[win_mask], filter_threshold)
    top_dest = _filter_reactions(dest_rates[win_mask], filter_threshold)
    top_prod = (
        prod_rates[win_mask][top_prod].max().sort_values(ascending=False).index.tolist()
    )
    top_dest = (
        dest_rates[win_mask][top_dest].max().sort_values(ascending=False).index.tolist()
    )

    # Geometric swap reactions are correction terms for ice layer growth/shrinkage;
    # always include them in Panel B regardless of filter_threshold.
    _swap_types = ("SURFSWAP_GEOMETRIC", "BULKSWAP_GEOMETRIC")
    for col in prod_rates.columns:
        if any(s in col for s in _swap_types) and col not in top_prod:
            top_prod.append(col)
    for col in dest_rates.columns:
        if any(s in col for s in _swap_types) and col not in top_dest:
            top_dest.append(col)

    # Companion species for Panel A
    reactant_species: set[str] = set()
    for rxn in top_prod + top_dest:
        reactant_species.update(_reactants_from_rxn(rxn))

    companion = sorted(reactant_species & set(chemistry_df.columns))
    companion = [
        sp for sp in companion if sp != species and (chemistry_df[sp][pos_mask] > 0).any()
    ]

    # For ice species, always show the surface/bulk partner so the two ice
    # reservoirs can be compared directly, regardless of filter_threshold.
    _ice_partner = {"#": "@", "@": "#"}
    if species and species[0] in _ice_partner:
        ice_partner = _ice_partner[species[0]] + species[1:]
        if (
            ice_partner in chemistry_df.columns
            and (chemistry_df[ice_partner][pos_mask] > 0).any()
            and ice_partner not in companion
        ):
            companion.insert(0, ice_partner)

    companion = companion[:max_species_show]

    # Figure layout: abundances / rate constants / rates (height ratio 3:1:3)
    if fig is None:
        fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 1, height_ratios=[3, 1, 3])
    ax_a = fig.add_subplot(gs[0])
    ax_c = fig.add_subplot(gs[1])
    ax_b = fig.add_subplot(gs[2])

    chem_filtered = chemistry_df[pos_mask]
    draw_panel_abundances(
        ax_a,
        time,
        species,
        chem_filtered,
        companion,
        reactant_species=reactant_species,
        color_registry=_reg,
    )
    draw_panel_rates(
        ax_b, time, prod_rates, dest_rates, top_prod, top_dest, color_registry=_reg
    )
    draw_panel_rate_constants(
        ax_c, time, prod_k, dest_k, top_prod, top_dest, bar=True, color_registry=_reg
    )

    if output_path is not None:
        base = Path(output_path)
        fig.savefig(base.with_suffix(".pdf"), format="pdf")
        fig.savefig(base.with_suffix(".png"), format="png")

    return fig, ax_a, ax_b, ax_c
