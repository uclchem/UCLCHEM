"""Visualization utilities for UCLCHEM model outputs.

This module provides specialized plotting functions for visualizing
chemical abundances and reaction rates from UCLCHEM models.

**Key Functions:**

- :func:`plot_rate_summary` - Visualize top production/destruction reactions
- :func:`plot_species` - Plot species abundances over time
- :func:`create_abundance_plot` - Create publication-ready abundance plots

**Example Usage:**
    >>> import uclchem
    >>>
    >>> model = uclchem.model.Cloud({})
    >>> model.check_error()
    Model ran successfully
    >>>
    >>> physics_df, chemistry_df, rate_constants_df = model.get_dataframes(
    ...     with_rate_constants=True,
    ... )
    >>> # Making a plot of the abundances over time
    >>> fig, ax = uclchem.plot.create_abundance_plot(
    ...     model.get_joined_dataframes(), # need both "Time" and abundance columns in one dataframe
    ...     ["H", "$H", "H2O", "$H2O", "CH3OH", "$CH3OH"],
    ... )
    >>>
    >>> # Making a plot of the main formation and destruction reactions
    >>> # at a specific timepoint
    >>> network = uclchem.makerates.network.Network.from_csv()
    >>> dy, reaction_rates = uclchem.analysis.rate_constants_to_dy_and_rates(
    ...     physics_df,
    ...     chemistry_df,
    ...     rate_constants_df,
    ...     network=network,
    ... )
    >>> production_df, destruction_df = uclchem.analysis.get_production_and_destruction(
    ...     "H2O",
    ...     reaction_rates,
    ... )
    >>>
    >>> # Plot top 5 reactions at a specific timestep
    >>> uclchem.plot.plot_rate_summary(
    ...      production_df,
    ...      destruction_df,
    ...      step=50,
    ...      top_k_rates=5
    ...  ) # doctest: +SKIP

**Note:**

Most plotting functionality is available through the model objects themselves
via methods like :meth:`~uclchem.model.AbstractModel.create_abundance_plot`.

**See Also:**

- :mod:`uclchem.analysis` - Analysis tools that include plotting functions
- :mod:`uclchem.model` - Model classes with built-in plotting methods
"""

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd


def plot_rate_summary(
    production_df: pd.DataFrame,
    destruction_df: pd.DataFrame,
    step: int,
    xlabel: str = "Reaction rate (abundance wrt H / s)",
    top_k_rates: int = 5,
) -> list[plt.Axes]:
    """Create a summary of the top k production and destruction reactions.

    Args:
        production_df (pd.DataFrame): dataframe with reaction rates of
            formation reactions of species of interest
        destruction_df (pd.DataFrame): dataframe with reaction rates of
            destruction reactions of species of interest
        step (int): time index of dataframes to plot.
        xlabel (str): xlabel. Default: "Reaction rate (abundance wrt H / s)"
        top_k_rates (int): Plot top k formation and destruction reactions.
            Default: 5

    Returns:
        axs (list[plt.Axes]): axes of the plot

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


def create_abundance_plot(
    df: pd.DataFrame,
    species: list[str],
    figsize: tuple[int | float, int | float] = (16, 9),
    plot_file: str | Path | None = None,
    plot_kwargs: dict[str, Any] | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Create a plot of the abundance of a list of species through time.

    Args:
        df (pd.DataFrame): Pandas dataframe containing the UCLCHEM output, see
            ``uclchem.analysis.read_output_file``, ``uclchem.model.load_model`` or
            ``uclchem.model.Model.get_dataframes``.
        species (list[str]): list of strings containing species names.
            Using a $ instead of # or @ will plot the sum of surface and bulk abundances.
        figsize (tuple[int | float]): Size of figure, width by height in inches.
            Defaults to (16, 9).
        plot_file (str | Path | None): Path to file where figure will be saved.
            If None, figure is not saved. Defaults to None.
        plot_kwargs (dict[str, Any] | None): keyword arguments passed to ``ax.plot``.

    Returns:
        fig (plt.Figure): created Figure object
        ax (plt.Axes): created axis object

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


def plot_species(
    ax: plt.Axes,
    df: pd.DataFrame,
    species: list[str],
    legend: bool = True,
    plot_kwargs: dict[str, Any] | None = None,
) -> plt.Axes:
    """Plot the abundance of a list of species through time directly onto an axis.

    Args:
        ax (plt.Axes): An axis object to plot on
        df (pd.DataFrame): A dataframe created by
            ``uclchem.analysis.read_output_file``, ``uclchem.model.load_model`` or
            ``uclchem.model.Model.get_dataframes``.
        species (list[str]): A list of species names to be plotted.
            If species name starts with "$" instead of "#" or "@",
            plots the sum of surface and bulk abundances
        legend (bool): Whether to add a legend to the plot. Default = True.
        plot_kwargs (dict[str, Any] | None): keyword arguments passed to ``ax.plot``.

    Returns:
        ax (plt.Axes): Modified input axis is returned

    Raises:
        KeyError: if no ``"Time"`` column is present in ``df``.

    """
    if plot_kwargs is None:
        plot_kwargs = {}
    for species_name in species:
        linestyle = "solid"
        if species_name[0] == "$":
            abundances = df[species_name.replace("$", "#")]
            linestyle = "dashed"
            if species_name.replace("$", "@") in df.columns:
                abundances = abundances + df[species_name.replace("$", "@")]
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
