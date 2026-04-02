"""Visualization utilities for UCLCHEM model outputs.

This module provides specialized plotting functions for visualizing
chemical abundances and reaction rates from UCLCHEM models.

**Key Functions:**

- :func:`plot_rate_summary` - Visualize top production/destruction reactions
- :func:`plot_species` - Plot species abundances over time
- :func:`create_abundance_plot` - Create publication-ready abundance plots

**Example Usage:**

.. code-block:: python

    from uclchem import plot

    # Plot top 5 reactions at a specific timestep
    plot.plot_rate_summary(
        production_df,
        destruction_df,
        step=50,
        top_k_rates=5
    )

**Note:**

Most plotting functionality is available through the model objects themselves
via methods like :meth:`~uclchem.model.AbstractModel.create_abundance_plot`.

**See Also:**

- :mod:`uclchem.analysis` - Analysis tools that include plotting functions
- :mod:`uclchem.model` - Model classes with built-in plotting methods
"""

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
    figsize: tuple[int | float] = (16, 9),
    plot_file: str | None = None,
    **plot_kwargs: dict[str, Any],
) -> tuple[plt.Figure, plt.Axes]:
    """Create a plot of the abundance of a list of species through time.

    Args:
        df (pd.DataFrame): Pandas dataframe containing the UCLCHEM output, see `read_output_file`
        species (list[str]): list of strings containing species names.
            Using a $ instead of # or @ will plot the sum of surface and bulk abundances.
        figsize (tuple[int | float]): Size of figure, width by height in inches.
            Defaults to (16, 9).
        plot_file (str | None): Path to file where figure will be saved.
            If None, figure is not saved. Defaults to None.
        plot_kwargs (dict[str, Any]): keyword arguments passed to `ax.plot`.

    Returns:
        fig (plt.Figure): created Figure object
        ax (plt.Axes): created axis object

    """
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
    **plot_kwargs: dict[str, Any],
) -> plt.Axes:
    """Plot the abundance of a list of species through time directly onto an axis.

    Args:
        ax (plt.Axes): An axis object to plot on
        df (pd.DataFrame): A dataframe created by `read_output_file`
        species (list[str]): A list of species names to be plotted.
            If species name starts with "$" instead of "#" or "@",
            plots the sum of surface and bulk abundances
        legend (bool): Whether to add a legend to the plot. Default = True.
        plot_kwargs (dict[str, Any]): keyword arguments passed to `ax.plot`.

    Returns:
        ax (plt.Axes): Modified input axis is returned

    Raises:
        KeyError: if no "Time" column is present in `df`.

    """
    for specIndx, specName in enumerate(species):
        linestyle = "solid"
        if specName[0] == "$":
            abundances = df[specName.replace("$", "#")]
            linestyle = "dashed"
            if specName.replace("$", "@") in df.columns:
                abundances = abundances + df[specName.replace("$", "@")]
        else:
            abundances = df[specName]
        plot_kwargs["linestyle"] = linestyle
        plot_kwargs["label"] = specName
        # Support legacy code that use either "age" or "Time" as the time variable
        if "age" in df.columns:
            timecolumn = "age"
        elif "Time" in df.columns:
            timecolumn = "Time"
        else:
            raise KeyError("No time variable in dataframe")
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
