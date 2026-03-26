"""Visualization utilities for UCLCHEM model outputs.

This module provides specialized plotting functions for visualizing
chemical abundances and reaction rates from UCLCHEM models.

**Key Functions:**

- :func:`plot_rate_summary` - Visualize top production/destruction reactions

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
