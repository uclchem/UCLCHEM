"""UCLCHEM Plotting Module

Visualization utilities for UCLCHEM model outputs.

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


def plot_rate_summary(prod, dest, step, xlabel="rate", top_k_rates=5):
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, top_k_rates))
    prod.iloc[step].sort_values(ascending=False)[:top_k_rates].plot.barh(
        ax=ax[0], title="production", logx=True
    )
    dest.iloc[step].sort_values(ascending=False)[:top_k_rates].plot.barh(
        ax=ax[1], title="destruction"
    )
    ax[1].set_xlabel(xlabel)
