"""Visualization utilities for UCLCHEM model outputs.

This module provides specialized plotting functions for visualizing
chemical abundances and reaction rates from UCLCHEM models.

**Key Functions:**

- :func:`plot_rate_summary` - Visualize top production/destruction reactions
- :func:`plot_species` - Plot species abundances over time
- :func:`create_abundance_plot` - Create publication-ready abundance plots
- :func:`draw_panel_abundances` - Draw species abundance panel (Panel A)
- :func:`draw_panel_rates` - Draw production/destruction rate panel (Panel B)
- :func:`draw_panel_rate_constants` - Draw mean rate-constant bar chart (Panel C)
- :func:`plot_rates_deepdive` - Compose all three panels into a deepdive figure

**Example Usage:**
    >>> import uclchem
    >>>
    >>> model = uclchem.model.Cloud({})
    >>> model.check_error()
    Model ran successfully
    >>>
    >>> physics_df, chemistry_df, rate_constants_df = model.get_split_dataframes(
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

from . import style as style
from .compositions import create_abundance_plot, plot_rate_summary, plot_rates_deepdive
from .panels import (
    draw_panel_abundances,
    draw_panel_rate_constants,
    draw_panel_rates,
    plot_species,
)

__all__ = [
    "create_abundance_plot",
    "draw_panel_abundances",
    "draw_panel_rate_constants",
    "draw_panel_rates",
    "plot_rate_summary",
    "plot_rates_deepdive",
    "plot_species",
]
