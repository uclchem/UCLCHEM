"""Example plotting script for UCLCHEM.

Runs the default static Cloud model and demonstrates every plotting routine
in ``uclchem.plot``, explaining the purpose of each one.

Run from the repository root::

    python scripts/example_plotting_script.py

Output figures are written to ``examples/plots/``.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt

import uclchem
import uclchem.plot

if __name__ == "__main__":
    # ---------------------------------------------------------------------------
    # Setup: run the model and collect dataframes
    # ---------------------------------------------------------------------------

    out_dir = Path("examples/plots")
    out_dir.mkdir(parents=True, exist_ok=True)

    model = uclchem.model.Cloud({})
    model.check_error()

    # Joined dataframe is convenient for abundance-vs-time plots.
    joined_df = model.get_joined_dataframes()

    # Separate dataframes are required for reaction-rate analysis.
    physics_df, chemistry_df, rate_constants_df = model.get_dataframes(
        with_rate_constants=True, joined=False
    )

    # A loaded network is needed to convert rate constants to actual rates.
    network = uclchem.makerates.network.Network.from_csv()

    SPECIES_OF_INTEREST = "H2O"
    COMPANION_SPECIES = ["H2O", "#H2O", "$H2O", "CO", "#CO", "CH3OH", "#CH3OH"]
    TIMESTEP = 50  # index used for single-step snapshots

    # ---------------------------------------------------------------------------
    # 1. create_abundance_plot
    #    Purpose: the fastest way to get an overview of species abundances through
    #    time.  Handles the figure/axis boilerplate for you and returns both so
    #    you can keep tweaking if needed.
    # ---------------------------------------------------------------------------

    fig, ax = uclchem.plot.create_abundance_plot(
        joined_df,
        COMPANION_SPECIES,
        figsize=(10, 5),
        plot_file=out_dir / "01_abundance_overview.png",
    )
    ax.set_title("1 · create_abundance_plot — quick abundance overview")
    fig.savefig(out_dir / "01_abundance_overview.png")
    plt.close(fig)

    # ---------------------------------------------------------------------------
    # 2. plot_species (panel-level)
    #    Purpose: draws species lines onto an axis *you* own, so you can compose
    #    it freely — here we overlay the gas temperature on a twin y-axis to
    #    correlate chemistry with the physical evolution.
    # ---------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(10, 5), tight_layout=True)
    ax = uclchem.plot.plot_species(ax, joined_df, COMPANION_SPECIES)
    ax.set_xscale("log")
    ax.set_xlabel("Time / years")
    ax.set_ylabel("Abundance (w.r.t. H)")
    ax.set_title("2 · plot_species — abundances + physical parameter overlay")

    ax_twin = ax.twinx()
    ax_twin.plot(joined_df["Time"], joined_df["gasTemp"], color="black", lw=1.2, alpha=0.5, linestyle=":")
    ax_twin.set_ylabel("Gas temperature / K")

    fig.savefig(out_dir / "02_species_with_temperature.png")
    plt.close(fig)

    # ---------------------------------------------------------------------------
    # 3. plot_rate_summary
    #    Purpose: instant snapshot of the top-k formation and destruction
    #    reactions at a single timestep.  Good for a quick sanity check before
    #    committing to a full deepdive.
    # ---------------------------------------------------------------------------

    _, rates = uclchem.analysis.rate_constants_to_dy_and_rates(
        physics_df, chemistry_df, rate_constants_df, network=network
    )
    production_df, destruction_df = uclchem.analysis.get_production_and_destruction(
        SPECIES_OF_INTEREST, rates
    )

    axs = uclchem.plot.plot_rate_summary(
        production_df, destruction_df, step=TIMESTEP, top_k_rates=6
    )
    axs[0].set_title(
        f"3 · plot_rate_summary — top reactions for {SPECIES_OF_INTEREST} at step {TIMESTEP}"
    )
    axs[0].figure.tight_layout()
    axs[0].figure.savefig(out_dir / "03_rate_summary.png")
    plt.close(axs[0].figure)

    # ---------------------------------------------------------------------------
    # 4. plot_rates_deepdive
    #    Purpose: the full three-panel diagnostic figure.
    #      Panel A — abundances of the target species and its key reactant partners
    #      Panel B — individual production (solid) and destruction (dashed) rates
    #                through time, plus totals and an "Other" aggregate
    #      Panel C — mean rate constants k for the top reactions, color-matched
    #                to Panel B so you can immediately see which reactions dominate
    #    Use this when you want to understand *why* an abundance evolves the way
    #    it does over the full time range.
    # ---------------------------------------------------------------------------

    fig, ax_a, ax_b, ax_c = uclchem.plot.plot_rates_deepdive(
        SPECIES_OF_INTEREST,
        physics_df,
        chemistry_df,
        rate_constants_df,
        network=network,
        filter_window=(1e3, 1e7),
        figsize=(9, 13),
    )
    fig.suptitle(
        f"4 · plot_rates_deepdive — full chemical diagnostic for {SPECIES_OF_INTEREST}",
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(out_dir / "04_deepdive.png", bbox_inches="tight")
    plt.close(fig)

    # ---------------------------------------------------------------------------
    # 5. draw_panel_* (individual panels)
    #    Purpose: side-by-side deepdive comparison of two species in a single
    #    figure.  The shared color_registry ensures that any species or reaction
    #    that appears in both columns gets the same color, making cross-species
    #    comparisons visually consistent.
    # ---------------------------------------------------------------------------

    COMPARE_SPECIES = ["H2O", "CO"]
    color_registry: dict[str, str] = {}

    fig = plt.figure(figsize=(16, 13))
    fig.suptitle(
        "5 · draw_panel_* — side-by-side deepdive using individual panel functions",
        y=1.01,
    )

    for col, sp in enumerate(COMPARE_SPECIES):
        prod_rates_full, dest_rates_full = uclchem.analysis.get_production_and_destruction(sp, rates)
        prod_k_full, dest_k_full = uclchem.analysis.get_production_and_destruction(sp, rate_constants_df)

        pos_mask = physics_df["Time"] > 0
        time = physics_df["Time"][pos_mask]
        prod_rates = prod_rates_full[pos_mask]
        dest_rates = dest_rates_full[pos_mask].abs()
        prod_k = prod_k_full[pos_mask]
        dest_k = dest_k_full[pos_mask].abs()

        gs = fig.add_gridspec(3, 2, height_ratios=[3, 1, 3], hspace=0.05, wspace=0.35)

        ax_a = fig.add_subplot(gs[0, col])
        ax_c = fig.add_subplot(gs[1, col])
        ax_b = fig.add_subplot(gs[2, col])

        companion = [
            s for s in chemistry_df.columns
            if s != sp and (chemistry_df[s][pos_mask] > 0).any()
        ][:8]

        uclchem.plot.draw_panel_abundances(
            ax_a, time, sp, chemistry_df[pos_mask],
            companion, color_registry=color_registry,
        )
        uclchem.plot.draw_panel_rates(
            ax_b, time, prod_rates, dest_rates,
            color_registry=color_registry,
        )
        uclchem.plot.draw_panel_rate_constants(
            ax_c, time, prod_k, dest_k, color_registry=color_registry,
        )
        ax_a.set_title(sp, fontsize=11, fontweight="bold")

    fig.savefig(out_dir / "05_panel_comparison.png", bbox_inches="tight")
    plt.close(fig)

    print(f"All figures written to {out_dir.resolve()}/")
