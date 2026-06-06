"""Plot the UCLCHEM tests ran by `run_uclchem_tests.py`.

Demonstration of plotfunctions. called from main UCLCHEM directory.
It reads full UCLCHEM output and saves a plot of the abudances of select species.
"""

from pathlib import Path

import matplotlib.pyplot as plt

import uclchem

plot_types = {
    "simple_molecules": [
        "H",
        "$H",
        "H2",
        "$H2",
        "H2O",
        "$H2O",
        "CO",
        "$CO",
        "CH3OH",
        "$CH3OH",
    ],
    "charge": ["E-", "C+", "H2O+", "H+", "HE+", "HCO+"],
    "Silicon-bearing": [
        "SI",
        "$SI",
        # "SIO",
        # "$SIO",
        "SIH",
        "$SIH",
        "SIH2",
        "$SIH2",
        "SIH3",
        "$SIH3",
        "SIH4",
        "$SIH4",
    ],
}

if __name__ == "__main__":
    save_file = Path("examples/test-output/models.h5")
    model_names = {
        "phase1": "Collapsing Cloud",
        "phase2": "Hot Core",
        "static": "Static Cloud",
    }

    # Load test-output models from HDF5 and extract DataFrames
    test_data = {}
    for model in ["phase1", "phase2", "static"]:
        loaded = uclchem.model.load_model(file=str(save_file), name=model)
        test_data[model] = loaded.get_dataframes()
        for spec in ["#SI", "@SI"]:
            if spec not in test_data[model]:
                test_data[model][spec] = 1.0e-30

    # Load reference example-output from legacy .dat files
    example_data = {}
    print_elemental_conservation = True
    for model in ["phase1", "phase2", "static"]:
        example_data[model] = uclchem.analysis.read_output_file(
            f"examples/example-output/{model}-full.dat"
        )
        for spec in ["#SI", "@SI"]:
            if spec not in example_data[model]:
                example_data[model][spec] = 1.0e-30
        if print_elemental_conservation:
            print(f"Testing element conservation for {model}")
            print("Printing fractional change in total abundance")
            conservation = uclchem.analysis.check_element_conservation(
                example_data[model],
                element_list=["H", "N", "C", "O", "SI", "S"],
            )
            print(conservation)
            if model == "static":
                print_elemental_conservation = False

    for plot_type in plot_types:
        fig, axes = plt.subplots(3, 3, figsize=(16, 12), tight_layout=True)
        axes = axes.flatten()
        i = 0

        speciesNames = plot_types[plot_type]

        for folder, data_dict in [
            ("example-output/", example_data),
            ("test-output/", test_data),
        ]:
            for model in ["phase1", "phase2", "static"]:
                axis = axes[i]
                data = data_dict[model]

                axis = uclchem.plot.plot_species(axis, data, speciesNames, legend=False)
                if folder == "test-output/":
                    axis.set_prop_cycle(None)
                    axis = uclchem.plot.plot_species(
                        axis,
                        example_data[model],
                        speciesNames,
                        alpha=0.5,
                        legend=False,
                    )
                if plot_type == "charge":
                    ions = [s for s in list(data.columns) if "+" in s]
                    charge_conservation = data.loc[:, "E-"] - data.loc[:, ions].sum(
                        axis=1
                    )

                axis.set(
                    xscale="log",
                    ylim=(1e-15, 1),
                    xlim=(1, 6e6),
                    yscale="log",
                    ylabel=r"Fractional abundance $x_i (n_i/n_{\mathrm{H,nuclei}}) [\mathrm{cm}^{-3}]$",
                    xlabel="Time [yr]",
                )
                axis.grid()
                if i == 0:
                    axis.legend(bbox_to_anchor=(-0.25, 0.5), loc="center right")
                if model == "Stage 1":
                    axis.set(xlim=(1e3, 6e6))
                axis.set_title(model_names[model])
                if plot_type == "charge":
                    axis.set_title(
                        axis.get_title()
                        + f" (Charge conservation: {charge_conservation.mean():.2e})"
                    )
                i += 1

        axes[0].text(
            0.02,
            0.98,
            "Example Row",
            horizontalalignment="left",
            verticalalignment="top",
            transform=axes[0].transAxes,
        )
        axes[3].text(
            0.02,
            0.98,
            "Your Row",
            horizontalalignment="left",
            verticalalignment="top",
            transform=axes[3].transAxes,
        )

        # Add third row for temperature and density plots (test-output only)
        for j, model in enumerate(["phase1", "phase2", "static"]):
            axis = axes[6 + j]
            data = test_data[model]

            ax1 = axis
            ax2 = ax1.twinx()

            color = "tab:red"
            ax1.set_xlabel("Time [yr]")
            ax1.set_ylabel("Temperature [K]", color=color)
            line1 = ax1.plot(
                data["Time"], data["gasTemp"], color=color, label="Gas Temperature"
            )
            ax1.tick_params(axis="y", labelcolor=color)
            ax1.set_xscale("log")
            ax1.set_yscale("log")
            ax1.set_xlim(1, 6e6)
            ax1.grid(True, alpha=0.3)

            color = "tab:blue"
            ax2.set_ylabel("Density [cm⁻³]", color=color)
            line2 = ax2.plot(data["Time"], data["Density"], color=color, label="Density")
            ax2.tick_params(axis="y", labelcolor=color)
            ax2.set_yscale("log")

            ax1.set_title(f"{model_names[model]} - Temperature & Density")

            lines = line1 + line2
            labels = [line.get_label() for line in lines]
            if j == 0:
                ax1.legend(lines, labels, bbox_to_anchor=(-0.25, 0.5), loc="center right")

        axes[6].text(
            0.02,
            0.98,
            "Temperature & Density",
            horizontalalignment="left",
            verticalalignment="top",
            transform=axes[6].transAxes,
        )

        fig.savefig(f"examples/comparisons_{plot_type}.png", dpi=300)

    # Separate elemental conservation diagnostic plot (test-output)
    fig_ec, axes_ec = plt.subplots(1, 3, figsize=(15, 4), tight_layout=True)
    elements = ["H", "N", "C", "O", "SI", "S"]
    for j, model in enumerate(["phase1", "phase2", "static"]):
        ax = axes_ec[j]
        data = test_data[model]

        for element in elements:
            total = uclchem.analysis.total_element_abundance(element, data)
            ax.plot(data["Time"], total / total.iloc[0], label=element)

        ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.set_xscale("log")
        ax.set_xlim(1, 6e6)
        ax.set_xlabel("Time [yr]")
        ax.set_ylabel("Total abundance / initial")
        ax.set_title(model_names[model])
        ax.grid(True, alpha=0.3)
        if j == 0:
            ax.legend(bbox_to_anchor=(-0.25, 0.5), loc="center right")

    fig_ec.suptitle("Elemental Conservation (test-output)")
    fig_ec.savefig("examples/elemental_conservation.png", dpi=300)
