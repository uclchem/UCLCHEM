# demonstration of plotfunctions. called from main UCLCHEM directory
# it reads full UCLCHEM output and saves a plot of the abudances of select species

import matplotlib.pyplot as plt

import uclchem

plot_types = {"simple_molecules": ["H", "$H", "H2", "$H2", "H2O", "$H2O", "CO", "$CO", "CH3OH", "$CH3OH"], 
              "charge": ["E-", "C+", "H2O+", "H+", "HE+", "HCO+"]}

print_elemental_conservation = True

for plot_type in plot_types:

    fig, axes = plt.subplots(3, 3, figsize=(16, 12), tight_layout=True)
    axes = axes.flatten()
    i = 0

    model_names = {
        "phase1": "Collapsing Cloud",
        "phase2": "Hot Core",
        "static": "Static Cloud",
    }
    model_data = {}
    speciesNames = plot_types[plot_type]
    for folder in ["example-output/", "test-output/"]:
        for model in ["phase1", "phase2", "static"]:
            axis = axes[i]
            # call read_uclchem.
            model_data[folder + model] = uclchem.analysis.read_output_file(
                "examples/" + folder + model + "-full.dat"
            )
            for spec in ["#SI", "@SI"]:
                if spec not in model_data[folder + model]:
                    model_data[folder + model][spec] = 1.0e-30
            # demonstrate checking element conservation
            if folder == "example-output/" and print_elemental_conservation:
                print(f"Testing element conservation for {model}")
                print("Printing fractional change in total abundance")
                conservation = uclchem.analysis.check_element_conservation(
                    model_data[folder + model]
                )
                print(conservation)
                # Only plot the elmental conservation once
                if model == "static":
                    print_elemental_conservation = False
            
            # plot species and save to test.png, alternatively send dens instead of time.
            axis = uclchem.analysis.plot_species(
                axis, model_data[folder + model], speciesNames, legend=False
            )
            if folder == "test-output/":
                axis.set_prop_cycle(None)
                axis = uclchem.analysis.plot_species(
                    axis, model_data["example-output/" + model], speciesNames, alpha=0.5, legend=False
                )
            if plot_type == "charge":
                ions = [s for s in list(model_data[folder+model].columns) if "+" in s]
                charge_conservation = model_data[folder+model].loc[:, "E-"] - model_data[folder+model].loc[:, ions].sum(axis=1)
            # plot species returns the axis so we can further edit
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
                axis.set_title(axis.get_title() + " (Charge conservation: {:.2e})".format(charge_conservation.mean()))
            i = i + 1
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
        axis = axes[6 + j]  # Third row starts at index 6
        
        # Get test-output data for this model
        data = model_data["test-output/" + model]
        
        # Create double y-axis plot
        ax1 = axis
        ax2 = ax1.twinx()
        
        # Plot temperature on left y-axis
        color = 'tab:red'
        ax1.set_xlabel('Time [yr]')
        ax1.set_ylabel('Temperature [K]', color=color)
        line1 = ax1.plot(data['Time'], data['gasTemp'], color=color,
                         label='Gas Temperature')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1, 6e6)
        ax1.grid(True, alpha=0.3)
        
        # Plot density on right y-axis
        color = 'tab:blue'
        ax2.set_ylabel('Density [cm⁻³]', color=color)
        line2 = ax2.plot(data['Time'], data['Density'], color=color,
                         label='Density')
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_yscale('log')
        
        # Set title
        ax1.set_title(f"{model_names[model]} - Temperature & Density")
        
        # Add legend combining both axes
        lines = line1 + line2
        labels = [line.get_label() for line in lines]
        if j == 0:  # Only add legend to first plot in row
            ax1.legend(lines, labels, bbox_to_anchor=(-0.25, 0.5),
                       loc="center right")
    
    # Add row label for third row
    axes[6].text(
        0.02,
        0.98,
        "Temperature & Density",
        horizontalalignment="left",
        verticalalignment="top",
        transform=axes[6].transAxes,
    )
    
    fig.savefig(f"examples/comparisons_{plot_type}.png", dpi=300)
