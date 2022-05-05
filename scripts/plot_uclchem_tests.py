# demonstration of plotfunctions. called from main UCLCHEM directory
# it reads full UCLCHEM output and saves a plot of the abudances of select species

import uclchem
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 3, figsize=(16, 9), tight_layout=True)
axes = axes.flatten()
i = 0

model_names = {"phase1": "Collapsing Cloud", "phase2": "Hot Core", "static": "Static Cloud"}
for folder in ["example-output/", "test-output/"]:
    for model in ["phase1", "phase2", "static"]:
        axis = axes[i]
        # pick species, any number is fine
        speciesNames = ["H", "H2", "$H", "H2O", "$H2O", "CO", "$CO", "$CH3OH", "CH3OH"]

        # call read_uclchem.
        model_data = uclchem.analysis.read_output_file("examples/" + folder + model + "-full.dat")
        for spec in ["#SI", "@SI"]:
            if spec not in model_data:
                model_data[spec] = 1.0e-30
        # demonstrate checking element conservation
        if folder == "example-output/":
            print(f"Testing element conservation for {model}")
            print("Printing fractional change in total abundance")
            conservation = uclchem.analysis.check_element_conservation(model_data)
            print(conservation)

        # plot species and save to test.png, alternatively send dens instead of time.
        axis = uclchem.analysis.plot_species(axis, model_data, speciesNames)

        # plot species returns the axis so we can further edit
        axis.set(xscale="log", ylim=(1e-15, 1), xlim=(1, 6e6), yscale="log")
        axis.legend()
        if model == "phase1":
            axis.set(xlim=(1e3, 6e6))
        axis.set_title(model_names[model])
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
fig.savefig("examples/example-comparisons.png", dpi=300)
