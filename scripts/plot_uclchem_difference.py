# demonstration of plotfunctions. called from main UCLCHEM directory
# it reads full UCLCHEM output and saves a plot of the abudances of select species

import uclchem
import matplotlib.pyplot as plt
import pandas as pd

fig, axes = plt.subplots(2, 3, figsize=(16, 9), tight_layout=True, sharey="row")
axes = axes.flatten()
i = 0

model_names = {
    "phase1": "Collapsing Cloud",
    "phase2": "Hot Core",
    "static": "Static Cloud",
}
# for folder in ["example-output/", "test-output/"]:
for model in ["phase1", "phase2", "static"]:
    # pick species, any number is fine
    speciesNames = ["H", "H2", "$H", "H2O", "$H2O", "CO", "$CO", "$CH3OH", "CH3OH"]

    # call read_uclchem.
    model_data1 = uclchem.analysis.read_output_file(
        "examples/" + "example-output/" + model + "-full.dat"
    )
    model_data2 = uclchem.analysis.read_output_file(
        "examples/" + "test-output/" + model + "-full.dat"
    )

    for spec in ["#SI", "@SI"]:
        if spec not in model_data1:
            model_data1[spec] = 1.0e-30
        if spec not in model_data2:
            model_data2[spec] = 1.0e-30
    # demonstrate checking element conservation
    # if folder == "example-output/":
    #     print(f"Testing element conservation for {model}")
    #     print("Printing fractional change in total abundance")
    #     conservation = uclchem.analysis.check_element_conservation(model_data)
    #     print(conservation)
    diff_data = pd.DataFrame()
    reldiff_data = pd.DataFrame()
    for column in model_data1.columns:
        if column in ["Time", "av", "gasTemp", "point", "radfield", "zeta"]:
            if column == "Time":
                print(model_data1["Time"] - model_data2["Time"])
                # assert sum(abs(model_data1["Time"] - model_data2["Time"])) == 0.0
            diff_data[column] = model_data2[column]
            reldiff_data[column] = model_data2[column]
        else:
            diff_data[column] = abs(model_data1[column] - model_data2[column])
            reldiff_data[column] = (
                abs(model_data1[column] - model_data2[column]) / model_data2[column]
            )

    # plot species and save to test.png, alternatively send dens instead of time.
    axis = uclchem.analysis.plot_species(axes[i], diff_data, speciesNames)

    # plot species returns the axis so we can further edit
    axis.set(xscale="log", ylim=(1e-15, 1), xlim=(1, 6e6), yscale="log")
    axis.legend()
    if model == "phase1":
        axis.set(xlim=(1e3, 6e6))
    axis.set_title(model_names[model])
    axis.set_title(f"Plot number {i}")
    axis = uclchem.analysis.plot_species(axes[i + 3], reldiff_data, speciesNames)
    axis.set_title(f"Plot number {i+3}")
    i = i + 1
axes[0].set_ylabel("Absolute difference")
axes[3].set_ylabel("Relative difference")
fig.savefig("examples/example-difference.png", dpi=300)
