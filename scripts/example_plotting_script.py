"""
A script largely intended for people unfamiliar with Python.
If you run run_uclchem_tests.py, you'll produce several example model outputs.
This script then uses those outputs to demonstrate some simple plotting.
"""

import uclchem
import matplotlib.pyplot as plt

# pick species, any number is fine
species_list = ["#H2O", "#CO", "#CH3OH", "CO", "CH3OH", "H2O"]
input_file = "examples/example-output/phase1-full.dat"
input_file = "output/test.csv"


# call read_uclchem.
model_data = uclchem.analysis.read_output_file(input_file)

# create_abundance_plot will return pyplot figure and axis objects where the axis
# contains a plot of the species abundance through time for all species in species_list
# optionally, save it to plot_file
fig, ax = uclchem.analysis.create_abundance_plot(
    model_data, species_list, plot_file="examples/example_plot.png"
)


# alternatively, we may already have an axis we'd like to plot to
# in which case, plot_species() may be more helpful

fig, ax = plt.subplots()
ax = uclchem.analysis.plot_species(ax, model_data, species_list)

# the returned object lets us make some edits
# lets plot the temperature on a second y axis
ax2 = ax.twinx()
ax2.plot(model_data["Time"], model_data["gasTemp"], color="black")
ax2.set(ylabel="Temperature / K")
# and make some slight adjustments to the plot before saving again
ax.set(xscale="log", xlim=(1, 1e7), ylim=(9e-31, 5e-4))
ax.set_title("This is a Test Plot")
fig.savefig("examples/improved_example_plot.png")


# plot_species allows us to do more complex things such as subplots
fig, axes = plt.subplots(2, 2, figsize=(16, 9))
axes = axes.flatten()  # turn [2,2] array into [4]

for i, species in enumerate(
    [["CO", "#CO"], ["H2O", "#H2O"], ["CH3OH", "#CH3OH"], ["CO2", "#CO2"]]
):
    axes[i] = uclchem.plot_species(axes[i], model_data, species)
fig.savefig("examples/multiplot_example.png")
