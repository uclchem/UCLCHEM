# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Heating and Cooling in UCLCHEM
#
# This notebook demonstheating how to use UCLCHEM with heating and cooling mechanisms enabled. UCLCHEM includes a comprehensive set of heating and cooling processes that can significantly affect the temperature evolution and chemistry in astrophysical environments.
#
# ## Overview
#
# UCLCHEM can track gas temperature changes due to:
#
# **Heating processes:**
# - Photoelectric heating from PAHs and dust grains
# - H₂ formation heating
# - H₂ photodissociation heating
# - Cosmic ray heating
# - Chemical reaction enthalpy (optional)
# - X-ray heating
# - Turbulent heating
# - And more...
#
# **Cooling processes:**
# - Line cooling (CO, H₂O, OH, etc.)
# - Continuum cooling from dust
# - Atomic fine structure cooling (C⁺, O, etc.)
# - Recombination cooling
# - And more...
#
# By default, heating is enabled (`heatingFlag=True`), but you can turn it off or customize which processes are included.

# %%
import uclchem
import matplotlib.pyplot as plt

# %% [markdown]
# ## Example 1: Basic Model with Heating Enabled
#
# Let's start with a simple cloud model with heating and cooling enabled (the default behavior). We'll model a molecular cloud core at constant density.

# %%
# Define parameters for a molecular cloud with heating enabled
param_dict_heating = {
    "initialDens": 1e4,  # Initial density (cm^-3)
    "initialTemp": 10.0,  # Initial temperature (K)
    "finalTime": 1.0e6,  # Final time (years)
    "baseAv": 10.0,  # Visual extinction
    "rout": 0.1,  # Outer radius (pc)
    "freefall": False,  # Keep constant density
    "heatingFlag": True,  # Enable heating (default)
}

# Run the model using the Cloud class
cloud = uclchem.model.Cloud(param_dict=param_dict_heating)

# Extract data from the model object
physics, abundances = cloud.get_dataframes(joined=False)
_, _, rates = cloud.get_dataframes(joined=False, with_rates=True)
_, _, heating = cloud.get_dataframes(joined=False, with_heating=True)
start_abund = cloud.next_starting_chemistry_array
flag = 0 if cloud.has_attr("_data") else -1

print(f"Model completed with flag: {flag}")
if flag < 0:
    print(f"Error: Model failed to complete")
else:
    print("Success! Model completed.")
    print(
        f"\nTemperature range: {physics['gasTemp'].min():.2f} - {physics['gasTemp'].max():.2f} K"
    )
    print(
        f"Density range: {physics['Density'].min():.2e} - {physics['Density'].max():.2e} cm^-3"
    )

# %% [markdown]
# ### Visualizing Temperature Evolution
#
# The `physics` DataFrame contains physical properties including gas temperature. Let's plot how temperature evolves over time.

# %%
# Plot temperature evolution
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(physics["Time"], physics["gasTemp"], linewidth=2, color="darkred")
ax.set_xlabel("Time (years)", fontsize=12)
ax.set_ylabel("Gas Temperature (K)", fontsize=12)
ax.set_xscale("log")
ax.set_title("Gas Temperature Evolution with Heating/Cooling", fontsize=14)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Example 2: Comparing Heating On vs Off
#
# To see the impact of heating/cooling processes, let's run two models side by side: one with heating enabled and one without.

# %%
# Model WITHOUT heating
param_dict_no_heating = param_dict_heating.copy()
param_dict_no_heating["heatingFlag"] = False

cloud_no_heat = uclchem.model.Cloud(param_dict=param_dict_no_heating)

# Extract data
physics_no_heat, abundances_no_heat = cloud_no_heat.get_dataframes(joined=False)
_, _, rates_no_heat = cloud_no_heat.get_dataframes(joined=False, with_rates=True)
_, _, heating_no_heat = cloud_no_heat.get_dataframes(joined=False, with_heating=True)
flag_no_heat = 0 if cloud_no_heat.has_attr("_data") else -1

print(f"No heating model completed with flag: {flag_no_heat}")

# Compare the two runs
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Temperature comparison
axes[0].plot(physics["Time"], physics["gasTemp"], label="Heating ON", linewidth=2)
axes[0].plot(
    physics_no_heat["Time"],
    physics_no_heat["gasTemp"],
    label="Heating OFF",
    linewidth=2,
    linestyle="--",
)
axes[0].set_xlabel("Time (years)", fontsize=12)
axes[0].set_ylabel("Gas Temperature (K)", fontsize=12)
axes[0].set_xscale("log")
axes[0].set_title("Temperature Evolution", fontsize=14)
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Compare a key chemical species (CO)
axes[1].plot(physics["Time"], abundances["H2O"], label="Heating ON", linewidth=2)
axes[1].plot(
    physics_no_heat["Time"],
    abundances_no_heat["H2O"],
    label="Heating OFF",
    linewidth=2,
    linestyle="--",
)
axes[1].set_xlabel("Time (years)", fontsize=12)
axes[1].set_ylabel("H$_2$O Abundance", fontsize=12)
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_title("CO Abundance Evolution", fontsize=14)
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print("\nFinal temperatures:")
print(f"  With heating: {physics['gasTemp'].iloc[-1]:.2f} K")
print(f"  Without heating: {physics_no_heat['gasTemp'].iloc[-1]:.2f} K")

# %% [markdown]
# ## Example 3: Accessing Heating and Cooling Rates
#
# UCLCHEM can return detailed heating and cooling rates when you use `return_dataframe=True` with the rates arrays. The `heating` DataFrame includes individual contributions from each heating and cooling process.

# %%
# The heating DataFrame contains heating and cooling information
print("Available heating/cooling columns:")
print(heating.columns.tolist()[:20], "...")  # Show first 20 columns

# Extract heating and cooling columns
heating_cols = [
    col for col in heating.columns if "heating" in col.lower() and col != "Time"
]
cooling_cols = [col for col in heating.columns if "cooling" in col.lower()]

print(f"\nFound {len(heating_cols)} heating processes")
print(f"Found {len(cooling_cols)} cooling processes")

# Show the heating processes
print("\nHeating processes:")
for col in heating_cols[:10]:  # Show first 10
    print(f"  - {col}")

# %%
heating

# %% [markdown]
# ### Plotting Individual Heating and Cooling Contributions
#
# Let's visualize the dominant heating and cooling processes over time.

# %%
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Plot heating processes
time = physics["Time"]
for col in heating_cols:
    # Only plot processes with significant contribution
    max_val = heating[col].abs().max()
    if max_val > 1e-30:  # Filter out negligible contributions
        label = col.replace("_heating", "").replace("_", " ")
        axes[0].plot(time, heating[col], label=label, linewidth=2, alpha=0.7)

axes[0].set_xlabel("Time (years)", fontsize=12)
axes[0].set_ylabel("Heating Rate (erg cm⁻³ s⁻¹)", fontsize=12)
axes[0].set_xscale("log")
axes[0].set_yscale("symlog", linthresh=1e-30)
axes[0].set_title("Heating Processes", fontsize=14)
axes[0].legend(fontsize=8, loc="best")
axes[0].grid(True, alpha=0.3)

# Plot cooling processes (select major ones to avoid cluttering)
for col in cooling_cols:
    max_val = heating[col].abs().max()
    if max_val > 1e-40:
        label = col.replace("_cooling", "").replace("_", " ")
        axes[1].plot(time, heating[col], label=label, linewidth=2, alpha=0.7)

axes[1].set_xlabel("Time (years)", fontsize=12)
axes[1].set_ylabel("Cooling Rate (erg cm⁻³ s⁻¹)", fontsize=12)
axes[1].set_xscale("log")
axes[1].set_yscale("symlog", linthresh=1e-30)
axes[1].set_title("Major Cooling Processes", fontsize=14)
axes[1].legend(fontsize=8, loc="best")
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Example 4: Writing Heating Rates to a File
#
# You can also save heating and cooling rates directly to a file during the model run using the `heatingFile` parameter. This requires not running any other cells to memory before!

# %%
# You can uncomment the cel below to run a model to disk.

# # Run model with heating file output
# param_dict_with_file = param_dict_heating.copy()
# param_dict_with_file["outputFile"] = "../examples/test-output/heating_demo.dat"
# param_dict_with_file["heatingFile"] = "../examples/test-output/heating_rates.dat"


# result = uclchem.model.cloud(
#     param_dict=param_dict_with_file, out_species=["CO", "H2O", "OH"]
# )

# print(f"Model completed with flag: {result[0]}")
# print(f"Final CO abundance: {result[1]:.2e}")
# print(f"Final H2O abundance: {result[2]:.2e}")
# print(f"Final OH abundance: {result[3]:.2e}")

# # Read the heating file
# if os.path.exists("../examples/test-output/heating_rates.dat"):
#     print("\nHeating rates file created successfully!")
#     heating_df = pd.read_csv("../examples/test-output/heating_rates.dat")
#     print(f"File contains {len(heating_df)} time steps")
#     print(f"Columns: {heating_df.columns.tolist()[:10]}...")  # Show first 10 columns

# %% [markdown]
# ## Example 5: Different Physical Conditions
#
# Let's explore how heating and cooling behave under different physical conditions. We'll compare a dense, cold core with a less dense, warmer environment.

# %%
# Dense, cold core (like Example 1)
param_cold_core = {
    "initialDens": 1e5,
    "initialTemp": 10.0,
    "finalTime": 1.0e6,
    "baseAv": 20.0,
    "rout": 0.05,
    "freefall": False,
    "heatingFlag": True,
}

# Less dense, warmer cloud
param_warm_cloud = {
    "initialDens": 1e3,
    "initialTemp": 30.0,
    "finalTime": 1.0e6,
    "baseAv": 5.0,
    "rout": 0.5,
    "freefall": False,
    "heatingFlag": True,
}

# Run both models
print("Running cold core model...")
cloud_cold = uclchem.model.Cloud(param_dict=param_cold_core)
phys_cold, abund_cold = cloud_cold.get_dataframes(joined=False)
_, _, rates_cold = cloud_cold.get_dataframes(joined=False, with_rates=True)
_, _, heating_cold = cloud_cold.get_dataframes(joined=False, with_heating=True)
flag_cold = 0 if cloud_cold.has_attr("_data") else -1

print("Running warm cloud model...")
cloud_warm = uclchem.model.Cloud(param_dict=param_warm_cloud)
phys_warm, abund_warm = cloud_warm.get_dataframes(joined=False)
_, _, rates_warm = cloud_warm.get_dataframes(joined=False, with_rates=True)
_, _, heating_warm = cloud_warm.get_dataframes(joined=False, with_heating=True)
flag_warm = 0 if cloud_warm.has_attr("_data") else -1

print(f"\nCold core flag: {flag_cold}, Warm cloud flag: {flag_warm}")

# Compare temperature evolution
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(
    phys_cold["Time"],
    phys_cold["gasTemp"],
    label=f"Dense Core (n={param_cold_core['initialDens']:.0e})",
    linewidth=2,
)
axes[0].plot(
    phys_warm["Time"],
    phys_warm["gasTemp"],
    label=f"Diffuse Cloud (n={param_warm_cloud['initialDens']:.0e})",
    linewidth=2,
)
axes[0].set_xlabel("Time (years)", fontsize=12)
axes[0].set_ylabel("Gas Temperature (K)", fontsize=12)
axes[0].set_xscale("log")
axes[0].set_title("Temperature Evolution", fontsize=14)
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Compare CO abundance
axes[1].plot(phys_cold["Time"], abund_cold["CO"], label="Dense Core", linewidth=2)
axes[1].plot(phys_warm["Time"], abund_warm["CO"], label="Diffuse Cloud", linewidth=2)
axes[1].set_xlabel("Time (years)", fontsize=12)
axes[1].set_ylabel("CO Abundance", fontsize=12)
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_title("CO Abundance", fontsize=14)
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print("\nFinal temperatures:")
print(f"  Cold core: {phys_cold['gasTemp'].iloc[-1]:.2f} K")
print(f"  Warm cloud: {phys_warm['gasTemp'].iloc[-1]:.2f} K")

# %% [markdown]
# ## Summary of Key Heating and Cooling Processes
#
# UCLCHEM includes the following heating and cooling mechanisms:
#
# ### Major Heating Processes:
# 1. **Photoelectric heating** - from PAHs and dust grains (Bakes & Tielens 1994)
# 2. **H₂ formation heating** - energy released during H₂ formation on grains
# 3. **H₂ photodissociation heating** - kinetic energy from photodissociation
# 4. **Cosmic ray heating** - direct and indirect (via ionization)
# 5. **Chemical reaction enthalpy** - exothermic reactions (optional, requires network configuration)
# 6. **X-ray heating** - in environments with X-ray sources
# 7. **Turbulent heating** - dissipation of turbulent energy
# 8. **Viscous heating** - for shock models
#
# ### Major Cooling Processes:
# 1. **Molecular line cooling** - CO, H₂O, OH, and other molecules
# 2. **Atomic fine structure cooling** - C⁺, O, C, etc.
# 3. **Dust continuum cooling** - gas-grain collisional cooling
# 4. **H₂ line cooling** - rovibrational transitions
# 5. **Recombination cooling** - electronic to kinetic energy
# 6. **Ly-α cooling** - Lyman-alpha photon emission
#
# ### Key Parameters:
# - `heatingFlag`: Enable/disable all heating processes (default: `True`)
# - `heatingFile`: Path to write detailed heating/cooling rates
# - `initialTemp`: Starting gas temperature
# - `baseAv`: Visual extinction (affects radiation field attenuation)
# - `zeta`: Cosmic ray ionization rate
#
# For a complete list of processes and their equations, see the `HEATING_COOLING_SUMMARY.md` documentation file.

# %% [markdown]
# ## Advanced: Chemical Reaction Enthalpy
#
# UCLCHEM can optionally include the enthalpy changes from chemical reactions as a heating/cooling source. This requires configuring the chemical network at build time using MakeRates.
#
# **Note:** This feature requires rebuilding UCLCHEM with specific settings in the `user_settings.yaml` file. The main settings are:
#
# - `add_delta_enthalpy`: Can be `False`, `GAS`, or `ALL`
#   - `False`: No reaction enthalpy (default)
#   - `GAS`: Include enthalpy from gas-phase reactions only
#   - `ALL`: Include enthalpy from all reactions (gas + surface)
#
# To enable reaction enthalpy:
# ```yaml
# # In Makerates/user_settings.yaml
# add_delta_enthalpy: GAS
# ```
#
# Then rebuild:
# ```bash
# cd Makerates
# python MakeRates.py
# cd ..
# pip install .
# ```
#
# For more details, see the `heating_cooling_benchmarks/` directory which contains examples with different enthalpy configurations.

# %% [markdown]
# ## Tips for Using Heating and Cooling
#
# 1. **Always check convergence**: Use `return_dataframe=True` and inspect temperature evolution to ensure physical results
#
# 2. **Monitor element conservation**: Use `uclchem.analysis.check_element_conservation()` to verify model accuracy
#
# 3. **Save heating rates for analysis**: Use `heatingFile` parameter to save detailed heating/cooling contributions
#
# 4. **Choose appropriate tolerances**: If you see unphysical temperature spikes, try adjusting `reltol` and `abstol` parameters
#
# 5. **Consider your environment**: Different astrophysical environments have different dominant heating/cooling processes:
#    - **Dense cores**: Line cooling (CO, H₂O) dominates
#    - **Diffuse clouds**: Photoelectric heating and fine-structure cooling
#    - **Shock regions**: Viscous/turbulent heating important
#    - **PDRs**: Strong FUV field means photoelectric heating dominates
#
# 6. **Radiation field matters**: The `baseAv` parameter controls UV attenuation, which significantly affects photoelectric heating and photodissociation rates
#
# ## Next Steps
#
# - Check out [7_heating_cooling_settings.ipynb](7_heating_cooling_settings.ipynb) for more advanced configuration options
# - See `HEATING_COOLING_SUMMARY.md` for detailed equations and references
# - Explore the `heating_cooling_benchmarks/` directory for systematic comparisons
