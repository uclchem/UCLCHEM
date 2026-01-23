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
# # Advanced Settings with GeneralSettings
#
# This notebook demonstrates how to use the `advanced.GeneralSettings` class to modify UCLCHEM parameters instead of using parameter dictionaries. This provides a more object-oriented interface with better autocomplete, type checking, and state management.

# %%
import uclchem
from uclchem import advanced
import os

# Ensure output directory exists
if not os.path.exists("output_6"):
    os.makedirs("output_6")

# %% [markdown]
# ## Using GeneralSettings Instead of param_dict
#
# The traditional way to configure UCLCHEM models is with parameter dictionaries that are passed when running the UCLCHEM model. The `GeneralSettings` class provides an alternative that directly modifies the Fortran module state before running a model. This has the benefit that you can change things globally, and access things that are not exposted through parameter dictionaries.
#
# **Key differences:**
# - Direct access to Fortran variables with autocomplete support
# - Tracks which settings have been modified
# - Provides context managers for temporary changes
# - Settings persist across model runs (not reset automatically)

# %%
# Create settings interface
settings = advanced.GeneralSettings()

# View current value
print(f"Current initial density: {settings.defaultparameters.initialdens.get()}")
print(f"Current initial temperature: {settings.defaultparameters.initialtemp.get()}")

# %% [markdown]
# ## Setting Parameters Directly
#
# Instead of creating a `param_dict`, we can set parameters directly:

# %%
# Set model parameters
settings.defaultparameters.initialdens = 1e4
settings.defaultparameters.initialtemp = 10.0
settings.defaultparameters.finaltime = 1.0e6
settings.defaultparameters.rout = 0.1
settings.defaultparameters.baseav = 1.0
settings.defaultparameters.freefall = False
settings.defaultparameters.endatfinaldensity = False

# Verify changes
print(f"New initial density: {settings.defaultparameters.initialdens.get()}")

# %% [markdown]
# ## Running a Model with GeneralSettings
#
# When using `GeneralSettings`, the Fortran code uses the values we've set directly. However, **file paths (outputFile, abundSaveFile) should still be set via `param_dict`** since these are handled specially by the model wrapper and file parsers.

# %%
# Ensure output directory exists
if not os.path.exists("output_6"):
    os.makedirs("output_6")

# Note: Output file paths should be set via param_dict, not GeneralSettings
# This is because file paths are handled specially by the model wrapper
param_dict = {
    "outputFile": "output_6/baseline.dat",
}

# Run model with param_dict for file I/O, but using GeneralSettings for other parameters
out_species = ["SO", "CO"]
cloud = uclchem.model.Cloud(param_dict=param_dict, out_species=out_species)
print(f"Model completed successfully")

# %% [markdown]
# ## Viewing Modified Settings
#
# You can check which settings have been changed from their defaults:

# %%
# Show all modified settings
settings.print_all_edited()

# %% [markdown]
# ## Using Context Managers for Temporary Changes
#
# The `temporary_changes()` context manager is useful for running models with temporary parameter modifications without affecting the global state:

# %%
# Set a baseline density
settings.defaultparameters.initialdens = 1e4
print(f"Baseline density: {settings.defaultparameters.initialdens.get()}")

# Run model with temporary higher density
with settings.temporary_changes():
    settings.defaultparameters.initialdens = 1e5
    print(f"Inside context: {settings.defaultparameters.initialdens.get()}")

    # Use param_dict for file paths
    param_dict_high = {"outputFile": "output_6/high_density.dat"}
    cloud_high = uclchem.model.Cloud(
        param_dict=param_dict_high, out_species=out_species
    )
    print(f"High density model completed successfully")

# Settings automatically restored
print(f"After context: {settings.defaultparameters.initialdens.get()}")

# %% [markdown]
# ## Comparing Results
#
# Let's load and compare the two model runs:

# %%
result_baseline = uclchem.analysis.read_output_file("output_6/baseline.dat")
result_high_dens = uclchem.analysis.read_output_file("output_6/high_density.dat")

print(f"Baseline final CO abundance: {result_baseline['CO'].iloc[-1]:.2e}")
print(f"High density final CO abundance: {result_high_dens['CO'].iloc[-1]:.2e}")

# %% [markdown]
# ## Resetting Settings
#
# You can reset individual settings or all settings to defaults:

# %%
# Reset a single setting
print(f"Before reset: {settings.defaultparameters.initialdens.get()}")
settings.defaultparameters.initialdens.reset()
print(f"After reset: {settings.defaultparameters.initialdens.get()}")

# %% [markdown]
# ## Searching for Settings
#
# You can search across all modules for settings matching a pattern:

# %%
# Search for settings related to "temp"
temp_settings = settings.search("temp")
print(f"Found {len(temp_settings)} settings related to 'temp':")
for name, setting in list(temp_settings.items())[:5]:
    print(f"  {name}: {setting.get()}")

# %% [markdown]
# ## Best Practices
#
# **When to use GeneralSettings:**
# - You want to run multiple models with the same base configuration
# - You need to modify settings that aren't easily accessible via `param_dict`
# - You want autocomplete support in your IDE
# - You're building a GUI or interactive application
#
# **When to use param_dict:**
# - One-off model runs
# - Grid computations where each run has different parameters
# - **File paths (outputFile, abundSaveFile)** - always use param_dict for these
# - You want parameters explicitly documented in the function call
#
# **Important Notes:**
# - ⚠️ Settings are **global** and persist across model runs
# - ⚠️ **Not thread-safe** - don't use with multiprocessing
# - ⚠️ **File I/O paths** should always be set via `param_dict`, not GeneralSettings
# - Always reset or use context managers when running multiple models

# %% [markdown]
# ## Summary
#
# The `GeneralSettings` class provides a powerful alternative to parameter dictionaries for configuring UCLCHEM models. Key features:
#
# - Direct access to all Fortran module variables
# - Edit tracking and reset capabilities
# - Context managers for temporary changes  
# - Search functionality across all settings
# - Better IDE autocomplete support
#
# See the documentation for more details on the `advanced` module.
