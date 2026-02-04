# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.5'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # DVODE Solver Statistics
# 
# UCLCHEM uses the DVODE ODE integrator to solve the chemical network at each timestep.
# This notebook demonstrates how to access per-timestep solver diagnostics, which are
# useful for understanding solver performance, identifying stiff timesteps, and diagnosing
# convergence issues.

# %%
import uclchem
import matplotlib.pyplot as plt
import numpy as np

# %% [markdown]
# ## Running a Model and Accessing Stats
# 
# Every model automatically collects solver statistics. Access them via `stats_array`
# (raw numpy array) or `get_dataframes(with_stats=True)` (pandas DataFrame).

# %%
params = {
    "endAtFinalDensity": False,
    "freefall": False,
    "initialDens": 1e5,
    "initialTemp": 10.0,
    "finalTime": 3e6,
    "writetimestep": True,

}

model = uclchem.model.Cloud(param_dict=params)
print(f"Model completed with flag: {model.success_flag}")

# %% [markdown]
# ## Viewing Stats as a DataFrame
# 
# The `get_dataframes` method supports a `with_stats` flag. When joined, the solver
# statistics are appended as extra columns alongside the physics and chemistry data.

# %%
df = model.get_dataframes(with_stats=True)
print(f"DataFrame shape: {df.shape}")
print(f"\nSolver stat columns: {uclchem.constants.DVODE_STAT_NAMES}")
df[["Time"] + uclchem.constants.DVODE_STAT_NAMES].head(10)

# %% [markdown]
# ## Plotting Solver Effort Over Time
# 
# The number of function evaluations (NFE) and steps (NST) per timestep show
# how hard the solver works at each point in the simulation. Spikes indicate
# chemically stiff regions.

# %%
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

active = df[df["Time"] > 0]

axes[0, 0].plot(active["Time"], active["NST"])
axes[0, 0].set_xlabel("Time (years)")
axes[0, 0].set_ylabel("Steps (NST)")
axes[0, 0].set_xscale("log")
axes[0, 0].set_title("ODE Steps per Timestep")

axes[0, 1].plot(active["Time"], active["NFE"])
axes[0, 1].set_xlabel("Time (years)")
axes[0, 1].set_ylabel("Function Evaluations (NFE)")
axes[0, 1].set_xscale("log")
axes[0, 1].set_title("RHS Evaluations per Timestep")

axes[1, 0].plot(active["Time"], active["NJE"])
axes[1, 0].set_xlabel("Time (years)")
axes[1, 0].set_ylabel("Jacobian Evaluations (NJE)")
axes[1, 0].set_xscale("log")
axes[1, 0].set_title("Jacobian Evaluations per Timestep")

axes[1, 1].plot(active["Time"], active["CPU_TIME"])
axes[1, 1].set_xlabel("Time (years)")
axes[1, 1].set_ylabel("CPU Time (s)")
axes[1, 1].set_xscale("log")
axes[1, 1].set_title("CPU Time per Timestep")

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Convergence Diagnostics
# 
# NCFN (convergence failures) and NETF (error test failures) indicate solver
# difficulties. Non-zero values suggest the chemistry is particularly stiff.

# %%
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(active["Time"], active["NCFN"], label="Convergence Failures (NCFN)")
ax.plot(active["Time"], active["NETF"], label="Error Test Failures (NETF)")
ax.set_xlabel("Time (years)")
ax.set_ylabel("Count")
ax.set_xscale("log")
ax.set_title("Solver Convergence Diagnostics")
ax.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Accessing Stats via the Raw Array
# 
# The raw stats array has shape `(timesteps, points, 18)` and is available
# directly on the model object.

# %%
print(f"Stats array shape: {model.stats_array.shape}")
print(f"Total function evaluations: {model.stats_array[:, 0, 6].sum():.0f}")
print(f"Total CPU time: {model.stats_array[:, 0, 17].sum():.3f} s")

# %% [markdown]
# ## Comparing Solver Effort Across Models
#
# Here we compare against a model with a temperature of 50K.

# %%
params["initialTemp"] = 50.0

model_50k = uclchem.model.Cloud(param_dict=params)
df_50k = model_50k.get_dataframes(with_stats=True)

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(active["Time"], active["NFE"], label="10K Static")
ax.plot(df_50k["Time"], df_50k["NFE"], label="50K static")
ax.set_xlabel("Time (years)")
ax.set_ylabel("Function Evaluations (NFE)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_title("Solver Effort: Static vs Freefall")
ax.legend()
plt.tight_layout()
plt.show()
