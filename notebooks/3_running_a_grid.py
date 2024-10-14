# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: UCLCHEM 3.4.0 Release Candidate
#     language: python
#     name: uclchem_rc3.4.0
# ---

# # Running a Grid
#
# A common task is to run UCLCHEM over a grid of parameter combinations. This notebook sets up a simple approach to doing so for regular grids.

import uclchem
import numpy as np
import pandas as pd
import os
from joblib import Parallel, delayed


# ## A Simple Grid
# ### Define Parameter Space
# First, we define our parameter space. We do this by using numpy and pandas to produce a table of all possible combinations of some parameters of interest.

# +
# This part can be substituted with any choice of grid
# here we just vary the density, temperature and zeta
temperatures = np.linspace(10, 50, 3)
densities = np.logspace(4, 6, 3)
zetas = np.logspace(1, 3, 3)

# meshgrid will give all combinations, then we shape into columns and put into a table
parameterSpace = np.asarray(np.meshgrid(temperatures, densities, zetas)).reshape(3, -1)
model_table = pd.DataFrame(parameterSpace.T, columns=["temperature", "density", "zeta"])

# keep track of where each model output will be saved and make sure that folder exists
model_table["outputFile"] = model_table.apply(
    lambda row: f"../grid_folder/{row.temperature}_{row.density}_{row.zeta}.csv", axis=1
)
print(f"{model_table.shape[0]} models to run")
if not os.path.exists("../grid_folder"):
    os.makedirs("../grid_folder")
# -

model_table.head()


# ### Set up the model
# Next, we need a function that will run our model. We write a quick function that takes a row from our dataframe and uses it to populate a parameter dictionary for UCLCHEM and then run a cloud model. We can then map our dataframe to that function.


def run_model(row):
    # basic set of parameters we'll use for this grid.
    ParameterDictionary = {
        "endatfinaldensity": False,
        "freefall": False,
        "initialDens": row.density,
        "initialTemp": row.temperature,
        "zeta": row.zeta,
        "outputFile": row.outputFile,
        "finalTime": 1.0e6,
        "baseAv": 10,
    }
    result = uclchem.model.cloud(param_dict=ParameterDictionary)
    return result[0]  # just the integer error code


# ### Run Grid
#
# #### The Simple Way
# We can use pandas apply to simply pass each row to our helper function in turn. This will take some time since we're running the models one by one. I'll use the `head` function just to run five rows as an example here.

# %%time
result = model_table.head().apply(run_model, axis=1)

# #### The Fast Way
# Alternatively, we can use multiprocessing to run the models in parallel. That will allow us to run many models simulataneously and make use of all the cores available on our machine.

# %%time
results = Parallel(n_jobs=4, verbose=100)(
    delayed(run_model)(row) for idx, row in model_table.iterrows()
)


# ## Checking Your Grid
# After running, we should do two things. First, let's add `results` to our dataframe as a new column. Positive results mean a successful UCLCHEM run and negative ones are unsuccessful. Then we can run each model through `check_element_conservation` to check the integration was successful. We'll use both these things to flag models that failed in some way.


def element_check(output_file):
    df = uclchem.analysis.read_output_file(output_file)
    # get conservation values
    conserves = uclchem.analysis.check_element_conservation(df)
    # check if any error is greater than 1%
    return all([float(x[:-1]) < 1 for x in conserves.values()])


model_table["run_result"] = results
model_table["elements_conserved"] = model_table["outputFile"].map(element_check)
# check both conditions are met
model_table["Successful"] = (model_table.run_result >= 0) & (
    model_table.elements_conserved
)

model_table.head()

# ## Complex Grid
#
# The above was straightforward enough but what about a modelling a grid of shocks? Not only do we want to loop over relevant parameters, we also need to run a few preliminary models to give ourselves starting abundances. We'll start by defining two helper functions, one to run our preliminary cloud and one to run the shock.
#
# Let's further imagine that we want to obtain the abundances of several species at the end of the model. We can use the `out_species` parameter to specify which species we want to track and return them to our dataframe.

# +
out_species = ["CO", "H2O", "CH3OH"]


def run_prelim(density):
    # basic set of parameters we'll use for this grid.
    ParameterDictionary = {
        "endatfinaldensity": True,
        "freefall": True,
        "initialDens": 1e2,
        "finalDens": density,
        "initialTemp": 10.0,
        "abundSaveFile": f"../grid_folder/starts/{density:.0f}.csv",
        "baseAv": 1,
    }
    result = uclchem.model.cloud(param_dict=ParameterDictionary)
    return result


def run_model(row):
    # basic set of parameters we'll use for this grid.
    ParameterDictionary = {
        "endatfinaldensity": False,
        "freefall": False,
        "initialDens": row.density,
        "initialTemp": 10.0,
        "outputFile": row.outputFile,
        "abundLoadFile": f"../grid_folder/starts/{row.density:.0f}.csv",
        "finalTime": 1.0e5,
        "abstol_factor": 1e-18,
        "reltol": 1e-12,
        "baseAv": 1,
    }
    result = uclchem.model.cshock(
        row.shock_velocity, param_dict=ParameterDictionary, out_species=out_species
    )
    # First check UCLCHEM's result flag to seeif it ran succesfully, if it is return the abundances
    if result[0] == 0:
        return result[:]
    # if not, return NaNs because model failed
    else:
        return [np.nan] * len(out_species)


# -

# Then we define our parameter space again. We'll create two folders, one to store a set of initial abundances for each starting density in our model and another to store our shock outputs.

# +
# This part can be substituted with any choice of grid
# here we just combine various initial and final densities into an easily iterable array
shock_velocities = np.linspace(10, 50, 3)
densities = np.logspace(4, 6, 3)

parameterSpace = np.asarray(np.meshgrid(shock_velocities, densities)).reshape(2, -1)
model_table = pd.DataFrame(parameterSpace.T, columns=["shock_velocity", "density"])
model_table["outputFile"] = model_table.apply(
    lambda row: f"../grid_folder/shocks/{row.shock_velocity}_{row.density}.csv", axis=1
)
print(f"{model_table.shape[0]} models to run")

for folder in ["starts", "shocks"]:
    if not os.path.exists(f"../grid_folder/{folder}"):
        os.makedirs(f"../grid_folder/{folder}")
# -

# We can then run our preliminary models followed by our science models. The science models will return the abundances at the final time step of each run so we can unpack those directly to our dataframe.

results = Parallel(n_jobs=4, verbose=100)(
    delayed(run_prelim)(dens) for dens in densities
)

results = Parallel(n_jobs=4, verbose=100)(
    delayed(run_model)(row) for idx, row in model_table.iterrows()
)

model_table[["Result", "Dissipation Time"] + out_species] = results

model_table

# ## Summary
#
# There are many ways to run grids of models and users will naturally develop their own methods. This notebook is just a simple example of how to run UCLCHEM for many parameter combinations whilst producing a useful output (the model_table) to keep track of all the combinations that were run. In a real script, we'd save the model file to csv at the end.
#
# For much larger grids, it's recommended that you find a way to make your script robust to failure. Over a huge range of parameters, it is quite likely UCLCHEM will hit integration trouble for at least a few parameter combinations. Very occasionally, UCLCHEM will get caught in a loop where it fails to integrate and cannot adjust its strategy to manage it. This isn't a problem for small grids as the model can be stopped and the tolerances adjusted. However, for very large grids, you may end up locking all threads as they each get stuck on a different model. The best solution we've found for this case is to add a check so that models in your dataframe are skipped if their file already exists, this allows you to stop and restart the grid script as needed.
#
