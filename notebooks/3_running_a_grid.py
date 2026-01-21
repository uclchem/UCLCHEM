# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Running a Grid
#
# A common task is to run UCLCHEM over a grid of parameter combinations. This notebook shows how to use the built-in GridModels class to doing so for regular grids.

import uclchem
import numpy as np


# ## A Simple Grid
# ### Define Parameter Space
# First, we define our parameter space. We do this by using numpy arrays or functions to produce a table of all possible combinations of some parameters of interest.

ParameterDictionary = {
    "param_dict": {
        "endatfinaldensity": False,
        "freefall": False,
        "initialDens": np.logspace(4, 6, 3),
        "initialTemp": np.linspace(10, 50, 3),
        "zeta": np.logspace(1, 3, 3),
        "finalTime": 1.0e6,
        "baseAv": 10,
        "abstol_min": 1e-22,
    }
}
grid = uclchem.model.GridModels(
    model_type="Cloud",
    full_parameters=ParameterDictionary,
    max_workers=18,
    grid_file="./GridRun.h5",
    model_name_prefix="",
    delay_run=True,
)
grid.flat_grids

# 27 models to run

# # Run Grid
#
# Using the GridModels class, we can take advantage of setting up our varying parameter space with numpy arrays. We can define the type of model to run, how many concurrent workers should be used, where to store the output models as well as the naming convention of the models in the stored .h5 file. By default, instantiating a uclchem.model.GridModels object will call the ```run()``` method, but for the sake of demonstration, we have delayed the run using ```delay_run=True``` in the call to instantiate the object ```grid```. To run the models, we now call the ```run()``` method

grid.run()

# # Checking Your Grid
# After running, we can do two things. First, we can validate element conservation using the ```check_conservation``` method. Second, we can check the status of the individual models by inspecting the list of models that were run.

grid.check_conservation()
grid.models



# Each model is stored in the "./GridRun.h5" under the name listed in the key 'Model' for each of the entries of ```grid.models```. Each of these models can be loaded individually using the ```load_model``` function. If we wanted to load the model '9' with initialDens=10000, initialTemp=30 and zeta=10 we could do the following.

cloud = uclchem.model.load_model("./GridRun.h5", "9")

# Now ```cloud``` behaves as any model object should, allowing us to perform the same analyses and plotting as done in previous notebooks.

# # Complex Grid
# The above was straightforward enough but what about a modelling a grid of shocks? Not only do we want to loop over relevant parameters, we also need to run preliminary models to give ourselves starting abundances. We can do this by taking advantage of the ```SequentialModel``` class.
#
# This class can be used in isolation, in the following way.

# +
sequential_model_params = {
    "Cloud": {
        "param_dict": {
            "endAtFinalDensity": False,  # stop at finalTime
            "freefall": True,  # increase density in freefall
            "initialDens": 1e2,  # starting density
            "finalDens": 1e6,  # final density
            "initialTemp": 10.0,  # temperature of gas
            "finalTime": 6.0e6,  # final time
            "rout": 0.1,  # radius of cloud in pc
            "baseAv": 1.0,  # visual extinction at cloud edge.
        }
    },
    "CShock": {
        "param_dict": {
            "endAtFinalDensity": False,
            "freefall": False,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "abstol_factor": 1e-14,
            "abstol_min": 1e-20,
            "reltol": 1e-6,
            "baseAv": 1,
        },
        "shock_vel": 10,
    },
}

SequentialCShock = uclchem.model.SequentialModel(
    sequenced_model_parameters=sequential_model_params,
    parameters_to_match=["finalDens"],
)
# -

# Running this in a grid, works the same way as it does for normal models, except we can pass the additional ```parameters_to_match``` through the ```full_parameters``` dictionary.

# +
models_to_run = {
    "Cloud": {
        "param_dict": {
            "endAtFinalDensity": False,  # stop at finalTime
            "freefall": True,  # increase density in freefall
            "initialDens": 1e2,  # starting density
            "finalDens": np.logspace(4, 6, 3),  # final density
            "initialTemp": 10.0,  # temperature of gas
            "finalTime": 6.0e6,  # final time
            "rout": 0.1,  # radius of cloud in pc
            "baseAv": 1.0,  # visual extinction at cloud edge.
        }
    },
    "CShock": {
        "param_dict": {
            "endAtFinalDensity": False,
            "freefall": False,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "abstol_factor": 1e-14,
            "abstol_min": 1e-20,
            "reltol": 1e-6,
            "baseAv": 1,
        },
        "shock_vel": np.linspace(10, 50, 3),
    },
    "parameters_to_match": ["finalDens"],
}

complex_grid = uclchem.model.GridModels(
    model_type="SequentialModel",
    full_parameters=models_to_run,
    max_workers=10,
    grid_file="./ComplexGrid.h5",
)
# -

complex_grid.models

# In the case of SequentialModels being run in a grid, each individual model is saved using the naming convention of "<Model name>_<Model Class>_<Order in Sequence>" so the 6th model in our grid would consist of both "6_Cloud_0" and "6_CShock_1".

# ## Summary
#
# There are many ways to run grids of models and users will naturally develop their own methods. This notebook is just a simple example of how to run UCLCHEM for many parameter combinations whilst producing a useful output (the model_table) to keep track of all the combinations that were run. In a real script, we'd save the model file to csv at the end.
#
# For much larger grids, it's recommended that you find a way to make your script robust to failure. Over a huge range of parameters, it is quite likely UCLCHEM will hit integration trouble for at least a few parameter combinations. Very occasionally, UCLCHEM will get caught in a loop where it fails to integrate and cannot adjust its strategy to manage it. This isn't a problem for small grids as the model can be stopped and the tolerances adjusted. However, for very large grids, you may end up locking all threads as they each get stuck on a different model. The best solution we've found for this case is to add a check so that models in your dataframe are skipped if their file already exists, this allows you to stop and restart the grid script as needed.
#
