# %% [markdown]
# # Running a Grid
#
# A common task is to run UCLCHEM over a grid of parameter combinations. This notebook sets up a simple approach to doing so for regular grids.

# %%
import uclchem
import numpy as np
import pandas as pd
from multiprocessing import Pool
import os


# ## A Simple Grid
# ### Define Parameter Space
# First, we define our parameter space. We do this by using numpy and pandas to produce a table of all possible combinations of some parameters of interest.
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

# %% [markdown]
# ### Set up the model
# Next, we need a function that will run our model. We write a quick function that takes a row from our dataframe and uses it to populate a parameter dictionary for UCLCHEM and then run a cloud model. We can then map our dataframe to that function.


# %%
def run_model(row):
    i, row = row  # pandas iterrows actually come as tuples with the row number
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
    return result


with Pool(processes=6) as pool:
    results = pool.map(run_model, model_table.iterrows())

# ## Checking Your Grid
# After running, we should do two things. First, let's add `results` to our dataframe as a new column. Positive results mean a successful UCLCHEM run and negative ones are unsuccessful. Then we can run each model through `check_element_conservation` to check the integration was successful. We'll use both these things to flag models that failed in some way.


# %%
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

model_table.to_csv("my_grid.csv")
