# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import uclchem
from uclchem.makerates.network import Network, LoadedNetwork
from uclchem.makerates.species import Species
from uclchem.makerates.reaction import Reaction
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
# load species
species = pd.read_csv("../src/uclchem/species.csv")
species = [Species(list(spec)) for idx, spec in species.iterrows()]

# %%
# load reactions
reactions = pd.read_csv("../src/uclchem/reactions.csv")
# fix annoying electron freezeout
reactions.iloc[527, 3] = ""
reactions = [Reaction(list(reac)) for idx, reac in reactions.iterrows()]

# %%
reactions[0]

# %%
network = LoadedNetwork(species, reactions)

# %%
# print("Running test models...")
# # set a parameter dictionary for static model
# outSpecies = ["OH", "OCS", "CO", "CS", "CH3OH"]
# params = {
#     "endAtFinalDensity": False,
#     "freefall": False,
#     "writeStep": 1,
#     "initialDens": 1e4,
#     "initialTemp": 10.0,
#     "finalDens": 1e5,
#     "finalTime": 5.0e6,
#     "outputFile": "examples/test-output/static-full.dat",
#     "abundSaveFile": "examples/test-output/startstatic.dat",
#     "rateFile": "static_rates.csv",
# }
# uclchem.model.cloud(param_dict=params, out_species=outSpecies, return_rates=True)

# # change to collapsing phase1 params
# params["freefall"] = True
# params["endAtFinalDensity"] = True
# params["initialDens"] = 1e2
# params["abundSaveFile"] = "examples/test-output/startcollapse.dat"
# params["outputFile"] = "examples/test-output/phase1-full.dat"
# params["columnFile"] = "examples/test-output/phase1-column.dat"
# params["rateFile"] = "stage1_rates.csv"
# uclchem.model.cloud(param_dict=params, out_species=outSpecies, return_rates=True)

# # finally, run phase 2 from the phase 1 model.
# params["initialDens"] = 1e5
# params["freezeFactor"] = 0.0
# params["thermdesorb"] = True
# params["endAtFinalDensity"] = False
# params["freefall"] = False
# params["finalTime"] = 1e6
# params.pop("abundSaveFile")
# params["abundLoadFile"] = "examples/test-output/startcollapse.dat"
# params["outputFile"] = "examples/test-output/phase2-full.dat"
# params.pop("columnFile")
# params["rateFile"] = "stage2_rates.csv"
# uclchem.model.hot_core(
#     3, 300.0, param_dict=params, out_species=outSpecies, return_rates=True
# )

# %%
print("Running test models...")
# set a parameter dictionary for static model
outSpecies = ["OH", "OCS", "CO", "CS", "CH3OH"]
params = {
    "endAtFinalDensity": False,
    "freefall": False,
    "writeStep": 1,
    "initialDens": 1e4,
    "initialTemp": 10.0,
    "finalDens": 1e5,
    "finalTime": 5.0e6,
}
static_phys, static_chem, static_rates, _, succesflag = uclchem.model.cloud(
    param_dict=params, out_species=outSpecies, return_dataframe=True, return_rates=True
)

# change to collapsing phase1 params
params["freefall"] = True
params["endAtFinalDensity"] = True
params["initialDens"] = 1e2
phase1_phys, phase1_chem, phase1_rates, start_abunds, succesflag = uclchem.model.cloud(
    param_dict=params, out_species=outSpecies, return_dataframe=True, return_rates=True
)

# finally, run phase 2 from the phase 1 model.
params["initialDens"] = 1e5
params["freezeFactor"] = 0.0
params["thermdesorb"] = True
params["endAtFinalDensity"] = False
params["freefall"] = False
params["finalTime"] = 1e6
params["rateFile"] = "stage2_rates.csv"
phase2_phys, phase2_chem, phase2_rates, _, succesflag = uclchem.model.hot_core(
    3,
    300.0,
    param_dict=params,
    out_species=outSpecies,
    return_dataframe=True,
    return_rates=True,
    starting_chemistry=start_abunds,
)

# %%
# TODO: fix the string that is attached to reactions.
phase2_rates

# %%
from uclchem.analysis import postprocess_rates_to_dy

# %%
dy, flux = postprocess_rates_to_dy(
    phase2_phys, phase2_chem, phase2_rates, network=network
)
