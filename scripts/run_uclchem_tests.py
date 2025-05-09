import os
from time import perf_counter

import uclchem

if not os.path.exists("examples/test-output/"):
    os.makedirs("examples/test-output/")

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
    "outputFile": "examples/test-output/static-full.dat",
    "abundSaveFile": "examples/test-output/startstatic.dat",
}

start = perf_counter()
uclchem.model.cloud(param_dict=params, out_species=outSpecies)
stop = perf_counter()
print(f"Static model in {stop-start:.1f} seconds")

# change to collapsing phase1 params
params["freefall"] = True
params["endAtFinalDensity"] = True
params["initialDens"] = 1e2
params["abundSaveFile"] = "examples/test-output/startcollapse.dat"
params["outputFile"] = "examples/test-output/phase1-full.dat"
params["columnFile"] = "examples/test-output/phase1-column.dat"
start = perf_counter()
uclchem.model.cloud(param_dict=params, out_species=outSpecies)
stop = perf_counter()
print(f"Phase 1 in {stop-start:.1f} seconds")

# finally, run phase 2 from the phase 1 model.
params["initialDens"] = 1e5
params["freezeFactor"] = 0.0
params["thermdesorb"] = True
params["endAtFinalDensity"] = False
params["freefall"] = False
params["finalTime"] = 1e6
params.pop("abundSaveFile")
params["abundLoadFile"] = "examples/test-output/startcollapse.dat"
params["outputFile"] = "examples/test-output/phase2-full.dat"
params.pop("columnFile")
start = perf_counter()
uclchem.model.hot_core(3, 300.0, param_dict=params, out_species=outSpecies)
stop = perf_counter()
print(f"Phase 2 in {stop-start:.1f} seconds")
