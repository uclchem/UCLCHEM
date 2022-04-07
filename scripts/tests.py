import uclchem
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

print("ODE element conservation")
print("------------------------")

species=np.loadtxt("src/species.csv",usecols=[0],delimiter=",",
                   dtype=str,comments=None)
result=uclchem.tests.test_ode_conservation(species)
print("Total rates of change:")
for key,value in result.items():
    print(f"{key} {value:.2e}")
assert result["H"]<1e-15, f"H not conserved with total rate of change {result['H']:.2e}"
print("Elements are conserved.")
print("------------------------")


print("Running test models...")
#set a parameter dictionary for static model
outSpecies=['OCS','CO','CS','CH3OH']
params = {"switch": 0, "freefall": 0, "writeStep": 1,
			"initialDens": 1e4, "initialTemp":10.0,
			"finalDens":1e5, "finalTime":5.0e6,
			"outputFile":"examples/test-output/static-full-test.dat",
			"abundSaveFile":"examples/test-output/startstatic.dat"
}

start=perf_counter()
uclchem.cloud(params,outSpecies)
stop=perf_counter()
print(f"Static model in {stop-start:.1f} seconds")

#change to collapsing phase1 params
params["freefall"]=1
params["switch"]=1
params["initialDens"]=1e2
params["abundSaveFile"]="examples/test-output/startcollapse.dat"
params["outputFile"]="examples/test-output/phase1-full-test.dat"
params["columnFile"]="examples/test-output/phase1-column.dat"
start=perf_counter()
uclchem.cloud(params,outSpecies)
stop=perf_counter()
print(f"Phase 1 in {stop-start:.1f} seconds")

#finally, run phase 2 from the phase 1 model.
params["initialDens"]=1e5
params["fr"]=0.0
params["thermdesorb"]=1
params["switch"]=0
params["freefall"]=0
params["finalTime"]=1e6
params.pop("abundSaveFile")
params["abundLoadFile"]="examples/test-output/startcollapse.dat"
params["outputFile"]="examples/test-output/phase2-full-test.dat"
params.pop("columnFile")
start=perf_counter()
uclchem.hot_core(3,300.0,params,outSpecies)
stop=perf_counter()
print(f"Phase 2 in {stop-start:.1f} seconds")