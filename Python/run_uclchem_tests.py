import uclchem
from time import perf_counter
from multiprocessing import Pool

#set a parameter dictionary for phase 1 collapse model
params = {"phase": 1, "switch": 0, "collapse": 0, "readAbunds": 0, "writeStep": 1,
			"outSpecies": 'OCS CO CS CH3OH', "initialDens": 1e4, "initialTemp":10.0,
			"finalDens":1e5, "finalTime":5.0e6,
			"outputFile":"examples/test-output/static-full.dat",
			"abundFile":"examples/test-output/startstatic.dat"
}

start=perf_counter()
uclchem.run_model(params.copy())
stop=perf_counter()
print(f"Static model in {stop-start:.1f} seconds")

#change to static cloud params
params["collapse"]=1
params["switch"]=1
params["initialDens"]=1e2
params["abundFile"]="examples/test-output/startcollapse.dat"
params["outputFile"]="examples/test-output/phase1-full.dat"
params["columnFile"]="examples/test-output/phase1-column.dat"
start=perf_counter()
uclchem.run_model(params.copy())
stop=perf_counter()
print(f"Phase 1 in {stop-start:.1f} seconds")

#finally, run phase 2 from the phase 1 model.
params["initialDens"]=1e5
params["tempindx"]=3
params["fr"]=0.0
#params["thermdesorb"]=0
params["readAbunds"]=1
params["phase"]=2
params["switch"]=0
params["finalTime"]=1e6
params["abundFile"]="examples/test-output/startcollapse.dat"
params["outputFile"]="examples/test-output/phase2-full.dat"
params.pop("columnFile")
start=perf_counter()
uclchem.run_model(params.copy())
stop=perf_counter()
print(f"Phase 2 in {stop-start:.1f} seconds")