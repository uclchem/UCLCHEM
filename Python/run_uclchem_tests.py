import uclchem
import numpy as np
import pandas as pd
from multiprocessing import Pool

#set a parameter dictionary for phase 1 collapse model
params = {"phase": 1, "switch": 1, "collapse": 1, "readAbunds": 0, "writeStep": 1,
			"outSpecies": 'SO CO', "initialDens": 1e2, "initialTemp":10.0,
			"finalDens":1e5, "finalTime":5.0e6,
			"outputFile":"examples/test-output/phase1-full.dat",
			"abundFile":"examples/test-output/startcollapse.dat"}

uclchem.run_model(params.copy())

#change to static cloud params
params["collapse"]=0
params["switch"]=0
params["initialDens"]=1e4
params["abundFile"]="examples/test-output/startstatic.dat"

params["outputFile"]="examples/test-output/static-full.dat"
uclchem.run_model(params.copy())

#finally, run phase 2 from the phase 1 model.
params["initialDens"]=1e5
params["tempindx"]=3
params["fr"]=0.0
#params["thermdesorb"]=0
params["readAbunds"]=1
params["phase"]=2
params["finalTime"]=1e6
params["abundFile"]="examples/test-output/startcollapse.dat"
params["outputFile"]="examples/test-output/phase2-full.dat"
uclchem.run_model(params.copy())