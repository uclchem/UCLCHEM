#Marcus Keil and Jon Holdship 13/03/2020
#Examples of a simple grid of models run in parallel
from __future__ import print_function
import uclchem
import numpy as np
from multiprocessing import Pool
import time

#uclchem general takes a dictionary of parameters where outSpecies=number of outspecies
#and a string outspeciesin which is a delimited list of species
#this wrapper helps Pool to work and also converts a delimited list in the dictionary into a number
#and send list separately as required.
def run_uclchem(param_dict):
    outSpecies = (param_dict['outSpecies'])
    param_dict['outSpecies'] = len(outSpecies.split())
    uclchem.general(dictionary=param_dict, outspeciesin=outSpecies)
    return None


#basic set of parameters we'll use for this grid. 
ParameterDictionary = {"phase": 1, "switch": 1, "collapse": 1, "readAbunds": 0, "writeStep": 1,
                       "outSpecies": 'OCS CO CS CH3OH', "initialTemp": 10.0}


# This part can be substituted with any choice of grid
# here we just combine various initial and final densities into an easily iterable array
initialDensGrid = np.linspace(100, 1000, 10)
finalDensGrid = np.linspace(10000, 100000, 10)
parameterSpace = np.asarray(np.meshgrid(initialDensGrid, finalDensGrid)).reshape(2, -1)

#then we loop through parameters, update parameter dictionary and run
#however, to use Pool we just store dictionaries for now
models=[]

#we'll number our models so store the parameters we used in a list
with open("examples/grid/model_list.dat","w") as f:
    f.write("modelNo initDens final Dens\n")
    for k in range(np.shape(parameterSpace)[1]):
        paramDict=ParameterDictionary.copy()
        paramDict["initialDens"] = parameterSpace[0,k]
        paramDict["finalDens"] = parameterSpace[1,k]
        paramDict["abundFile"]= "examples/grid/startcollapseModel_"+str(k)+".dat"
        paramDict["outputFile"]= "examples/grid/phase1-fullModel_"+str(k)+".dat"
        paramDict["columnFile"]= "examples/grid/phase1-columnModel_"+str(k)+".dat"
        models.append(paramDict)
        f.write(f"{k} {parameterSpace[0,k]} {parameterSpace[1,k]}\n")

#use pool.map to run each dictionary throuh our helper function
start=time.time()
pool=Pool(4)
pool.map(run_uclchem,models)
pool.close()
pool.join()
end=time.time()
end=(end-start)/60.0
print(f"grid in {end} minutes")

