#Marcus Keil and Jon Holdship 13/03/2020
#Examples of a simple grid of models run in parallel
from __future__ import print_function
import uclchem
import numpy as np
import pandas as pd
from multiprocessing import Pool
import time

#uclchem general takes a dictionary of parameters where outSpecies=number of outspecies
#and a string outspeciesin which is a delimited list of species
#this wrapper helps Pool to work and also converts a delimited list in the dictionary into a number
#and send list separately as required.
def run_uclchem(param_dict):
    outSpecies = (param_dict['outSpecies'])
    param_dict['outSpecies'] = len(outSpecies.split())
    
    #this will run uclchem to final time step and return an array of abundances
    #arrays cannot be variably sized with f2py so you need to set length of output array in src/wrap.f90
    abunds=uclchem.wrap.run_model_for_abundances(dictionary=param_dict, outspeciesin=outSpecies)
    
    #altenatively, you can run uclchem to file as you normally would
    #just put columnfile and/or output file in the dictionary
    #uclchem.wrap.run_model_to_file(dictionary, outSpeciesIn)
    return abunds


#basic set of parameters we'll use for this grid. 
ParameterDictionary = {"phase": 1, "switch": 0, "collapse": 0, "readAbunds": 0, "writeStep": 1,
                       "outSpecies": 'SO H3O+', "initialDens": 1e4,"finalTime":2.0e6,"baseAv":10}


# This part can be substituted with any choice of grid
# here we just combine various initial and final densities into an easily iterable array
initialTempGrid = np.linspace(50, 300, 6)
crgrid = np.logspace(1, 5, 5)
parameterSpace = np.asarray(np.meshgrid(initialTempGrid, crgrid)).reshape(2, -1)

print(parameterSpace.shape)

#then we loop through parameters, update parameter dictionary and run
#however, to use Pool we just store dictionaries for now
models=[]

#we'll number our models so store the parameters we used in a list
for k in range(np.shape(parameterSpace)[1]):
    paramDict=ParameterDictionary.copy()
    paramDict["initialTemp"] = parameterSpace[0,k]
    paramDict["zeta"] = parameterSpace[1,k]
    models.append(paramDict)

#use pool.map to run each dictionary throuh our helper function
start=time.time()
pool=Pool(10)
result=pool.map(run_uclchem,models)
result=np.asarray(result)
pool.close()
pool.join()
df=pd.DataFrame({"Temperature":parameterSpace[0,:],"CR":parameterSpace[1,:],
                "SO":result[:,0],"H3O+":result[:,1]})

df.to_csv("grid.csv",index=False)
end=time.time()
end=(end-start)/60.0
print(f"grid in {end} minutes")

