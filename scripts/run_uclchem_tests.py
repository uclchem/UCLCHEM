#Jon Holdship 23/11/2020
#Run the 3 test models and store their outputs for the plotting script
from sys import path
path.insert(0,'./')
path.insert(0,'../')
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
    abunds=uclchem.wrap.run_model_to_file(dictionary=param_dict, outspeciesin=outSpecies)
    
    #altenatively, you can run uclchem to file as you normally would
    #just put columnfile and/or output file in the dictionary
    #uclchem.wrap.run_model_to_file(dictionary, outSpeciesIn)
    return abunds


#set a parameter dictionary for phase 1 collapse model
params = {"phase": 1, "switch": 1, "collapse": 1, "readAbunds": 0, "writeStep": 1,
               "outSpecies": 'SO CO', "initialDens": 1e2, "initialTemp":10.0,
               "finalDens":1e5, "finalTime":5.0e6,
               "outputFile":"examples/test-output/phase1-full.dat",
               "abundFile":"examples/test-output/startcollapse.dat"}

run_uclchem(params.copy())

#change to static cloud params
params["collapse"]=0
params["switch"]=0
params["initialDens"]=1e4
params["abundFile"]="examples/test-output/startstatic.dat"

params["outputFile"]="examples/test-output/static-full.dat"
run_uclchem(params.copy())

#finally, run phase 2 from the phase 1 model.
params["initialDens"]=1e5
params["tempindx"]=3
params["fr"]=0.0
params["readAbunds"]=1
params["phase"]=2
params["finalTime"]=1e6
params["abundFile"]="examples/test-output/startcollapse.dat"
params["outputFile"]="examples/test-output/phase2-full.dat"
run_uclchem(params.copy())
