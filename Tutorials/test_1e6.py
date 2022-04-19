from regex import P
import uclchem

# set a parameter dictionary for phase 1 collapse model
param_dict = {
    "switch": 0,#stop at finalTime
    "freefall": 1,#increase density in freefall
    "initialDens": 1e2, #starting density
    "finalDens":1e6, #final density
    "initialTemp": 10.0,#temperature of gas
    "finalTime": 6.0e6, #final time
    "rout":0.1, #radius of cloud in pc
    "baseAv":1.0, #visual extinction at cloud edge.
    "abundSaveFile": "../examples/test-output/startcollapse.dat",#save final abundances to file
}
#change other bits of input to set up phase 2
param_dict["initialDens"]=1e6
param_dict["finalTime"]=1e6
param_dict["switch"]=0
param_dict["freefall"]=0
param_dict["abundLoadFile"]=param_dict.pop("abundSaveFile") #this is still set to startcollapse.dat from phase 1 so remove it or change it.
param_dict["outputFile"]="../examples/test-output/phase2-full.dat"

result=uclchem.model.hot_core(temp_indx=3,max_temperature=300.0,param_dict=param_dict)