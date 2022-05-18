
import uclchem
# set a parameter dictionary for static model
outSpecies = ["OCS", "CO", "CS", "CH3OH"]
params = {
    "endAtFinalDensity": False,
    "freefall": False,
    "writeStep": 1,
    "initialDens": 1e4,
    "initialTemp": 10.0,
    "finalDens": 1e5,
    "finalTime": 5.0e6,
    "outputFile": "test-output/static-full.dat",
    "abundSaveFile": "test-output/startstatic.dat",
}

uclchem.model.cloud(param_dict=params, out_species=outSpecies)