import uclchem
import numpy as np
import matplotlib.pyplot as plt

################################################
# User Inputs Go Here
################################################

species_name = "CH3OH"
result_file = "examples/test-output/phase2-full.dat"
reaction_file = "src/reactions.csv"
species_file = "src/species.csv"
output="analysis.dat"

################################################
uclchem.analysis.analysis(species_name, result_file, output)