#! /usr/bin/python
#
####################################################################################################
# 				MakeRates
# 		Current version by Jonathan Holdship & Antonios Makrymallis. Original by Tom Bell.
# 		MakeRates reads in lists of species and reactions and produces the network files needed
# 		by UCLCHEM to run. It also performs basic cleaning and sanity checks on the network.
#
####################################################################################################
import src.io_functions as io 
from src.network import *
import os
import yaml

param_list = [
    "species_file",
    "database_reaction_file",
    "database_reaction_type",
    "custom_reaction_file",
    "custom_reaction_type",
    "three_phase",
]

with open("user_settings.yaml", "r") as f:
    user_params = yaml.safe_load(f)
for param in param_list:
    try:
        print(f"{param} : {user_params[param]}")
    except:
        raise KeyError(f"{param} not found in user_settings.yaml")

try:
    user_output_dir = user_params["output_directory"]
    if not os.path.exists(user_output_dir):
        os.makedirs(user_output_dir)
except:
    user_output_dir = None

#################################################################################################


print("\n################################################")
print("Reading and checking input")
print("################################################\n")

# Read user inputs
species_list, user_defined_bulk = io.read_species_file(user_params["species_file"])
reactions1, dropped_reactions = io.read_reaction_file(
    user_params["database_reaction_file"], species_list, user_params["database_reaction_type"]
)
reactions2, dropped_reactions = io.read_reaction_file(user_params["custom_reaction_file"], species_list, user_params["custom_reaction_type"])

# Create Network
network = Network(species=species_list, reactions=reactions1 + reactions2, three_phase=user_params["three_phase"], user_defined_bulk=user_defined_bulk)

io.output_drops(dropped_reactions, user_output_dir)

# check network to see if there are potential problems
print("Checking Network")
network.check_network()


print("\n################################################")
print("Checks complete, writing output files")
print("################################################\n")

io.write_outputs(network, user_output_dir)

io.pickle_network(network, "test.pl")

# debug:

print(network.__dict__)

# end debug



ngrain = len([x for x in species_list if x.is_surface_species()])


print(f"Total number of species = {len(network.species_list)}")
print(f"Number of surface species = {ngrain}")
print(f"Number of reactions = {len(network.reaction_list)}")
