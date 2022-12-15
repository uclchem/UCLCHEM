
import uclchem.makerates.io_functions as io 
from uclchem.makerates import Network
import os
import yaml
import logging
from logging import Logger



param_list = [
    "species_file",
    "database_reaction_file",
    "database_reaction_type",
    "custom_reaction_file",
    "custom_reaction_type",
    "three_phase",
]


def run_makerates(configuration_file="user_settings.yaml", verbosity="INFO", write_files=True):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=verbosity)
    
    with open(configuration_file, "r") as f:
        user_params = yaml.safe_load(f)
    
    for param in param_list:
        try:
            # Check if we can access the needed parameter:
            user_params[param]
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

    print("\n################################################\n"+
                "Reading and checking input\n" + 
                "################################################\n")

    # Read user inputs
    species_list, user_defined_bulk = io.read_species_file(user_params["species_file"])
    reactions1, dropped_reactions1 = io.read_reaction_file(
        user_params["database_reaction_file"], species_list, user_params["database_reaction_type"]
    )
    reactions2, dropped_reactions2 = io.read_reaction_file(user_params["custom_reaction_file"], species_list, user_params["custom_reaction_type"])

    # Create Network
    network = Network(species=species_list, reactions=reactions1 + reactions2, three_phase=user_params["three_phase"], user_defined_bulk=user_defined_bulk)

    io.output_drops(dropped_reactions1 + dropped_reactions2, user_output_dir)

    # check network to see if there are potential problems
    print("Checking Network")
    network.check_network()

    print("\n################################################\n"+
                "Checks complete, writing output files\n"+
                "################################################\n")
    if write_files:
        io.write_outputs(network, user_output_dir)        

    ngrain = len([x for x in species_list if x.is_surface_species()])
    print(f"Total number of species = {len(network.species_list)}")
    print(f"Number of surface species = {ngrain}")
    print(f"Number of reactions = {len(network.reaction_list)}")
    # return the network such that the object can be reused in code/notebooks
    return network
