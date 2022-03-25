#! /usr/bin/python
#
####################################################################################################
# 				MakeRates
# 		Current version by Jonathan Holdship & Antonios Makrymallis. Original by Tom Bell.
# 		MakeRates reads in lists of species and reactions and produces the network files needed
# 		by UCLCHEM to run. It also performs basic cleaning and sanity checks on the network.
#
####################################################################################################
from src.io_functions import *
from src.network import *
import os


###################################################################################################
# USER OPTIONS
###################################################################################################

database_reaction_file = "inputFiles/umist12-ucledit.csv"
database_reaction_type = "UMIST"
# database_reaction_file = "inputFiles/kida.uva.2014.dat"
# database_reaction_type="KIDA"

custom_reaction_file = "inputFiles/default_grain_network.csv"


species_file = "inputFiles/default_species.csv"

three_phase = True


#################################################################################################
if not os.path.exists("outputFiles"):
    os.makedirs("outputFiles")

print("\n################################################")
print("Reading and checking input")
print("################################################\n")

# Read user inputs
species_list = read_species_file(species_file)
reactions1, dropped_reactions = read_reaction_file(
    database_reaction_file, species_list, database_reaction_type
)
reactions2, dropped_reactions = read_reaction_file(custom_reaction_file, species_list, "UCL")

# Create Network
network = Network(species=species_list, reactions=reactions1 + reactions2, three_phase=three_phase)


# Print dropped reactions from grain file or write if many
if len(dropped_reactions) < 6:
    print("Reactions dropped from grain file:\n")
    for reaction in dropped_reactions:
        print(reaction)
else:
    print("\nReactions dropped from grain file written to outputFiles/dropped_reactions.csv\n")
    f = open("outputFiles/dropped_reactions.csv", "w")
    writer = csv.writer(
        f, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL, lineterminator="\n"
    )
    for reaction in dropped_reactions:
        writer.writerow(reaction)
    f.close()

# check network to see if there are potential problems
print("Checking Network")
network.check_network()


print("\n################################################")
print("Checks complete, writing output files")
print("################################################\n")

# Create the species file
print("\nWriting final species file...")
filename = "outputFiles/species.csv"
write_species(filename, network.species_list)
print("\tFinal Species File:", filename)

# Create the reaction file
print("Writing final reaction file...")
filename = "outputFiles/reactions.csv"
write_reactions(filename, network.reaction_list)
print("\tFinal Reaction File:", filename)

# Write the ODEs in the appropriate language format
print("Writing system of ODEs in F95 format...")
filename = "outputFiles/odes.f90"
write_odes_f90(filename, network.species_list, network.reaction_list, three_phase)
# filename = 'outputFiles/jacobian.f90'
# write_jacobian(filename,speciesList)
print("\tFinal ODE file:", filename)

print("Writing Network File...")
filename = ""
write_network_file(network)
print("\tFinal Network file:", filename)

ngrain = len([x for x in species_list if x.is_surface_species()])


print(f"Total number of species = {len(network.species_list)}")
print(f"Number of surface species = {ngrain}")
print(f"Number of reactions = {len(network.reaction_list)}")
