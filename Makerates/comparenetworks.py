# enter two reaction files to see all reactions that are one file but not the other

from uclchem.makerates.io_functions import read_species_file, read_reaction_file


reactions1 = "inputFiles/default_grain_network.csv"
reactions2 = "inputFiles/freeze_only_grain_network.csv"
speciesFile = "inputFiles/default_species.csv"

# differences are only relevant insofar as the missing reactions contain your species
speciesList = read_species_file(speciesFile)

print("\nReading reactions")
reactions2, drops = read_reaction_file(reactions2, speciesList, "UCL")
reactions1, drops = read_reaction_file(reactions1, speciesList, "UCL")

print("\nReactions from file 1 not in file 2")
for reaction1 in reactions1:
    match = False
    for reaction2 in reactions2:
        if reaction2 == reaction1:
            match = True
            break
    if not match:
        print(reaction1)

print("Reactions from file 2 not in file 1")
for reaction1 in reactions2:
    match = False
    for reaction2 in reactions1:
        if reaction2 == reaction1:
            match = True
            break
    if not match:
        print(reaction1)

print("Reactions with different coefficients")
for reaction1 in reactions1:
    for reation2 in reactions2:
        if reaction1 == reaction2:
            if reaction1.alpha != reaction2.alpha:
                print(reaction1)
                print(f"alpha 1 = {reaction1.alpha}, alpha 2 = {reaction2.alpha}")
            if reaction1.beta != reaction2.beta:
                print(reaction1)
                print(f"beta 1 = {reaction1.beta}, beta 2 = {reaction2.beta}")
            if reaction1.gamma != reaction2.gamma:
                print(reaction1)
                print(f"gamma 1 = {reaction1.gamma}, gamma 2 = {reaction2.gamma}")
