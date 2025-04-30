# enter two reaction files to see all reactions that are one file but not the other

from uclchem.makerates.io_functions import read_species_file, read_reaction_file


reactionsFile1 = "inputFiles/default_grain_network.csv"
# reactionsFile2 = "inputFiles/freeze_only_grain_network.csv"
reactionsFile2 = "inputFiles/excited_grain_network.csv"
speciesFile = "inputFiles/default_species.csv"

# differences are only relevant insofar as the missing reactions contain your species
speciesList, _ = read_species_file(speciesFile)
print(f"I found {len(speciesList)} species in the file {speciesFile}")

print("\nReading reactions")
reactions1, drops = read_reaction_file(reactionsFile1, speciesList, "UCL")
print(
    f"I found {len(reactions1)} reactions in the file {reactionsFile1}, I dropped {len(drops)} reactions."
)
# If you need to see which reactions are dropped:
# print("\n".join([str(drop) for drop in drops]))
reactions2, drops = read_reaction_file(reactionsFile2, speciesList, "UCL")
print(
    f"I found {len(reactions2)} reactions in the file {reactionsFile2}, I dropped {len(drops)} reactions."
)


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
            if reaction1.get_alpha() != reaction2.get_alpha():
                print(reaction1)
                print(f"alpha 1 = {reaction1.alpha}, alpha 2 = {reaction2.alpha}")
            if reaction1.get_beta() != reaction2.get_beta():
                print(reaction1)
                print(f"beta 1 = {reaction1.beta}, beta 2 = {reaction2.beta}")
            if reaction1.get_gamma() != reaction2.get_gamma():
                print(reaction1)
                print(f"gamma 1 = {reaction1.gamma}, gamma 2 = {reaction2.gamma}")
