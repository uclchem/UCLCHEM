# enter two species files to see all species that are one file but not the other

from uclchem.makerates.io_functions import read_species_file


speciesFile = "inputFiles/default_species.csv"
speciesFile2 = "inputFiles/excited_species.csv"

species_list, _ = read_species_file(speciesFile)
species_list2, _ = read_species_file(speciesFile2)

print(f"File {speciesFile} contains {len(species_list)} species.")
print(f"File {speciesFile2} contains {len(species_list2)} species.")

print(f"Species in {speciesFile2} but not {speciesFile}")
missing = [x for x in species_list2 if x not in species_list]
print(len(missing))
for species in missing:
    print(species.name)

print(f"Species in {speciesFile} but not {speciesFile2}")
missing = [x for x in species_list if x not in species_list2]
print(len(missing))

for species in missing:
    print(species.name)
