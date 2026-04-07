"""Compare two species files.

Helper functions to enter two species files to see
all species that are in one file but not the other
"""

from uclchem.makerates.io_functions import read_species_file

if __name__ == "__main__":
    species_file = "inputFiles/default_species.csv"
    species_file2 = "inputFiles/excited_species.csv"

    species_list, _ = read_species_file(species_file)
    species_list2, _ = read_species_file(species_file2)

    print(f"File {species_file} contains {len(species_list)} species.")
    print(f"File {species_file2} contains {len(species_list2)} species.")

    print(f"Species in {species_file2} but not {species_file}")
    missing = [x for x in species_list2 if x not in species_list]
    print(len(missing))
    for species in missing:
        print(species.name)

    print(f"Species in {species_file} but not {species_file2}")
    missing = [x for x in species_list if x not in species_list2]
    print(len(missing))

    for species in missing:
        print(species.name)
