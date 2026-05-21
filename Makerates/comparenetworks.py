"""Compare two reaction files.

Helper functions to enter two reaction files to see
all reactions that are in one file but not the other

"""

from uclchem.makerates.io_functions import read_reaction_file, read_species_file

if __name__ == "__main__":
    reactions_file1 = "inputFiles/default_grain_network.csv"
    reactions_file2 = "inputFiles/excited_grain_network.csv"
    species_file = "inputFiles/default_species.csv"

    # differences are only relevant insofar as the missing reactions contain your species
    species_list, _ = read_species_file(species_file)
    print(f"I found {len(species_list)} species in the file {species_file}")

    print("\nReading reactions")
    reactions1, drops = read_reaction_file(reactions_file1, species_list, "UCL")
    print(
        f"I found {len(reactions1)} reactions in the file {reactions_file1}, I dropped {len(drops)} reactions."
    )
    # If you need to see which reactions are dropped:
    # print("\n".join([str(drop) for drop in drops])) # noqa: ERA001
    reactions2, drops = read_reaction_file(reactions_file2, species_list, "UCL")
    print(
        f"I found {len(reactions2)} reactions in the file {reactions_file2}, I dropped {len(drops)} reactions."
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
        for reaction2 in reactions2:
            if reaction1 == reaction2:
                if reaction1.get_alpha() != reaction2.get_alpha():
                    print(reaction1)
                    print(
                        f"alpha 1 = {reaction1.get_alpha()}, alpha 2 = {reaction2.get_alpha()}"
                    )
                if reaction1.get_beta() != reaction2.get_beta():
                    print(reaction1)
                    print(
                        f"beta 1 = {reaction1.get_beta()}, beta 2 = {reaction2.get_beta()}"
                    )
                if reaction1.get_gamma() != reaction2.get_gamma():
                    print(reaction1)
                    print(
                        f"gamma 1 = {reaction1.get_gamma()}, gamma 2 = {reaction2.get_gamma()}"
                    )
