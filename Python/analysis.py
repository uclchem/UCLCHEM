import uclchem
import numpy as np
import matplotlib.pyplot as plt

################################################
# User Inputs Go Here
################################################

speciesName = "CH3OH"
result_file = "examples/test-output/phase2-full.dat"
reaction_file = "src/reactions.csv"
species_file = "src/species.csv"
output="analysis.dat"

################################################

result_df = uclchem.read_output_file(result_file)
species = np.loadtxt(
    species_file, usecols=[0], dtype=str, skiprows=0, unpack=True, delimiter=",", comments="%"
)
reactions = np.loadtxt(
    reaction_file, dtype=str, skiprows=0, delimiter=",", usecols=[0, 1, 2, 3, 4, 5, 6], comments="%"
)

fortran_reac_indxs = [i + 1 for i, reaction in enumerate(reactions) if speciesName in reaction]
reac_indxs = [i for i, reaction in enumerate(reactions) if speciesName in reaction]

old_key_reactions = []
with open(output, "w") as f:
    for i, row in result_df.iterrows():

        param_dict = uclchem.param_dict_from_output(row)
        rates = uclchem.get_species_rates(param_dict, row[species], fortran_reac_indxs)
        change_reacs, changes = uclchem.get_rates_of_change(
            rates, reactions[reac_indxs], species, speciesName, row
        )


        change_reacs = uclchem.format_reactions(change_reacs)

        total_formation, total_destruct, key_reactions, key_changes = uclchem.remove_slow_reactions(
            changes, change_reacs
        )

        if old_key_reactions != key_reactions:
            old_key_reactions = key_reactions[:]
            uclchem.write_analysis(
                f,row["Time"], total_formation, total_destruct, key_reactions, key_changes
            )
