"""Analyze the formation and destruction reactions of a species.

Deprecated.
"""

import uclchem

if __name__ == "__main__":
    ################################################
    # User Inputs Go Here
    ################################################

    species_name = "#CO"
    result_file = "examples/test-output/phase1-full.dat"
    output = "analysis.dat"

    ################################################
    uclchem.analysis.analysis(species_name, result_file, output)
