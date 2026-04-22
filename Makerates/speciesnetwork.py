"""Simple python script that reads a Makerates produced reaction file
and prints every reaction that forms or destroys a chosen species.
"""

import csv
import pathlib

if __name__ == "__main__":
    species = "HCN"

    with pathlib.Path("../src/uclchem/reactions.csv").open() as f:
        reader = csv.reader(f, delimiter=",", quotechar="|")

        forms = "Species Formed in: \n"
        dest = "Species Destroyed in: \n"
        for row in reader:
            if len(row) > 1:
                if species in {row[0], row[1], row[2]}:
                    dest = (
                        dest
                        + row[0]
                        + " + "
                        + row[1]
                        + " --> "
                        + row[3]
                        + " + "
                        + row[4]
                        + " + "
                        + row[5]
                        + "\n"
                    )
                elif species in {row[3], row[4], row[5]}:
                    forms = (
                        forms
                        + row[0]
                        + " + "
                        + row[1]
                        + " --> "
                        + row[3]
                        + " + "
                        + row[4]
                        + " + "
                        + row[5]
                        + "\n"
                    )

        print(forms)
        print(dest)
