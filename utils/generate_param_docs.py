"""Simple script to parse defaultparameters.f90 and generate a markdown
file that can be copied over to the website docs. Run it with
python generate_param_docs.py default_parameter_fortran_file output_markdown_file
"""

import ast
from sys import argv

param_file = argv[1]
output_file = argv[2]
constants_file = argv[3]
default_param_dictionary = {}

with open(param_file, "r") as f:
    with open(output_file, "w") as output:
        for line in f.readlines():
            if line.startswith("!"):
                if not line.startswith("!!"):
                    output.write(line.replace("!", ""))
            elif line.startswith(
                ("USE", "MODULE", "IMPLICIT", "END"),
            ):
                # Do not read the placeholder function that was introduced for f2py
                continue
            elif line.startswith("CONTAINS"):
                break
            else:
                if "=" in line:
                    new_line = line.split("=")
                    temp_line = new_line[0].split("::")
                    key = temp_line[-1].strip()
                    type_of_value = temp_line[0].strip()
                    new_line = new_line[1].split("!")
                    value = new_line[0]
                    description = new_line[1]
                    line = (
                        "|"
                        + key
                        + "|"
                        + new_line[0]
                        + "|"
                        + new_line[1].strip()
                        + "|\n"
                    )
                    output.write(line)
                    if "REAL" in type_of_value:
                        default_param_dictionary[key.lower()] = float(
                            value.replace("d", "e").strip()
                        )
                    elif "LOGICAL" in type_of_value:
                        default_param_dictionary[key.lower()] = bool(value[1:-1])
                    elif "CHARACTER" in type_of_value:
                        if '"' in value:
                            value = value[value.find('"') + 1 : value.rfind('"')]
                        elif "'" in value:
                            value = value[value.find("'") + 1 : value.rfind("'")]
                        if value == "":
                            value = None

                        default_param_dictionary[key.lower()] = value
                    elif "INTEGER" in type_of_value:
                        default_param_dictionary[key.lower()] = int(value)

    # Read constants and potentially modify them to update default_para_dictionary
    with open(constants_file, "r") as constants_r:
        lines = constants_r.readlines()

    default_dict_found = False
    for i, line in enumerate(lines):
        if default_dict_found:
            break
        elif "default_param_dictionary" in line:
            default_dict_found = True
            if ast.literal_eval(line.split("=")[-1]) != default_param_dictionary:
                lines[i] = f"default_param_dictionary={default_param_dictionary}\n"
            break

    if not default_dict_found:
        lines += [f"default_param_dictionary = {default_param_dictionary}\n"]

    with open(constants_file, "w") as constants_w:
        constants_w.writelines(lines)
