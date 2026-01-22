"""Generate parameter documentation from f2py uclchemwrap module.

This script introspects the compiled uclchemwrap.defaultparameters module
to generate comprehensive parameter documentation. It replaces the old
Fortran source parsing approach with direct inspection of the f2py-generated
wrapper, ensuring documentation stays in sync with the actual compiled code.

Usage:
    python generate_param_docs_from_wrap.py <output_markdown_file>
"""

import sys
from typing import Any, Dict, Tuple

import numpy as np


def get_parameter_info() -> Dict[str, Tuple[Any, str, str]]:
    """Extract all parameters from uclchemwrap.defaultparameters.

    Returns:
        Dictionary mapping parameter names to tuples of (value, type_str, description)
        The description is extracted from Fortran comments when available.
    """
    try:
        import uclchemwrap

        dp = uclchemwrap.defaultparameters
    except ImportError:
        print("ERROR: Cannot import uclchemwrap")
        print("Make sure UCLCHEM is installed with: pip install -e .")
        sys.exit(1)
    except AttributeError:
        print("ERROR: uclchemwrap.defaultparameters module not found")
        print("Make sure UCLCHEM is properly compiled.")
        sys.exit(1)

    # Parameter descriptions from Fortran source comments
    # These are maintained here as f2py doesn't preserve Fortran comments
    descriptions = {
        # Physical Variables
        "initialtemp": "Initial gas temperature in Kelvin for all gas parcels in model.",
        "initialdens": "Initial gas density in H nuclei per cm$^{-3}$ for all gas parcels in model.",
        "finaldens": "Final gas density achieved via freefall.",
        "currenttime": "Time at start of model in years.",
        "finaltime": "Time to stop model in years, if not using `endAtFinalDensity`.",
        "radfield": "Interstellar radiation field in Habing units.",
        "zeta": "Cosmic ray ionisation rate as multiple of $1.3 \\times 10^{-17}$ s$^{-1}$.",
        "rout": "Outer radius of cloud being modelled in pc.",
        "rin": "Minimum radial distance from cloud centre to consider.",
        "baseav": "Extinction at cloud edge, Av of a parcel at rout.",
        "points": "Number of gas parcels equally spaced between rin to rout to consider.",
        "bm0": "Magnetic parameter [microgauss]: B0 = bm0*sqrt(initialDens).",
        # Behavioural Controls
        "freezefactor": "Modify freeze out rate of gas parcels by this factor.",
        "endatfinaldensity": "Choose to end model at final density, otherwise end at final time.",
        "freefall": "Controls whether model density increases following freefall equation.",
        "freefallfactor": "Modify freefall rate by factor, usually to slow it.",
        "desorb": "Toggles all non-thermal desorption processes on or off.",
        "h2desorb": "Individually toggle non-thermal desorption due to H2 formation.",
        "crdesorb": "Individually toggle non-thermal desorption due to cosmic rays.",
        "uvdesorb": "Individually toggle non-thermal desorption due to UV photons.",
        "thermdesorb": "Toggle continuous thermal desorption.",
        "instantsublimation": "Toggle instantaneous sublimation of the ices at t=0.",
        "cosmicrayattenuation": "Use column density to attenuate cosmic ray ionisation rate following Padovani et al. 2018.",
        "ionmodel": "L/H model for cosmic ray attenuation (Padovani et al. 2018).",
        "improvedh2crpdissociation": "Use H2 CRP dissociation rate from Padovani et al. 2018b.",
        "heatingflag": "If True, heating is applied to the gas parcels.",
        "enforcechargeconservation": "Enforce charge conservation by tracking charged ions.",
        # Input and Output
        "outputfile": "File to write full output of UCLCHEM (physical parameters and all abundances at every time step).",
        "columnfile": "File to write specific species abundances (see outSpecies).",
        "rateconstantfile": "File to write rate 'constants' at each timestep.",
        "ratesfile": "File to write reaction rates (flux) at each timestep.",
        "heatingfile": "File to write heating and cooling rates at each timestep.",
        "writestep": "Writing to columnFile only happens every writeStep timesteps.",
        "abundsavefile": "File to save final abundances at end of model for use in future models.",
        "abundloadfile": "File to load initial abundances from (created through abundSaveFile).",
        "coolantdatadir": "Directory where the collisional rate files are stored.",
        # Initial Abundances
        "metallicity": "Scale the abundances of all elements heavier than He by this factor.",
        "ion": "Sets how much elemental C is initially atomic (0=all atomic/1=50:50/2=fully ionized).",
        "fh": "Fraction of H initially in atomic H (rest goes to H2). Total H abundance is always 1.",
        "fhe": "Total elemental abundance of He.",
        "fc": "Total elemental abundance of C.",
        "fo": "Total elemental abundance of O.",
        "fn": "Total elemental abundance of N.",
        "fs": "Total elemental abundance of S.",
        "fmg": "Total elemental abundance of Mg.",
        "fsi": "Total elemental abundance of Si.",
        "fcl": "Total elemental abundance of Cl.",
        "fp": "Total elemental abundance of P.",
        "ffe": "Total elemental abundance of Fe.",
        "ff": "Total elemental abundance of F.",
        "fd": "Total elemental abundance of D.",
        "fli": "Total elemental abundance of Li.",
        "fna": "Total elemental abundance of Na.",
        "fpah": "Total initial abundance of PAHs.",
        "f15n": "Total initial abundance of $^{15}$N.",
        "f13c": "Total initial abundance of $^{13}$C.",
        "f18o": "Total initial abundance of $^{18}$O.",
        # Integration Controls
        "reltol": "Relative tolerance for ODE integration (see integration troubleshooting).",
        "abstol_factor": "Absolute tolerance is calculated by multiplying species abundance by this factor.",
        "abstol_min": "Minimum value absolute tolerances can take.",
        "mxstep": "Maximum steps allowed in integration before warning is thrown.",
        # Advanced Parameters
        "ebmaxh2": "Maximum binding energy of species desorbed by H2 formation.",
        "ebmaxcr": "Maximum binding energy of species desorbed by cosmic ray ionisation.",
        "ebmaxuvcr": "Maximum binding energy of species desorbed by UV photons.",
        "epsilon": "Number of molecules desorbed per H2 formation.",
        "uv_yield": "Number of molecules desorbed per UV photon (extrapolated from Oberg et al. 2009).",
        "phi": "Number of molecules desorbed per cosmic ray ionisation.",
        "uvcreff": "Ratio of CR induced UV photons to ISRF UV photons.",
        "omega": "Dust grain albedo.",
    }

    params = {}

    # Get all attributes from defaultparameters module
    for name in dir(dp):
        if name.startswith("_"):
            continue

        try:
            value = getattr(dp, name)

            # Skip functions and modules
            if callable(value):
                continue

            # Get type information
            if isinstance(value, np.ndarray):
                type_str = f"array[{value.dtype}]"
            elif isinstance(value, (bool, np.bool_)):
                type_str = "bool"
            elif isinstance(value, (int, np.integer)):
                type_str = "int"
            elif isinstance(value, (float, np.floating)):
                type_str = "float"
            elif isinstance(value, str):
                type_str = "str"
            elif isinstance(value, bytes):
                type_str = "str"
                value = value.decode("utf-8").strip()
            else:
                type_str = type(value).__name__

            # Get description
            description = descriptions.get(name.lower(), "")

            params[name.lower()] = (value, type_str, description)

        except Exception as e:
            print(f"Warning: Could not process parameter '{name}': {e}")
            continue

    return params


def categorize_parameters(params: Dict[str, Tuple[Any, str, str]]) -> Dict[str, list]:
    """Organize parameters into logical categories.

    Args:
        params: Dictionary of parameter info

    Returns:
        Dictionary mapping category names to lists of parameter names
    """
    categories = {
        "Physical Variables": [
            "initialtemp",
            "initialdens",
            "finaldens",
            "currenttime",
            "finaltime",
            "radfield",
            "zeta",
            "rout",
            "rin",
            "baseav",
            "points",
            "bm0",
        ],
        "Behavioural Controls": [
            "freezefactor",
            "endatfinaldensity",
            "freefall",
            "freefallfactor",
            "desorb",
            "h2desorb",
            "crdesorb",
            "uvdesorb",
            "thermdesorb",
            "instantsublimation",
            "cosmicrayattenuation",
            "ionmodel",
            "improvedh2crpdissociation",
            "heatingflag",
            "enforcechargeconservation",
        ],
        "Input and Output": [
            "outputfile",
            "columnfile",
            "rateconstantfile",
            "ratesfile",
            "heatingfile",
            "writestep",
            "abundsavefile",
            "abundloadfile",
            "coolantdatadir",
        ],
        "Initial Abundances": [
            "metallicity",
            "ion",
            "fh",
            "fhe",
            "fc",
            "fo",
            "fn",
            "fs",
            "fmg",
            "fsi",
            "fcl",
            "fp",
            "ffe",
            "ff",
            "fd",
            "fli",
            "fna",
            "fpah",
            "f15n",
            "f13c",
            "f18o",
        ],
        "Integration Controls": ["reltol", "abstol_factor", "abstol_min", "mxstep"],
        "Advanced Parameters": [
            "ebmaxh2",
            "ebmaxcr",
            "ebmaxuvcr",
            "epsilon",
            "uv_yield",
            "phi",
            "uvcreff",
            "omega",
        ],
    }

    return categories


def format_value(value: Any) -> str:
    """Format a parameter value for display."""
    if isinstance(value, (bool, np.bool_)):
        return ".True." if value else ".False."
    elif isinstance(value, (int, np.integer)):
        # Check if it's a boolean disguised as int (Fortran LOGICAL)
        if value in (0, 1):
            return ".True." if value == 1 else ".False."
        return str(value)
    elif isinstance(value, (float, np.floating)):
        # Use scientific notation for very small/large numbers
        if abs(value) < 0.001 or abs(value) > 10000:
            return f"{value:.2e}"
        else:
            return f"{value:.3g}"
    elif isinstance(value, (bytes, np.bytes_)):
        # Handle Fortran strings (bytes)
        decoded = value.decode("utf-8").strip()
        return '""' if not decoded else f'"{decoded}"'
    elif isinstance(value, str):
        return f'"{value}"' if value else '""'
    else:
        return str(value)


def generate_markdown(params: Dict[str, Tuple[Any, str, str]], output_file: str):
    """Generate markdown documentation from parameter information.

    Args:
        params: Dictionary of parameter info
        output_file: Path to output markdown file
    """
    categories = categorize_parameters(params)

    with open(output_file, "w") as f:
        f.write("# UCLCHEM Parameters\n\n")
        f.write(
            "*Auto-generated from compiled uclchemwrap.defaultparameters module*\n\n"
        )
        f.write(
            "UCLCHEM will default to these values unless they are overridden by the user. "
        )
        f.write(
            "Users can override these by adding the variable name to the `param_dict` "
        )
        f.write(
            "argument of any UCLCHEM model function. Parameter names are case-insensitive.\n\n"
        )

        for category_name, param_names in categories.items():
            f.write(f"## {category_name}\n\n")

            # Check if we actually have params in this category
            available_params = [
                (name, params[name]) for name in param_names if name in params
            ]

            if not available_params:
                f.write("*No parameters in this category*\n\n")
                continue

            # Write table header
            f.write("| Parameter | Default Value | Description |\n")
            f.write("| --------- | ------------- | ----------- |\n")

            # Write parameters
            for name, (value, type_str, description) in available_params:
                value_str = format_value(value)
                f.write(f"| `{name}` | {value_str} | {description} |\n")

            f.write("\n")

        # Add footer note
        f.write("## Notes\n\n")
        f.write("- Parameter names are **case-insensitive** in `param_dict`\n")
        f.write("- Elemental abundances use the heavily depleted case from ")
        f.write(
            "[Jenkins et al. 2009](https://ui.adsabs.harvard.edu/abs/2009ApJ...700.1299J)\n"
        )
        f.write(
            "- Total H abundance is always 1.0 by definition (abundances are relative to H nuclei)\n"
        )
        f.write(
            '- Empty string values (`""`) indicate parameters that should be set by the user\n\n'
        )


def main():
    """Main entry point."""
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_markdown_file>")
        sys.exit(1)

    output_file = sys.argv[1]

    print("Extracting parameters from uclchemwrap.defaultparameters...")
    params = get_parameter_info()
    print(f"Found {len(params)} parameters")

    print(f"Generating markdown documentation to {output_file}...")
    generate_markdown(params, output_file)

    print("Done!")


if __name__ == "__main__":
    main()
