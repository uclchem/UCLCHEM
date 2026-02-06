"""UCLCHEM Utility Functions

Helper functions and utilities for UCLCHEM operations.

This module provides utility functions for:
- Error handling and reporting
- Physics calculations (shock dissipation times)
- Parameter validation
- File path management

**Key Functions:**

- :func:`check_error` - Convert UCLCHEM error codes to messages
- :func:`cshock_dissipation_time` - Calculate C-shock dissipation timescale

**Example Usage:**

.. code-block:: python

    import uclchem.utils as utils

    # Check error from model run
    success_flag = cloud.success_flag
    if success_flag < 0:
        error_msg = utils.check_error(success_flag)
        print(f"Model failed: {error_msg}")

    # Calculate shock timescale
    t_diss = utils.cshock_dissipation_time(
        shock_vel=50.0,  # km/s
        initial_dens=1e4  # cm^-3
    )
    print(f"Dissipation time: {t_diss:.1e} years")

**Error Codes:**

UCLCHEM model functions return negative integer error codes on failure:

- ``-1``: Parameter read failed (misspelled parameter)
- ``-2``: Physics initialization failed (invalid parameters)
- ``-3``: Chemistry initialization failed
- ``-4``: Integrator error (DVODE failed)

Use :func:`check_error` to get human-readable error messages.

**See Also:**

- :mod:`uclchem.model` - Model classes that use these utilities
"""

from os import path

import pandas as pd

_ROOT = path.dirname(path.abspath(__file__))


def cshock_dissipation_time(shock_vel: float, initial_dens: float) -> float:
    """A simple function used to calculate the dissipation time of a C-type shock.
    Use to obtain a useful timescale for your C-shock model runs. Velocity of
    ions and neutrals equalizes at dissipation time and full cooling takes a few dissipation times.

    Args:
        shock_vel (float): Velocity of the shock in km/s
        initial_dens (float): Preshock density of the gas in cm$^{-3}$

    Returns:
        float: The dissipation time of the shock in years
    """
    pc = 3.086e18  # parsec in cgs
    SECONDS_PER_YEAR = 3.15569e7
    dlength = 12.0 * pc * shock_vel / initial_dens
    return (dlength * 1.0e-5 / shock_vel) / SECONDS_PER_YEAR


def check_error(error_code: int) -> str:
    """Converts the UCLCHEM integer result flag to a simple messaging explaining what went wrong"

    Args:
        error_code (int): Error code returned by UCLCHEM models, the first element of the results list.

    Returns:
        str: Error message
    """
    errors = {
        -1: "Parameter read failed. Likely due to a mispelled parameter name, compare your dictionary to the parameters docs.",
        -2: "Physics intiialization failed. Often due to user chosing unacceptable parameters such as hot core masses or collapse modes that don't exist. Check the docs for your model function.",
        -3: "Chemistry initialization failed",  # this doesn't exist yet
        -4: "Unrecoverable integrator error, DVODE failed to integrate the ODEs in a way that UCLCHEM could not fix. Run UCLCHEM tests to check your network works at all then try to see if bad parameter combination is at play.",
        -5: "Too many integrator fails. DVODE failed to integrate the ODE and UCLCHEM repeatedly altered settings to try to make it pass but tried too many times without success so code aborted to stop infinite loop.",
        -6: "The model was stopped because there are not enough time points allocated in the time array. Increase the number of time points in the time array in constants.py and try again.",
    }
    try:
        return errors[error_code]
    except KeyError:
        raise ValueError(f"Unknown error code: {error_code}")


def get_species_table() -> pd.DataFrame:
    """A simple function to load the list of species in the UCLCHEM network into a pandas dataframe.

    Returns:
        pandas.DataFrame: A dataframe containing the species names and their details
    """

    species_list = pd.read_csv(path.join(_ROOT, "species.csv"))
    return species_list


def get_species() -> list[str]:
    """A simple function to load the list of species present in the UCLCHEM network

    Returns:
        list[str] : A list of species names
    """

    species_list = pd.read_csv(path.join(_ROOT, "species.csv")).iloc[:, 0].tolist()
    return species_list


def get_reaction_table() -> pd.DataFrame:
    """A function to load the reaction table from the UCLCHEM network into a pandas dataframe.

    Returns:
        pandas.DataFrame: A dataframe containing the reactions and their rates
    """

    reactions = pd.read_csv(path.join(_ROOT, "reactions.csv"))
    return reactions


def find_number_of_consecutive_digits(string: str, start: int) -> int:
    """Determine the number of consecutive digits in a string.

    Args:
        string (str): the string
        start (int): the starting index

    Returns:
        num_digits (int): the number of consecutive digits in the string
            starting from "start".

    Examples:
        >> find_number_of_consecutive_digits("Hello123", 0) -> 0,
        >> find_number_of_consecutive_digits("Hello123", 5) -> 3,
        >> find_number_of_consecutive_digits("Hello123", 6) -> 2,
        >> find_number_of_consecutive_digits("He1llo23", 2) -> 1,

    """
    num_digits = 0
    while start + num_digits < len(string) and string[start + num_digits].isdigit():
        num_digits += 1
    return num_digits
