from os import path

import pandas as pd

_ROOT = path.dirname(path.abspath(__file__))


def cshock_dissipation_time(shock_vel, initial_dens):
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


def check_error(error_code):
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


def get_species_table():
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


def get_reaction_table():
    """A function to load the reaction table from the UCLCHEM network into a pandas dataframe.

    Returns:
        pandas.DataFrame: A dataframe containing the reactions and their rates
    """

    reactions = pd.read_csv(path.join(_ROOT, "reactions.csv"))
    return reactions
