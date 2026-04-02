"""Helper functions and utilities for UCLCHEM operations.

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

import enum
from pathlib import Path
from typing import Self

import pandas as pd

UCLCHEM_ROOT_DIR: Path = Path(__file__).parent.resolve().absolute()


def cshock_dissipation_time(shock_vel: float, initial_dens: float) -> float:
    """Calculate the dissipation time of a C-type shock.
    Use to obtain a useful timescale for your C-shock model runs.
    Velocity of ions and neutrals equalizes at dissipation time and
    full cooling takes a few dissipation times.

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


def get_species_table() -> pd.DataFrame:
    """Load the list of species in the UCLCHEM network into a pandas dataframe.

    Returns:
        pandas.DataFrame: A dataframe containing the species names and their details

    """
    species_list = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv")
    return species_list


def get_species() -> list[str]:
    """Load the list of species present in the UCLCHEM network.

    Returns:
        list[str] : A list of species names

    """
    species_list = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv").iloc[:, 0].tolist()
    return species_list


def get_reaction_table() -> pd.DataFrame:
    """Load the reaction table from the UCLCHEM network into a pandas dataframe.

    Returns:
        pandas.DataFrame: A dataframe containing the reactions and their rates

    """
    reactions = pd.read_csv(UCLCHEM_ROOT_DIR / "reactions.csv")
    return reactions


def find_number_of_consecutive_digits(string: str, start: int) -> int:
    """Determine the number of consecutive digits in a string, starting
    from some index `start`.

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


@enum.verify(enum.UNIQUE)
class SuccessFlag(enum.IntEnum):
    """SuccessFlag indicates whether the model failed or ran successfully,
    and if it failed how.

    """

    @staticmethod
    def _generate_next_value_(
            name: str, start: int, count: int, last_values: list[int | tuple[int, str]]  # noqa: ARG004
    ) -> int:
        last_value = last_values[-1]
        if isinstance(last_value, tuple):
            return last_value[0] - 1
        return last_value - 1

    def __new__(cls, value: int, docstring: str) -> Self:  # noqa: D102
        member = int.__new__(cls, value)

        member._value_ = value
        member.__doc__ = docstring

        return member

    SUCCESS = 0, "Model ran successfully"
    PARAMETER_READ_ERROR = enum.auto(), "Parameter read failed."
    PHYSICS_INIT_ERROR = enum.auto(), "Physics initialization failed."
    CHEM_INIT_ERROR = enum.auto(), "Chemistry initialization failed."
    INT_UNRECOVERABLE_ERROR = enum.auto(), "Unrecoverable integrator error occured."
    INT_TOO_MANY_FAILS_ERROR = enum.auto(), "Too many integrator fails occured."
    NOT_ENOUGH_TIMEPOINTS_ERROR = enum.auto(), "Not enough time points allocated in the time array."
    PHYSICS_UPDATE_ERROR = enum.auto(), "Error updating physics during integration."
    SOLVER_STATS_OVERFLOW_ERROR = enum.auto(), "Solver statistics array overflowed."
    COOLANT_FILE_ERROR = enum.auto(), "Coolant data file could not be openend."
    COOLANT_DATA_ERROR = enum.auto(), "Coolant data file has invalid format."
    COOLANT_FREQ_TOL_ERROR = enum.auto(), "Coolant frequency tolerance exceeded."
    COOLANT_POP_TOL_ERROR = enum.auto(), "LTE population sum tolerance exceeded."
    COOLANT_SOLVER_ERROR = enum.auto(), "Coolant solver numerical error occured."
    COOLANT_CONFIG_ERROR = enum.auto(), "Coolant configuration error occured."
    NEGATIVE_ABUNDANCE_ERROR = enum.auto(), "A negative abundance was detected."

    def check_error(self, only_error: bool = False, raise_on_error: bool = True) -> str:
        """Converts the UCLCHEM integer result flag to a message explaining what went wrong.

        Args:
            only_error (bool): If True, skip printing if the model ran successfully, and only
                error out if it did not. Default = False.
            raise_on_error (bool): If True, raises RuntimeError if the `self` is not
                `SuccessFlag.SUCCESS`. If False, returns the message string. Default = True.

        Returns:
            str: Error message | None

        Raises:
            RuntimeError: If `raise_on_error` is True.

        """
        if self == SuccessFlag.SUCCESS:
            if only_error:
                return None
            return "Model ran successfully"
        error_msg_dict = {
            SuccessFlag.PARAMETER_READ_ERROR: "Parameter read failed. Likely due to a misspelled parameter name, compare your dictionary to the parameters docs.",
            SuccessFlag.PHYSICS_INIT_ERROR: "Physics initialization failed. Often due to user choosing unacceptable parameters such as hot core masses or collapse modes that don't exist. Check the docs for your model function.",
            SuccessFlag.CHEM_INIT_ERROR: "Chemistry initialization failed.",
            SuccessFlag.INT_UNRECOVERABLE_ERROR: "Unrecoverable integrator error, DVODE failed to integrate the ODEs in a way that UCLCHEM could not fix. Run UCLCHEM tests to check your network works at all then try to see if bad parameter combination is at play.",
            SuccessFlag.INT_TOO_MANY_FAILS_ERROR: "Too many integrator fails. DVODE failed to integrate the ODE and UCLCHEM repeatedly altered settings to try to make it pass but tried too many times without success so code aborted to stop infinite loop.",
            SuccessFlag.NOT_ENOUGH_TIMEPOINTS_ERROR: "The model was stopped because there are not enough time points allocated in the time array. Increase the number of time points in the time array in constants.py and try again.",
            SuccessFlag.PHYSICS_UPDATE_ERROR: "Physics update error during integration.",
            SuccessFlag.SOLVER_STATS_OVERFLOW_ERROR: "Solver statistics array overflow.",
            SuccessFlag.COOLANT_FILE_ERROR: "Coolant data file could not be opened. Check that coolant data files exist in the expected directory.",
            SuccessFlag.COOLANT_DATA_ERROR: "Coolant data file has invalid format (bad NLEVEL, invalid partner ID, or too many temperature values).",
            SuccessFlag.COOLANT_FREQ_TOL_ERROR: "Frequency tolerance exceeded: the deviation between energy-level-computed and LAMDA-file frequencies exceeds freq_rel_tol. Increase freq_rel_tol or check your coolant data files.",
            SuccessFlag.COOLANT_POP_TOL_ERROR: "LTE population sum tolerance exceeded: level populations do not sum to total density within pop_rel_tol.",
            SuccessFlag.COOLANT_SOLVER_ERROR: "Coolant solver numerical error (NaN in matrix, singular matrix, or negative populations). The statistical equilibrium solver failed for a coolant species.",
            SuccessFlag.COOLANT_CONFIG_ERROR: "Coolant configuration error: parent species not found in network, or unphysical abundance detected.",
            SuccessFlag.NEGATIVE_ABUNDANCE_ERROR: "Negative abundance detected.",
        }

        msg = error_msg_dict[self]
        if raise_on_error:
            raise RuntimeError(f"UCLCHEM error (code {self.name}, {self.value}): {msg}")
        return msg
