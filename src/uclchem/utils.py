"""Helper functions and utilities for UCLCHEM operations.

This module provides utility functions for:
- Error handling and reporting
- Physics calculations (shock dissipation times)
- Parameter validation
- File path management

**Key Functions:**

- :meth:`SuccessFlag.check_error` - Convert UCLCHEM error codes to messages
- :func:`cshock_dissipation_time` - Calculate C-shock dissipation timescale

**Example Usage:**

    >>> import uclchem
    >>>
    >>> model = uclchem.model.Cloud({})
    >>> success_flag = model.success_flag
    >>>
    >>> # Check error from model run
    >>> success_flag.check_error()
    Model ran successfully
    >>>
    >>> # Only print if an error occurred
    >>> success_flag.check_error(only_error=True)

    >>> # Calculate shock timescale
    >>> t_diss = uclchem.utils.cshock_dissipation_time(
    ...     shock_vel=50.0,  # km/s
    ...     initial_dens=1e4  # cm^-3
    ... )
    >>> print(f"Dissipation time: {t_diss:.1e} years")
    Dissipation time: ... years

**Error Codes:**

UCLCHEM model functions return :class:`SuccessFlag` instances:

- ``-1``: Parameter read failed (misspelled parameter)
- ``-2``: Physics initialization failed (invalid parameters)
- ``-3``: Chemistry initialization failed
- ``-4``: Integrator error (DVODE failed)

and more...

Use :meth:`SuccessFlag.check_error` to get human-readable error messages.

**See Also:**

- :mod:`uclchem.model` - Model classes that use these utilities
"""

import enum
import logging
import sys
from io import TextIOWrapper
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Self, TypeAlias

if TYPE_CHECKING:
    from uclchem.model import Collapse

import numpy as np
import pandas as pd

from uclchem.constants import CENTIMETERS_PER_PARSEC, SECONDS_PER_YEAR

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
    dlength = 12.0 * CENTIMETERS_PER_PARSEC * shock_vel / initial_dens
    return (dlength * 1.0e-5 / shock_vel) / SECONDS_PER_YEAR


def get_species_table() -> pd.DataFrame:
    """Load the list of species in the UCLCHEM network into a pandas dataframe.

    Returns:
        species (pd.DataFrame): A dataframe containing the species names and their details

    """
    species = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv")
    return species


def get_species() -> list[str]:
    """Load the list of species present in the UCLCHEM network.

    Returns:
        species_list (list[str]): A list of species names

    """
    species_list = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv").iloc[:, 0].tolist()
    return species_list


def get_reaction_table() -> pd.DataFrame:
    """Load the reaction table from the UCLCHEM network into a pandas dataframe.

    Returns:
        reactions (pd.DataFrame): A dataframe containing the reactions and their rates

    """
    reactions = pd.read_csv(UCLCHEM_ROOT_DIR / "reactions.csv")
    return reactions


def find_number_of_consecutive_digits(string: str, start: int) -> int:
    """Determine the number of consecutive digits in a string, starting
    from some index ``start``.

    Args:
        string (str): the string
        start (int): the starting index

    Returns:
        num_digits (int): the number of consecutive digits in the string
            starting from ``start``.

    Examples:
        >>> find_number_of_consecutive_digits("Hello123", 0)
        0
        >>> find_number_of_consecutive_digits("Hello123", 5)
        3
        >>> find_number_of_consecutive_digits("Hello123", 6)
        2
        >>> find_number_of_consecutive_digits("He1llo23", 2)
        1

    """
    num_digits = 0
    while start + num_digits < len(string) and string[start + num_digits].isdigit():
        num_digits += 1
    return num_digits


# ---------------------------------------------------------------------------
# Collapse radial velocity — Priestley et al. 2018
# ---------------------------------------------------------------------------

# Physical constants matching collapse.f90
_PC = 3.086e18  # parsec in cm
_MH = 1.6736e-24  # hydrogen mass in g
_KB = 1.38e-16  # Boltzmann constant in erg/K
_G = 6.67e-8  # gravitational constant in cgs
_RHO0_FILAMENT = 2.2e4  # reference density for filament/ambipolar (cm^-3)
_TWO_PI_G = 2.0 * np.pi * _G


def _filament_units():
    """Return (unitr_pc, unitt_yr) for filament (mode 3) collapse."""
    two_pi_g_rho0_mh = _TWO_PI_G * _RHO0_FILAMENT * _MH
    cs = np.sqrt(_KB * 10.0 / (2.0 * _MH))  # sound speed at 10 K
    unitr = cs * two_pi_g_rho0_mh ** (-0.5) / _PC  # in pc
    unitt = two_pi_g_rho0_mh ** (-0.5) / SECONDS_PER_YEAR  # in yr
    return unitr, unitt


def _rminfit(t_yr: float, mode: int) -> float:
    """Fit to time evolution of the radius of minimum velocity.

    Returns:
        Radius of minimum velocity (pc for mode 3, normalized units for mode 4).
    """
    if mode == 3:
        _, unitt = _filament_units()
        tnew = t_yr / unitt
        if tnew == 0.0:
            return 7.2
        elif np.log(tnew) < 1.6:
            return -1.149 * tnew + 7.2
        elif np.log(tnew) < 1.674:
            return -9.2 * np.log(tnew) + 16.25
        else:
            return -22.0 * np.log(tnew) + 37.65
    else:  # mode 4
        t6 = 1e-6 * t_yr
        if t6 <= 10.2:
            return -0.0039 * t6 + 0.49
        elif t6 <= 15.1:
            return -0.0306 * (t6 - 10.2) + 0.45
        else:
            return -0.282 * (t6 - 15.1) + 0.3


def _vminfit(t_yr: float, mode: int) -> float:
    """Fit to time evolution of minimum velocity (dimensionless units).

    Returns:
        Minimum velocity in dimensionless units.
    """
    if mode == 3:
        _, unitt = _filament_units()
        tnew = t_yr / unitt
        if tnew == 0.0:
            return 0.0
        elif np.log(tnew) < 1.6:
            return 0.0891 * tnew
        elif np.log(tnew) < 1.674:
            return 5.5 * np.log(tnew) - 8.37
        else:
            return 18.9 * np.log(tnew) - 30.8
    else:  # mode 4
        t6 = 1e-6 * t_yr
        return 3.44 * (16.138 - t6) ** (-0.35) - 0.7


def _avfit(t_yr: float, mode: int) -> float:
    """Fit to velocity a-parameter (mode 4) or velocity at r=0.5 (mode 3).

    Returns:
        Velocity a-parameter (mode 4) or velocity at r=0.5 (mode 3).
    """
    if mode == 3:
        _, unitt = _filament_units()
        tnew = t_yr / unitt
        if tnew == 0.0:
            return 0.4
        elif np.log(tnew) < 1.6:
            return 0.0101 * tnew + 0.4
        elif np.log(tnew) < 1.674:
            return 0.695 * np.log(tnew) - 0.663
        else:
            return 2.69 * np.log(tnew) - 4.0
    else:  # mode 4
        t6 = 1e-6 * t_yr
        if t6 <= 10.2:
            return 0.143 * t6
        else:
            return 0.217 * (t6 - 10.2) + 1.46


def _vrfit(r_pc: float, rmin: float, vmin: float, av: float, mode: int) -> float:
    """Radial velocity fit in cm/s (Priestley et al. 2018).

    Modes 3 (filament) and 4 (ambipolar) only.

    Returns:
        Radial velocity in cm/s.
    """
    if mode == 3:
        unitr, _ = _filament_units()
        cs = np.sqrt(_KB * 10.0 / (2.0 * _MH))
        new_r = r_pc / unitr - rmin
        if new_r < 0.0:
            vr = vmin * ((new_r / rmin) ** 2 - 1.0)
        else:
            vr = vmin * (np.exp(-2.0 * av * new_r) - 2.0 * np.exp(-av * new_r))
        return cs * vr
    else:  # mode 4
        rmid = 0.5
        r75 = r_pc / 0.75
        new_r = r75 - rmin
        if r75 < rmin:
            vr = vmin * ((new_r / rmin) ** 2 - 1.0)
        elif r75 <= rmid:
            vr = (vmin - av) * (new_r / (rmid - rmin)) ** 0.3 - vmin
        else:
            vr = av / (1.0 - rmid) * (r75 - rmid) - av
        return 1e3 * vr  # convert from 1e-2 km/s to cm/s


def collapse_radial_velocity(model: "Collapse", point: int = 0) -> pd.Series:
    """Return the radial velocity (cm/s) for a parcel of a Collapse model.

    For filament (mode 3) and ambipolar (mode 4) collapse modes, uses the
    analytical radial-velocity fit functions from Priestley et al. (2018).

    For BE1.1 and BE4 modes (1 & 2), the radius is tracked via mass-conservation
    integration in Fortran, not a velocity fit. The radial velocity is therefore
    a finite-difference approximation of parcel_radius — it is NOT the relationship
    used to generate the model and should be treated as an estimate only.

    Args:
        model: A successfully run :class:`~uclchem.model.Collapse` instance.
        point: Parcel index (0-based). Defaults to 0.

    Returns:
        pd.Series: Radial velocity in cm/s, indexed by time in years.
                   Negative values indicate infall.

    Raises:
        TypeError: If ``model`` is not a Collapse model instance.
    """
    from uclchem.model import Collapse

    if not isinstance(model, Collapse):
        msg = f"model must be a Collapse instance, got {type(model).__name__}"
        raise TypeError(msg)

    df: pd.DataFrame = model.get_dataframes(point=point)  # type: ignore[assignment]
    t_yr = df["Time"]
    r_pc = df["parcel_radius"]
    mode = model.collapse  # integer 1-4

    if mode in (3, 4):
        vr = np.array(
            [
                _vrfit(r, _rminfit(t, mode), _vminfit(t, mode), _avfit(t, mode), mode)
                for t, r in zip(t_yr, r_pc)
            ]
        )
    else:
        # BE-sphere modes: approximate via finite differences of parcel_radius.
        # This is NOT the relationship used to generate the model.
        t_s = t_yr * SECONDS_PER_YEAR
        r_cm = r_pc * _PC
        vr = np.gradient(r_cm, t_s)

    return pd.Series(vr, index=t_yr, name="radial_velocity_cm_s")


@enum.verify(enum.UNIQUE)
class SuccessFlag(enum.IntEnum):
    """SuccessFlag indicates whether the model failed or ran successfully,
    and if it failed how.

    """

    def __new__(cls, value: int, docstring: str = "") -> Self:  # noqa: D102
        member = int.__new__(cls, value)

        member._value_ = value
        member.__doc__ = docstring

        return member

    # Zen line two: Explicit is better than implicit.
    SUCCESS = 0, "Model ran successfully"
    PARAMETER_READ_ERROR = -1, "Parameter read failed."
    PHYSICS_INIT_ERROR = -2, "Physics initialization failed."
    CHEM_INIT_ERROR = -3, "Chemistry initialization failed."
    INT_UNRECOVERABLE_ERROR = -4, "Unrecoverable integrator error occurred."
    INT_TOO_MANY_FAILS_ERROR = -5, "Too many integrator fails occurred."
    NOT_ENOUGH_TIMEPOINTS_ERROR = (
        -6,
        "Not enough time points allocated in the time array.",
    )
    PHYSICS_UPDATE_ERROR = -7, "Error updating physics during integration."
    SOLVER_STATS_OVERFLOW_ERROR = -8, "Solver statistics array overflowed."
    COOLANT_FILE_ERROR = -9, "Coolant data file could not be opened."
    COOLANT_DATA_ERROR = -10, "Coolant data file has invalid format."
    COOLANT_FREQ_TOL_ERROR = -11, "Coolant frequency tolerance exceeded."
    COOLANT_POP_TOL_ERROR = -12, "LTE population sum tolerance exceeded."
    COOLANT_SOLVER_ERROR = -13, "Coolant solver numerical error occurred."
    COOLANT_CONFIG_ERROR = -14, "Coolant configuration error occurred."
    NEGATIVE_ABUNDANCE_ERROR = -15, "A negative abundance was detected."

    def check_error(
        self, only_error: bool = False, raise_on_error: bool = True
    ) -> str | None:
        """Converts the UCLCHEM integer result flag to a message explaining what went wrong.

        Args:
            only_error (bool): If True, skip printing if the model ran successfully, and only
                error out if it did not. Default = False.
            raise_on_error (bool): If True, raises RuntimeError if the ``self`` is not
                ``SuccessFlag.SUCCESS``. If False, returns the message string. Default = True.

        Returns:
            str | None: Error message. Returns None if no error was found.

        Raises:
            RuntimeError: If ``raise_on_error`` is True.

        """
        if self == SuccessFlag.SUCCESS:
            if not only_error:
                print("Model ran successfully")
            return None
        error_msg_dict = {
            SuccessFlag.PARAMETER_READ_ERROR: "Parameter read failed. Likely due to a misspelled parameter name, compare your dictionary to the parameters docs.",
            SuccessFlag.PHYSICS_INIT_ERROR: "Physics initialization failed. Often due to user choosing unacceptable parameters such as hot core masses or collapse modes that don't exist. Check the docs for your model function.",
            SuccessFlag.CHEM_INIT_ERROR: "Chemistry initialization failed.",
            SuccessFlag.INT_UNRECOVERABLE_ERROR: "Unrecoverable integrator error, DVODE failed to integrate the ODEs in a way that UCLCHEM could not fix. Run UCLCHEM tests to check your network works at all then try to see if bad parameter combination is at play.",
            SuccessFlag.INT_TOO_MANY_FAILS_ERROR: "Too many integrator fails. DVODE failed to integrate the ODE and UCLCHEM repeatedly altered settings to try to make it pass but tried too many times without success so code aborted to stop infinite loop.",
            SuccessFlag.NOT_ENOUGH_TIMEPOINTS_ERROR: "The model was stopped because there are not enough time points allocated in the time array. Increase the number of time points in the time array in constants.py and try again.",
            SuccessFlag.PHYSICS_UPDATE_ERROR: "Physics update error during integration.",
            SuccessFlag.SOLVER_STATS_OVERFLOW_ERROR: "Solver statistics array overflows, add more timepoints",
            SuccessFlag.COOLANT_FILE_ERROR: "Coolant data file could not be opened. Check that coolant data files exist in the expected directory.",
            SuccessFlag.COOLANT_DATA_ERROR: "Coolant data file has invalid format (bad NLEVEL, invalid partner ID, or too many temperature values).",
            SuccessFlag.COOLANT_FREQ_TOL_ERROR: "Frequency tolerance exceeded: the deviation between energy-level-computed and LAMDA-file frequencies exceeds freq_rel_tol. Increase freq_rel_tol or check your coolant data files.",
            SuccessFlag.COOLANT_POP_TOL_ERROR: "LTE population sum tolerance exceeded: level populations do not sum to total density within pop_rel_tol.",
            SuccessFlag.COOLANT_SOLVER_ERROR: "Coolant solver numerical error (NaN in matrix, singular matrix, or negative populations). The statistical equilibrium solver failed for a coolant species.",
            SuccessFlag.COOLANT_CONFIG_ERROR: "Coolant configuration error: parent species not found in network, or unphysical abundance detected.",
            SuccessFlag.NEGATIVE_ABUNDANCE_ERROR: "Negative abundance detected. That exceeds solver tolerances, consider adjusting the negative_abundance_tol parameter in param_dict to a larger magnitude",
        }

        msg = error_msg_dict[self]
        if raise_on_error:
            msg = f"UCLCHEM error (code {self.name}, {self.value}): {msg}"
            raise RuntimeError(msg)
        return msg


def check_expected_type(
    variable: Any, expected_type: type[Any], name: str | None = None
) -> None:
    """Check that the type of a variable matches the expected type.

    Args:
        variable (Any): variable to check type of.
        expected_type (Type[Any]): expected type.
        name (str | None): Name of variable. If None, no name information will be printed.
            Defaults to None.

    Raises:
        TypeError: If ``variable`` is not an instance of ``expected_type``.

    """
    if isinstance(variable, expected_type):
        return

    if name is not None:
        msg = f"{name} was supposed to be type {expected_type} but got type {type(variable)}"
    else:
        msg = f"Expected type {expected_type} but got type {type(variable)}"
    raise TypeError(msg)


ArrayLike: TypeAlias = list | pd.Series | np.ndarray


def configure_logging(
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] | int = "INFO",
    stream: TextIOWrapper | str | Path | None = sys.stdout,  # type: ignore[assignment]
) -> None:
    """Configure logging of UCLCHEM.

    Args:
        level (Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'] | int): Level of logs to be
            logged. Case insensitive. Default = "INFO".
        stream (TextIOWrapper | str | Path | None): stream to write to. If string or Path, write to file
            in append mode. If None, do not write logs at all. Default = sys.stdout (write to stdout).

    Raises:
        TypeError: If stream is not an instance of TextIOWrapper, str, Path, or None.

    Examples:
        >>> # Get INFO or higher messages in stdout.
        >>> configure_logging()

        >>> # Also get DEBUG messages in stdout
        >>> configure_logging(level="DEBUG")

        >>> # Write WARNING or higher messages to "uclchem.log"
        >>> configure_logging(level="WARNING", stream="uclchem.log")

        >>> # Suppress all logging messages
        >>> configure_logging(stream = None)

    """
    if isinstance(level, str):
        level = level.upper()  # type: ignore[assignment]

    handler: logging.Handler
    if stream is None:
        handler = logging.NullHandler()
    elif isinstance(stream, TextIOWrapper):
        handler = logging.StreamHandler(stream=stream)
    elif isinstance(stream, str | Path):
        handler = logging.FileHandler(filename=stream)
    else:
        msg = f"stream should be None, or type TextIOWrapper (such as sys.stdout), string or Path, but got type {type(stream)}"
        raise TypeError(msg)

    logger = logging.getLogger("uclchem")
    logger.propagate = False  # Do not propagate to the root logger
    logger.setLevel(level)

    handler.setLevel(level)
    formatter = logging.Formatter(fmt="%(levelname)-8s %(message)s")
    handler.setFormatter(formatter)

    logger.addHandler(handler)

    logger.debug(
        f"Logging configured with level {logging.getLevelName(logger.getEffectiveLevel())}, with handler {handler}"
    )
