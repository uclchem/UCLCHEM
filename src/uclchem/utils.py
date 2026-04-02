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

from pathlib import Path

import numpy as np
import pandas as pd

UCLCHEM_ROOT_DIR: Path = Path(__file__).parent.resolve().absolute()


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


def check_error(error_code: int, raise_on_error: bool = True) -> str:
    """Converts the UCLCHEM integer result flag to a message explaining what went wrong.

    Args:
        error_code (int): Error code returned by UCLCHEM models, the first element of the results list.
        raise_on_error (bool): If True (default), raises RuntimeError. If False, returns the message string.

    Returns:
        str: Error message

    Raises:
        RuntimeError: If raise_on_error is True and error_code is recognized.
    """
    errors = {
        -1: "Parameter read failed. Likely due to a misspelled parameter name, compare your dictionary to the parameters docs.",
        -2: "Physics initialization failed. Often due to user choosing unacceptable parameters such as hot core masses or collapse modes that don't exist. Check the docs for your model function.",
        -3: "Chemistry initialization failed.",
        -4: "Unrecoverable integrator error, DVODE failed to integrate the ODEs in a way that UCLCHEM could not fix. Run UCLCHEM tests to check your network works at all then try to see if bad parameter combination is at play.",
        -5: "Too many integrator fails. DVODE failed to integrate the ODE and UCLCHEM repeatedly altered settings to try to make it pass but tried too many times without success so code aborted to stop infinite loop.",
        -6: "The model was stopped because there are not enough time points allocated in the time array. Increase the number of time points in the time array in constants.py and try again.",
        -7: "Physics update error during integration.",
        -8: "Solver statistics array overflow.",
        -9: "Coolant data file could not be opened. Check that coolant data files exist in the expected directory.",
        -10: "Coolant data file has invalid format (bad NLEVEL, invalid partner ID, or too many temperature values).",
        -11: "Frequency tolerance exceeded: the deviation between energy-level-computed and LAMDA-file frequencies exceeds freq_rel_tol. Increase freq_rel_tol or check your coolant data files.",
        -12: "LTE population sum tolerance exceeded: level populations do not sum to total density within pop_rel_tol.",
        -13: "Coolant solver numerical error (NaN in matrix, singular matrix, or negative populations). The statistical equilibrium solver failed for a coolant species.",
        -14: "Coolant configuration error: parent species not found in network, or unphysical abundance detected.",
    }
    msg = errors.get(error_code, f"Unknown error code: {error_code}")
    if raise_on_error:
        raise RuntimeError(f"UCLCHEM error (code {error_code}): {msg}")
    return msg


def get_species_table() -> pd.DataFrame:
    """A simple function to load the list of species in the UCLCHEM network into a pandas dataframe.

    Returns:
        pandas.DataFrame: A dataframe containing the species names and their details
    """

    species_list = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv")
    return species_list


def get_species() -> list[str]:
    """A simple function to load the list of species present in the UCLCHEM network

    Returns:
        list[str] : A list of species names
    """

    species_list = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv").iloc[:, 0].tolist()
    return species_list


def get_reaction_table() -> pd.DataFrame:
    """A function to load the reaction table from the UCLCHEM network into a pandas dataframe.

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


# ---------------------------------------------------------------------------
# Collapse radial velocity — Priestley et al. 2018
# ---------------------------------------------------------------------------

# Physical constants matching collapse.f90
_PC = 3.086e18           # parsec in cm
_MH = 1.6736e-24         # hydrogen mass in g
_KB = 1.38e-16           # Boltzmann constant in erg/K
_G = 6.67e-8             # gravitational constant in cgs
_SECONDS_PER_YEAR = 3.15569e7
_RHO0_FILAMENT = 2.2e4   # reference density for filament/ambipolar (cm^-3)
_TWO_PI_G = 2.0 * np.pi * _G


def _filament_units():
    """Return (unitr_pc, unitt_yr) for filament (mode 3) collapse."""
    two_pi_g_rho0_mh = _TWO_PI_G * _RHO0_FILAMENT * _MH
    cs = np.sqrt(_KB * 10.0 / (2.0 * _MH))          # sound speed at 10 K
    unitr = cs * two_pi_g_rho0_mh ** (-0.5) / _PC   # in pc
    unitt = two_pi_g_rho0_mh ** (-0.5) / _SECONDS_PER_YEAR  # in yr
    return unitr, unitt


def _rminfit(t_yr: float, mode: int) -> float:
    """Fit to time evolution of the radius of minimum velocity."""
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
    """Fit to time evolution of minimum velocity (dimensionless units)."""
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
    """Fit to velocity a-parameter (mode 4) or velocity at r=0.5 (mode 3)."""
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


def collapse_radial_velocity(model, point: int = 0) -> pd.Series:
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
        pd.Series: Radial velocity in cm s⁻¹, indexed by time in years.
                   Negative values indicate infall.

    Raises:
        TypeError: If *model* is not a Collapse model instance.
    """
    from uclchem.model import Collapse

    if not isinstance(model, Collapse):
        raise TypeError(
            f"model must be a Collapse instance, got {type(model).__name__}"
        )

    df = model.get_dataframes(point=point)
    t_yr = df["Time"].values
    r_pc = df["parcel_radius"].values
    mode = model.collapse  # integer 1-4

    if mode in (3, 4):
        vr = np.array([
            _vrfit(r, _rminfit(t, mode), _vminfit(t, mode), _avfit(t, mode), mode)
            for t, r in zip(t_yr, r_pc)
        ])
    else:
        # BE-sphere modes: approximate via finite differences of parcel_radius.
        # This is NOT the relationship used to generate the model.
        t_s = t_yr * _SECONDS_PER_YEAR
        r_cm = r_pc * _PC
        vr = np.gradient(r_cm, t_s)

    return pd.Series(vr, index=t_yr, name="radial_velocity_cm_s")
