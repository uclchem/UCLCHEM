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
    ...     initial_dens=1e4,  # cm^-3
    ... )
    >>> print(f"Dissipation time: {t_diss:.1e} years")  # doctest: +ELLIPSIS
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
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Self

if TYPE_CHECKING:
    from uclchem.model import Collapse

import numpy as np
import pandas as pd

UCLCHEM_ROOT_DIR: Path = Path(__file__).parent.resolve().absolute()
"""UCLCHEM root directory"""


MISSING_VALUE_INTEGER: int = -1
"""Integer to indicate a missing value"""

MISSING_VALUE_FLOAT: float = -1.0
"""Float to indicate a missing value"""

NO_REACTANT_OR_PRODUCT: int = 9999
"""Integer to indicate that there is no reactant or product."""

DTYPE_MAPPING = MappingProxyType(
    {
        "fp128": np.float128,
        "fp64": np.float64,
        "fp32": np.float32,
        "fp16": np.float16,
    }
)
"""Mapping from shorthand strings to numpy dtypes."""


def get_dtype(
    dtype: str | type[np.dtype] | np.dtype | np.typing.DTypeLike,
) -> np.typing.DTypeLike:
    """Convert a dtype shorthand string or numpy dtype to a numpy dtype.

    Parameters
    ----------
    dtype : str | type[np.dtype] | np.dtype | np.typing.DTypeLike
        Either a shorthand string ("fp128", "fp64", "fp32", "fp16") or a numpy dtype/type.

    Returns
    -------
    np.typing.DTypeLike
        The corresponding numpy dtype or type.

    Raises
    ------
    ValueError
        If `dtype` is a string not in :data:`DTYPE_MAPPING`.

    """
    if isinstance(dtype, str):
        if dtype not in DTYPE_MAPPING:
            msg = f"Unknown dtype shorthand: {dtype!r}"
            raise ValueError(msg)
        return np.dtype(DTYPE_MAPPING[dtype])
    return np.dtype(dtype)


def would_overflow(value: float, dtype: np.dtype) -> bool:
    """Whether `value` would overflow (become ``np.inf``) if cast to `dtype`.

    Parameters
    ----------
    value : float
        Value to check
    dtype : np.dtype
        Dtype to cast to

    Returns
    -------
    bool
        Whether `value` would be ``np.inf`` if cast to `dtype`.

    """
    return np.abs(value) > np.finfo(dtype).max


def configure_logging(
    level: str = "WARNING",
    stream: str | Any | None = None,
) -> None:
    """Configure the ``uclchem`` logger.

    Parameters
    ----------
    level : str
        Logging level name (e.g. ``"DEBUG"``, ``"WARNING"``). Default ``"WARNING"``.
    stream : str | Any | None
        Where to send log output.
        ``None`` suppresses all output (no handlers, propagation disabled).
        A string is treated as a file path and a ``FileHandler`` is added.
        Any other value (e.g. ``sys.stdout``) is wrapped in a ``StreamHandler``.
        Defaults to ``None``.

    """
    uclchem_logger = logging.getLogger("uclchem")
    for handler in uclchem_logger.handlers[:]:
        uclchem_logger.removeHandler(handler)

    level_value = getattr(logging, level) if isinstance(level, str) else level
    uclchem_logger.setLevel(level_value)

    if stream is None:
        uclchem_logger.propagate = False
        return

    if isinstance(stream, str):
        handler = logging.FileHandler(stream)
    else:
        handler = logging.StreamHandler(stream)

    handler.setLevel(level_value)
    uclchem_logger.addHandler(handler)
    uclchem_logger.propagate = True


# Pandas default NA strings minus element symbols that collide ("NA" = sodium, "#NA")
_NAN_STRINGS: list[str] = [
    "",
    "#N/A",
    "#N/A N/A",
    "-1.#IND",
    "-1.#QNAN",
    "-NaN",
    "-nan",
    "1.#IND",
    "1.#QNAN",
    "<NA>",
    "N/A",
    "NULL",
    "NaN",
    "None",
    "n/a",
    "nan",
    "null",
]


def get_species_table(file: str | Path | None = None) -> pd.DataFrame:
    """Load the list of species in the UCLCHEM network into a pandas dataframe.

    Parameters
    ----------
    file : str | Path | None
        Path to the species CSV file. If None, uses
        ``UCLCHEM_ROOT_DIR / "species.csv"``. Default = None.

    Returns
    -------
    species : pd.DataFrame
        A dataframe containing the species names and their details

    """
    if file is None:
        file = UCLCHEM_ROOT_DIR / "species.csv"

    species = pd.read_csv(file, na_values=_NAN_STRINGS, keep_default_na=False)
    return species


def get_species(file: str | Path | None = None) -> list[str]:
    """Load the list of species present in the UCLCHEM network.

    Parameters
    ----------
    file : str | Path | None
        Path to the species CSV file. If None, uses
        ``UCLCHEM_ROOT_DIR / "species.csv"``. Default = None.

    Returns
    -------
    species_list : list[str]
        A list of species names

    """
    species_list = get_species_table(file=file)["NAME"].tolist()
    return species_list


def get_reaction_table(file: str | Path | None = None) -> pd.DataFrame:
    """Load the reaction table from the UCLCHEM network into a pandas dataframe.

    Parameters
    ----------
    file : str | Path | None
        Path to the reactions CSV file. If None, uses
        ``UCLCHEM_ROOT_DIR / "reactions.csv"``. Default = None.

    Returns
    -------
    reactions : pd.DataFrame
        A dataframe containing the reactions and their rates

    """
    if file is None:
        file = UCLCHEM_ROOT_DIR / "reactions.csv"

    reactions = pd.read_csv(file, na_values=_NAN_STRINGS, keep_default_na=False)
    return reactions


@enum.unique
class PrestellarCoreMass(enum.IntEnum):
    """Index which indicates which index to use for heating in ``hotcore.f90``.

    The different values correspond to different solar masses of the prestellar core.

    """

    M_SUN_1 = 1
    M_SUN_5 = 2
    M_SUN_10 = 3
    M_SUN_15 = 4
    M_SUN_25 = 5
    M_SUN_60 = 6


# ---------------------------------------------------------------------------
# Collapse radial velocity — Priestley et al. 2018
# ---------------------------------------------------------------------------

# Physical constants matching collapse.f90
_PC = 3.086e18  # parsec in cm
_MH = 1.6736e-24  # hydrogen mass in g
_KB = 1.38e-16  # Boltzmann constant in erg/K
_G = 6.67e-8  # gravitational constant in cgs
_SECONDS_PER_YEAR = 3.15569e7
_RHO0_FILAMENT = 2.2e4  # reference density for filament/ambipolar (cm^-3)
_TWO_PI_G = 2.0 * np.pi * _G


def cshock_dissipation_time(shock_vel: float, initial_dens: float) -> float:
    """Calculate the dissipation time of a C-type shock.

    Use to obtain a useful timescale for your C-shock model runs.
    Velocity of ions and neutrals equalizes at dissipation time and
    full cooling takes a few dissipation times.

    Parameters
    ----------
    shock_vel : float
        Velocity of the shock in km/s
    initial_dens : float
        Preshock density of the gas in cm$^{-3}$

    Returns
    -------
    float
        The dissipation time of the shock in years

    """
    dlength = 12.0 * _PC * shock_vel / initial_dens
    return (dlength * 1.0e-5 / shock_vel) / _SECONDS_PER_YEAR


@enum.unique
class CollapseMode(enum.IntEnum):
    """Collapse mode."""

    BE1_1 = 1
    BE4 = 2
    FILAMENT = 3
    AMBIPOLAR = 4


def get_collapse_mode(mode: str | int | CollapseMode) -> CollapseMode:
    """Get the CollapseMode instance.

    Integers are converted to the CollapseMode, CollapseMode instances are returned
    unchanged, and strings are converted according to the names of the different modes.

    Parameters
    ----------
    mode : str | int | CollapseMode
        Collapse mode.

    Returns
    -------
    CollapseMode
        collapse mode Enum

    Raises
    ------
    TypeError
        If ``mode`` is not a string, integer or :class:`CollapseMode`.
    ValueError
        If ``mode`` is a string, but not one of
        ``["BE1.1", "BE4", "filament", "ambipolar"]`` (case insensitive).

    """
    if isinstance(mode, CollapseMode):
        return mode
    if isinstance(mode, int):
        return CollapseMode(mode)
    if not isinstance(mode, str):
        msg = f"Expected 'mode' to be type string, int or CollapseMode, but got type {type(mode)}"
        raise TypeError(msg)
    mode = mode.lower()
    if mode == "be1.1":
        return CollapseMode.BE1_1
    if mode == "be4":
        return CollapseMode.BE4
    if mode == "filament":
        return CollapseMode.FILAMENT
    if mode == "ambipolar":
        return CollapseMode.AMBIPOLAR
    msg = f"If 'mode' is a string, it should be one of ['BE1.1', 'BE4', 'filament', 'ambipolar'], but got {mode}"
    raise ValueError(msg)


def _filament_units():
    """Return (unitr_pc, unitt_yr) for filament (mode 3) collapse.

    Returns
    -------
    tuple[float, float]
        Unit conversion factors ``(unitr_pc, unitt_yr)``.

    """
    two_pi_g_rho0_mh = _TWO_PI_G * _RHO0_FILAMENT * _MH
    cs = np.sqrt(_KB * 10.0 / (2.0 * _MH))  # sound speed at 10 K
    unitr = cs * two_pi_g_rho0_mh ** (-0.5) / _PC  # in pc
    unitt = two_pi_g_rho0_mh ** (-0.5) / _SECONDS_PER_YEAR  # in yr
    return unitr, unitt


def _rminfit(t_yr: float, mode: CollapseMode) -> float:
    """Fit to time evolution of the radius of minimum velocity.

    Parameters
    ----------
    t_yr : float
        Time in years.
    mode : CollapseMode
        Collapse mode. One of ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    Returns
    -------
    float
        Radius of minimum velocity (pc for mode 3, normalized units for mode 4).

    Raises
    ------
    ValueError
        If ``mode`` is not one of
        ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    """
    if mode == CollapseMode.FILAMENT:
        _, unitt = _filament_units()
        tnew = t_yr / unitt
        if tnew == 0:
            return 7.2
        if np.log(tnew) < 1.6:  # ruff: ignore[magic-value-comparison]
            return -1.149 * tnew + 7.2
        if np.log(tnew) < 1.674:  # ruff: ignore[magic-value-comparison]
            return -9.2 * np.log(tnew) + 16.25
        return -22.0 * np.log(tnew) + 37.65
    if mode == CollapseMode.AMBIPOLAR:
        t6 = 1e-6 * t_yr
        if t6 <= 10.2:  # ruff: ignore[magic-value-comparison]
            return -0.0039 * t6 + 0.49
        if t6 <= 15.1:  # ruff: ignore[magic-value-comparison]
            return -0.0306 * (t6 - 10.2) + 0.45
        return -0.282 * (t6 - 15.1) + 0.3
    msg = "mode was not one of CollapseMode.FILAMENT or CollapseMode.AMBIPOLAR."
    raise ValueError(msg)


def _vminfit(t_yr: float, mode: CollapseMode) -> float:
    """Fit to time evolution of minimum velocity (dimensionless units).

    Parameters
    ----------
    t_yr : float
        Time in years.
    mode : CollapseMode
        Collapse mode. One of ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    Returns
    -------
    float
        Minimum velocity in dimensionless units.

    Raises
    ------
    ValueError
        If ``mode`` is not one of
        ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    """
    if mode == CollapseMode.FILAMENT:
        _, unitt = _filament_units()
        tnew = t_yr / unitt
        if tnew == 0:
            return 0.0
        if np.log(tnew) < 1.6:  # ruff: ignore[magic-value-comparison]
            return 0.0891 * tnew
        if np.log(tnew) < 1.674:  # ruff: ignore[magic-value-comparison]
            return 5.5 * np.log(tnew) - 8.37
        return 18.9 * np.log(tnew) - 30.8
    if mode == CollapseMode.AMBIPOLAR:
        t6 = 1e-6 * t_yr
        return 3.44 * (16.138 - t6) ** (-0.35) - 0.7
    msg = "mode was not one of CollapseMode.FILAMENT or CollapseMode.AMBIPOLAR."
    raise ValueError(msg)


def _avfit(t_yr: float, mode: CollapseMode) -> float:
    """Fit to velocity a-parameter (mode 4) or velocity at r=0.5 (mode 3).

    Parameters
    ----------
    t_yr : float
        Time in years.
    mode : CollapseMode
        Collapse mode. One of ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    Returns
    -------
    float
        Velocity a-parameter (mode 4) or velocity at r=0.5 (mode 3).

    Raises
    ------
    ValueError
        If ``mode`` is not one of
        ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    """
    if mode == CollapseMode.FILAMENT:
        _, unitt = _filament_units()
        tnew = t_yr / unitt
        if tnew == 0:
            return 0.4
        if np.log(tnew) < 1.6:  # ruff: ignore[magic-value-comparison]
            return 0.0101 * tnew + 0.4
        if np.log(tnew) < 1.674:  # ruff: ignore[magic-value-comparison]
            return 0.695 * np.log(tnew) - 0.663
        return 2.69 * np.log(tnew) - 4.0
    if mode == CollapseMode.AMBIPOLAR:
        t6 = 1e-6 * t_yr
        if t6 <= 10.2:  # ruff: ignore[magic-value-comparison]
            return 0.143 * t6
        return 0.217 * (t6 - 10.2) + 1.46
    msg = "mode was not one of CollapseMode.FILAMENT or CollapseMode.AMBIPOLAR."
    raise ValueError(msg)


def _vrfit(r_pc: float, rmin: float, vmin: float, av: float, mode: CollapseMode) -> float:
    """Radial velocity fit in cm/s (Priestley et al. 2018).

    Modes 3 (filament) and 4 (ambipolar) only.

    Parameters
    ----------
    r_pc : float
        Radius in pc.
    rmin : float
        Radius of minimum velocity in the same units as the fit.
    vmin : float
        Minimum velocity in dimensionless units.
    av : float
        Velocity amplitude parameter for the outer profile fit.
    mode : CollapseMode
        Collapse mode. One of ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    Returns
    -------
    float
        Radial velocity in cm/s.

    Raises
    ------
    ValueError
        If ``mode`` is not one of
        ``CollapseMode.FILAMENT`` or ``CollapseMode.AMBIPOLAR``.

    """
    if mode == CollapseMode.FILAMENT:
        unitr, _ = _filament_units()
        cs = np.sqrt(_KB * 10.0 / (2.0 * _MH))
        new_r = r_pc / unitr - rmin
        if new_r < 0.0:
            vr = vmin * ((new_r / rmin) ** 2 - 1.0)
        else:
            vr = vmin * (np.exp(-2.0 * av * new_r) - 2.0 * np.exp(-av * new_r))
        return cs * vr
    if mode == CollapseMode.AMBIPOLAR:
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
    msg = "mode was not one of CollapseMode.FILAMENT or CollapseMode.AMBIPOLAR."
    raise ValueError(msg)


def collapse_radial_velocity(model: "Collapse", point: int = 0) -> pd.Series:
    """Return the radial velocity (cm/s) for a parcel of a Collapse model.

    For filament (mode 3) and ambipolar (mode 4) collapse modes, uses the
    analytical radial-velocity fit functions from Priestley et al. (2018).

    For BE1.1 and BE4 modes (1 & 2), the radius is tracked via mass-conservation
    integration in Fortran, not a velocity fit. The radial velocity is therefore
    a finite-difference approximation of parcel_radius — it is NOT the relationship
    used to generate the model and should be treated as an estimate only.

    Parameters
    ----------
    model : Collapse
        A successfully run :class:`~uclchem.model.Collapse` instance.
    point : int
        Parcel index (0-based). Defaults to 0.

    Returns
    -------
    pd.Series
        Radial velocity in cm s⁻¹, indexed by time in years.
        Negative values indicate infall.

    Raises
    ------
    TypeError
        If *model* is not a Collapse model instance.
    TypeError
        If ``model.collapse`` is not an instance of class:`CollapseMode`.

    """
    from uclchem.model import (  # ruff: ignore[import-outside-top-level] — avoid circular import at module level
        Collapse,
    )

    if not isinstance(model, Collapse):
        msg = f"model must be a Collapse instance, got {type(model).__name__}"
        raise TypeError(msg)

    df_result = model.get_dataframes(point=point)
    # get_dataframes with joined=True returns a single DataFrame (the return type
    # annotation is overly broad; cast to narrow for mypy)
    df: pd.DataFrame = df_result  # type: ignore[assignment]
    t_yr: np.ndarray = np.asarray(df["Time"].values)
    r_pc: np.ndarray = np.asarray(df["parcel_radius"].values)
    mode = model.collapse  # CollapseMode

    if not isinstance(mode, CollapseMode):
        msg = f"Expected 'model.collapse' to be an instance of CollapseMode, but got type {type(mode)}"
        raise TypeError(msg)

    if mode in {CollapseMode.FILAMENT, CollapseMode.AMBIPOLAR}:
        vr = np.array(
            [
                _vrfit(r, _rminfit(t, mode), _vminfit(t, mode), _avfit(t, mode), mode)
                for t, r in zip(t_yr, r_pc, strict=False)
            ]
        )
    else:
        # BE-sphere modes: approximate via finite differences of parcel_radius.
        # This is NOT the relationship used to generate the model.
        t_s = t_yr * _SECONDS_PER_YEAR
        r_cm = r_pc * _PC
        vr = np.gradient(r_cm, t_s)

    return pd.Series(vr, index=t_yr, name="radial_velocity_cm_s")


@enum.verify(enum.UNIQUE)
class SuccessFlag(enum.IntEnum):
    """SuccessFlag indicates whether the model failed or ran successfully,.

    and if it failed how.

    """

    def __new__(cls, value: int, docstring: str) -> Self:  # ruff: ignore[undocumented-public-method]
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
    CONSERVATION_ERROR = -16, "Runtime element conservation tolerance exceeded."
    ZERO_INNER_RADIUS_ERROR = -17, "rin must be > 0 when enable_radiative_transfer=True."

    def check_error(
        self, *, only_error: bool = False, raise_on_error: bool = True
    ) -> str | None:
        """Convert the UCLCHEM integer result flag to a message explaining what went wrong.

        Parameters
        ----------
        only_error : bool
            If True, skip printing if the model ran successfully, and only
            error out if it did not. Default = False.
        raise_on_error : bool
            If True, raises RuntimeError if the `self` is not
            `SuccessFlag.SUCCESS`. If False, returns the message string. Default = True.

        Returns
        -------
        str | None
            Error message, or ``None`` if the model ran successfully.

        Raises
        ------
        RuntimeError
            If `raise_on_error` is True.

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
            SuccessFlag.CONSERVATION_ERROR: "Runtime conservation check failed: an element changed by more than runtime_conservation_tolerance. This usually indicates a severe integrator error or network bug. Set runtime_conservation_tolerance to a negative value to disable the check.",
            SuccessFlag.ZERO_INNER_RADIUS_ERROR: "rin must be > 0 when enable_radiative_transfer=True with multiple points. The innermost parcel would be placed at r=0, causing division by zero in G0_internal_at_r. Set rin to a small positive value (e.g. rout/100).",
        }

        msg = error_msg_dict[self]
        if raise_on_error:
            error_msg = f"UCLCHEM error (code {self.name}, {self.value}): {msg}"
            raise RuntimeError(error_msg)
        return msg


# -------------------------------------------------------------------------------
# 1D hot core helper routines:
# -------------------------------------------------------------------------------

L_SUN: float = 3.828e26  # W — IAU 2015 nominal solar luminosity


class TempMode(enum.Enum):
    """Interpolation scheme for :func:`get_protostellar_Teff`."""

    BINS = "bins"  # discrete stepwise lookup,   1 to 1e6  L_sun
    SMOOTH = "smooth"  # global power-law fit,       1 to 1e6  L_sun
    SMOOTH_LOW = "smooth_low"  # power-law fit to bins 1 to 3,  1 to 1e4  L_sun
    SMOOTH_HIGH = "smooth_high"  # power-law fit to bins 4 to 7,  1e3 to 1e6 L_sun
    SMOOTH_CUSTOM = "smooth_custom"  # user-provided fit parameters;


# (L_lower [L_sun], T_eff [K]) — ordered high to low; first match wins
_BINS: tuple[tuple[float, float], ...] = (
    (5.0e5, 40_000.0),
    (1.0e5, 30_000.0),
    (1.0e4, 25_000.0),
    (1.0e3, 20_000.0),
    (1.0e2, 10_000.0),
    (1.0e1, 8_000.0),
    (1.0e0, 6_000.0),
)

# Least-squares power-law fits in log-log space on bin geometric midpoints.
# Fit form: T_eff = T0 * (L_star / L_sun) ** alpha
_FIT_SMOOTH = (4_788.72, 0.156034)  # all 7 bins,  valid 1 to 1e6  L_sun
_FIT_SMOOTH_LOW = (5_337.78, 0.110924)  # bins 1 to 3,    valid 1 to 1e4  L_sun
_FIT_SMOOTH_HIGH = (7_343.74, 0.120552)  # bins 4 to 7,    valid 1e3 to 1e6 L_sun

# Valid luminosity ranges [L_sun] per mode
_RANGES: dict[TempMode, tuple[float, float]] = {
    TempMode.BINS: (1.0, 1.0e6),
    TempMode.SMOOTH: (1.0, 1.0e6),
    TempMode.SMOOTH_LOW: (1.0, 1.0e4),
    TempMode.SMOOTH_HIGH: (1.0e3, 1.0e6),
    TempMode.SMOOTH_CUSTOM: (
        0.0,
        float("inf"),
    ),  # no limits; user must ensure fit is valid for their L_star
}


def get_protostellar_Teff(
    L_star: float,
    *,
    mode: TempMode = TempMode.BINS,
    custom_T0: float | None = None,
    custom_alpha: float | None = None,
    L_star_unit: str = "L_sun",
) -> float:
    """Return the effective temperature [K] for a protostellar object.

    Parameters
    ----------
    L_star : float
        Stellar luminosity. Valid range depends on *mode*; see :data:`_RANGES`.
    mode : TempMode
        Interpolation scheme:

        - :attr:`TempMode.SMOOTH` (default) — global power-law fit to all
          seven bins, valid 1 to 1e6 L_sun.
        - :attr:`TempMode.BINS` — discrete stepwise lookup, valid 1 to 1e6 L_sun.
        - :attr:`TempMode.SMOOTH_LOW` — power-law fit to bins 1 to 3,
          valid 1 to 1e4 L_sun.
        - :attr:`TempMode.SMOOTH_HIGH` — power-law fit to bins 4 to 7,
          valid 1e3 to 1e6 L_sun.

        Defaults to ``TempMode.BINS``.
    custom_T0 : float | None
        Temperature constant T0 in the power-law fit T = T0 * L^alpha, used only when
        mode is ``TempMode.SMOOTH_CUSTOM``. Defaults to ``None``.
    custom_alpha : float | None
        Power-law exponent alpha in T = T0 * L^alpha, used only when mode is
        ``TempMode.SMOOTH_CUSTOM``. Defaults to ``None``.
    L_star_unit : str
        Unit of ``L_star``. Either ``'L_sun'`` or ``'W'``. Defaults to ``'L_sun'``.

    Returns
    -------
    float
        Effective temperature in Kelvin.

    Raises
    ------
    TypeError
        If ``L_star`` is not a real number.
    ValueError
        If ``L_star`` is outside the valid range for the chosen mode, or if custom
        fit parameters are provided but mode is not ``TempMode.SMOOTH_CUSTOM``.
    RuntimeError
        If no bin matches ``L_star`` in ``TempMode.BINS`` mode (indicates a bug).

    Examples
    --------
    >>> get_protostellar_Teff(1.0, mode=TempMode.SMOOTH)    # 1 L_sun, smooth
    4788.72
    >>> get_protostellar_Teff(1000.0, mode=TempMode.BINS)   # 1000 L_sun, bins
    20000.0

    """
    if not isinstance(L_star, (int, float)):
        msg = f"L_star must be a real number, got {type(L_star).__name__!r}"
        raise TypeError(msg)

    if (
        custom_T0 is not None or custom_alpha is not None
    ) and mode != TempMode.SMOOTH_CUSTOM:
        msg = "Custom fit parameters provided but mode is not TempMode.SMOOTH_CUSTOM."
        raise ValueError(msg)

    if L_star_unit not in {"L_sun", "W"}:
        msg = "L_star_unit must be 'L_sun' or 'W'."
        raise ValueError(msg)

    # Work in L_sun units throughout
    l_star_lsun: float = L_star / L_SUN if L_star_unit == "W" else L_star
    L_min, L_max = _RANGES[mode]

    if not (L_min <= l_star_lsun <= L_max):
        msg = (
            f"L_star = {l_star_lsun:.3e} W ({l_star_lsun:.3e} L_sun) is outside "
            f"the valid range [{L_min:.0e}, {L_max:.0e}] L_sun for {mode.value!r}."
        )
        raise ValueError(msg)

    # Return binned values
    if mode is TempMode.BINS:
        for L_lower, T_eff in _BINS:
            if l_star_lsun >= L_lower:
                return T_eff
        msg = f"No bin matched L_star={l_star_lsun!r} — this is a bug."
        raise RuntimeError(msg)

    # Return power-law fit values
    if mode is TempMode.SMOOTH_CUSTOM:
        if custom_T0 is None or custom_alpha is None:
            msg = "mode is TempMode.SMOOTH_CUSTOM but custom_T0 or custom_alpha is None."
            raise ValueError(msg)
        return custom_T0 * (l_star_lsun**custom_alpha)

    T0, alpha = {
        TempMode.SMOOTH: _FIT_SMOOTH,
        TempMode.SMOOTH_LOW: _FIT_SMOOTH_LOW,
        TempMode.SMOOTH_HIGH: _FIT_SMOOTH_HIGH,
    }[mode]
    return T0 * (l_star_lsun**alpha)


def get_protostellar_model_index(L_star: float) -> int:
    """Obtain the right model index for protostellar core models in 1D mode based on the.

    stellar luminosity

    Parameters
    ----------
    L_star : float
        Stellar luminosity in units of solar luminosity (L_sun)

    Returns
    -------
    int
        The model index for the protostellar core model.

    Raises
    ------
    ValueError
        If the stellar luminosity is below the minimum valid value.

    """
    if L_star >= 1.0e6:  # ruff: ignore[magic-value-comparison]
        return 6
    if L_star >= 1.0e5:  # ruff: ignore[magic-value-comparison]
        return 5
    if L_star >= 1.0e4:  # ruff: ignore[magic-value-comparison]
        return 4
    if L_star >= 1.0e3:  # ruff: ignore[magic-value-comparison]
        return 3
    if L_star >= 1.0e2:  # ruff: ignore[magic-value-comparison]
        return 2
    if L_star >= 1.0e1:  # ruff: ignore[magic-value-comparison]
        return 1
    if L_star >= 1.0e0:
        return 0
    msg = f"L_star = {L_star:.3e} W is below the minimum of 1 L_sun = {L_SUN:.3e} W."
    raise ValueError(msg)


def remove_keys_with_none(dct: dict[str, Any]) -> None:
    """Remove keys with None values from the dictionary, modified in-place.

    Parameters
    ----------
    dct : dict[str, Any]
        Dictionary, modified in-place

    Examples
    --------
    >>> dct = {'key_a': 'a', 'key_b': None}
    >>> remove_keys_with_none(dct)
    >>> dct
    {'key_a': 'a'}

    """
    keys_to_delete = [k for k, v in dct.items() if v is None]
    for k in keys_to_delete:
        del dct[k]


def get_lowercase_copy(dct: dict[str, Any]) -> dict[str, Any]:
    """Get a copy of a dictionary where the keys have been replaced by lowercase keys.

    Parameters
    ----------
    dct : dict[str, Any]
        Dictionary to get lowercase copy of.

    Returns
    -------
    lowercase_dct : dict[str, Any]
        Dictionary with same values as `dct`, but lowercase keys.

    Raises
    ------
    ValueError
        If duplicate keys are found after lowercasing, e.g. the original dictionary
        has key ``"a"`` and ``"A"``.

    Examples
    --------
    >>> get_lowercase_copy({'KEY_A': 'A', 'kEy_B': 'b'})
    {'key_a': 'A', 'key_b': 'b'}

    """
    new_dct = {}
    for k, v in dct.items():
        if k.lower() in new_dct:
            msg = f"Duplicate lower case key {k} is already in the dict, stopping"
            raise ValueError(msg)
        new_dct[k.lower()] = v
    return new_dct


def pad_to_length(array: np.ndarray, length: int, **kwargs) -> np.ndarray:
    """Pad an array to a specific length.

    Parameters
    ----------
    array : np.ndarray
        Array to pad
    length : int
        Padded length
    **kwargs : dict[str, Any]
        Keyword arguments passed to :func:`np.pad`.

    Returns
    -------
    np.ndarray
        Array padded to a length of `length`.

    Raises
    ------
    ValueError
        If the length of `array` is larger than `length`.

    Examples
    --------
    >>> pad_to_length(np.array([1, 2]), 4)
    array([1, 2, 0, 0])
    >>> # Takes same dtype as input array (by default)
    >>> pad_to_length(np.array([1., 2.]), 4)
    array([1., 2., 0., 0.])
    >>> # Custom padding values can be specified by passing a kwarg
    >>> pad_to_length(np.array([1., 2.]), 4, mode='symmetric')
    array([1., 2., 2., 1.])

    """
    if len(array) > length:
        msg = f"Tried padding array, but already had length ({len(array)}) larger than desired padded length ({length})"
        raise ValueError(msg)
    return np.pad(array, (0, length - len(array)), **kwargs)
