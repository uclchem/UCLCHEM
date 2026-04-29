"""Functional API for UCLCHEM models.

This module provides a functional interface to UCLCHEM's physics models, designed for backward
compatibility with UCLCHEM v2 while leveraging the modern object-oriented architecture.

The functional API offers convenience functions that wrap the model classes (:mod:`uclchem.model`)
with a simpler calling convention. Each function follows a consistent pattern:

- **Input**: parameter dictionary, output species list, and return format flags
- **Output**: Either file output (default) or in-memory arrays/DataFrames
- **Modes**: Choose between disk I/O or memory-based workflows

Available Models
----------------
- :func:`cloud` - Static cloud models
- :func:`collapse` - Collapsing cloud models (BE, filament, ambipolar)
- :func:`prestellar_core` (alias :func:`hot_core`) - Prestellar core with temperature evolution
- :func:`cshock` - C-type shock models
- :func:`jshock` - J-type shock models

Usage Patterns
--------------

**Disk mode (default)** - Results written to files:

    >>> import uclchem
    >>>
    >>> success_flag, out_species_abundances = uclchem.functional.cloud(
    ...     param_dict={'initialDens': 1e4, 'outputFile': 'cloud.dat'},
    ...     out_species=['CO', 'H2O']
    ... )
    >>>
    >>> success_flag.check_error()
    Model ran successfully
    >>>
    >>> # Print abundances of CO and H2O at the end of the model
    >>> print(out_species_abundances)
    [...]

**Memory mode** - Results returned as arrays:

    >>> phys, chem, rate_constants, heat, final_abun, success_flag = uclchem.functional.cloud(
    ...     param_dict={'initialDens': 1e4},
    ...     out_species=['CO', 'H2O'],
    ...     return_array=True,
    ...     return_rate_constants=True,
    ...     return_heating=True
    ... )

**DataFrame mode** - Results as pd.DataFrames:

    >>> (
    ...     phys_df,
    ...     chem_df,
    ...     rate_constants_df,
    ...     heat_df,
    ...     final_abun,
    ...     success_flag,
    ... ) = uclchem.functional.cloud(
    ...     param_dict={'initialDens': 1e4},
    ...     out_species=['CO'],
    ...     return_dataframe=True,
    ...     return_rate_constants=True
    ... )

.. note::
   You cannot mix file I/O parameters (``outputFile``, ``abundSaveFile``) with memory return
   modes (``return_array``, ``return_dataframe``). Choose one approach per run.

.. tip::
   For interactive work and complex analysis, the object-oriented API (:mod:`uclchem.model`)
   offers more flexibility. Use the functional API when you need backward compatibility or
   prefer a simpler interface.

Seealso
-------
- :mod:`uclchem.model` - Object-oriented model classes
- :mod:`uclchem.utils` - Utility functions including ``check_error()``
- :doc:`/tutorials/index` - Interactive tutorials for both APIs
"""

from typing import Any, Literal

import numpy as np
import pandas as pd

from uclchem.constants import TIMEPOINTS, default_elements_to_check
from uclchem.model import AbstractModel, Cloud, Collapse, CShock, JShock, PrestellarCore


def _validate_functional_api_params_(
    param_dict: dict | None,
    return_array: bool,
    return_dataframe: bool,
    return_rate_constants: bool,
    return_heating: bool,
    starting_chemistry: np.ndarray | None,  # noqa: ARG001
    return_stats: bool = False,
) -> None:
    """Validate functional API specific constraints.

    Checks that return_* parameters are not mixed with file parameters.

    Parameters
    ----------
    param_dict : dict | None
        The parameter dictionary
    return_array : bool
        Whether arrays are being returned
    return_dataframe : bool
        Whether DataFrames are being returned
    return_rate_constants : bool
        Whether rate constants are being returned
    return_heating : bool
        Whether heating arrays are being returned
    starting_chemistry : np.ndarray | None
        Starting chemistry array if provided
    return_stats : bool
        Whether DVODE statistics should be returned. Default = False.

    Raises
    ------
    RuntimeError
        If file parameters are mixed with memory return parameters
    """
    # Determine if this is a memory return request (user wants data returned, not written)
    memory_return_requested = (
        return_array
        or return_dataframe
        or return_rate_constants
        or return_heating
        or return_stats
    )

    # Check file parameter mixing with memory return parameters
    # (can't both return data AND write to files - user must choose one)
    if memory_return_requested:
        file_params = ["outputFile", "abundSaveFile", "abundLoadFile", "columnFile"]
        if param_dict is not None and any(k in param_dict for k in file_params):
            msg = (
                "return_array or return_dataframe cannot be used if any output or input file is specified. "
                "These parameters are mutually exclusive: use either file I/O (outputFile, abundSaveFile, etc.) "
                "OR in-memory returns (return_array, return_dataframe), but not both."
            )
            raise RuntimeError(msg)


def _functional_return_(
    model_object: AbstractModel,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
) -> tuple[Any, ...]:
    """Take in the object that was modeled and return the values based on the specified booleans.

    Parameters
    ----------
    model_object : AbstractModel
        model_object of a class that inherited from AbstractModel,
        from which the results should be returned.
    return_array : bool
        A boolean on whether a np.array should be returned to a user.
        If both return_array and return_dataframe are false, the function will return
        the success_flag, dissipation_time if the model_object has that attribute,
        and the final abundances of the out_species. Default = False.
    return_dataframe : bool
        Whether a pd.DataFrame should be returned to a user.
        If both return_array and return_dataframe are false, the function will return
        the success_flag, dissipation_time if the model_object has that attribute,
        and the final abundances of the out_species. Default = False.
    return_rate_constants : bool
        A boolean on whether the reaction rate constants
        should be returned to a user. Default = False.
    return_heating : bool
        A boolean on whether the heating/cooling rates
        should be returned to a user. Default = False.
    return_stats : bool
        Whether DVODE statistics should be returned. Default = False.

    Returns
    -------
    tuple[Any, ...]
        tuple with entries depending on which bools are True.
    """
    result: list[Any]
    if return_dataframe:
        # If multiple spatial points are present, return DataFrames concatenated across points
        points = model_object._param_dict.get("points", 1)
        rate_constants_df = None
        heating_df = None
        stats_df = None

        if points > 1:
            physics_list = []
            chemistry_list = []
            rate_constants_list = []
            heating_list = []
            stats_list = []
            for pt in range(points):
                res = model_object.get_dataframes(
                    point=pt,
                    with_rate_constants=return_rate_constants,
                    with_heating=return_heating,
                    with_stats=return_stats,
                )
                phys = res[0].copy()
                chem = res[1].copy()
                phys["Point"] = pt + 1
                chem["Point"] = pt + 1
                physics_list.append(phys)
                chemistry_list.append(chem)
                idx = 2
                if return_rate_constants and len(res) > idx:
                    rate_constants_list.append(res[idx].assign(Point=pt + 1))
                    idx += 1
                if return_heating and len(res) > idx:
                    heating_list.append(res[idx].assign(Point=pt + 1))
                    idx += 1
                if return_stats and len(res) > idx:
                    stats_list.append(res[idx].assign(Point=pt + 1))

            phys_df = pd.concat(physics_list, ignore_index=True)
            chem_df = pd.concat(chemistry_list, ignore_index=True)
            rate_constants_df = (
                pd.concat(rate_constants_list, ignore_index=True)
                if rate_constants_list
                else None
            )
            heating_df = (
                pd.concat(heating_list, ignore_index=True) if heating_list else None
            )
            stats_df = pd.concat(stats_list, ignore_index=True) if stats_list else None
        else:
            # Single point: behave as before but include a Point column
            result_dfs = model_object.get_dataframes(
                with_rate_constants=return_rate_constants,
                with_heating=return_heating,
                with_stats=return_stats,
            )
            phys_df = result_dfs[0]  # type: ignore[assignment]
            chem_df = result_dfs[1]  # type: ignore[assignment]
            idx = 2
            if return_rate_constants and len(result_dfs) > idx:
                rate_constants_df = result_dfs[idx]  # type: ignore[assignment]
                idx += 1
            if return_heating and len(result_dfs) > idx:
                heating_df = result_dfs[idx]  # type: ignore[assignment]
                idx += 1
            if return_stats and len(result_dfs) > idx:
                stats_df = result_dfs[idx]  # type: ignore[assignment]
            phys_df["Point"] = 1
            chem_df["Point"] = 1

        # Build result tuple - stats only included when requested (for backward compatibility)
        result = [phys_df, chem_df, rate_constants_df, heating_df]
        if return_stats:
            result.append(stats_df)
        if hasattr(model_object, "dissipation_time"):
            result.append(model_object.dissipation_time)
        result.append(model_object.next_starting_chemistry_array)
        result.append(model_object.success_flag)
        return tuple(result)
    elif return_array:
        # Build result tuple - stats only included when requested (for backward compatibility)
        result = [
            model_object.physics_array,
            model_object.chemical_abun_array,
            model_object.rate_constants_array if return_rate_constants else None,
            model_object.heat_array if return_heating else None,
        ]
        if return_stats:
            result.append(model_object.stats_array)
        if hasattr(model_object, "dissipation_time"):
            result.append(model_object.dissipation_time)
        result.append(model_object.next_starting_chemistry_array)
        result.append(model_object.success_flag)
        return tuple(result)
    # Disk mode with file output
    # FIXED format: [success_flag, abundances] OR [success_flag, dissipation_time, abundances]
    elif hasattr(model_object, "dissipation_time"):
        return_value: tuple[Any, ...] = (
            model_object.success_flag,
            model_object.dissipation_time,
        )

        return return_value + (model_object.out_species_abundances_array,)
    else:
        return_value = (model_object.success_flag,)
        return return_value + (model_object.out_species_abundances_array,)


def cloud(
    param_dict: dict | None = None,
    out_species: list[str] | None = None,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray | None = None,
    timepoints: int = TIMEPOINTS,
) -> tuple[Any, ...]:
    """Run cloud model from UCLCHEM.

    Parameters
    ----------
    param_dict : dict | None
        A dictionary of parameters where keys are any of the variables in
        ``defaultparameters.f90`` and values are value for current run. Default = None.
    out_species : list[str] | None
        A list of species for which final abundance will be returned.
        If None, no abundances will be returned.
        Defaults to ``uclchem.constants.default_elements_to_check``.
    return_array : bool
        A boolean on whether a np.array should be returned to a user.
        If both return_array and return_dataframe are false,
        this function will default to writing outputs to a file. Default = False.
    return_dataframe : bool
        A boolean on whether a pd.DataFrame should be returned to a user.
        If both return_array and return_dataframe are false,
        this function will default to writing outputs to a file. Default = False.
    return_rate_constants : bool
        A boolean on whether the reaction rate constants should be
        returned to a user. Default = False.
    return_heating : bool
        A boolean on whether the heating/cooling arrays should
        be returned to a user. Default = False.
    return_stats : bool
        Whether DVODE statistics should be returned. Default = False.
    starting_chemistry : np.ndarray | None
        Array containing the starting chemical abundances
        needed by uclchem. If None, use default initial abundanes. Default = None.
    timepoints : int
        Integer value of how many timesteps should be calculated before
        aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.

    Returns
    -------
    tuple[Any, ...]
        tuple with entries depending on which bools are True.
    """
    if out_species is None:
        out_species = default_elements_to_check

    # Validate functional API constraints
    _validate_functional_api_params_(
        param_dict=param_dict,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        starting_chemistry=starting_chemistry,
        return_stats=return_stats,
    )

    model_object = Cloud(
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints,
    )

    return _functional_return_(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def collapse(
    collapse: Literal["BE1.1", "BE4", "filament", "ambipolar"],
    param_dict: dict | None = None,
    out_species: list[str] | None = None,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray | None = None,
    timepoints: int = TIMEPOINTS,
) -> tuple[Any, ...]:
    """Run collapse model from UCLCHEM based on Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract).

    Parameters
    ----------
    collapse : Literal['BE1.1', 'BE4', 'filament', 'ambipolar']
        A string containing the collapse
        type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'.
    param_dict : dict | None
        A dictionary of parameters where keys are any of the variables in
        ``defaultparameters.f90`` and values are value for current run. Default = None.
    out_species : list[str] | None
        A list of species for which final abundance will be returned.
        If None, no abundances will be returned.
        Defaults to ``uclchem.constants.default_elements_to_check``.
    return_array : bool
        A boolean on whether a np.array should be returned to a user.
        If both return_array and return_dataframe are false,
        this function will default to writing outputs to a file. Default = False.
    return_dataframe : bool
        A boolean on whether a pd.DataFrame should be returned to a user.
        If both return_array and return_dataframe are false,
        this function will default to writing outputs to a file. Default = False.
    return_rate_constants : bool
        A boolean on whether the reaction rate constants should
        be returned to a user. Default = False.
    return_heating : bool
        A boolean on whether the heating/cooling arrays should be
        returned to a user. Default = False.
    return_stats : bool
        Whether DVODE statistics should be returned. Default = False.
    starting_chemistry : np.ndarray | None
        Array containing the starting chemical abundances
        needed by UCLCHEM. If None, uses default initial abundances. Default = None.
    timepoints : int
        Integer value of how many timesteps should be calculated before aborting
        the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``

    Returns
    -------
    tuple[Any, ...]
        tuple with entries depending on which bools are True.
    """
    if out_species is None:
        out_species = default_elements_to_check

    _validate_functional_api_params_(
        param_dict,
        return_array,
        return_dataframe,
        return_rate_constants,
        return_heating,
        starting_chemistry,
        return_stats=return_stats,
    )

    model_object = Collapse(
        collapse=collapse,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints,
    )

    return _functional_return_(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def prestellar_core(
    temp_index: int = 1,
    max_temperature: float = 300.0,
    param_dict: dict | None = None,
    out_species: list[str] | None = None,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray | None = None,
    timepoints: int = TIMEPOINTS,
) -> tuple[Any, ...]:
    """Run prestellar core model from UCLCHEM, based on Viti et al. 2004 and Collings et al. 2004.

    This model type was previously known as hot core.

    Parameters
    ----------
    temp_index : int
        Used to select the mass of prestellar core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        Default = 1.
    max_temperature : float
        Value at which gas temperature will stop increasing.
        Default = 300.0.
    param_dict : dict | None
        A dictionary of parameters where keys are any of the variables
        in ``defaultparameters.f90`` and values are value for current run.
        Default = None.
    out_species : list[str] | None
        A list of species for which final abundance will be returned.
        If None, no abundances will be returned.
        Defaults to ``uclchem.constants.default_elements_to_check``.
    return_array : bool
        A boolean on whether a np.array should be returned to a user.
        If both return_array and return_dataframe are false,
        this function will default to writing outputs to a file. Default = False.
    return_dataframe : bool
        A boolean on whether a pd.DataFrame should be returned to a user.
        If both return_array and return_dataframe are false,
        this function will default to writing outputs to a file. Default = False.
    return_rate_constants : bool
        A boolean on whether the reaction rate constants
        should be returned to a user. Default = False.
    return_heating : bool
        A boolean on whether the heating/cooling arrays should be
        returned to a user. Default = False.
    return_stats : bool
        Whether DVODE statistics should be returned. Default = False.
    starting_chemistry : np.ndarray | None
        Array containing the starting chemical abundances
        needed by uclchem. If None, use default initial abundances. Default = None.
    timepoints : int
        Integer value of how many timesteps should be calculated before aborting
        the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.

    Returns
    -------
    tuple[Any, ...]
        tuple with entries depending on which bools are True.
    """
    if out_species is None:
        out_species = default_elements_to_check

    _validate_functional_api_params_(
        param_dict,
        return_array,
        return_dataframe,
        return_rate_constants,
        return_heating,
        starting_chemistry,
        return_stats=return_stats,
    )

    model_object = PrestellarCore(
        temp_index=temp_index,
        max_temperature=max_temperature,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints,
    )

    return _functional_return_(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def cshock(
    shock_vel: float,
    timestep_factor: float = 0.01,
    minimum_temperature: float = 0.0,
    param_dict: dict | None = None,
    out_species: list[str] | None = default_elements_to_check,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray | None = None,
    timepoints: int = TIMEPOINTS,
) -> tuple[Any, ...]:
    """Run C-type shock model from UCLCHEM.

    Parameters
    ----------
    shock_vel : float
        Velocity of the shock in km/s
    timestep_factor : float
        Whilst the time is less than 2 times the dissipation time of shock,
        timestep is timestep_factor*dissipation time. Essentially controls
        how well resolved the shock is in your model. Defaults to 0.01.
    minimum_temperature : float
        Minimum post-shock temperature. Defaults to 0.0 (no minimum).
        The shocked gas typically cools to ``initialTemp`` if this is not set.
    param_dict : dict | None
        A dictionary of parameters where keys are any of the variables
        in ``defaultparameters.f90`` and values are value for current run.
        Default = None.
    out_species : list[str] | None
        A list of species for which final
        abundance will be returned. If None, no abundances will be returned.
        Default = ``uclchem.constants.default_elements_to_check``.
    return_array : bool
        Whether a np.array should be returned.
        If both return_array and return_dataframe are false,
        this function will default to writing outputs to a file. Default = False.
    return_dataframe : bool
        Whether a pd.DataFrame should be returned.
        If both return_array and return_dataframe are False,
        this function will default to writing outputs to a file. Default = False.
    return_rate_constants : bool
        Whether the reaction rate constants should be returned.
        Default = False.
    return_heating : bool
        Whether the heating/cooling arrays should be returned to a user.
        Default = False.
    return_stats : bool
        Whether DVODE statistics should be returned. Default = False.
    starting_chemistry : np.ndarray | None
        np.ndarray containing the starting chemical abundances
        needed by UCLCHEM. If None, use default initial abundances. Default = None.
    timepoints : int
        Integer value of how many timesteps should be calculated before
        aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.

    Returns
    -------
    tuple[Any, ...]
        tuple with entries depending on which bools are True.
    """
    if out_species is None:
        out_species = default_elements_to_check

    _validate_functional_api_params_(
        param_dict,
        return_array,
        return_dataframe,
        return_rate_constants,
        return_heating,
        starting_chemistry,
        return_stats=return_stats,
    )

    model_object = CShock(
        shock_vel=shock_vel,
        timestep_factor=timestep_factor,
        minimum_temperature=minimum_temperature,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints,
    )

    return _functional_return_(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def jshock(
    shock_vel: float,
    param_dict: dict | None = None,
    out_species: list[str] | None = None,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray | None = None,
    timepoints: int = TIMEPOINTS,
) -> tuple[Any, ...]:
    """Run J-type shock model from UCLCHEM.

    Parameters
    ----------
    shock_vel : float
        Velocity of the shock
    param_dict : dict | None
        A dictionary of parameters where keys are any of the variables in
        ``defaultparameters.f90`` and values are value for current run.
    out_species : list[str] | None
        A list of species for which final abundance will be returned.
        If None, no abundances will be returned.
        Defaults to ``uclchem.constants.default_elements_to_check``.
    return_array : bool
        A boolean on whether a np.array should be returned to a user.
        If both return_array and return_dataframe are false, this function will default
        to writing outputs to a file
    return_dataframe : bool
        A boolean on whether a pd.DataFrame should be
        returned to a user. If both return_array and return_dataframe are False,
        this function will default to writing outputs to a file. Default = False.
    return_rate_constants : bool
        A boolean on whether the reaction rate constants
        should be returned to a user. Default = False.
    return_heating : bool
        A boolean on whether the heating/cooling arrays
        should be returned to a user. Default = False.
    return_stats : bool
        Whether to return DVODE stats. Default = False.
    starting_chemistry : np.ndarray | None
        np.ndarray containing the starting
        chemical abundances needed by UCLCHEM. Default = None.
    timepoints : int
        Integer value of how many timesteps should be calculated
        before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

    Returns
    -------
    tuple[Any, ...]
        tuple with entries depending on which bools are True.
    """
    if out_species is None:
        out_species = default_elements_to_check

    _validate_functional_api_params_(
        param_dict,
        return_array,
        return_dataframe,
        return_rate_constants,
        return_heating,
        starting_chemistry,
        return_stats=return_stats,
    )

    model_object = JShock(
        shock_vel=shock_vel,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints,
    )

    return _functional_return_(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


hot_core = prestellar_core  # Alias for backward compatibility
