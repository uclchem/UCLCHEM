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
    >>> success_flag, abundance_CO, abundance_H2O = uclchem.functional.cloud(
    ...     param_dict={'initialDens': 1e4, 'outputFile': 'cloud.dat'},
    ...     out_species=['CO', 'H2O']
    ... )
    >>>
    >>> success_flag.check_error()
    Model ran successfully

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

See Also:
--------
- :mod:`uclchem.model` - Object-oriented model classes
- :mod:`uclchem.utils` - Utility functions including ``check_error()``
- :doc:`/tutorials/index` - Interactive tutorials for both APIs

"""

import numpy as np
import pandas as pd

from uclchem.constants import TIMEPOINTS, default_elements_to_check
from uclchem.model import AbstractModel, Cloud, Collapse, CShock, JShock, PrestellarCore


def __validate_functional_api_params__(
    param_dict: dict,
    return_array: bool,
    return_dataframe: bool,
    return_rate_constants: bool,
    return_heating: bool,
    starting_chemistry: np.ndarray,  # noqa: ARG001
    return_stats: bool = False,
) -> None:
    """Validate functional API specific constraints.
    Checks that return_* parameters are not mixed with file parameters.

    Args:
        param_dict: The parameter dictionary
        return_array: Whether arrays are being returned
        return_dataframe: Whether DataFrames are being returned
        return_rate_constants: Whether rate constants are being returned
        return_heating: Whether heating arrays are being returned
        starting_chemistry: Starting chemistry array if provided
        return_stats (bool): Whether DVODE statistics should be returned. Default = False.

    Raises:
        RuntimeError: If file parameters are mixed with memory return parameters

    Note:
        The system always uses the memory interface internally. File writing
        is controlled by the presence of outputFile, abundSaveFile, etc.
        This validation ensures users don't request both data return AND file writing.

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
            raise RuntimeError(
                "return_array or return_dataframe cannot be used if any output or input file is specified. "
                "These parameters are mutually exclusive: use either file I/O (outputFile, abundSaveFile, etc.) "
                "OR in-memory returns (return_array, return_dataframe), but not both."
            )


def __functional_return__(
    model_object: AbstractModel,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
) -> tuple:
    """Return function that takes in the object that was modelled and returns the values
    based on the specified booleans.

    Args:
        model_object (AbstractModel): model_object of a class that inherited from AbstractModel,
            from which the results should be returned.
        return_array (bool): A boolean on whether a np.array should be returned to a user.
            If both return_array and return_dataframe are false, the function will return
            the success_flag, dissipation_time if the model_object has that attribute,
            and the final abundances of the out_species.
        return_dataframe (bool): Whether a pd.DataFrame should be returned to a user.
            If both return_array and return_dataframe are false, the function will return
            the success_flag, dissipation_time if the model_object has that attribute,
            and the final abundances of the out_species.
        return_rate_constants (bool): A boolean on whether the reaction rate constants
            should be returned to a user.
        return_heating (bool): A boolean on whether the heating/cooling rates
            should be returned to a user.
        return_stats (bool): Whether DVODE statistics should be returned. Default = False.

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model
                failed to run. Can be passed to `uclchem.utils.check_error()` to see more details.
                If the model succeeded, and the model_object has the dissipation_time attribute
                the second element is the dissipation time.
                Further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (np.ndarray): array containing the physical outputs
                for each written timestep
            - chemicalAbunArray (np.ndarray): array containing the chemical abundances
                for each written timestep
            - rateConstantsArray (np.ndarray): array containing reaction rate constants
                (if return_rate_constants=True)
            - heatArray (np.ndarray): array containing heating/cooling rates
                (if return_heating=True)
            - dissipation_time (float): dissipation time in years
                (if model_object contains the dissipation_time attribute)
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last
                timestep in the format uclchem needs in order to perform an additional
                run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pd.DataFrame): DataFrame containing the physical outputs
                for each written timestep
            - chemicalDF (pd.DataFrame): DataFrame containing the chemical abundances
                for each written timestep
            - rate_constants_df (pd.DataFrame): DataFrame containing reaction rate constants
                (if return_rate_constants=True)
            - heatingDF (pd.DataFrame): DataFrame containing heating/cooling rates
                (if return_heating=True)
            - dissipation_time (float): dissipation time in years
                (if model_object contains the dissipation_time attribute)
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last timestep
                in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.

    """
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
                    joined=False,
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
                joined=False,
                with_rate_constants=return_rate_constants,
                with_heating=return_heating,
                with_stats=return_stats,
            )
            phys_df = result_dfs[0]
            chem_df = result_dfs[1]
            idx = 2
            if return_rate_constants and len(result_dfs) > idx:
                rate_constants_df = result_dfs[idx]
                idx += 1
            if return_heating and len(result_dfs) > idx:
                heating_df = result_dfs[idx]
                idx += 1
            if return_stats and len(result_dfs) > idx:
                stats_df = result_dfs[idx]
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
    else:
        # Disk mode with file output
        # FIXED format: [success_flag, abundances] OR [success_flag, dissipation_time, abundances]
        if hasattr(model_object, "dissipation_time"):
            return (
                model_object.success_flag,
                model_object.dissipation_time,
            ) + tuple(model_object.out_species_abundances_array)
        else:
            return (model_object.success_flag,) + tuple(
                model_object.out_species_abundances_array
            )


def __cloud__(
    param_dict: dict | None = None,
    out_species: list[str] | None = None,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray | None = None,
    timepoints: int = TIMEPOINTS,
):
    """Run cloud model from UCLCHEM.

    Args:
        param_dict (dict): A dictionary of parameters where keys are any of the variables in
            `defaultparameters.f90` and values are value for current run.
        out_species (list): A list of species for which final abundance will be returned.
            If None, no abundances will be returned.
            Defaults to `uclchem.constants.default_elements_to_check`.
        return_array (bool): A boolean on whether a np.array should be returned to a user.
            If both return_array and return_dataframe are false,
            this function will default to writing outputs to a file. Default = False.
        return_dataframe (bool): A boolean on whether a pd.DataFrame should be returned to a user.
            If both return_array and return_dataframe are false,
            this function will default to writing outputs to a file. Default = False.
        return_rate_constants (bool): A boolean on whether the reaction rate constants should be
            returned to a user.
        return_heating (bool): A boolean on whether the heating/cooling arrays should
            be returned to a user.
        return_stats (bool): Whether DVODE statistics should be returned. Default = False.
        starting_chemistry (np.ndarray): Array containing the starting chemical abundances
            needed by uclchem.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed
                to run. Can be passed to `uclchem.utils.check_error()` to see more details.
                If the `out_species` parameter is provided, the remaining elements of this list
                will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (np.ndarray): array containing the physical outputs for each written timestep
            - chemicalAbunArray (np.ndarray): array containing the chemical abundances for
                each written timestep
            - ratesArray (np.ndarray | None): array containing reaction rate constants for each timestep
                (if return_rate_constants=True)
            - heatArray (np.ndarray | None): array containing heating/cooling terms for each timestep
                (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last timestep
                in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pd.DataFrame): DataFrame containing the physical outputs
                for each written timestep
            - chemicalDF (pd.DataFrame): DataFrame containing the chemical abundances
                for each written timestep
            - ratesConstantsDF (pd.DataFrame or None): DataFrame containing reaction rate constants
                for each timestep (if return_rate_constants=True)
            - heatingDF (pd.DataFrame or None): DataFrame containing heating/cooling terms
                for each timestep (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances
                of the last timestep in the format uclchem needs in order to perform an additional
                run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.

    """
    if out_species is None:
        out_species = default_elements_to_check

    # Validate functional API constraints
    __validate_functional_api_params__(
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

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def __collapse__(
    collapse: str,
    physics_output: str,
    param_dict: dict | None = None,
    out_species: list[str] | None = None,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray | None = None,
    timepoints: int = TIMEPOINTS,
):
    """Run collapse model from UCLCHEM based on Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract).

    Args:
        collapse (str): A string containing the collapse type, options are
            'BE1.1', 'BE4', 'filament', or 'ambipolar'
        physics_output (str): Filename to store physics output, only relevant for
            'filament' and 'ambipolar' collapses. If None, no physics output will be saved.
        param_dict (dict): A dictionary of parameters where keys are any of the variables in
            `defaultparameters.f90` and values are value for current run.
        out_species (list): A list of species for which final abundance will be returned.
            If None, no abundances will be returned.
            Defaults to `uclchem.constants.default_elements_to_check`.
        return_array (bool): A boolean on whether a np.array should be returned to a user.
            If both return_array and return_dataframe are false,
            this function will default to writing outputs to a file. Default = False.
        return_dataframe (bool): A boolean on whether a pd.DataFrame should be returned to a user.
            If both return_array and return_dataframe are false,
            this function will default to writing outputs to a file. Default = False.
        return_rate_constants (bool): A boolean on whether the reaction rate constants should
            be returned to a user. Default = False.
        return_heating (bool): A boolean on whether the heating/cooling arrays should be
            returned to a user. Default = False.
        return_stats (bool): Whether DVODE statistics should be returned. Default = False.
        starting_chemistry (np.ndarray): Array containing the starting chemical abundances
            needed by UCLCHEM.
        timepoints (int): Integer value of how many timesteps should be calculated before aborting
            the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model
                failed to run. Can be passed to `uclchem.utils.check_error()` to see more details.
                If the `out_species` parameter is provided, the remaining elements of this list
                will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (np.ndarray): array containing the physical outputs for each written timestep
            - chemicalAbunArray (np.ndarray): array containing the chemical abundances
                for each written timestep
            - ratesArray (np.ndarray | None): array containing reaction rate constants for each timestep
                (if return_rate_constants=True)
            - heatArray (np.ndarray | None): array containing heating/cooling terms for each timestep
                (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last timestep
                in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pd.DataFrame): DataFrame containing the physical outputs
                for each written timestep
            - chemicalDF (pd.DataFrame): DataFrame containing the chemical abundances
                for each written timestep
            - rateConstantsDF (pd.DataFrame or None): DataFrame containing reaction rate constants
                for each timestep (if return_rate_constants=True)
            - heatingDF (pd.DataFrame or None): DataFrame containing heating/cooling terms
                for each timestep (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last timestep
                in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.

    """
    if out_species is None:
        out_species = default_elements_to_check

    __validate_functional_api_params__(
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
        physics_output=physics_output,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints,
    )

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def __prestellar_core__(
    temp_indx: int = 1,
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
):
    """Run prestellar core model from UCLCHEM, based on Viti et al. 2004 and Collings et al. 2004.
    This model type was previously known as hot core.

    Args:
        temp_indx (int): Used to select the mass of prestellar core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        max_temperature (float): Value at which gas temperature will stop increasing.
        param_dict (dict): A dictionary of parameters where keys are any of the variables
            in `defaultparameters.f90` and values are value for current run.
        out_species (list): A list of species for which final abundance will be returned.
            If None, no abundances will be returned.
            Defaults to `uclchem.constants.default_elements_to_check`.
        return_array (bool): A boolean on whether a np.array should be returned to a user.
            If both return_array and return_dataframe are false,
            this function will default to writing outputs to a file
        return_dataframe (bool): A boolean on whether a pd.DataFrame should be returned to a user.
            If both return_array and return_dataframe are false,
            this function will default to writing outputs to a file
        return_rate_constants (bool): A boolean on whether the reaction rate constants
            should be returned to a user. Default = False.
        return_heating (bool): A boolean on whether the heating/cooling arrays should be
            returned to a user.
        return_stats (bool): Whether DVODE statistics should be returned. Default = False.
        starting_chemistry (np.ndarray): Array containing the starting chemical abundances
            needed by uclchem
        timepoints (int): Integer value of how many timesteps should be calculated before aborting
            the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed
                to run and can be sent to `uclchem.utils.check_error()` to see more details.
                If the `out_species` parametere is provided, the remaining elements of this list
                will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (np.ndarray): array containing the physical outputs for each written timestep
            - chemicalAbunArray (np.ndarray): array containing the chemical abundances for each
                written timestep
            - rateConstantsArray (np.ndarray | None): array containing reaction rate constants
                for each timestep (if return_rate_constants=True)
            - heatArray (np.ndarray | None): array containing heating/cooling terms for each timestep
                (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last timestep
                in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pd.DataFrame): DataFrame containing the physical outputs for each
                written timestep
            - chemicalDF (pd.DataFrame): DataFrame containing the chemical abundances for each
                written timestep
            - rateConsantsDF (pd.DataFrame or None): DataFrame containing reaction rate constants
                for each timestep (if return_rate_constants=True)
            - heatingDF (pd.DataFrame | None): DataFrame containing heating/cooling terms for
                each timestep (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last timestep
                in the format uclchem needs in order to perform an additional run
                after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.

    """
    if out_species is None:
        out_species = default_elements_to_check

    __validate_functional_api_params__(
        param_dict,
        return_array,
        return_dataframe,
        return_rate_constants,
        return_heating,
        starting_chemistry,
        return_stats=return_stats,
    )

    model_object = PrestellarCore(
        temp_indx=temp_indx,
        max_temperature=max_temperature,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints,
    )

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def __cshock__(
    shock_vel: float,
    timestep_factor: float = 0.01,
    minimum_temperature: float = 0.0,
    param_dict: dict = None,
    out_species: list[str] | None = default_elements_to_check,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.array = None,
    timepoints: int = TIMEPOINTS,
):
    """Run C-type shock model from UCLCHEM.

    Args:
        shock_vel (float): Velocity of the shock in km/s
        timestep_factor (float): Whilst the time is less than 2 times the dissipation time of shock,
            timestep is timestep_factor*dissipation time. Essentially controls
            how well resolved the shock is in your model. Defaults to 0.01.
        minimum_temperature (float): Minimum post-shock temperature. Defaults to 0.0 (no minimum).
            The shocked gas typically cools to `initialTemp` if this is not set.
        param_dict (dict): A dictionary of parameters where keys are any of the variables
            in `defaultparameters.f90` and values are value for current run.
        out_species (list[str | None]): A list of species for which final
            abundance will be returned. If None, no abundances will be returned.
            Default = `uclchem.constants.default_elements_to_check`.
        return_array (bool): Whether a np.array should be returned.
            If both return_array and return_dataframe are false,
            this function will default to writing outputs to a file. Default = False.
        return_dataframe (bool): Whether a pd.DataFrame should be returned.
            If both return_array and return_dataframe are False,
            this function will default to writing outputs to a file. Default = False.
        return_rate_constants (bool): Whether the reaction rate constants should be returned.
        return_heating (bool): Whether the heating/cooling arrays should be returned to a user.
        return_stats (bool): Whether DVODE statistics should be returned. Default = False.
        starting_chemistry (np.ndarray): np.array containing the starting chemical abundances needed
            by UCLCHEM.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model
                failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
                If the model succeeded, the second element is the dissipation time and further
                elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (np.ndarray): array containing the physical outputs for each written
                timestep
            - chemicalAbunArray (np.ndarray): array containing the chemical abundances for each
                written timestep
            - rateConstantsArray (np.ndarray | None): array containing reaction rate constants for
                each timestep (if return_rate_constants=True)
            - heatArray (np.ndarray | None): array containing heating/cooling terms for each timestep
                (if return_heating=True)
            - disspation_time (float): dissipation time in years
            - abundanceStart (np.ndarray): array containing the chemical abundances of the last timestep
                in the format uclchem needs in order to perform an additional run
                after the initial model
            - success_flag (integer): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pd.DataFrame): DataFrame containing the physical outputs for
                each written timestep
            - chemicalDF (pd.DataFrame): DataFrame containing the chemical abundances for
                each written timestep
            - rateConstantsDF (pd.DataFrame or None): DataFrame containing reaction rate
                constants for each timestep (if return_rate_constants=True)
            - heatingDF (pd.DataFrame or None): DataFrame containing heating/cooling terms
                for each timestep (if return_heating=True)
            - disspation_time (float): dissipation time in years
            - abundanceStart (np.array): array containing the chemical abundances of the last
                timestep in the format uclchem needs in order to perform an additional run
                after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.

    """
    if out_species is None:
        out_species = default_elements_to_check

    __validate_functional_api_params__(
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

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


def __jshock__(
    shock_vel: float,
    param_dict: dict | None = None,
    out_species: list[str] | None = default_elements_to_check,
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rate_constants: bool = False,
    return_heating: bool = False,
    return_stats: bool = False,
    starting_chemistry: np.ndarray = None,
    timepoints: int = TIMEPOINTS,
):
    """Run J-type shock model from UCLCHEM.

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict | None): A dictionary of parameters where keys are any of the variables in
            `defaultparameters.f90` and values are value for current run.
        out_species (list | None): A list of species for which final abundance will be returned.
            If None, no abundances will be returned.
            Defaults to `uclchem.constants.default_elements_to_check`.
        return_array (bool): A boolean on whether a np.array should be returned to a user.
            If both return_array and return_dataframe are false, this function will default
            to writing outputs to a file
        return_dataframe (bool): A boolean on whether a pd.DataFrame should be
            returned to a user. If both return_array and return_dataframe are False,
            this function will default to writing outputs to a file. Default = False.
        return_rate_constants (bool): A boolean on whether the reaction rate constants
            should be returned to a user. Default = False.
        return_heating (bool): A boolean on whether the heating/cooling arrays
            should be returned to a user. Default = False.
        return_stats (bool): Whether to return DVODE stats. Default = False.
        starting_chemistry (np.ndarray | None): np.ndarray containing the starting
            chemical abundances needed by UCLCHEM. Default = None.
        timepoints (int): Integer value of how many timesteps should be calculated
            before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative
                if the model failed to run. Can be passed to `uclchem.utils.check_error()`
                to see more details. If the model succeeded, the second element
                is the dissipation time and further elements are the abundances
                of all species in `out_species`.
        if return_array is True:
            - physicsArray (np.ndarray): array containing the physical outputs
                for each written timestep
            - chemicalAbunArray (np.ndarray): array containing the chemical abundances
                for each written timestep
            - rateConstantsArray (np.ndarray | None): array containing reaction rate constants
                for each timestep (if return_rate_constants=True)
            - heatArray (np.ndarray | None): array containing heating/cooling terms for
                each timestep (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances of
                the last timestep in the format uclchem needs in order to perform
                an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pd.DataFrame): DataFrame containing the physical outputs for
                each written timestep
            - chemicalDF (pd.DataFrame): DataFrame containing the chemical abundances
                for each written timestep
            - rateConstantsDF (pd.DataFrame or None): DataFrame containing reaction rate
                constants for each timestep (if return_rate_constants=True)
            - heatingDF (pd.DataFrame or None): DataFrame containing heating/cooling
                terms for each timestep (if return_heating=True)
            - abundanceStart (np.ndarray): array containing the chemical abundances of
                the last timestep in the format uclchem needs in order to perform
                an additional run after the initial model
            - success_flag (int): which is negative if the model failed to run.
                Can be passed to `uclchem.utils.check_error()` to see more details.

    """
    if out_species is None:
        out_species = default_elements_to_check

    __validate_functional_api_params__(
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

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rate_constants=return_rate_constants,
        return_heating=return_heating,
        return_stats=return_stats,
    )


# Expose the functional API functions at module level
cloud = __cloud__
collapse = __collapse__
prestellar_core = __prestellar_core__
hot_core = __prestellar_core__  # Alias for backward compatibility
cshock = __cshock__
jshock = __jshock__
