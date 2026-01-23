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

**Disk mode (default)** - Results written to files::

    import uclchem

    result = uclchem.functional.cloud(
        param_dict={'initialDens': 1e4, 'outputFile': 'cloud.dat'},
        out_species=['CO', 'H2O']
    )
    # Returns: (success_flag, abundance_CO, abundance_H2O)

**Memory mode** - Results returned as arrays::

    phys, chem, rates, heat, final_abun, flag = uclchem.functional.cloud(
        param_dict={'initialDens': 1e4},
        out_species=['CO', 'H2O'],
        return_array=True,
        return_rates=True,
        return_heating=True
    )

**DataFrame mode** - Results as pandas DataFrames::

    phys_df, chem_df, rates_df, heat_df, final_abun, flag = uclchem.functional.cloud(
        param_dict={'initialDens': 1e4},
        out_species=['CO'],
        return_dataframe=True,
        return_rates=True
    )

.. note::
   You cannot mix file I/O parameters (``outputFile``, ``abundSaveFile``) with memory return
   modes (``return_array``, ``return_dataframe``). Choose one approach per run.

.. tip::
   For interactive work and complex analysis, the object-oriented API (:mod:`uclchem.model`)
   offers more flexibility. Use the functional API when you need backward compatibility or
   prefer a simpler interface.

See Also
--------
- :mod:`uclchem.model` - Object-oriented model classes
- :mod:`uclchem.utils` - Utility functions including ``check_error()``
- :doc:`/tutorials/index` - Interactive tutorials for both APIs
"""

import numpy as np

from uclchem.constants import TIMEPOINTS
from uclchem.model import (
    AbstractModel,
    Cloud,
    Collapse,
    CShock,
    JShock,
    PrestellarCore,
)


def __validate_functional_api_params__(
    param_dict: dict,
    return_array: bool,
    return_dataframe: bool,
    return_rates: bool,
    return_heating: bool,
    starting_chemistry: np.ndarray,
):
    """
    Validate functional API specific constraints.
    Checks that return_* parameters are not mixed with file parameters.

    Args:
        param_dict: The parameter dictionary
        return_array: Whether arrays are being returned
        return_dataframe: Whether DataFrames are being returned
        return_rates: Whether rates are being returned
        return_heating: Whether heating arrays are being returned
        starting_chemistry: Starting chemistry array if provided

    Raises:
        RuntimeError: If file parameters are mixed with memory return parameters

    Note:
        The system always uses the memory interface internally. File writing
        is controlled by the presence of outputFile, abundSaveFile, etc.
        This validation ensures users don't request both data return AND file writing.
    """
    # Determine if this is a memory return request (user wants data returned, not written)
    memory_return_requested = (
        return_array or return_dataframe or return_rates or return_heating
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
    return_rates: bool = False,
    return_heating: bool = False,
):
    """
    return function that takes in the object that was modelled and returns the values based on the specified booleans.

    Args:
        model_object: model_object of a class that inherited from AbstractModel, from which the results should be returned.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, the function will return
            the success_flag, dissipation_time if the model_object has that attribute, and the final abundances of the out_species.
        return_dataframe: A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, the function will return
            the success_flag, dissipation_time if the model_object has that attribute, and the final abundances of the out_species.
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        return_heating (bool, optional): A boolean on whether the heating/cooling rates should be returned to a user.
    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, and the model_object has the dissipation_time attribute the second element is the dissipation time. Further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - ratesArray (array): array containing reaction rates (if return_rates=True)
            - heatArray (array): array containing heating/cooling rates (if return_heating=True)
            - dissipation_time (float): dissipation time in years (if model_object contains the dissipation_time attribute)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - ratesDF (pandas.DataFrame): DataFrame containing reaction rates (if return_rates=True)
            - heatingDF (pandas.DataFrame): DataFrame containing heating/cooling rates (if return_heating=True)
            - dissipation_time (float): dissipation time in years (if model_object contains the dissipation_time attribute)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.

    """
    if return_dataframe:
        result_dfs = model_object.get_dataframes(
            joined=False, with_rates=return_rates, with_heating=return_heating
        )
        phys_df = result_dfs[0]
        chem_df = result_dfs[1]
        rates_df = None
        heating_df = None

        # Extract rates and heating based on what was requested
        idx = 2
        if return_rates and len(result_dfs) > idx:
            rates_df = result_dfs[idx]
            idx += 1
        if return_heating and len(result_dfs) > idx:
            heating_df = result_dfs[idx]

        # FIXED format: [phys, chem, rates/None, heating/None, abundances/dissipation, flag]
        # For models with dissipation_time: [phys, chem, rates/None, heating/None, dissipation, abundances, flag]
        if hasattr(model_object, "dissipation_time"):
            return (
                phys_df,
                chem_df,
                rates_df,
                heating_df,
                model_object.dissipation_time,
                model_object.next_starting_chemistry_array,
                model_object.success_flag,
            )
        else:
            return (
                phys_df,
                chem_df,
                rates_df,
                heating_df,
                model_object.next_starting_chemistry_array,
                model_object.success_flag,
            )
    elif return_array:
        # FIXED format: [phys, chem, rates/None, heating/None, abundances/dissipation, flag]
        # For models with dissipation_time: [phys, chem, rates/None, heating/None, dissipation, abundances, flag]
        if hasattr(model_object, "dissipation_time"):
            return (
                model_object.physics_array,
                model_object.chemical_abun_array,
                model_object.rates_array if return_rates else None,
                model_object.heat_array if return_heating else None,
                model_object.dissipation_time,
                model_object.next_starting_chemistry_array,
                model_object.success_flag,
            )
        else:
            return (
                model_object.physics_array,
                model_object.chemical_abun_array,
                model_object.rates_array if return_rates else None,
                model_object.heat_array if return_heating else None,
                model_object.next_starting_chemistry_array,
                model_object.success_flag,
            )
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
    param_dict: dict = None,
    out_species: list = ["H", "N", "C", "O"],
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rates: bool = False,
    return_heating: bool = False,
    starting_chemistry: np.array = None,
    timepoints: int = TIMEPOINTS,
):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file.
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file.
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        return_heating (bool, optional): A boolean on whether the heating/cooling arrays should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem.
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS
    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - ratesArray (array or None): array containing reaction rates for each timestep (if return_rates=True)
            - heatArray (array or None): array containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - ratesDF (pandas.DataFrame or None): DataFrame containing reaction rates for each timestep (if return_rates=True)
            - heatingDF (pandas.DataFrame or None): DataFrame containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    # Validate functional API constraints
    __validate_functional_api_params__(
        param_dict=param_dict,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rates=return_rates,
        return_heating=return_heating,
        starting_chemistry=starting_chemistry,
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
        return_rates=return_rates,
        return_heating=return_heating,
    )


def __collapse__(
    collapse: str,
    physics_output: str,
    param_dict: dict = None,
    out_species: list = ["H", "N", "C", "O"],
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rates: bool = False,
    return_heating: bool = False,
    starting_chemistry: np.array = None,
    timepoints: int = TIMEPOINTS,
):
    """Run collapse model from UCLCHEM based on Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)

    Args:
        collapse (str): A string containing the collapse type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'
        physics_output(str): Filename to store physics output, only relevant for 'filament' and 'ambipolar' collapses. If None, no physics output will be saved.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        return_heating (bool, optional): A boolean on whether the heating/cooling arrays should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - ratesArray (array or None): array containing reaction rates for each timestep (if return_rates=True)
            - heatArray (array or None): array containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - ratesDF (pandas.DataFrame or None): DataFrame containing reaction rates for each timestep (if return_rates=True)
            - heatingDF (pandas.DataFrame or None): DataFrame containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    __validate_functional_api_params__(
        param_dict,
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
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
        return_rates=return_rates,
        return_heating=return_heating,
    )


def __prestellar_core__(
    temp_indx: int,
    max_temperature: float,
    param_dict: dict = None,
    out_species: list = ["H", "N", "C", "O"],
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rates: bool = False,
    return_heating: bool = False,
    starting_chemistry: np.array = None,
    timepoints: int = TIMEPOINTS,
):
    """Run prestellar core model from UCLCHEM, based on Viti et al. 2004 and Collings et al. 2004. This model type was previously known as hot core

    Args:
        temp_indx (int): Used to select the mass of prestellar core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        max_temperature (float): Value at which gas temperature will stop increasing.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        return_heating (bool, optional): A boolean on whether the heating/cooling arrays should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - ratesArray (array or None): array containing reaction rates for each timestep (if return_rates=True)
            - heatArray (array or None): array containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - ratesDF (pandas.DataFrame or None): DataFrame containing reaction rates for each timestep (if return_rates=True)
            - heatingDF (pandas.DataFrame or None): DataFrame containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    __validate_functional_api_params__(
        param_dict,
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
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
        return_rates=return_rates,
        return_heating=return_heating,
    )


def __cshock__(
    shock_vel: float,
    timestep_factor: float = 0.01,
    minimum_temperature: float = 0.0,
    param_dict: dict = None,
    out_species: list = ["H", "N", "C", "O"],
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rates: bool = False,
    return_heating: bool = False,
    starting_chemistry: np.array = None,
    timepoints: int = TIMEPOINTS,
):
    """Run C-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock in km/s
        timestep_factor (float, optional): Whilst the time is less than 2 times the dissipation time of shock, timestep is timestep_factor*dissipation time. Essentially controls
        how well resolved the shock is in your model. Defaults to 0.01.
        minimum_temperature (float, optional): Minimum post-shock temperature. Defaults to 0.0 (no minimum). The shocked gas typically cools to `initialTemp` if this is not set.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        return_heating (bool, optional): A boolean on whether the heating/cooling arrays should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, the second element is the dissipation time and further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - ratesArray (array or None): array containing reaction rates for each timestep (if return_rates=True)
            - heatArray (array or None): array containing heating/cooling terms for each timestep (if return_heating=True)
            - disspation_time (float): dissipation time in years
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - ratesDF (pandas.DataFrame or None): DataFrame containing reaction rates for each timestep (if return_rates=True)
            - heatingDF (pandas.DataFrame or None): DataFrame containing heating/cooling terms for each timestep (if return_heating=True)
            - disspation_time (float): dissipation time in years
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    __validate_functional_api_params__(
        param_dict,
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
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
        return_rates=return_rates,
        return_heating=return_heating,
    )


def __jshock__(
    shock_vel: float,
    param_dict: dict = None,
    out_species: list = ["H", "N", "C", "O"],
    return_array: bool = False,
    return_dataframe: bool = False,
    return_rates: bool = False,
    return_heating: bool = False,
    starting_chemistry: np.array = None,
    timepoints: int = TIMEPOINTS,
):
    """Run J-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        return_heating (bool, optional): A boolean on whether the heating/cooling arrays should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

    Returns:if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, the second element is the dissipation time and further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - ratesArray (array or None): array containing reaction rates for each timestep (if return_rates=True)
            - heatArray (array or None): array containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - ratesDF (pandas.DataFrame or None): DataFrame containing reaction rates for each timestep (if return_rates=True)
            - heatingDF (pandas.DataFrame or None): DataFrame containing heating/cooling terms for each timestep (if return_heating=True)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.

    """
    __validate_functional_api_params__(
        param_dict,
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
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
        return_rates=return_rates,
        return_heating=return_heating,
    )


# Expose the functional API functions at module level
cloud = __cloud__
collapse = __collapse__
prestellar_core = __prestellar_core__
hot_core = __prestellar_core__  # Alias for backward compatibility
cshock = __cshock__
jshock = __jshock__
