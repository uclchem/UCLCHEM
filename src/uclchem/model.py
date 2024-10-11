from pathlib import Path

import numpy as np
import pandas as pd

from uclchem.constants import (
    N_PHYSICAL_PARAMETERS,
    PHYSICAL_PARAMETERS,
    TIMEPOINTS,
    n_species,
)

from .uclchemwrap import uclchemwrap as wrap


def _reform_inputs(param_dict, out_species):
    """Copies param_dict so as not to modify user's dictionary. Then reformats out_species from pythonic list
    to a string of space separated names for Fortran.
    """
    if param_dict is None:
        param_dict = {}
    else:
        # lower case (and conveniently copy so we don't edit) the user's dictionary
        # this is key to UCLCHEM's "case insensitivity"
        new_param_dict = {}
        for k, v in param_dict.items():
            assert (
                k.lower() not in new_param_dict
            ), f"Lower case key {k} is already in the dict, stopping"
            if isinstance(v, Path):
                v = str(v)
            new_param_dict[k.lower()] = v
        param_dict = new_param_dict.copy()
        del new_param_dict
    if out_species is not None:
        n_out = len(out_species)
        param_dict["outspecies"] = n_out
        out_species = " ".join(out_species)
    else:
        out_species = ""
        n_out = 0
    return n_out, param_dict, out_species


def _format_output(n_out, abunds, success_flag):
    if success_flag < 0 or n_out == 0:
        abunds = []
    else:
        abunds = list(abunds[:n_out])
    return [success_flag] + abunds


def _create_fortranarray(param_dict, nPhysParam, timepoints=TIMEPOINTS):
    # fencepost problem, need to add 1 to timepoints to account for the 0th timestep
    physicsArray = np.zeros(
        shape=(timepoints + 1, param_dict["points"], nPhysParam),
        dtype=np.float64,
        order="F",
    )
    chemicalAbunArray = np.zeros(
        shape=(timepoints + 1, param_dict["points"], n_species),
        dtype=np.float64,
        order="F",
    )
    return physicsArray, chemicalAbunArray


def _return_array_checks(params):
    if any([key.endswith("File") for key in list(params.keys())]):
        raise RuntimeError(
            "return_array or return_dataframe cannot be used if any output of input file is specified."
        )


def _array_clean(
    physicalParameterArray, chemicalAbundanceArray, specname, n_physics_params
):
    """Clean the array

    Args:
        physicalParameterArray (np.ndarray): Array with the UCLCHEM physical parameters
        chemicalAbundanceArray (np.ndarray): Array with the output chemical abundances
        specname (np.ndarray): Numpy array with the names of all the species
        nPhysParam (int): The number of physical parameters you are interested in.

    Returns:
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray:
        The physical parameters, the abundances, overtime species names, the last
        nonzero physical parameters, the last abundances respectively.
    """
    specname_new = specname.astype(str)
    specname_new = np.array([x.strip() for x in specname_new if x != ""])

    # Find the first element with all the zeros
    last_timestep_index = physicalParameterArray[:, 0, 0].nonzero()[0][-1]
    # Get the arrays for only the simulated timesteps (not the zero padded ones)
    physicsArray = physicalParameterArray[: last_timestep_index + 1, :, :]
    chemArray = chemicalAbundanceArray[: last_timestep_index + 1, :, :]
    # Get the last arrays simulated, easy for starting another model.
    abundanceStart = chemicalAbundanceArray[last_timestep_index, 0, :]
    return physicsArray, chemArray, specname_new, abundanceStart


def outputArrays_to_DataFrame(
    physicalParameterArray, chemicalAbundanceArray, specname, physParameter
):
    """Convert the output arrays to a pandas dataframe

    Args:
        physicalParameterArray (np.array): Array with the output physical parameters
        chemicalAbundanceArray (np.array): Array with the output chemical abundances
        specname (list): List with the names of all the species
        physParameter (list): Array with all the physical parameter names

    Returns:
        _type_: _description_
    """
    # Create a physical parameter dataframe
    physics_df = pd.DataFrame(
        physicalParameterArray[:, 0, : len(physParameter)],
        index=None,
        columns=physParameter,
    )
    # Create a abundances dataframe.
    chemistry_df = pd.DataFrame(
        chemicalAbundanceArray[:, 0, :], index=None, columns=specname
    )
    return physics_df, chemistry_df


def cloud(
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    give_start_abund = starting_chemistry is not None
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    # Check to make sure no output files are specified, if so, halt the execution.
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=timepoints
    )
    _, _, abunds, specname, success_flag = wrap.cloud(
        dictionary=param_dict,
        outspeciesin=out_species,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        returnarray=return_array or return_dataframe,
        givestartabund=give_start_abund,
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        abundancestart=starting_chemistry,
    )
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, N_PHYSICAL_PARAMETERS
        )
    if return_dataframe:
        # Mode that returns dataframes
        physics_df, abundances_df = outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            PHYSICAL_PARAMETERS[:N_PHYSICAL_PARAMETERS],
        )
        return physics_df, abundances_df, abundanceStart, success_flag
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def collapse(
    collapse,
    physics_output,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run collapse model from UCLCHEM based on Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)

    Args:
        collapse (str): A string containing the collapse type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'
        physics_output(str): Filename to store physics output, only relevant for 'filament' and 'ambipolar' collapses. If None, no physics output will be saved.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    collapse_dict = {"BE1.1": 1, "BE4": 2, "filament": 3, "ambipolar": 4}
    try:
        collapse = collapse_dict[collapse]
    except KeyError:
        raise ValueError(
            "collapse must be one of 'BE1.1', 'BE4', 'filament', or 'ambipolar'"
        )
    write_physics = physics_output is not None
    if not write_physics:
        physics_output = ""
    give_start_abund = starting_chemistry is not None

    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, len(PHYSICAL_PARAMETERS), timepoints=timepoints
    )
    abunds, specname, success_flag = wrap.collapse(
        collapseIn=collapse,
        collapseFileIn=physics_output,
        writeOut=write_physics,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        givesstartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        abundanceStart=starting_chemistry,
    )
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, N_PHYSICAL_PARAMETERS
        )
    if return_dataframe:
        physicsDF, chemicalDF = outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            PHYSICAL_PARAMETERS[:N_PHYSICAL_PARAMETERS],
        )
        return physicsDF, chemicalDF, abundanceStart, success_flag
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def hot_core(
    temp_indx,
    max_temperature,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run hot core model from UCLCHEM, based on Viti et al. 2004 and Collings et al. 2004

    Args:
        temp_indx (int): Used to select the mass of hot core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        max_temperature (float): Value at which gas temperature will stop increasing.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    if return_array or return_dataframe:
        # Check to make sure no output files are specified, if so, halt the execution.
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, len(PHYSICAL_PARAMETERS), timepoints=timepoints
    )

    give_start_abund = starting_chemistry is not None

    abunds, specname, success_flag = wrap.hot_core(
        temp_indx=temp_indx,
        max_temp=max_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        abundancestart=starting_chemistry,
    )
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, N_PHYSICAL_PARAMETERS
        )
    if return_dataframe:
        physicsDF, chemicalDF = outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            PHYSICAL_PARAMETERS[:N_PHYSICAL_PARAMETERS],
        )
        return physicsDF, chemicalDF, abundanceStart, success_flag
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def cshock(
    shock_vel,
    timestep_factor=0.01,
    minimum_temperature=0.0,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run C-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock in km/s
        timestep_factor (float, optional): Whilst the time is less than 2 times the dissipation time of shock, timestep is timestep_factor*dissipation time. Essentially controls
        how well resolved the shock is in your model. Defaults to 0.01.
        minimum_temperature (float, optional): Minimum post-shock temperature. Defaults to 0.0 (no minimum). The shocked gas typically cools to `initialTemp` if this is not set.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, the second element is the dissipation time and further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - disspation_time (float): dissipation time in years
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - disspation_time (float): dissipation time in years
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=timepoints
    )
    give_start_abund = starting_chemistry is not None
    abunds, disspation_time, specname, success_flag = wrap.cshock(
        shock_vel=shock_vel,
        timestep_factor=timestep_factor,
        minimum_temperature=minimum_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        abundancestart=starting_chemistry,
    )
    if success_flag < 0:
        disspation_time = None
        abunds = []
    else:
        abunds = list(abunds[:n_out])
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, N_PHYSICAL_PARAMETERS
        )
    if return_dataframe:
        physicsDF, chemicalDF = outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            PHYSICAL_PARAMETERS[:N_PHYSICAL_PARAMETERS],
        )
        return (
            physicsDF,
            chemicalDF,
            disspation_time,
            abundanceStart,
            success_flag,
        )
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            disspation_time,
            abundanceStart,
            success_flag,
        )
    else:
        if success_flag < 0:
            disspation_time = None
            abunds = []
        else:
            abunds = list(abunds[:n_out])
        result = [success_flag, disspation_time] + abunds
        return result


def jshock(
    shock_vel,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run J-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, the second element is the dissipation time and further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.

    """
    if "points" not in param_dict:
        param_dict["points"] = 1
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=timepoints
    )
    give_start_abund = starting_chemistry is not None
    abunds, specname, success_flag = wrap.jshock(
        shock_vel=shock_vel,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        abundancestart=starting_chemistry,
    )
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, N_PHYSICAL_PARAMETERS
        )
    if return_dataframe:
        physicsDF, chemicalDF = outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            PHYSICAL_PARAMETERS[:N_PHYSICAL_PARAMETERS],
        )
        return physicsDF, chemicalDF, abundanceStart, success_flag
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def postprocess(
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    starting_chemistry=None,
    time_array=None,
    density_array=None,
    gas_temperature_array=None,
    dust_temperature_array=None,
    zeta_array=None,
    radfield_array=None,
    coldens_H_array=None,
    coldens_H2_array=None,
    coldens_CO_array=None,
    coldens_C_array=None,
):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    # Assure that every array is cast to a fortran array
    postprocess_arrays = dict(
        timegrid=time_array,
        densgrid=density_array,
        gastempgrid=gas_temperature_array,
        dusttempgrid=dust_temperature_array,
        radfieldgrid=radfield_array,
        zetagrid=zeta_array,
        nhgrid=coldens_H_array,
        nh2grid=coldens_H2_array,
        ncogrid=coldens_CO_array,
        ncgrid=coldens_C_array,
    )
    for key, array in postprocess_arrays.items():
        if array is not None:
            # Convert single values into arrays that can be used
            if isinstance(array, float):
                array = np.ones(shape=time_array.shape) * array
            # Assure lengths are correct
            assert len(array) == len(time_array), "All arrays must be the same length"
            # Ensure Fortran memory
            array = np.asfortranarray(array, dtype=np.float64)
            postprocess_arrays[key] = array
    give_start_abund = starting_chemistry is not None
    if not give_start_abund:
        starting_chemistry = np.zeros(
            shape=(n_species),
            dtype=np.float64,
            order="F",
        )
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    # Check to make sure no output files are specified, if so, halt the execution.
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=len(time_array)
    )
    abunds, specname, success_flag = wrap.postprocess(
        dictionary=param_dict,
        outspeciesin=out_species,
        timepoints=len(time_array),
        gridpoints=param_dict["points"],
        returnarray=return_array or return_dataframe,
        givestartabund=give_start_abund,
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        abundancestart=starting_chemistry,
        usecoldens=coldens_H_array is not None,
        **postprocess_arrays,
    )
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, N_PHYSICAL_PARAMETERS
        )
    if return_dataframe:
        # Mode that returns dataframes
        physics_df, abundances_df = outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            PHYSICAL_PARAMETERS[:N_PHYSICAL_PARAMETERS],
        )
        return physics_df, abundances_df, abundanceStart, success_flag
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)
        return _format_output(n_out, abunds, success_flag)
