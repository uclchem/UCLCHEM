import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import uclchemwrap
from uclchemwrap import uclchemwrap as wrap

from uclchem.constants import (
    N_PHYSICAL_PARAMETERS,
    PHYSICAL_PARAMETERS,
    TIMEPOINTS,
    n_heating_terms,
    n_reactions,
    n_species,
)

OUTPUT_MODE = ""


def set_collisional_rates_directory():
    coolant_directory = (
        os.path.dirname(os.path.abspath(__file__)) + "/data/collisional_rates/"
    )
    # Provide the correct path to the coolant files:
    assert len(coolant_directory) < 256, (
        "Coolant directory path is too long, please shorten it. Path is "
        + coolant_directory
    )
    try:
        uclchemwrap.defaultparameters.coolantdatadir = coolant_directory
        assert (
            str(np.char.decode(uclchemwrap.defaultparameters.coolantdatadir)).strip()
            == coolant_directory
        ), "Coolant directory path is not set correctly, please check the path."
    except AttributeError:
        logging.warning(
            "Cannot set the coolant directory path, please set 'coolantDataDir' correctly at runtime."
        )


def reaction_line_formatter(line):
    reactants = list(filter(lambda x: not str(x).lower().endswith("nan"), line[0:3]))
    products = list(filter(lambda x: not str(x).lower().endswith("nan"), line[3:7]))
    return " + ".join(reactants) + " -> " + " + ".join(products)


class ReactionNamesStore:
    def __init__(self):
        self.reaction_names = None

    def __call__(self):
        # Only load the reactions once, after that use the cached version
        if self.reaction_names is None:
            reactions = pd.read_csv(
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)), "reactions.csv"
                )
            )
            # format the reactions:
            self.reaction_names = [
                reaction_line_formatter(line) for idx, line in reactions.iterrows()
            ]
        return self.reaction_names


# Before doing anything else, set the right collision rate directory and get the reaction names.
set_collisional_rates_directory()
get_reaction_names = ReactionNamesStore()


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
            assert k.lower() not in new_param_dict, (
                f"Lower case key {k} is already in the dict, stopping"
            )
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


def _get_standard_array_specs(return_rates=False, return_heating=False):
    """Get standard array specifications for UCLCHEM models

    Args:
        return_rates (bool): Whether to include rates array
        return_heating (bool): Whether to include heating array

    Returns:
        dict: Array specifications
    """
    specs = {
        "physicsarray": {
            "third_dim": N_PHYSICAL_PARAMETERS,
            "dtype": "float64",
            "dummy": False,
        },
        "chemicalabunarray": {
            "third_dim": n_species,
            "dtype": "float64",
            "dummy": False,
        },
    }
    if return_rates:
        specs["ratesarray"] = {
            "third_dim": n_reactions,
            "dtype": "float64",
            "dummy": False,
        }
    if return_heating:
        specs["heatarray"] = {
            "third_dim": n_heating_terms,
            "dtype": "float64",
            "dummy": False,
        }
    return specs


def _create_arrays(param_dict, array_specs, timepoints=TIMEPOINTS):
    """Create Fortran arrays based on specifications

    Args:
        param_dict (dict): Parameter dictionary containing 'points'
        array_specs (dict): Dictionary specifying arrays to create. Format:
            {
                'array_name': {
                    'third_dim': int,  # Size of third dimension
                    'dtype': str,      # 'float64' or 'int64'
                    'dummy': bool   # Whether to create a dummy array rather than a real one.
                }
            }
        timepoints (int): Number of timepoints

    Returns:
        dict: Dictionary containing the created arrays
    """
    arrays = {}
    points = param_dict.get("points", 1)
    dtype_map = {"float64": np.float64, "int64": np.int64}
    for array_name, spec in array_specs.items():
        dtype = dtype_map.get(spec["dtype"], np.float64)
        # If the array is a dummy, use 2 timepoints; This cause UCLCHEM to
        # ignore the array; opting to write to file instead.
        actual_time_points = timepoints + 1 if not spec.get("dummy", True) else 2
        # Get the array:
        arrays[array_name] = np.zeros(
            shape=(actual_time_points, points, spec["third_dim"]),
            dtype=dtype,
            order="F",
        )
    return arrays


def pre_flight_checklist(
    return_array,
    return_dataframe,
    return_rates,
    return_heating,
    starting_chemistry=None,
    user_params={},
):
    global OUTPUT_MODE
    """Function that ensures that we aren't mixing in memory and write to disk mode.

    Args:
        return_array (bool): Whether to return arrays
        return_dataframe (bool): whether to return dataframes
        starting_chemistry (np.array): Starting chemistry array
        return rates (bool): Whether to return reaction rates
        return_heating (bool): Whether to return heating and cooling rates
        user_params (dict): The user parameter that has to be specified

    Raises:
        RuntimeError: If anything causes undefined behaviour during UCLCHEM runtime.
    """
    if starting_chemistry is not None:
        assert return_array or return_dataframe, (
            "starting_chemistry can only be used with return_array or return_dataframe set to True;\n"
            "Instead specify 'abundLoadFile' in the param_dict to load starting abundances from a file."
        )
    # Check that we aren't mixing in memory and write to disk mode
    if return_array or return_dataframe or return_rates or return_heating:
        file_keys = [k for k in user_params.keys() if k.lower().endswith("file")]
        if file_keys:
            raise RuntimeError(
                "return_array or return_dataframe cannot be used if any "
                "output of input file is specified.\n"
                + f"Offending keys: {', '.join(file_keys)}"
            )
        if return_rates:
            assert return_array or return_dataframe, (
                "return_rates and return_heating can only be used with return_array or return_dataframe set to True; "
            )
        # Check we didn't run a disk based model before:
        if OUTPUT_MODE:
            assert OUTPUT_MODE == "memory", (
                f"Cannot run an in memory based model after running a disk based one in the same session. Found prior run with: {OUTPUT_MODE}"
            )
        OUTPUT_MODE = "memory"
    else:
        # Ensure we never run a disk based model after running an in memory one
        if OUTPUT_MODE:
            assert OUTPUT_MODE == "disk", (
                "Cannot run a disk based model after running an in memory one in the same session."
            )
        OUTPUT_MODE = "disk"


def _array_clean(
    physicalParameterArray,
    chemicalAbundanceArray,
    specname,
    ratesArray=None,
    heatArray=None,
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
    # Also clean the rates array, only if we have it:
    if ratesArray is not None:
        ratesArray = ratesArray[: last_timestep_index + 1, :, :]
    # Also clean the heating array, only if we have it:
    if heatArray is not None:
        heatArray = heatArray[: last_timestep_index + 1, :, :]
    # Get the last arrays simulated, easy for starting another model.
    abundanceStart = chemicalAbundanceArray[last_timestep_index, 0, :]
    return physicsArray, chemArray, specname_new, ratesArray, heatArray, abundanceStart


def outputArrays_to_DataFrame(
    physicalParameterArray, chemicalAbundanceArray, specname, ratesArray, heatArray=None
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
        physicalParameterArray[:, 0, :N_PHYSICAL_PARAMETERS],
        index=None,
        columns=PHYSICAL_PARAMETERS,
    )
    # Create a abundances dataframe.
    chemistry_df = pd.DataFrame(
        chemicalAbundanceArray[:, 0, :], index=None, columns=specname
    )
    if ratesArray is not None:
        # Create a rates dataframe.
        rates_df = pd.DataFrame(
            ratesArray[:, 0, :], index=None, columns=get_reaction_names()
        )
    else:
        rates_df = None

    if heatArray is not None:
        # Create a heating dataframe.
        heating_columns = [
            "Time",
            "Atomic Cooling",
            "Collisionally Induced Emission",
            "Compton Scattering Cooling",
            "Continuum Emission Cooling",
            "Line Cooling 1",
            "Line Cooling 2",
            "Line Cooling 3",
            "Line Cooling 4",
            "Line Cooling 5",
            "Photoelectric Heating",
            "H2Formation Heating",
            "FUVPumping Heating",
            "Photodissociation Heating",
            "CIonization Heating",
            "Cosmic Ray Heating",
            "Turbulent Heating",
            "Gas-Grain Collisions",
            "Chemical Heating",
        ]
        heating_df = pd.DataFrame(
            heatArray[:, 0, :], index=None, columns=heating_columns
        )
    else:
        heating_df = None

    return physics_df, chemistry_df, rates_df, heating_df


def cloud(
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    return_heating=False,
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
    pre_flight_checklist(
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
        param_dict,
    )

    # Create all arrays using the generalized approach
    array_specs = _get_standard_array_specs(
        return_rates=return_rates, return_heating=return_heating
    )
    arrays = _create_arrays(param_dict, array_specs, timepoints)
    physicsArray = arrays["physicsarray"]
    chemicalAbunArray = arrays["chemicalabunarray"]
    ratesArray = arrays.get("ratesarray")
    heatArray = arrays.get("heatarray")
    _, _, _, _, abunds, specname, success_flag = wrap.cloud(
        dictionary=param_dict,
        outspeciesin=out_species,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        heatarray=heatArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    # Overwrite the heatArray with None if its not used:
    if not return_heating or not (return_array or return_dataframe):
        heatArray = None
    if return_array or return_dataframe:
        (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        ) = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray, heatArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
        ) + (
            abundanceStart,
            success_flag,
        )
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            heatArray,
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
    return_rates=False,
    return_heating=False,
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
    pre_flight_checklist(
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
        param_dict,
    )
    array_specs = _get_standard_array_specs(
        return_rates=return_rates, return_heating=return_heating
    )
    arrays = _create_arrays(param_dict, array_specs, timepoints)
    physicsArray = arrays["physicsarray"]
    chemicalAbunArray = arrays["chemicalabunarray"]
    ratesArray = arrays.get("ratesarray")
    heatArray = arrays.get("heatarray")
    _, _, _, _, abunds, specname, success_flag = wrap.collapse(
        collapsein=collapse,
        collapsefilein=physics_output,
        writeout=write_physics,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        heatarray=heatArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    # Overwrite the heatArray with None if its not used:
    if not return_heating or not (return_array or return_dataframe):
        heatArray = None
    if return_array or return_dataframe:
        (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        ) = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray, heatArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
        ) + (
            abundanceStart,
            success_flag,
        )
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            heatArray,
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
    return_rates=False,
    return_heating=False,
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
    pre_flight_checklist(
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
        param_dict,
    )

    # Create all arrays using the generalized approach
    array_specs = _get_standard_array_specs(
        return_rates=return_rates, return_heating=return_heating
    )
    arrays = _create_arrays(param_dict, array_specs, timepoints)
    physicsArray = arrays["physicsarray"]
    chemicalAbunArray = arrays["chemicalabunarray"]
    ratesArray = arrays.get("ratesarray")
    heatArray = arrays.get("heatarray")
    give_start_abund = starting_chemistry is not None
    _, _, _, _, abunds, specname, success_flag = wrap.hot_core(
        temp_indx=temp_indx,
        max_temp=max_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        ratesarray=ratesArray,
        chemicalabunarray=chemicalAbunArray,
        heatarray=heatArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    # Overwrite the heatArray with None if its not used:
    if not return_heating or not (return_array or return_dataframe):
        heatArray = None
    if return_array or return_dataframe:
        (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        ) = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray, heatArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
        ) + (
            abundanceStart,
            success_flag,
        )
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            heatArray,
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
    return_rates=False,
    return_heating=False,
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
    pre_flight_checklist(
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
        param_dict,
    )

    # Create all arrays using the generalized approach
    array_specs = _get_standard_array_specs(
        return_rates=return_rates, return_heating=return_heating
    )
    arrays = _create_arrays(param_dict, array_specs, timepoints)
    physicsArray = arrays["physicsarray"]
    chemicalAbunArray = arrays["chemicalabunarray"]
    ratesArray = arrays.get("ratesarray")
    heatArray = arrays.get("heatarray")
    give_start_abund = starting_chemistry is not None
    _, _, _, _, abunds, disspation_time, specname, success_flag = wrap.cshock(
        shock_vel=shock_vel,
        timestep_factor=timestep_factor,
        minimum_temperature=minimum_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        ratesarray=ratesArray,
        chemicalabunarray=chemicalAbunArray,
        heatarray=heatArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates:
        ratesArray = None
    # Overwrite the heatArray with None if its not used:
    if not return_heating or not (return_array or return_dataframe):
        heatArray = None
    if success_flag < 0:
        disspation_time = None
        abunds = []
    else:
        abunds = list(abunds[:n_out])
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    if return_array or return_dataframe:
        (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        ) = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray, heatArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
        ) + (
            disspation_time,
            abundanceStart,
            success_flag,
        )
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            heatArray,
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
        (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        ) = [success_flag, disspation_time] + abunds
        return (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        )


def jshock(
    shock_vel,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    return_heating=False,
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
    pre_flight_checklist(
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
        param_dict,
    )

    # Create all arrays using the generalized approach
    array_specs = _get_standard_array_specs(
        return_rates=return_rates, return_heating=return_heating
    )
    arrays = _create_arrays(param_dict, array_specs, timepoints)
    physicsArray = arrays["physicsarray"]
    chemicalAbunArray = arrays["chemicalabunarray"]
    ratesArray = arrays.get("ratesarray")
    heatArray = arrays.get("heatarray")
    give_start_abund = starting_chemistry is not None

    _, _, _, _, abunds, specname, success_flag = wrap.jshock(
        shock_vel=shock_vel,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        heatarray=heatArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    # Overwrite the heatArray with None if its not used:
    if not return_heating or not (return_array or return_dataframe):
        heatArray = None
    if return_array or return_dataframe:
        (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        ) = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray, heatArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
        ) + (
            abundanceStart,
            success_flag,
        )
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            heatArray,
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
    return_rates=False,
    return_heating=False,
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
    pre_flight_checklist(
        return_array,
        return_dataframe,
        return_rates,
        return_heating,
        starting_chemistry,
        param_dict,
    )
    # Create all arrays using the generalized approach
    array_specs = _get_standard_array_specs(
        return_rates=return_rates, return_heating=return_heating
    )
    arrays = _create_arrays(param_dict, array_specs, len(time_array))
    physicsArray = arrays["physicsarray"]
    chemicalAbunArray = arrays["chemicalabunarray"]
    ratesArray = arrays.get("ratesarray")
    heatArray = arrays.get("heatarray")
    _, _, _, _, abunds, specname, success_flag = wrap.postprocess(
        dictionary=param_dict,
        outspeciesin=out_species,
        timepoints=len(time_array),
        gridpoints=param_dict["points"],
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        heatarray=heatArray,
        abundancestart=starting_chemistry,
        usecoldens=coldens_H_array is not None,
        **postprocess_arrays,
    )
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    # Overwrite the heatArray with None if its not used:
    if not return_heating or not (return_array or return_dataframe):
        heatArray = None
    if return_array or return_dataframe:
        (
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
            abundanceStart,
        ) = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray, heatArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
            heatArray,
        ) + (
            abundanceStart,
            success_flag,
        )
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            heatArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)
