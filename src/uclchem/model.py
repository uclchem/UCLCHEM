try:
    from .uclchemwrap import uclchemwrap as wrap
except ImportError as error:
    print("No UCLCHEM module, run ``make python'' in src/")
    print("Utility and plotting functions available but UCLCHEM based functions will fail")
    print("Import error was:")
    print(error)
    print("\n\n")


def _reform_inputs(param_dict, out_species):
    """Copies param_dict so as not to modify user's dictionary. Then reformats out_species from pythonic list
    to a string of space separated names for Fortran.
    """
    if param_dict is None:
        param_dict = {}
    else:
        # lower case (and conveniently copy so we don't edit) the user's dictionary
        # this is key to UCLCHEM's "case insensitivity"
        param_dict = {k.lower(): v for k, v in param_dict.items()}
    if out_species is not None:
        n_out = len(out_species)
        param_dict["outspecies"] = n_out
        out_species = " ".join(out_species)
    else:
        out_species = ""
        n_out = 0
    return n_out, param_dict, out_species


def cloud(param_dict=None, out_species=None):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.

    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
    """
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    abunds, success_flag = wrap.cloud(dictionary=param_dict, outspeciesin=out_species)
    if success_flag < 0 or n_out == 0:
        return success_flag
    else:
        return abunds[:n_out]


def collapse(collapse, physics_output, param_dict=None, out_species=None):
    """Run collapse model from UCLCHEM based on Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)

    Args:
        collapse (str): A string containing the collapse type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'
        physics_output(str): Filename to store physics output, only relevant for 'filament' and 'ambipolar' collapses. If None, no physics output will be saved.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.

    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
    """
    collapse_dict = {"BE1.1": 1, "BE4": 2, "filament": 3, "ambipolar": 4}
    try:
        collapse = collapse_dict[collapse]
    except:
        raise ValueError("collapse must be one of 'BE1.1', 'BE4', 'filament', or 'ambipolar'")
    write_physics = physics_output is not None
    if not write_physics:
        physics_output = ""

    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    abunds, success_flag = wrap.collapse(
        collapse, physics_output, write_physics, dictionary=param_dict, outspeciesin=out_species
    )
    if success_flag < 0 or n_out == 0:
        return success_flag
    else:
        return abunds[:n_out]


def hot_core(temp_indx, max_temperature, param_dict=None, out_species=None):
    """Run hot core model from UCLCHEM, based on Viti et al. 2004 and Collings et al. 2004

    Args:
        temp_indx (int): Used to select the mass of hot core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        max_temperature (float): Value at which gas temperature will stop increasing.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.

    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
    """
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    abunds, success_flag = wrap.hot_core(
        temp_indx=temp_indx,
        max_temp=max_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
    )
    if success_flag < 0 or n_out == 0:
        return success_flag
    else:
        return abunds[:n_out]


def cshock(shock_vel, timestep_factor=0.01, param_dict=None, out_species=None):
    """Run C-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock in km/s
        timestep_factor (float, optional): Whilst the time is less than 2 times the dissipation time of shock, timestep is timestep_factor*dissipation time. Essentially controls
        how well resolved the shock is in your model. Defaults to 0.01.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
        float: The dissipation time of the shock in years
    """
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)

    abunds, disspation_time, success_flag = wrap.cshock(
        shock_vel,
        timestep_factor=timestep_factor,
        dictionary=param_dict,
        outspeciesin=out_species,
    )
    result = abunds[:n_out] if n_out > 0 else success_flag
    if success_flag < 0:
        disspation_time = None
    return result, disspation_time


def jshock(shock_vel, param_dict=None, out_species=None):
    """Run J-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
    """
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)

    abunds, success_flag = wrap.jshock(
        shock_vel,
        dictionary=param_dict,
        outspeciesin=out_species,
    )
    if success_flag < 0 or n_out == 0:
        return success_flag
    else:
        return abunds[:n_out]
