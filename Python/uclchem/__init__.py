try:
    from .uclchem import wrap
except ImportError as error:
    print("No UCLCHEM module, run ``make python'' in src/")
    print("Utility and plotting functions available but UCLCHEM based functions will fail")
    print("Import error was:")
    print(error)
    print("\n\n")
from . import analysis

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import color_palette


def cloud(param_dict, out_species=None):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.

    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
    """
    param_dict = param_dict.copy()
    if out_species is not None:
        n_out = len(out_species)
        param_dict["outSpecies"] = n_out
        out_species = " ".join(out_species)
    else:
        out_species = ""
        n_out = 0

    abunds, success_flag = wrap.cloud(dictionary=param_dict, outspeciesin=out_species)
    if success_flag < 0 or n_out == 0:
        return success_flag
    else:
        return abunds[: n_out]


def hot_core(temp_indx, max_temperature, param_dict, out_species):
    """Run hot core model from UCLCHEM, following Viti et al. 2004.

    Args:
        temp_indx (int): Used to select the mass of hot core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        max_temperature (float): Value at which gas temperature will stop increasing.
        param_dict (dict): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.

    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
    """
    param_dict = param_dict.copy()
    if out_species is not None:
        n_out = len(out_species)
        param_dict["outSpecies"] = n_out
        out_species = " ".join(out_species)
    else:
        out_species = ""
        n_out = 0

    abunds, success_flag = wrap.hot_core(
        temp_indx=temp_indx,
        max_temp=max_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
    )
    if success_flag < 0 or n_out == 0:
        return success_flag
    else:
        return abunds[: n_out]


def cshock(shock_vel, param_dict, out_species=None):
    """Run C-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
        float: The dissipation time of the shock
    """
    param_dict = param_dict.copy()
    if out_species is not None:
        n_out = len(out_species)
        param_dict["outSpecies"] = n_out
        out_species = " ".join(out_species)
    else:
        out_species = ""
        n_out = 0

    abunds, disspation_time, success_flag = wrap.cshock(
        shock_vel,
        dictionary=param_dict,
        outspeciesin=out_species,
    )
    if success_flag < 0 or n_out == 0:
        disspation_time=None
        return success_flag, disspation_time
    else:
        return  abunds[: n_out], disspation_time
    

def jshock(shock_vel, param_dict, out_species=None):
    """Run J-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
    Returns:
        int,list: A integer which is negative if the model failed to run, or a list of abundances of all species in `outSpecies`
    """
    param_dict = param_dict.copy()
    if out_species is not None:
        n_out = len(out_species)
        param_dict["outSpecies"] = n_out
        out_species = " ".join(out_species)
    else:
        out_species = ""
        n_out = 0
        
    print(out_species)

    print(param_dict)
    abunds, success_flag = wrap.jshock(
        shock_vel,
        dictionary=param_dict,
        outspeciesin=out_species,
    )
    if success_flag < 0 or n_out == 0:
        return success_flag
    else:
        return  abunds[: n_out]


def get_species_rates(param_dict, input_abundances, reac_indxs):
    """
    Get the rate of up to 500 reactions from UCLCHEM for a given set of parameters and abundances.
    Intended for use within the analysis script.
    :param param_dict:  A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
    :param input_abundances: Abundance of every species in network
    :param reac_indxs: Index of reactions of interest in the network's reaction list.

    :returns: (ndarray) Array containing the rate of every reaction specified by reac_indxs
    """
    input_abund = np.zeros(500)
    input_abund[: len(input_abundances)] = input_abundances
    rate_indxs = np.zeros(500)
    rate_indxs[: len(reac_indxs)] = reac_indxs
    rates, success_flag = wrap.get_rates(param_dict, input_abund, rate_indxs)
    if success_flag<0:
        raise RuntimeError("UCLCHEM failed to return rates for these parameters")
    return rates[: len(reac_indxs)]


def read_output_file(output_file):
    """
    Read the output of a UCLCHEM run created with the outputFile parameter into a pandas DataFrame

    :param output_file: - (str) path to file containing a full UCLCHEM output

    :return: (dataframe) A dataframe containing the abundances and physical parameters of the model at every time step.
    """
    f = open(output_file)
    f.readline()
    bits = f.readline().split()
    radfield = float(bits[1])
    zeta = float(bits[3])
    data = pd.read_csv(f)
    data["zeta"] = zeta
    data["radfield"] = radfield
    data.columns = data.columns.str.strip()
    return data


def create_abundance_plot(df, species, figsize=(16, 9), plot_file=None):
    """
    Produce a plot of the abundances of chosen species through time, returning the pyplot
    figure and axis objects

    :param df: A dataframe created by `read_output_file`
    :param species: A list of species names to be plotted
    :param plot_file: optional argument with path to file where the plot should be saved

    :return: fig (matplotlib figure) A figure object and ax (matplotlib axis) An axis object which contains the plot
    """
    fig, ax = plt.subplots(figsize=figsize, tight_layout=True)

    ax = plot_species(ax, df, species)
    ax.legend(loc=4, fontsize="small")

    ax.set_xlabel("Time / years")
    ax.set_ylabel("X$_{Species}$")

    ax.set_yscale("log")
    if plot_file is not None:
        fig.savefig(plot_file)
    return fig, ax


def plot_species(ax, df, species):
    """
    Plot the abundance of several species through time onto an existing pyplot axis

    :param ax: pyplot axis on which to plot
    :param df: A dataframe created by `read_output_file`
    :param species: A list of species names to be plotted

    :returns: ax (matplotlib ax) The input axis is returned
    """
    color_palette(n_colors=len(species))
    for specIndx, specName in enumerate(species):
        if specName[0] == "$":
            abundances = df[specName.replace("$", "#")]
            if specName.replace("$", "@") in df.columns:
                abundances = abundances + df[specName.replace("$", "@")]
        else:
            abundances = df[specName]
        ax.plot(df["Time"], abundances, label=specName, lw=2)
        ax.set(yscale="log")
        ax.legend()
    return ax


def param_dict_from_output(output_line):
    """
    Generate a parameter dictionary with enough variables to correctly estimate the rates of
    reactions.

    :param output_line: (pandas series) any row from the relevant UCLCHEM output
    """
    param_dict = {
        "initialDens": output_line["Density"],
        "initialTemp": output_line["gasTemp"],
        "zeta": output_line["zeta"],
        "radfield": output_line["radfield"],
        "baseAv": 0.0,
        "rout": output_line["av"] * (1.6e21) / output_line["Density"],
    }
    return param_dict
