import itertools
import os
import sys
import json
import time
import types
import signal
import warnings
from contextlib import contextmanager
from copy import deepcopy
from itertools import product

import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
from typing import Optional, Dict, List, AnyStr, Type, Literal, Iterator

# UCLCHEM related imports
from uclchemwrap import uclchemwrap as wrap
from uclchem.analysis import check_element_conservation, create_abundance_plot, plot_species
from uclchem.constants import (
    N_PHYSICAL_PARAMETERS,
    PHYSICAL_PARAMETERS,
    TIMEPOINTS,
    n_reactions,
    n_species,
    default_param_dictionary,
)
# /UCLCHEM related imports

# Multiprocessing imports
import threading
import multiprocessing as mp
from multiprocessing import shared_memory
from typing import Optional, Dict, Any, Tuple
import h5netcdf
from concurrent.futures import ProcessPoolExecutor, as_completed
# /Multiprocessing imports

# Global variables determining formats of write files
PHYSICAL_PARAMETERS_HEADER_FORMAT = "%10s"
# in the below variable, the outputs were chosen according to the spacing needed for
# "      Time,    Density,    gasTemp,   dustTemp,         Av,   radfield,       zeta,      point,"
PHYSICAL_PARAMETERS_VALUE_FORMAT = "%10.3E, %10.4E, %10.2f, %10.2f, %10.4E, %10.4E, %10.4E, %10i"
SPECNAME_HEADER_FORMAT = "%11s"
SPECNAME_VALUE_FORMAT = "%9.5E"
#

# Model registration is intended to prevent code injection during loading time.
REGISTRY: Dict[str, Type["AbstractModel"]] = {}

def register_model(cls: Type["AbstractModel"]):
    name = getattr(cls, "MODEL_NAME", cls.__name__)
    if name in REGISTRY and REGISTRY[name] is not cls:
        raise ValueError(f"Duplicate model registration for {name}")
    REGISTRY[name] = cls
    return cls
# /Global variables determining formats of write files

# Reaction and Species name retrieval classes to reduce file read repetition.
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
            reactions = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "reactions.csv"))
            # format the reactions:
            self.reaction_names = [reaction_line_formatter(line) for idx, line in reactions.iterrows()]
        return self.reaction_names

get_reaction_names = ReactionNamesStore()

class SpeciesNameStore:
    def __init__(self):
        self.species_names = None

    def __call__(self):
        # Only load the species once, after that use the cached version
        if self.species_names is None:
            species = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "species.csv"))
            self.species_names = species["NAME"].tolist()
        return self.species_names

get_species_names = SpeciesNameStore()
# /Reaction and Species name retrieval classes to reduce file read repetition.

# Universal model loader
def load_model(file: str, name: str = 'default', engine: str = "h5netcdf", debug: bool = False):
    """
    load_model bypasses __init__ in order to load a pre-existing model from a file.

    Args:
        file (str): Path to a file that contains previously run and stored models.
        name (str, optional): Name of the stored object, if none was provided `default` will have been used. Defaults to 'default'
        engine (str, optional): “netcdf4”, “h5netcdf” or “zarr”, depending on the engine to be used. Defaults to "h5netcdf".
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False.
    Returns:
        obj (object): Loaded object that inherited from AbstractModel and has the class of to the model found in the loaded file.
    """
    try:
        loaded_data = xr.open_dataset(filename_or_obj=file, group=name, engine=engine)
    except FileNotFoundError:
        print(f"Unable to find file {file}")
        raise FileNotFoundError
    except OSError:
        print(f"Could not find model with name: {name}, in file {file}.")
        raise OSError
    except ValueError:
        print(f"Engine {engine}, is incompatible with xr.open_dataset")
        raise ValueError
    model_class = json.loads(loaded_data["attributes_dict"].item())["model_type"]
    cls = REGISTRY.get(model_class)
    if cls is None:
        raise ValueError(f"Unrecognized model type '{model_class}'. Not in trusted registry.")
    return cls.load_from_dataset(model_ds=loaded_data, debug=debug)
# /Universal model loader

# Worker entry for parallel jobs
def _worker_entry(model_class: str, init_kwargs: dict, shm_descs: dict, result_queue: mp.Queue):
    cls = REGISTRY.get(model_class)
    if cls is None:
        raise ValueError(f"Unrecognized model type '{model_class}'. Not in trusted registry.")
    model = cls.worker_build(init_kwargs=init_kwargs, shm_desc=shm_descs)
    output = model.run_fortran()
    result_queue.put(output)
    return
# /Worker entry for parallel jobs

#TODO Add catch of ctrl+c or other aborts so that it saves model and a full output to files of year, month, day, time type.
class AbstractModel(ABC):
    """Base model class used for inheritance only

    The AbstractModel class serves as an abstract class from which other model classes can be built. It is not intended
    to be used as a standalone class for running UCLCHEM.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values if not provided.
        out_specie_list (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False.
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_specie_list: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: object = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        self._data = xr.Dataset()
        self._pickle_dict = {}
        """Initiates the model with all common factors found within models"""
#["_param_dict", "_data", "_debug", "_prev_handler", "_shm_desc", "_shm_handles", "_orig_sigint", "_proc_handle"]
        # Shared memory
        self.run_type = run_type
        self._shm_desc = {}
        self._shm_handles = {}
        self._proc_handle = None
        self.shared_memory_types = ["managed", "external"]
        self.separate_worker_types = ["managed"]
        # /Shared memory

        # Signal Interrupt
        self._was_interrupted = False
        self._orig_sigint = signal.getsignal(signal.SIGINT)
        # /Signal Interrupt

        self.model_type = str(self.__class__.__name__)
        self._param_dict = {}
        self.out_species_list = out_specie_list
        self.out_species = ""
        self.full_array = None
        self._debug = debug
        self.success_flag = None
        self.specname = get_species_names()

        self.n_out = 0 if read_file is None else None
        self.timepoints = timepoints if read_file is None else None
        self.was_read = False if read_file is None else True
        self.PHYSICAL_PARAMETERS = PHYSICAL_PARAMETERS if read_file is None else None

        self._reform_inputs(param_dict, self.out_species_list)
        if "points" not in self._param_dict:
            self._param_dict["points"] = 1

        self.outputFile = param_dict.pop("outputFile") if "outputFile" in param_dict else None
        self.abundSaveFile = param_dict.pop("abundSaveFile") if "abundSaveFile" in param_dict else None
        self.abundLoadFile = param_dict.pop("abundLoadFile") if "abundLoadFile" in param_dict else None

        self.starting_chemistry_array = None
        if previous_model is None and self.abundLoadFile is None:
            self._create_starting_array(starting_chemistry)
        elif self.abundLoadFile is not None:
            self.read_starting_chemistry_output_file()
        elif previous_model.has_attr('next_starting_chemistry_array'):
            self._create_starting_array(previous_model.next_starting_chemistry_array)

        self.give_start_abund = self.starting_chemistry_array is not None
        self.next_starting_chemistry_array = None

        self.physics_array = None
        self.chemical_abun_array = None
        self._create_fortran_array()
        self.rates_array = None
        self._create_rates_array()
        self.out_species_abundances_array = None

        if read_file is not None:
            self.read_output_file(read_file)
        return

    def __del__(self):
        if hasattr(self, "_shm_handles"):
            self._coordinator_unlink_memory()

    # Separate class building method(s)
    @classmethod
    def load_from_dataset(cls, model_ds: xr.Dataset, debug: bool = False):
        obj = cls.__new__(cls)
        obj._param_dict = json.loads(model_ds["_param_dict"].item())
        model_ds.__delitem__("_param_dict")
        obj._data = xr.Dataset()
        obj._data = model_ds.copy()
        model_ds.close()
        temp_attribute_dict = json.loads(obj._data["attributes_dict"].item())
        for k, v in temp_attribute_dict.items():
            obj._data[k] = v
        obj._data.__delitem__("attributes_dict")
        obj.debug = debug
        return obj

    @classmethod
    def worker_build(cls, init_kwargs, shm_desc):
        obj = cls.__new__(cls)
        for k, v in init_kwargs.items():
            object.__setattr__(obj, k, v)
        obj._reform_array_in_worker(shm_desc)
        return obj
    # /Separate class building methods

    # Class utility methods
    def __getattr__(self, key):
        if (key.startswith("_") and key != "_data"):
            return super().__getattribute__(key)
        elif key in super().__getattribute__("_data"):
            values = super().__getattribute__("_data")[key].values
            if np.shape(values) != ():
                if type(values) == tuple:
                    return values[1]
                else:
                    return values
            else:
                return self._data[key].item()
        else:
            raise AttributeError(f'{self.__class__.__name__} has no attribute of name: "{key}".')

    def __setattr__(self, key, value):
        if key.startswith("_"):
            super().__setattr__(key, value)
        else:
            if key in self._data:
                self._data.__delitem__(key)

            if np.ndim(value) == 3 and "_array" in key:
                self._data[key] = (["time_step", "point", key.replace('array', 'values')], value)
            elif np.ndim(value) == 2 and "_array" in key:
                self._data[key] = (["point", key], value)
            else:
                self._data[key] = value
        return

    def has_attr(self, key):
        """Method to check if the object has an attribute stored in self._data"""
        return key in self._data
    # /Class utility method

    # UCLCHEM utility and analysis wrappers
    def check_conservation(self,
                           element_list: list = None,
                           percent: bool = True
                           ):
        """Utility method to check conservation of the chemical abundances

        Args:
            element_list (list, optional): List of elements to check conservation for. Defaults to self.out_species_lists.
            percent (bool, optional): Flag on if percentage values should be used. Defaults to True.
        """
        if element_list is None:
            element_list = self.out_species_list

        if self._param_dict["points"] > 1:
            for i in range(self._param_dict["points"]):
                print(f"Element conservation report for point {i + 1} of {self._param_dict['points']}")
                print(check_element_conservation(self.get_dataframes(i), element_list, percent))
        else:
            print(f"Element conservation report")
            print(check_element_conservation(self.get_dataframes(0), element_list, percent))

    def check_error(self, only_error: bool = False):
        """
        Prints the error message of the model based on self.success_flag, this method was originally an uclchem.utils function.

        Args:
            only_error (bool, optional): Flag to inform check_error to only print a message if self.success_flag was
                not 0
        """
        if self.success_flag != 0 and self.success_flag is not None:
            errors = {
                -1: "Parameter read failed. Likely due to a mispelled parameter name, compare your dictionary to the parameters docs.",
                -2: "Physics intiialization failed. Often due to user chosing unacceptable parameters such as prestellar core masses or collapse modes that don't exist. Check the docs for your model function.",
                -3: "Chemistry initialization failed",  # this doesn't exist yet
                -4: "Unrecoverable integrator error, DVODE failed to integrate the ODEs in a way that UCLCHEM could not fix. Run UCLCHEM tests to check your network works at all then try to see if bad parameter combination is at play.",
                -5: "Too many integrator fails. DVODE failed to integrate the ODE and UCLCHEM repeatedly altered settings to try to make it pass but tried too many times without success so code aborted to stop infinite loop.",
                -6: "The model was stopped because there are not enough time points allocated in the time array. Increase the number of time points in the time array in constants.py and try again.",
            }
            try:
                print(f'{errors[self.success_flag]}')
            except KeyError:
                raise ValueError(f"Unknown error code: {self.success_flag}")
        elif self.success_flag == 0 and not only_error:
            print(f'Model ran successfully.')
        elif self.success_flag is None:
            print(f'Model has not been run.')
        return

    def create_abundance_plot(self,
                              species: list = None,
                              figsize: tuple[2] = (16, 9),
                              point: int = 0,
                              plot_file=None):
        """uclchem.analysis.create_abundance_plot wrapper method
        Args:
            element_list (list, optional): List of elements to check conservation for. Defaults to  self.out_species_list.
            figsize (tuple[2], optional): The figure size to use for matplotlib Defaults to (16, 9).
            point (int, optional): Integer referring to which point of the UCLCHEM model to use. Defaults to 0.
        Returns:
            fig,ax: matplotlib figure and axis objects
        """
        if species is None:
            species = self.out_species_list

        if point > self._param_dict["points"]:
            raise Exception("'point' must be less than number of modelled points.")
        return create_abundance_plot(self.get_dataframes(point), species, figsize, plot_file)

    def get_dataframes(self, point: int = 0, joined: bool = True, with_rates: bool = False):
        """Converts the model physics and chemical_abun arrays from numpy to pandas arrays.
        Args:
            point (int, optional): Integer referring to which point of the UCLCHEM model to return as pandas does not support higher
                than 2D data structures. Defaults to 0.
            joined (bool, optional): Flag on whether the returned pandas dataframe should be one, or if two dataframes should be
                returned. One physical, one chemical_abun dataframe. Defaults to True.
            with_rates (bool, optional): Flag on whether to include reaction rates in the dataframe, and/or as a separate
                dataframe depending on the value of `joined`. Defaults to False.
        Returns:
            return_df (pandas.DataFrame): Dataframe of the joined arrays for point 'point' if joined = True
            physics_df (pandas.DataFrame): Dataframe of the physical parameters for point 'point' if joined = False
            chemistry_df (pandas.DataFrame): Dataframe of the chemical abundances  for point 'point' if joined = False
            rates_df (pandas.DataFrame): Dataframe of the reaction rates  for point 'point' if joined = False and with_rates = True
        """
        # Create a physical parameter dataframe
        physics_df = pd.DataFrame(
            self.physics_array[:, point, : len(self.PHYSICAL_PARAMETERS)],
            index=None,
            columns=self.PHYSICAL_PARAMETERS,
        )
        # Create an abundances dataframe.
        chemistry_df = pd.DataFrame(
            self.chemical_abun_array[:, point, :], index=None, columns=self.specname
        )
        if self.rates_array is not None and with_rates:
            # Create a rates dataframe.
            rates_df = pd.DataFrame(
                self.rates_array[:, point, :], index=None, columns=get_reaction_names()
            )
        else:
            rates_df = None
        if joined:
            if rates_df is not None and with_rates:
                return_df = physics_df.join(chemistry_df.join(rates_df))
            else:
                return_df = physics_df.join(chemistry_df)
            return return_df
        else:
            if rates_df is not None and with_rates:
                return physics_df, chemistry_df, rates_df
            else:
                return physics_df, chemistry_df

    def plot_species(self, ax: plt.axes, species: list[str] = None, point: int = 0, legend: bool = True, **plot_kwargs):
        """uclchem.analysis.plot(species) wrapper method
        Args:
            ax (plt.axes):
            species (list, optional):
            point (int, optional):
            legend (bool, optional):
            plot_kwargs (dict, optional):
        """
        if species is None:
            species = self.out_species_list
        return plot_species(ax, self.get_dataframes(point), species, legend, **plot_kwargs)
    # /UCLCHEM utility and analysis wrappers

    # Methods to start run of model
    def run(self):
        """__run__ resets the Fortran arrays if the model was not read, allowing the arrays to be reused for new runs."""
        if self.was_read:
            raise RuntimeError("This model was read. It can not be run. ")
        #        self._prev_handler = signal.getsignal(signal.SIGINT)
        #        signal.signal(signal.SIGINT, self.__on_sigint__)
        if self._debug:
            print(f"About to run {self.__class__.__name__} model with the following options:")
            print(f'dictionary = {self._param_dict},')
            print(f'outspeciesin = {self.out_species},')
            print(f'timepoints = {self.timepoints},')
            print(f'gridpoints = {self._param_dict["points"]},')
            print(f'returnarray = True,')
            print(f'returnrates = True,')
            print(f'givestartabund = {self.give_start_abund},')
            print(f'physicsarray = {self.physics_array},')
            print(f'ratesarray = {self.rates_array},')
            print(f'chemicalabunarray = {self.chemical_abun_array},')
            print(f'abundancestart = {self.starting_chemistry_array},')

        def _handler(signum, frame):
            try:
                self.on_interrupt()  # your “final steps”
            finally:
                # Restore default and re-raise KeyboardInterrupt to stop execution
                signal.signal(signal.SIGINT, self._orig_sigint)
                raise KeyboardInterrupt
        if self.run_type in self.shared_memory_types:
            signal.signal(signal.SIGINT, _handler)

        if not self.run_type in self.separate_worker_types:
            output = self.run_fortran()
        elif self.run_type in self.separate_worker_types:
            init_kwargs = self._create_init_dict()
            ctx = mp.get_context("spawn")
            result_queue = ctx.Queue()
            self._proc_handle = ctx.Process(
                target=_worker_entry,
                args=(self.model_type, init_kwargs, self._shm_desc, result_queue),
                daemon=False,
            )
            self._proc_handle.start()
            output = result_queue.get()
            result_queue.close()
            self._proc_handle.join()
            self._proc_handle.close()
            self._proc_handle = None
        else:
            raise ValueError(f"run_type of {self.run_type} is not a valid value.")

        if self.run_type in self.shared_memory_types:
            signal.signal(signal.SIGINT, self._orig_sigint)

        if bool(self._shm_desc):
            self._coordinator_unlink_memory()

        for k, v in output.items():
            self.__setattr__(k, v)

        self._array_clean()
        self.check_error(only_error=True)
        if self.outputFile is not None:
            self.legacy_write_full()
        if self.abundSaveFile is not None:
            self.legacy_write_starting_chemistry()
        return

    @abstractmethod
    def run_fortran(self):
        raise NotImplementedError
    # /Methods to start run of model

    # Model saving
    def save_model(self, file: str, name: str = 'default', engine: str = "h5netcdf", single_precision: bool = False, overwrite: bool = False):
        """
        save_model saves a model to a file on disk. Multiple models can be saved into the same file if different names are used to store them.

        Args:
            file (str): Path to a file to store models.
            name (str, optional): Name to use for the group of the object. Defaults to 'default'
            engine (str, optional): “netcdf4”, “h5netcdf” or “zarr”, depending on the engine to be used. Defaults to "h5netcdf".
            overwrite (bool, optional): Boolean on whether to overwrite pre-existing models, or error out. Defaults to False
        """
        #TODO: Allow for toggling of saving float64 or float32 for the arrays
        if os.path.isfile(file):
            with xr.open_datatree(filename_or_obj=file, engine=engine, phony_dims='sort') as tree:
                if "/" + name in tree.groups:
                    if not overwrite:
                        warnings.warn(f"Model with name: `{name}` already exists in `{file}` but overwrite is set to False. Unable to save model.")
                        return
                tree.close()

        temp_attribute_dict = {}
        for v in self._data.variables:
            if "_array" not in v and v not in ["_orig_sigint"]:
                if np.shape(self._data[v].values) != ():
                    if type(self._data[v].values) == tuple:
                        temp_attribute_dict[v] = self._data[v].values[1].tolist()
                        self._data = self._data.drop_vars(v)
                    else:
                        temp_attribute_dict[v] = self._data[v].values.tolist()
                        self._data = self._data.drop_vars(v)
                else:
                    temp_attribute_dict[v] = self._data[v].item()
                    self._data = self._data.drop_vars(v)

        self._data["attributes_dict"] = xr.DataArray([json.dumps(temp_attribute_dict)])
        self._data["_param_dict"] = xr.DataArray([json.dumps(self._param_dict)])

        self._data.to_netcdf(file, group=name, engine=engine, mode='a')
        return
    # /Model saving

    # Model Passing through Pickling
    def pickle(self):
        for k, v in self._data.items():
            if np.shape(v.values) != ():
                if type(v.values) == tuple:
                    self._pickle_dict[k] = v.values[1]
                else:
                    self._pickle_dict[k] = v.values
            else:
                self._pickle_dict[k] = v.item()
        self._data = None
        return

    def un_pickle(self):
        self._data = xr.Dataset()
        for k, v in self._pickle_dict.items():
            if np.ndim(v) == 3 and "_array" in k:
                self._data[k] = (["time_step", "point", k.replace('array', 'values')], v)
            elif np.ndim(v) == 2 and "_array" in k:
                self._data[k] = (["point", k], v)
            else:
                self._data[k] = v
        return
    # /Model Passing through Pickling

    # Legacy in & output support
    def legacy_read_output_file(self, read_file: str, rates_load_file: str = None):
        '''Perform classic output file reading.
        Args:
            read_file (str): path to file containing a full UCLCHEM output
            rates_load_file (str, optional): path to file containing the reaction rates output from UCLCHEM. Defaults
                to None. #TODO Add the code to read the rates files.
        '''
        self.was_read = True
        columns = np.char.strip(np.loadtxt(read_file, delimiter=",", max_rows=1, dtype=str, comments='%'))
        array = np.loadtxt(read_file, delimiter=",", skiprows=1)
        point_index = np.where(columns == 'point')[0][0]
        self._param_dict["points"] = int(np.max(array[:, point_index]))
        if self._param_dict["points"] > 1:
            array = np.loadtxt(read_file, delimiter=",", skiprows=2)
        row_count = int(np.shape(array)[0] / self._param_dict["points"])

        self.PHYSICAL_PARAMETERS = [p for p in PHYSICAL_PARAMETERS if p in columns]
        self.specname = get_species_names()

        self.physics_array = np.empty((row_count, self._param_dict["points"], len(self.PHYSICAL_PARAMETERS) + 1))
        self.chemical_abun_array = np.empty((row_count, self._param_dict["points"], len(self.specname)))
        for p in range(self._param_dict["points"]):
            self.physics_array[:, p, :] = \
            array[np.where(array[:, point_index] == p + 1), :len(self.PHYSICAL_PARAMETERS) + 1][0]
            self.chemical_abun_array[:, p, :] = \
            array[np.where(array[:, point_index] == p + 1), (len(self.PHYSICAL_PARAMETERS) + 1):][0]
        self._array_clean()
        last_timestep_index = self.physics_array[:, 0, 0].nonzero()[0][-1]
        self.next_starting_chemistry_array = self.chemical_abun_array[last_timestep_index, :, :]
        return

    def legacy_read_starting_chemistry(self):
        """Method to read the starting chemistry from the self.abundLoadFile provided in _param_dict"""
        self._create_starting_array(np.loadtxt(self.abundLoadFile, delimiter=","))
        return

    def legacy_write_full(self):
        """Perform classic output file writing to the file self.outputFile provided in _param_dict"""
        phys = self.physics_array.reshape(-1, self.physics_array.shape[-1])
        chem = self.chemical_abun_array.reshape(-1, self.chemical_abun_array.shape[-1])
        full_array = np.append(phys, chem, axis=1)
        #TODO Move away from the magic numbers seen here.
        string_fmt_string = f'{", ".join([PHYSICAL_PARAMETERS_HEADER_FORMAT] * (len(self.PHYSICAL_PARAMETERS)))}, {", ".join([SPECNAME_HEADER_FORMAT] * len(self.specname))}'
        # Magic numbers here to match/improve the formatting of the classic version
        #TODO Move away from the magic numbers seen here.
        number_fmt_string = f'{PHYSICAL_PARAMETERS_VALUE_FORMAT}, {", ".join([SPECNAME_VALUE_FORMAT] * len(self.specname))}'
        columns = np.array([self.PHYSICAL_PARAMETERS[:-1].tolist() + ["point"] + self.specname.tolist()])
        np.savetxt(self.outputFile, columns, fmt=string_fmt_string)
        with open(self.outputFile, "ab") as f:
            np.savetxt(f, full_array, fmt=number_fmt_string)
        return

    def legacy_write_starting_chemistry(self):
        """Perform classic starting abundance file writing to the file self.abundSaveFile provided in _param_dict"""
        last_timestep_index = self.chemical_abun_array[:, 0, 0].nonzero()[0][-1]
        #TODO Move away from the magic numbers seen here.
        number_fmt_string = f' {", ".join(["%9.5E"] * len(self.specname))}'
        with open(self.abundSaveFile, "wb") as f:
            np.savetxt(f, self.chemical_abun_array[last_timestep_index, :, :], fmt=number_fmt_string)
        return
    # /Legacy in & output support

    # Cleaning of array & inptus
    def _array_clean(self):
        """Internal Method.
        Clean the arrays changed by UCLCHEM Fortran code.
        """
        # Find the first element with all the zeros
        if self._debug:
            print(f"in _array_clean: physics_array = {self.physics_array}")
            print(f"in _array_clean: physics_array type = {type(self.physics_array)}")
        last_timestep_index = self.physics_array[:, 0, 0].nonzero()[0][-1]
        # Get the arrays for only the simulated timesteps (not the zero padded ones)
        self._data = self._data.isel(time_step=slice(0, last_timestep_index+1))
        self.next_starting_chemistry_array = self.chemical_abun_array[last_timestep_index, :, :]
        return

    def _reform_inputs(self,
                       param_dict: dict,
                       out_species: list):
        """Internal Method.
        Copies param_dict so as not to modify user's dictionary. Then reformats out_species from pythonic list
        to a string of space separated names for Fortran.

        Args:
            param_dict (dict): Parameter dictionary passed by the user to the model.
            out_species (list): List of output species that are considered important for this model.
        """
        if param_dict is None:
            self._param_dict = default_param_dictionary
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
            self._param_dict = {**default_param_dictionary, **new_param_dict.copy()}
            del new_param_dict
        for k, v in default_param_dictionary.items():
            if v is None:
                del self._param_dict[k]
        if out_species is not None:
            self.n_out = len(out_species)
            self._param_dict["outspecies"] = self.n_out
            self.out_species = " ".join(out_species)
        else:
            self.out_species = ""
            self.n_out = 0
        return
    # /Cleaning of array & inptus

    # Creation of arrays
    def _create_fortran_array(self):
        """Internal Method.
        Creates Fortran compliant np.arrays that can be passed to the Fortran part of UCLCHEM.
        """
        # For shared memory:
        if self.run_type in self.shared_memory_types:
            self._shm_handles["physics_array"], self._shm_desc["physics_array"], self.physics_array = self._create_shared_memory_allocation(
                (self.timepoints + 1, self._param_dict["points"], N_PHYSICAL_PARAMETERS))
            self._shm_handles["chemical_abun_array"], self._shm_desc["chemical_abun_array"], self.chemical_abun_array = self._create_shared_memory_allocation(
                (self.timepoints + 1, self._param_dict["points"], n_species))
        else:
        # fencepost problem, need to add 1 to timepoints to account for the 0th timestep
            self.physics_array = np.zeros(
                shape=(self.timepoints + 1, self._param_dict["points"], N_PHYSICAL_PARAMETERS),
                dtype=np.float64,
                order="F",
            )
            self.chemical_abun_array = np.zeros(
                shape=(self.timepoints + 1, self._param_dict["points"], n_species),
                dtype=np.float64,
                order="F",
            )
        return

    def _create_rates_array(self):
        """Internal Method.
        Creates Fortran compliant np.array for rates that can be passed to the Fortran part of UCLCHEM.
        """
        # For shared memory:
        if self.run_type in self.shared_memory_types:
            self._shm_handles["rates_array"], self._shm_desc["rates_array"], self.rates_array = self._create_shared_memory_allocation(
                (self.timepoints + 1, self._param_dict["points"], n_reactions))
        else:
            self.rates_array = np.zeros(shape=(self.timepoints + 1, self._param_dict["points"], n_reactions),
                                       dtype=np.float64, order="F")
        return

    def _create_starting_array(self, starting_chemistry):
        if starting_chemistry is None:
            self.starting_chemistry_array = None
        else:
            if len(np.shape(starting_chemistry)) == 1:
                starting_chemistry = starting_chemistry[np.newaxis, :]
            # For shared memory:
            if self.run_type in self.shared_memory_types:
                self._shm_handles["starting_chemistry_array"], self._shm_desc["starting_chemistry_array"], self.starting_chemistry_array = self._create_shared_memory_allocation(
                    np.shape(starting_chemistry))
                np.copyto(self.starting_chemistry_array, starting_chemistry, casting='no')
            else:
                self.starting_chemistry_array = np.asfortranarray(starting_chemistry, dtype=np.float64)
        return
    # /Creation of arrays

    # Signal Interrupt Catch
    def on_interrupt(self, grid = False, model_name = None):
        if self._proc_handle is not None:
            self._proc_handle.terminate()
            self._proc_handle.join()
            self._proc_handle = None

        if bool(self._shm_desc):
            self._coordinator_unlink_memory()
        self._array_clean()

        error_time = datetime.now().strftime("%y_%m_%d_%H_%M")
        if not grid:
            if self.outputFile is None:
                self.outputFile = "./" + error_time + ".dat"
            elif "/" in self.outputFile:
                self.outputFile = (self.outputFile[:self.outputFile.rfind("/") + 1] +
                                   error_time +
                                   self.outputFile[self.outputFile.rfind("."):])
            else:
                self.outputFile = "./" + error_time + ".dat"
            self.legacy_write_full()

        self._was_interrupted = True
        self.save_model(
            file="./" + error_time + ".h5" if not grid else "./grid_interrupted_models.h5",
            name=model_name if model_name is not None else "interrupted",
            overwrite=True
        )
        return
    # /Signal Interrupt Catch

    # Shared memory handlers
    @staticmethod
    def _create_shared_memory_allocation(shape:tuple):
        nbytes = int(np.prod(shape) * np.dtype(np.float64).itemsize)
        shm = shared_memory.SharedMemory(create=True, size=nbytes)
        array = np.ndarray(shape, dtype=np.float64, buffer=shm.buf, order="F")
        array.fill(0.0)
        spec = {"name": shm.name, "shape": shape}
        return shm, spec, array

    def _reform_array_in_worker(self, shm_desc):
        object.__setattr__(self, '_shm_handles', {})
        for k, v in shm_desc.items():
            shm = shared_memory.SharedMemory(name=v["name"], create=False)
            object.__setattr__(
                self,
                k,
                np.ndarray(
                    shape = v["shape"],
                    dtype=np.float64,
                    buffer=shm.buf,
                    order="F"
                )
            )
            self._shm_handles[k] = shm
            del shm
        return

    def _worker_close_memory(self):
        for k, v in self._shm_desc.items():
            try:
                self._shm_handles[k].close()
            except:
                pass
            finally:
                del self._shm_handles[k]
        return

    def _coordinator_unlink_memory(self):
        for k, v in self._shm_desc.items():
            try:
                self.__setattr__(k, self.__getattr__(k).copy())
                self._shm_handles[k].close()
                self._shm_handles[k].unlink()
            except:
                pass
            finally:
                del self._shm_handles[k]
        del self._shm_desc
        self._shm_desc = {}
        del self._shm_handles
        self._shm_handles = {}
        return
    # /Shared memory handlers

    @abstractmethod
    def _create_init_dict(self):
        raise NotImplementedError

@register_model
class Cloud(AbstractModel):
    """Cloud model class inheriting from AbstractModel.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file,
            run_type
        )
        if self.run_type != "external":
            self.run()
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        if self._debug:
            print("got to run_fortran")
            print(f"using "
                  f"dictionary={self._param_dict},"
                  f"outspeciesin={self.out_species},"
                  f"timepoints={self.timepoints},"
                  f"gridpoints={self._param_dict['points']},"
                  f"returnarray={True},"
                  f"returnrates={True},"
                  f"givestartabund={self.give_start_abund},"
                  f"physicsarray={self.physics_array},"
                  f"chemicalabunarray={self.chemical_abun_array},"
                  f"ratesarray={self.rates_array},"
                  f"abundancestart={self.starting_chemistry_array},"
                  )
        _, _, _, out_species_abundances_array, _, success_flag = wrap.cloud(
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            ratesarray=self.rates_array,
            abundancestart=self.starting_chemistry_array if "starting_chemistry_array" in object.__getattribute__(self,"__dict__") else None,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = len(self.out_species_list) if self.out_species_list is not None else 0
            out_species_abundances_array = list(out_species_abundances_array[:out_species_length])
        return {"success_flag": success_flag, "out_species_abundances_array": out_species_abundances_array}

    def _create_init_dict(self):
        return {
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug
        }

@register_model
class Collapse(AbstractModel):
    """Collapse model class inheriting from AbstractModel.

    Args:
        collapse (str, optional):A string containing the collapse type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'.
            Defaults to 'BE1.1'.
        physics_output (str, optional): Filename to store physics output, only relevant for 'filament' and 'ambipolar' collapses.
            If None, no physics output will be saved.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 collapse: str = "BE1.1",
                 physics_output: str = None,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        if collapse not in ["filament", "ambipolar"] and physics_output is None:
            warnings.warn("`physics_output` was None but `collapse` was `filament` or `ambipolar`. No output file will be created.", UserWarning)
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file,
            run_type
        )
        if read_file is None:
            collapse_dict = {"BE1.1": 1, "BE4": 2, "filament": 3, "ambipolar": 4}
            try:
                self.collapse = collapse_dict[collapse]
            except KeyError:
                raise ValueError(
                    f"collapse must be in {collapse_dict.keys()}"
                )
            self.physics_output = physics_output
            self.write_physics = self.physics_output is not None
            if not self.write_physics:
                self.physics_output = ""
            if self.run_type != "external":
                self.run()
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, out_species_abundances_array, _, success_flag = wrap.collapse(
            collapseIn=self.collapse,
            collapseFileIn=self.physics_output,
            writeOut=self.write_physics,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            ratesarray=self.rates_array,
            abundancestart=self.starting_chemistry_array if "starting_chemistry_array" in object.__getattribute__(self,"__dict__") else None,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = len(self.out_species_list) if self.out_species_list is not None else 0
            out_species_abundances_array = list(out_species_abundances_array[:out_species_length])
        return {"success_flag": success_flag, "out_species_abundances_array": out_species_abundances_array}

    def _create_init_dict(self):
        return {
            "collapse": self.collapse,
            "physics_output": self.physics_output,
            "write_physics": self.write_physics,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug
        }

@register_model
class PrestellarCore(AbstractModel):
    """PrestellarCore model class inheriting from AbstractModel. This model type was previously known as hot core.

    Args:
        temp_indx (int, optional): Used to select the mass of the prestellar core from the following selection
            [1=1Msun, 2=5, 3=10, 4=15, 5=25,6=60]. Defaults to 1, which is 1 Msun
        max_temperature (float, optional): Value at which gas temperature will stop increasing. Defaults to 300.0.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 temp_indx: int = 1,
                 max_temperature: float = 300.0,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file,
            run_type
        )
        if read_file is None:
            if temp_indx is None or max_temperature is None:
                raise ("temp_indx and max_temperature must be specified if not reading from file.")
            self.temp_indx = temp_indx
            self.max_temperature = max_temperature
            if self.run_type != "external":
                self.run()
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, out_species_abundances_array, _, success_flag = wrap.hot_core(
            temp_indx=self.temp_indx,
            max_temp=self.max_temperature,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            ratesarray=self.rates_array,
            abundancestart=self.starting_chemistry_array if "starting_chemistry_array" in object.__getattribute__(self,"__dict__") else None,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = len(self.out_species_list) if self.out_species_list is not None else 0
            out_species_abundances_array = list(out_species_abundances_array[:out_species_length])
        return {"success_flag": success_flag, "out_species_abundances_array": out_species_abundances_array}

    def _create_init_dict(self):
        return {
            "temp_indx": self.temp_indx,
            "max_temperature": self.max_temperature,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug
        }

@register_model
class CShock(AbstractModel):
    """C-Shock model class inheriting from AbstractModel.

    Args:
        shock_vel (float, optional): Velocity of the shock in km/s. Defaults to 10.0.
        timestep_factor (float, optional): Whilst the time is less than 2 times the dissipation time of shock,
            timestep is timestep_factor*dissipation time. Essentially controls how well resolved the shock is
            in your model. Defaults to 0.01.
        minimum_temperature (float, optional): Minimum post-shock temperature. Defaults to 0.0 (no minimum). The
            shocked gas typically cools to `initialTemp` if this is not set.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 shock_vel: float = 10.0,
                 timestep_factor: float = 0.01,
                 minimum_temperature: float = 0.0,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file,
            run_type
        )
        if read_file is None:
            if shock_vel is None:
                raise ("shock_vel must be specified if not reading from file.")
            self.shock_vel = shock_vel
            self.timestep_factor = timestep_factor
            self.minimum_temperature = minimum_temperature
            self.dissipation_time = -1
            if self.run_type != "external":
                self.run()
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, out_species_abundances_array, dissipation_time, _, success_flag = wrap.cshock(
            shock_vel=self.shock_vel,
            timestep_factor=self.timestep_factor,
            minimum_temperature=self.minimum_temperature,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            ratesarray=self.rates_array,
            abundancestart=self.starting_chemistry_array,
        )
        if success_flag < 0:
            dissipation_time = None
            out_species_abundances_array = np.array([])
        else:
            out_species_length = len(self.out_species_list) if self.out_species_list is not None else 0
            out_species_abundances_array = list(out_species_abundances_array[:out_species_length])
        return {"success_flag": success_flag, "dissipation_time": dissipation_time, "out_species_abundances_array": out_species_abundances_array}

    def _create_init_dict(self):
        return {
            "shock_vel": self.shock_vel,
            "timestep_factor": self.timestep_factor,
            "minimum_temperature": self.minimum_temperature,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug
        }

@register_model
class JShock(AbstractModel):
    """J-Shock model class inheriting from AbstractModel.

    Args:
        shock_vel (float): Velocity of the shock. Defaults to 10.0.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 shock_vel: float = 10.0,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file,
            run_type
        )
        if read_file is None:
            if shock_vel is None:
                raise ("shock_vel must be specified if not reading from file.")
            self.shock_vel = shock_vel
            if self.run_type != "external":
                self.run()
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, out_species_abundances_array, _, success_flag = wrap.jshock(
            shock_vel=self.shock_vel,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            ratesarray=self.rates_array,
            abundancestart=self.starting_chemistry_array,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = len(self.out_species_list) if self.out_species_list is not None else 0
            out_species_abundances_array = list(out_species_abundances_array[:out_species_length])
        return {"success_flag": success_flag, "out_species_abundances_array": out_species_abundances_array}

    def _create_init_dict(self):
        return {
            "shock_vel": self.shock_vel,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug
        }

@register_model
class Postprocess(AbstractModel):
    """Postprocess represents a model class with additional controls. It inherits from AbstractModel.

    Postprocess allows for additional controls of the time, density, gas temperature, radiation field, cosmic ray
    ionisation rate, atomic and molecular Hydrogen, CO and C column densities through the use of arrays. Using these
    arrays allows for experimental model crafting beyond the standard models in other model classes.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        time_array (np.array, optional): Represents the time grid to be used for the model. This sets the target timesteps for
            which outputs will be stored.
        density_array (np.array, optional): Represents the value of the density at different timepoints found in time_array.
        gas_temperature_array (np.array, optional): Represents the value of the gas temperature at different timepoints found in time_array.
        dust_temperature_array (np.array, optional):Represents the value of the dust temperature at different timepoints found in time_array.
        zeta_array (np.array, optional): Represents the value of the cosmic ray ionisation rate at different timepoints found in time_array.
        radfield_array (np.array, optional): Represents the value of the UV radiation field at different timepoints found in time_array.
        coldens_H_array (np.array, optional): Represents the value of the column density of H at different timepoints found in time_array.
        coldens_H2_array (np.array, optional): Represents the value of the column density of H2 at different timepoints found in time_array.
        coldens_CO_array (np.array, optional): Represents the value of the column density of CO at different timepoints found in time_array.
        coldens_C_array (np.array, optional): Represents the value of the column density of C at different timepoints found in time_array.
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 time_array: np.array = None,
                 density_array: np.array = None,
                 gas_temperature_array: np.array = None,
                 dust_temperature_array: np.array = None,
                 zeta_array: np.array = None,
                 radfield_array: np.array = None,
                 coldens_H_array: np.array = None,
                 coldens_H2_array: np.array = None,
                 coldens_CO_array: np.array = None,
                 coldens_C_array: np.array = None,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            len(time_array),
            debug,
            read_file,
            run_type
        )
        if read_file is None and time_array is not None:
            self.postprocess_arrays = dict(
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
            for key, array in self.postprocess_arrays.items():
                if array is not None:
                    # Convert single values into arrays that can be used
                    if isinstance(array, float):
                        array = np.ones(shape=time_array.shape) * array
                    # Assure lengths are correct
                    assert len(array) == len(time_array), "All arrays must be the same length"
                    # Ensure Fortran memory
                    array = np.asfortranarray(array, dtype=np.float64)
                    self.postprocess_arrays[key] = array
            self.time_array = time_array
            self.coldens_H_array = coldens_H_array
            self.usecoldens = self.coldens_H_array is not None
            if not self.give_start_abund:
                self.starting_chemistry_array = np.zeros(
                    shape=n_species,
                    dtype=np.float64,
                    order="F",
                )
            if self.run_type != "external":
                self.run()
        elif time_array is None and read_file is None:
            raise ValueError(f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given.")
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, out_species_abundances_array, _, success_flag = wrap.postprocess(
            usecoldens=self.usecoldens,
            **self.postprocess_arrays,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            ratesarray=self.rates_array,
            abundancestart=self.starting_chemistry_array,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = len(self.out_species_list) if self.out_species_list is not None else 0
            out_species_abundances_array = list(out_species_abundances_array[:out_species_length])
        return {"success_flag": success_flag, "out_species_abundances_array": out_species_abundances_array}

    def _create_init_dict(self):
        return {
            "usecoldens": self.coldens_H_array is not None,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug,
            **self.postprocess_arrays,
        }

@register_model
class Model(AbstractModel):
    """Model, like Postprocess, represents a model class with additional controls. It inherits from AbstractModel.

    Model follows the same logic as Postprocess but without the coldens Arguments. It allows for additional controls of
    the time, density, gas temperature, radiation field, and cosmic ray ionisation rate through the use of arrays.
    Using these arrays allows for experimental model crafting beyond the standard models in other model classes.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        time_array (np.array, optional): Represents the time grid to be used for the model. This sets the target timesteps for
            which outputs will be stored. While listed as optional, this is only done to allow
        density_array (np.array, optional): Represents the value of the density at different timepoints found in time_array.
        gas_temperature_array (np.array, optional): Represents the value of the gas temperature at different timepoints found in time_array.
        dust_temperature_array (np.array, optional):Represents the value of the dust temperature at different timepoints found in time_array.
        zeta_array (np.array, optional): Represents the value of the cosmic ray ionisation rate at different timepoints found in time_array.
        radfield_array (np.array, optional): Represents the value of the UV radiation field at different timepoints found in time_array.
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 time_array: np.array = None,
                 density_array: np.array = None,
                 gas_temperature_array: np.array = None,
                 dust_temperature_array: np.array = None,
                 zeta_array: np.array = None,
                 radfield_array: np.array = None,
                 debug: bool = False,
                 read_file: str = None,
                 run_type: Literal["local", "managed", "external"] = "local"
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            len(time_array),
            debug,
            read_file,
            run_type
        )
        if read_file is None and time_array is not None:
            self.time_array = time_array
            self.postprocess_arrays = dict(
                timegrid=time_array,
                densgrid=density_array,
                gastempgrid=gas_temperature_array,
                dusttempgrid=dust_temperature_array,
                radfieldgrid=radfield_array,
                zetagrid=zeta_array
            )
            for key, array in self.postprocess_arrays.items():
                if array is not None:
                    # Convert single values into arrays that can be used
                    if isinstance(array, float):
                        array = np.ones(shape=time_array.shape) * array
                    # Assure lengths are correct
                    assert len(array) == len(time_array), "All arrays must be the same length"
                    # Ensure Fortran memory
                    array = np.asfortranarray(array, dtype=np.float64)
                    self.postprocess_arrays[key] = array
            if not self.give_start_abund:
                self.starting_chemistry_array = np.zeros(
                    shape=n_species,
                    dtype=np.float64,
                    order="F",
                )
            if self.run_type != "external":
                self.run()
        elif time_array is None and read_file is None:
            raise ValueError(f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given.")
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, out_species_abundances_array, _, success_flag = wrap.postprocess(
            usecoldens=False,
            **self.postprocess_arrays,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            ratesarray=self.rates_array,
            abundancestart=self.starting_chemistry_array,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = len(self.out_species_list) if self.out_species_list is not None else 0
            out_species_abundances_array = list(out_species_abundances_array[:out_species_length])
        return {"success_flag": success_flag, "out_species_abundances_array": out_species_abundances_array}

    def _create_init_dict(self):
        return {
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug,
            **self.postprocess_arrays,
        }

@register_model
class SequentialModel:
    def __init__(
            self,
            model_type_sequence: List[AnyStr],
            joined_parameter_sequence: List[Dict],
            parameters_to_match: List = None
    ):
        self.models = {}
        for model in model_type_sequence:
            assert(model != SequentialModel)

        for model_number in range(len(model_type_sequence)):
            cls = REGISTRY.get(model_type_sequence[model_number])
            if model_number > 0 and parameters_to_match is not None:
                joined_parameter_sequence[model_number]["_param_dict"] = \
                    {**joined_parameter_sequence[model_number - 1]["_param_dict"],
                     **joined_parameter_sequence[model_number]["_param_dict"]}
                print("Parameter matching to be implemented")
                for parameter in parameters_to_match:
                    if parameter == "finalDens":
                        joined_parameter_sequence[model_number]["_param_dict"]["initialDens"] = (
                            self.models[f"{model_number-1}"].physics_array[-1,0,1]
                        )
                        continue
                    elif parameter == "finalTemp":
                        joined_parameter_sequence[model_number]["_param_dict"]["initialTemp"] = (
                            self.models[f"{model_number-1}"].physics_array[-1, 0, 2]
                        )
                    else:
                        print(f"Parameter '{parameter}' has not been implemented for parameter matching")
                self.models[f"{model_number}"] = cls(**joined_parameter_sequence[model_number],
                                    previous_model=self.models[f"{model_number-1}"])
            if model_number > 0:
                joined_parameter_sequence[model_number]["_param_dict"] = \
                    {**joined_parameter_sequence[model_number-1]["_param_dict"],
                     **joined_parameter_sequence[model_number]["_param_dict"]}
                self.models[f"{model_number}"] = cls(**joined_parameter_sequence[model_number],
                                     previous_model=self.models[f"{model_number-1}"])
            else:
                self.models[f"{model_number}"] = cls(**joined_parameter_sequence[model_number])
        return

    def save_models(self, file_name, sequence_name: str = 'default', overwrite: bool = False):
        for k, v in self.models.items():
            v.save_model(file_name, f"{sequence_name}_{k}", overwrite)
        return


def run_grid_model(model_id, model_type, param_dict, pending_model):
    """Run a single model (called by compute workers)"""
    cls = REGISTRY.get(model_type)
    model_obj = cls(
        param_dict=param_dict,
        **pending_model,
        run_type="external"
    )
    model_obj.run()
    model_obj._coordinator_unlink_memory()
    model_obj.pickle()
    return (model_id, model_obj)


class GridModels:
    def __init__(self,
                 model_type: AnyStr,
                 full_parameters: Dict,
                 max_workers: int = 8,
                 grid_file: str = "./default_grid_out.h5"
                 ):
        self.model_type = model_type
        self.full_parameters = full_parameters
        self.max_workers = max_workers-1 if max_workers < int(os.cpu_count()/2) else int(os.cpu_count()/2)-1
        self.grid_file = grid_file if ".h5" in grid_file else grid_file+".h5"

        self._orig_sigint = signal.getsignal(signal.SIGINT)

        self.parameters_to_grid = {}
        for k, v in self.full_parameters.items():
            if type(v) == list:
                self.parameters_to_grid[k] = v
                self.parameters_to_grid[k] = np.array(v, dtype=object)
        grids = np.meshgrid(*self.parameters_to_grid.values(), indexing="xy")
        self.flat_grids = np.reshape(grids,
                                        (
                                            len(self.parameters_to_grid),
                                            int(np.prod(np.shape(grids)) / len(self.parameters_to_grid))
                                        )
                                     )
        self.models = {}
        self.run()
        return

    def run(self):
        signal.signal(signal.SIGINT, self._handler)
        pending = self.grid_iter(self.full_parameters, list(self.parameters_to_grid.keys()), self.flat_grids)
        with mp.Pool(self.max_workers) as pool:
            active = 0
            submitted = 0
            completed = 0
            total_models = len(self.flat_grids[0])

            def submit_next() -> bool:
                nonlocal active, submitted
                try:
                    pending_model = next(pending)
                except StopIteration:
                    return False

                active += 1
                submitted += 1
                model_id = pending_model.pop("id")

                param_dict = {}
                for key in list(pending_model.keys()):
                    if key.lower() in default_param_dictionary:
                        param_dict[key] = pending_model.pop(key)
                try:
                    pool.apply_async(
                        run_grid_model,
                        args=(model_id, self.model_type, param_dict, pending_model),
                        callback=lambda result: on_result(result),
                        error_callback=lambda exc, _model_id=model_id: on_error(exc, model_id)
                    )
                except Exception as e:
                    print(f"exception found in submit_next: {e}")
                return True

            def on_result(result):
                nonlocal active, completed
                active -= 1
                completed += 1
                model_id, model_object = result

                try:
                    save_name = f'grid_model_{model_id}'
                    print(f"saving model: '{save_name}'")
                    model_object.un_pickle()
                    model_object.save_model(
                        file=self.grid_file,
                        name=save_name,
                        engine="h5netcdf",
                        overwrite=True
                    )
                    # Update reference to just the name
                    self.models[model_id] = save_name
                except Exception as e:
                    print(f"Error saving model {model_id}: {e}")
                    import traceback
                    traceback.print_exc()

                if not submit_next() and active == 0:
                    print("No more models to submit")

            def on_error(_exc, _model_id):
                nonlocal active
                active -= 1
                # Optionally record/log the exception here.
                self.models[_model_id]._coordinator_unlink_memory()
                print(f"on_error: {_exc}; model_id: {_model_id}")
                if not submit_next() and active == 0:
                    print("No more models to submit")

            # Submit initial batch
            for _ in range(np.shape(self.flat_grids)[1]):
                if not submit_next():
                    break

            pool.close()
            pool.join()

        signal.signal(signal.SIGINT, self._orig_sigint)
        print(f"GridModels complete. Output: {self.grid_file}")

        #with ProcessPoolExecutor(max_workers=self.max_workers) as p:
        #    iteration = self.grid_iter(self.full_parameters, list(self.parameters_to_grid.keys()), self.flat_grids)
        #    futures = {p.submit(self.run_grid_model, x): x for x in itertools.islice(iteration, self.max_workers)}
        #    results = []
        #    for fut in as_completed(futures):
        #        results.append(fut.result())
        #        try:
        #            x = next(iteration)  # get next item
        #            futures[p.submit(self.run_grid_model, x)] = x
        #        except StopIteration:
        #            pass
        #with ProcessPoolExecutor(max_workers=self.max_workers) as p:
        #    models_list = p.map(self.run_grid_model, list(self.grid_iter(self.full_parameters, list(self.parameters_to_grid.keys()), self.flat_grids)))
        #for model in models_list:
        #    key = next(iter(model))
        #    self.models[key] = model[key]

    def on_interrupt(self):
        return

    def _handler(self, signum, frame):
        try:
            self.on_interrupt()  # your “final steps”
        finally:
            # Restore default and re-raise KeyboardInterrupt to stop execution
            signal.signal(signal.SIGINT, self._orig_sigint)
            raise KeyboardInterrupt

    @staticmethod
    def grid_iter(full_parameters: dict, param_keys: list, flattened_grids: np.ndarray) -> Iterator[Dict[str, Any]]:
        for i in range(len(flattened_grids[0])):
            combo = ()
            for j in range(np.shape(flattened_grids)[0]):
                combo += (flattened_grids[j][i],)
            grid_dict = {k: v for k, v in zip(param_keys, combo)}
            yield {**full_parameters, **grid_dict, "id": i}


'''
Functional Submodule. 

The following functions are wrappers of the above classes. The format of the wrappers has been chosen such that they follow
the legacy methods of running UCLCHEM both in inputs and in return values, while depending on the Class objects defined
above to facilitate ease of maintenance.
'''

def __functional_return__(model_object: AbstractModel,
                          return_array: bool = False,
                          return_dataframe: bool = False,
                          return_rates: bool = False):
    """
    return function that takes in the object that was modelled and returns the values based on the specified booleans.

    Args:
        model_object: model_object of a class that inherited from AbstractModel, from which the results should be returned.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, the function will return
            the success_flag, dissipation_time if the model_object has that attribute, and the final abundances of the out_species.
        return_dataframe: A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, the function will return
            the success_flag, dissipation_time if the model_object has that attribute, and the final abundances of the out_species.
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, and the model_object has the dissipation_time attribute the second element is the dissipation time. Further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - dissipation_time (float): dissipation time in years (if model_object contains the dissipation_time attribute)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - dissipation_time (float): dissipation time in years (if model_object contains the dissipation_time attribute)
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.

    """
    if return_dataframe:
        if return_rates:
            phys_df, chem_df, rates_df = model_object.get_dataframes(joined=False, with_rates=return_rates)
        else:
            phys_df, chem_df = model_object.get_dataframes(joined=False, with_rates=return_rates)
            rates_df = None

        if hasattr(model_object, 'dissipation_time'):
            return (
                phys_df,
                chem_df,
                rates_df,
                model_object.dissipation_time,
                model_object.next_starting_chemistry_array,
                model_object.success_flag,
            )
        else:
            return (
                phys_df,
                chem_df,
                rates_df,
                model_object.next_starting_chemistry_array,
                model_object.success_flag,
            )
    elif return_array:
        if hasattr(model_object, 'dissipation_time'):
            return (
                model_object.physics_array,
                model_object.chemical_abun_array,
                model_object.rates_array if return_rates else None,
                model_object.dissipation_time,
                model_object.next_starting_chemistry_array,
                model_object.success_flag
            )
        else:
            return (
                model_object.physics_array,
                model_object.chemical_abun_array,
                model_object.rates_array if return_rates else None,
                model_object.next_starting_chemistry_array,
                model_object.success_flag
            )
    else:
        if hasattr(model_object, 'dissipation_time'):
            return [model_object.success_flag, model_object.dissipation_time] + model_object.out_species_abundances_array
        else:
            return [model_object.success_flag] + model_object.out_species_abundances_array


def __cloud__(
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        return_array: bool = False,
        return_dataframe: bool = False,
        return_rates: bool = False,
        starting_chemistry: np.array = None,
        timepoints: int = TIMEPOINTS):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file.
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file.
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem.
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS
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

    model_object = Cloud(
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints
    )

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rates=return_rates
    )


def __collapse__(
        collapse: str,
        physics_output: str,
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        return_array: bool = False,
        return_dataframe: bool = False,
        return_rates: bool = False,
        starting_chemistry: np.array = None,
        timepoints: int = TIMEPOINTS):
    """Run collapse model from UCLCHEM based on Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)

    Args:
        collapse (str): A string containing the collapse type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'
        physics_output(str): Filename to store physics output, only relevant for 'filament' and 'ambipolar' collapses. If None, no physics output will be saved.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

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

    model_object = Collapse(
        collapse=collapse,
        physics_output=physics_output,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints
    )

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rates=return_rates
    )


def __prestellar_core__(
        temp_indx: int,
        max_temperature: float,
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        return_array: bool = False,
        return_dataframe: bool = False,
        return_rates: bool = False,
        starting_chemistry: np.array = None,
        timepoints: int = TIMEPOINTS):
    """Run prestellar core model from UCLCHEM, based on Viti et al. 2004 and Collings et al. 2004. This model type was previously known as hot core

    Args:
        temp_indx (int): Used to select the mass of prestellar core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        max_temperature (float): Value at which gas temperature will stop increasing.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

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

    model_object = PrestellarCore(
        temp_indx=temp_indx,
        max_temperature=max_temperature,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints
    )

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rates=return_rates
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
        starting_chemistry: np.array = None,
        timepoints: int = TIMEPOINTS):
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
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

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

    model_object = CShock(
        shock_vel=shock_vel,
        timestep_factor=timestep_factor,
        minimum_temperature=minimum_temperature,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints
    )

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rates=return_rates
    )


def __jshock__(
        shock_vel: float,
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        return_array: bool = False,
        return_dataframe: bool = False,
        return_rates: bool = False,
        starting_chemistry: np.array = None,
        timepoints: int = TIMEPOINTS):
    """Run J-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned. Defaults to ["H", "N", "C", "O"].
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_rates (bool, optional): A boolean on whether the reaction rates should be returned to a user.
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model. Defaults to uclchem.constants.TIMEPOINTS

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

    model_object = JShock(
        shock_vel=shock_vel,
        param_dict=param_dict,
        out_species=out_species,
        starting_chemistry=starting_chemistry,
        timepoints=timepoints
    )

    return __functional_return__(
        model_object=model_object,
        return_array=return_array,
        return_dataframe=return_dataframe,
        return_rates=return_rates
    )


functional = types.ModuleType(__name__ + ".functional")
functional.cloud = __cloud__
functional.collapse = __collapse__
functional.prestellar_core = __prestellar_core__
functional.cshock = __cshock__
functional.jshock = __jshock__
sys.modules[__name__ + ".functional"] = functional