import json
import logging

# /UCLCHEM related imports
# Multiprocessing imports
import multiprocessing as mp
import os
import signal
import warnings
from abc import ABC, abstractmethod
from datetime import datetime
from multiprocessing import shared_memory
from pathlib import Path
from typing import Any, AnyStr, Dict, Iterator, List, Literal, Type

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uclchemwrap
import xarray as xr

# UCLCHEM related imports
from uclchemwrap import uclchemwrap as wrap

from uclchem.analysis import (
    check_element_conservation,
    create_abundance_plot,
    plot_species,
)
from uclchem.constants import (
    N_PHYSICAL_PARAMETERS,
    PHYSICAL_PARAMETERS,
    TIMEPOINTS,
    default_param_dictionary,
    n_reactions,
    n_species,
)

# /Multiprocessing imports

# Global variables determining formats of write files
PHYSICAL_PARAMETERS_HEADER_FORMAT = "%10s"
# in the below variable, the outputs were chosen according to the spacing needed for
# "      Time,    Density,    gasTemp,   dustTemp,         Av,   radfield,       zeta,      point,"
PHYSICAL_PARAMETERS_VALUE_FORMAT = (
    "%10.3E, %10.4E, %10.2f, %10.2f, %10.4E, %10.4E, %10.4E, %10i"
)
SPECNAME_HEADER_FORMAT = "%11s"
SPECNAME_VALUE_FORMAT = "%9.5E"
#


# Set collisional rates directory for heating/cooling calculations
def set_collisional_rates_directory():
    # TODO: move this functionality into the advanced heating suite.
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
        actual_dir = str(
            np.char.decode(uclchemwrap.defaultparameters.coolantdatadir)
        ).strip()
        assert actual_dir == coolant_directory, (
            f"Coolant directory path not set correctly. "
            f"Expected: {coolant_directory} "
            f"Got: {actual_dir}"
        )
    except AttributeError:
        logging.warning(
            "Cannot set the coolant directory path, please set 'coolantDataDir' correctly at runtime."
        )


# Call set_collisional_rates_directory to initialize heating/cooling support
try:
    set_collisional_rates_directory()
except Exception as e:
    logging.warning(f"Could not set collisional rates directory: {e}")

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


get_reaction_names = ReactionNamesStore()

# Before doing anything else, set the right collision rate directory.
set_collisional_rates_directory()


class SpeciesNameStore:
    def __init__(self):
        self.species_names = None

    def __call__(self):
        # Only load the species once, after that use the cached version
        if self.species_names is None:
            species = pd.read_csv(
                os.path.join(os.path.dirname(os.path.abspath(__file__)), "species.csv")
            )
            self.species_names = species["NAME"].tolist()
        return self.species_names


get_species_names = SpeciesNameStore()
# /Reaction and Species name retrieval classes to reduce file read repetition.


# Universal model loader
def load_model(
    file: str, name: str = "default", engine: str = "h5netcdf", debug: bool = False
):
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
        raise ValueError(
            f"Unrecognized model type '{model_class}'. Not in trusted registry."
        )
    return cls.load_from_dataset(model_ds=loaded_data, debug=debug)


# /Universal model loader


# Worker entry for parallel jobs
def _worker_entry(
    model_class: str, init_kwargs: dict, shm_descs: dict, result_queue: mp.Queue
):
    cls = REGISTRY.get(model_class)
    if cls is None:
        raise ValueError(
            f"Unrecognized model type '{model_class}'. Not in trusted registry."
        )
    model = cls.worker_build(init_kwargs=init_kwargs, shm_desc=shm_descs)
    output = model.run_fortran()
    result_queue.put(output)
    return


# /Worker entry for parallel jobs


# Short compatibility helper for legacy parameter `endAtFinalDensity`
def _convert_legacy_stopping_param(param_dict: dict) -> dict:
    """Minimal conversion of legacy `endAtFinalDensity` to `parcelStoppingMode`.
    Rules (short and strict):
      - If both keys are present: raise RuntimeError
      - If `endAtFinalDensity` is present and points>1: raise RuntimeError
      - Otherwise convert True->1, False->2 and remove the old key
    
    Note: This function assumes param_dict is already a copy and is case-normalized (lowercase keys).
    """
    if param_dict is None:
        return param_dict
    
    has_old = 'endatfinaldensity' in param_dict
    has_new = 'parcelstoppingmode' in param_dict
    if has_old and has_new:
        raise RuntimeError("Cannot specify both 'endAtFinalDensity' and 'parcelStoppingMode'. Use 'parcelStoppingMode' only.")
    if has_old:
        points = param_dict.get('points', 1)
        if points > 1:
            raise RuntimeError("endAtFinalDensity is no longer supported for multi-point models (points > 1). Use 'parcelStoppingMode' instead.")
        old_val = param_dict.pop('endatfinaldensity')
        param_dict['parcelstoppingmode'] = 1 if old_val else 2
    return param_dict


# TODO Add catch of ctrl+c or other aborts so that it saves model and a full output to files of year, month, day, time type.
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

    def __init__(
        self,
        param_dict: dict = None,
        out_specie_list: list = ["H", "N", "C", "O"],
        starting_chemistry: np.ndarray = None,
        previous_model: object = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["local", "managed", "external"] = "local",
    ):
        self._data = xr.Dataset()
        self._pickle_dict = {}
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

        self.outputFile = (
            param_dict.pop("outputFile") if "outputFile" in param_dict else None
        )
        self.abundSaveFile = (
            param_dict.pop("abundSaveFile") if "abundSaveFile" in param_dict else None
        )
        self.abundLoadFile = (
            param_dict.pop("abundLoadFile") if "abundLoadFile" in param_dict else None
        )

        self.starting_chemistry_array = None
        if previous_model is None and self.abundLoadFile is None:
            self._create_starting_array(starting_chemistry)
        elif self.abundLoadFile is not None:
            self.legacy_read_starting_chemistry()
        elif previous_model.has_attr("next_starting_chemistry_array"):
            self._create_starting_array(previous_model.next_starting_chemistry_array)

        self.give_start_abund = self.starting_chemistry_array is not None
        self.next_starting_chemistry_array = None

        self.physics_array = None
        self.chemical_abun_array = None
        self._create_fortran_array()
        self.rates_array = None
        self._create_rates_array()
        self.heat_array = None
        self._create_heating_array()
        self.out_species_abundances_array = None

        if read_file is not None:
            self.read_output_file(read_file)
        return

    def __del__(self):
        if hasattr(self, "_shm_desc"):
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
        obj._coord_assign()
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
        if key.startswith("_") and key != "_data":
            return super().__getattribute__(key)
        elif key in super().__getattribute__("_data"):
            values = super().__getattribute__("_data")[key].values
            if np.shape(values) != ():
                if isinstance(values, tuple):
                    return values[1]
                else:
                    return values
            else:
                return self._data[key].item()
        else:
            raise AttributeError(
                f'{self.__class__.__name__} has no attribute of name: "{key}".'
            )

    def __setattr__(self, key, value):
        if key.startswith("_"):
            super().__setattr__(key, value)
        else:
            if key in self._data:
                self._data.__delitem__(key)

            if np.ndim(value) == 3 and "_array" in key:
                self._data[key] = (
                    ["time_step", "point", key.replace("array", "values")],
                    value,
                )
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
    def check_conservation(self, element_list: list = None, percent: bool = True):
        """Utility method to check conservation of the chemical abundances

        Args:
            element_list (list, optional): List of elements to check conservation for. Defaults to self.out_species_lists.
            percent (bool, optional): Flag on if percentage values should be used. Defaults to True.
        """
        if element_list is None:
            element_list = self.out_species_list

        if self._param_dict["points"] > 1:
            for i in range(self._param_dict["points"]):
                print(
                    f"Element conservation report for point {i + 1} of {self._param_dict['points']}"
                )
                print(
                    check_element_conservation(
                        self.get_dataframes(i), element_list, percent
                    )
                )
        else:
            print("Element conservation report")
            print(
                check_element_conservation(
                    self.get_dataframes(0), element_list, percent
                )
            )

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
                print(f"{errors[self.success_flag]}")
            except KeyError:
                raise ValueError(f"Unknown error code: {self.success_flag}")
        elif self.success_flag == 0 and not only_error:
            print("Model ran successfully.")
        elif self.success_flag is None:
            print("Model has not been run.")
        return

    def create_abundance_plot(
        self,
        species: list = None,
        figsize: tuple[2] = (16, 9),
        point: int = 0,
        plot_file=None,
    ):
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
        return create_abundance_plot(
            self.get_dataframes(point), species, figsize, plot_file
        )

    def get_dataframes(
        self,
        point: int = 0,
        joined: bool = True,
        with_rates: bool = False,
        with_heating: bool = False,
    ):
        """Converts the model physics and chemical_abun arrays from numpy to pandas arrays.
        Args:
            point (int, optional): Integer referring to which point of the UCLCHEM model to return as pandas does not support higher
                than 2D data structures. Defaults to 0.
            joined (bool, optional): Flag on whether the returned pandas dataframe should be one, or if two dataframes should be
                returned. One physical, one chemical_abun dataframe. Defaults to True.
            with_rates (bool, optional): Flag on whether to include reaction rates in the dataframe, and/or as a separate
                dataframe depending on the value of `joined`. Defaults to False.
            with_heating (bool, optional): Flag on whether to include heating/cooling rates in the dataframe, and/or as a separate
                dataframe depending on the value of `joined`. Defaults to False.
        Returns:
            return_df (pandas.DataFrame): Dataframe of the joined arrays for point 'point' if joined = True
            physics_df (pandas.DataFrame): Dataframe of the physical parameters for point 'point' if joined = False
            chemistry_df (pandas.DataFrame): Dataframe of the chemical abundances  for point 'point' if joined = False
            rates_df (pandas.DataFrame): Dataframe of the reaction rates  for point 'point' if joined = False and with_rates = True
            heating_df (pandas.DataFrame): Dataframe of the heating/cooling rates  for point 'point' if joined = False and with_heating = True
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

        if self.heat_array is not None and with_heating:
            # Create a heating dataframe dynamically using labels from heating.f90
            heating_columns = ["Time"]

            # Add cooling mechanism labels
            cooling_labels = [
                str(np.char.decode(label)).strip() + " Cooling"
                for label in uclchemwrap.heating.coolinglabels
            ]
            heating_columns.extend(cooling_labels)

            # Add line cooling labels with species names
            line_cooling_labels = [
                str(np.char.decode(label)).strip()
                for label in uclchemwrap.f2py_constants.coolantnames
            ]
            for label in line_cooling_labels:
                heating_columns.append(f"{label} Line Cooling")

            # Add heating mechanism labels
            heating_labels = [
                str(np.char.decode(label)).strip() + " Heating"
                for label in uclchemwrap.heating.heatinglabels
            ]
            heating_columns.extend(heating_labels)

            heating_columns.append("Chemical Heating")

            heating_df = pd.DataFrame(
                self.heat_array[:, point, :], index=None, columns=heating_columns
            )
        else:
            heating_df = None

        if joined:
            return_df = physics_df.join(chemistry_df)
            if with_rates and rates_df is not None:
                return_df = return_df.join(rates_df)
            if with_heating and heating_df is not None:
                return_df = return_df.join(heating_df)
            return return_df
        else:
            result = [physics_df, chemistry_df]
            if with_rates:
                result.append(rates_df)
            if with_heating:
                result.append(heating_df)
            return tuple(result)

    def plot_species(
        self,
        ax: plt.axes,
        species: list[str] = None,
        point: int = 0,
        legend: bool = True,
        **plot_kwargs,
    ):
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
        return plot_species(
            ax, self.get_dataframes(point), species, legend, **plot_kwargs
        )

    # /UCLCHEM utility and analysis wrappers

    # Methods to start run of model
    def run(self):
        """__run__ resets the Fortran arrays if the model was not read, allowing the arrays to be reused for new runs."""
        if self.was_read:
            raise RuntimeError("This model was read. It can not be run. ")
            self.physics_array = None
            self.chemical_abun_array = None
            self.ratesArray = None
            self.heatArray = None
            self._create_fortran_array()
            self._create_rates_array()
            self._create_heating_array()

        def _handler(signum, frame):
            try:
                self.on_interrupt()  # your “final steps”
            finally:
                # Restore default and re-raise KeyboardInterrupt to stop execution
                signal.signal(signal.SIGINT, self._orig_sigint)
                raise KeyboardInterrupt

        if self.run_type in self.shared_memory_types:
            signal.signal(signal.SIGINT, _handler)

        if self.run_type not in self.separate_worker_types:
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

        if hasattr(self, "_shm_handles"):
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
    def save_model(
        self,
        file: str,
        name: str = "default",
        engine: str = "h5netcdf",
        overwrite: bool = False,
    ):
        """
        save_model saves a model to a file on disk. Multiple models can be saved into the same file if different names are used to store them.

        Args:
            file (str): Path to a file to store models.
            name (str, optional): Name to use for the group of the object. Defaults to 'default'
            engine (str, optional): “netcdf4”, “h5netcdf” or “zarr”, depending on the engine to be used. Defaults to "h5netcdf".
            overwrite (bool, optional): Boolean on whether to overwrite pre-existing models, or error out. Defaults to False
        """
        # TODO: Allow for toggling of saving float64 or float32 for the arrays
        if os.path.isfile(file):
            with xr.open_datatree(
                filename_or_obj=file, engine=engine, phony_dims="sort"
            ) as tree:
                if "/" + name in tree.groups:
                    if not overwrite:
                        warnings.warn(
                            f"Model with name: `{name}` already exists in `{file}` but overwrite is set to False. Unable to save model."
                        )
                        return
                tree.close()

        temp_attribute_dict = {}
        for v in self._data.variables:
            if "_array" not in v and v not in ["_orig_sigint"]:
                if np.shape(self._data[v].values) != ():
                    if isinstance(self._data[v].values, tuple):
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

        self._data.to_netcdf(file, group=name, engine=engine, mode="a")
        return

    # /Model saving

    # Model Passing through Pickling
    def pickle(self):
        if self._data is not None and not bool(self._pickle_dict):
            for v in self._data.variables:
                if np.shape(self._data[v].values) != ():
                    if isinstance(self._data[v].values, tuple):
                        self._pickle_dict[v] = self._data[v].values[1].tolist()
                    else:
                        self._pickle_dict[v] = self._data[v].values.tolist()
                else:
                    self._pickle_dict[v] = self._data[v].item()
            self._data = None
        return self

    def un_pickle(self):
        if self._data is None and bool(self._pickle_dict):
            self._data = xr.Dataset()
            for k, v in self._pickle_dict.items():
                if np.ndim(v) == 3 and "_array" in k:
                    self._data[k] = (
                        ["time_step", "point", k.replace("array", "values")],
                        v,
                    )
                elif np.ndim(v) == 2 and "_array" in k:
                    self._data[k] = (["point", k], v)
                elif "_values" in k:
                    pass
                else:
                    self._data[k] = v
            self._pickle_dict = {}
        self._coord_assign()
        return self

    # /Model Passing through Pickling

    # Legacy in & output support
    def legacy_read_output_file(self, read_file: str, rates_load_file: str = None):
        """Perform classic output file reading.
        Args:
            read_file (str): path to file containing a full UCLCHEM output
            rates_load_file (str, optional): path to file containing the reaction rates output from UCLCHEM. Defaults
                to None. #TODO Add the code to read the rates files.
        """
        self.was_read = True
        columns = np.char.strip(
            np.loadtxt(read_file, delimiter=",", max_rows=1, dtype=str, comments="%")
        )
        array = np.loadtxt(read_file, delimiter=",", skiprows=1)
        point_index = np.where(columns == "point")[0][0]
        self._param_dict["points"] = int(np.max(array[:, point_index]))
        if self._param_dict["points"] > 1:
            array = np.loadtxt(read_file, delimiter=",", skiprows=2)
        row_count = int(np.shape(array)[0] / self._param_dict["points"])

        self.PHYSICAL_PARAMETERS = [p for p in PHYSICAL_PARAMETERS if p in columns]
        self.specname = get_species_names()

        self.physics_array = np.empty(
            (row_count, self._param_dict["points"], len(self.PHYSICAL_PARAMETERS) + 1)
        )
        self.chemical_abun_array = np.empty(
            (row_count, self._param_dict["points"], len(self.specname))
        )
        for p in range(self._param_dict["points"]):
            self.physics_array[:, p, :] = array[
                np.where(array[:, point_index] == p + 1),
                : len(self.PHYSICAL_PARAMETERS) + 1,
            ][0]
            self.chemical_abun_array[:, p, :] = array[
                np.where(array[:, point_index] == p + 1),
                (len(self.PHYSICAL_PARAMETERS) + 1) :,
            ][0]
        self._array_clean()
        nonzero_indices = self.physics_array[:, 0, 0].nonzero()[0]
        if len(nonzero_indices) > 0:
            last_timestep_index = nonzero_indices[-1]
            self.next_starting_chemistry_array = self.chemical_abun_array[
                last_timestep_index, :, :
            ]
        else:
            # Model failed immediately, no valid timesteps
            self.next_starting_chemistry_array = None
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
        # TODO Move away from the magic numbers seen here.
        string_fmt_string = f'{", ".join([PHYSICAL_PARAMETERS_HEADER_FORMAT] * (len(self.PHYSICAL_PARAMETERS)))}, {", ".join([SPECNAME_HEADER_FORMAT] * len(self.specname))}'
        # Magic numbers here to match/improve the formatting of the classic version
        # TODO Move away from the magic numbers seen here.
        number_fmt_string = f'{PHYSICAL_PARAMETERS_VALUE_FORMAT}, {", ".join([SPECNAME_VALUE_FORMAT] * len(self.specname))}'
        columns = np.array(
            [
                self.PHYSICAL_PARAMETERS[:-1].tolist()
                + ["point"]
                + self.specname.tolist()
            ]
        )
        np.savetxt(self.outputFile, columns, fmt=string_fmt_string)
        with open(self.outputFile, "ab") as f:
            np.savetxt(f, full_array, fmt=number_fmt_string)
        return

    def legacy_write_starting_chemistry(self):
        """Perform classic starting abundance file writing to the file self.abundSaveFile provided in _param_dict"""
        last_timestep_index = self.chemical_abun_array[:, 0, 0].nonzero()[0][-1]
        # TODO Move away from the magic numbers seen here.
        number_fmt_string = f' {", ".join(["%9.5E"] * len(self.specname))}'
        with open(self.abundSaveFile, "wb") as f:
            np.savetxt(
                f,
                self.chemical_abun_array[last_timestep_index, :, :],
                fmt=number_fmt_string,
            )
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

        nonzero_indices = self.physics_array[:, 0, 0].nonzero()[0]
        if len(nonzero_indices) == 0:
            # Model failed immediately, keep only the first row to indicate failure
            self._data = self._data.isel(time_step=slice(0, 1))
            self._coord_assign()
            self.next_starting_chemistry_array = self.chemical_abun_array[0, :, :]
            return

        last_timestep_index = nonzero_indices[-1]
        # Get the arrays for only the simulated timesteps (not the zero padded ones)
        self._data = self._data.isel(time_step=slice(0, last_timestep_index + 1))
        self._coord_assign()
        self.next_starting_chemistry_array = self.chemical_abun_array[
            last_timestep_index, :, :
        ]
        return

    def _coord_assign(self):
        self._data = self._data.assign_coords(
            {"physics_values": self.PHYSICAL_PARAMETERS}
        )
        self._data = self._data.assign_coords({"chemical_abun_values": self.specname})
        return

    def _reform_inputs(self, param_dict: dict, out_species: list):
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
            
            # Handle deprecated endAtFinalDensity parameter (after lowercasing)
            new_param_dict = _convert_legacy_stopping_param(new_param_dict)
            
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
            (
                self._shm_handles["physics_array"],
                self._shm_desc["physics_array"],
                self.physics_array,
            ) = self._create_shared_memory_allocation(
                (self.timepoints + 1, self._param_dict["points"], N_PHYSICAL_PARAMETERS)
            )
            (
                self._shm_handles["chemical_abun_array"],
                self._shm_desc["chemical_abun_array"],
                self.chemical_abun_array,
            ) = self._create_shared_memory_allocation(
                (self.timepoints + 1, self._param_dict["points"], n_species)
            )
        else:
            # fencepost problem, need to add 1 to timepoints to account for the 0th timestep
            self.physics_array = np.zeros(
                shape=(
                    self.timepoints + 1,
                    self._param_dict["points"],
                    N_PHYSICAL_PARAMETERS,
                ),
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
            (
                self._shm_handles["rates_array"],
                self._shm_desc["rates_array"],
                self.rates_array,
            ) = self._create_shared_memory_allocation(
                (self.timepoints + 1, self._param_dict["points"], n_reactions)
            )
        else:
            self.rates_array = np.zeros(
                shape=(self.timepoints + 1, self._param_dict["points"], n_reactions),
                dtype=np.float64,
                order="F",
            )
        return

    def _create_heating_array(self):
        """Internal Method.
        Creates Fortran compliant np.array for heating/cooling rates that can be passed to the Fortran part of UCLCHEM.
        """
        try:
            heating_array_size = (
                2
                + uclchemwrap.heating.ncooling
                + uclchemwrap.heating.nheating
                + uclchemwrap.f2py_constants.ncoolants
            )
            # For shared memory:
            if self.run_type in self.shared_memory_types:
                (
                    self._shm_handles["heat_array"],
                    self._shm_desc["heat_array"],
                    self.heat_array,
                ) = self._create_shared_memory_allocation(
                    (
                        self.timepoints + 1,
                        self._param_dict["points"],
                        heating_array_size,
                    )
                )
            else:
                self.heat_array = np.zeros(
                    shape=(
                        self.timepoints + 1,
                        self._param_dict["points"],
                        heating_array_size,
                    ),
                    dtype=np.float64,
                    order="F",
                )
        except AttributeError:
            # Heating module not available, likely compiled without heating support
            logging.debug("Heating module not available in uclchemwrap")
            self.heat_array = None
        return

    def _create_starting_array(self, starting_chemistry):
        if starting_chemistry is None:
            self.starting_chemistry_array = None
        else:
            if len(np.shape(starting_chemistry)) == 1:
                starting_chemistry = starting_chemistry[np.newaxis, :]
            # For shared memory:
            if self.run_type in self.shared_memory_types:
                (
                    self._shm_handles["starting_chemistry_array"],
                    self._shm_desc["starting_chemistry_array"],
                    self.starting_chemistry_array,
                ) = self._create_shared_memory_allocation(np.shape(starting_chemistry))
                np.copyto(
                    self.starting_chemistry_array, starting_chemistry, casting="no"
                )
            else:
                self.starting_chemistry_array = np.asfortranarray(
                    starting_chemistry, dtype=np.float64
                )
        return

    # /Creation of arrays

    # Signal Interrupt Catch
    def on_interrupt(self, grid=False, model_name=None):
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
                self.outputFile = (
                    self.outputFile[: self.outputFile.rfind("/") + 1]
                    + error_time
                    + self.outputFile[self.outputFile.rfind(".") :]
                )
            else:
                self.outputFile = "./" + error_time + ".dat"
            self.legacy_write_full()

        self._was_interrupted = True
        self.save_model(
            file="./" + error_time + ".h5"
            if not grid
            else "./grid_interrupted_models.h5",
            name=model_name if model_name is not None else "interrupted",
            overwrite=True,
        )
        return

    # /Signal Interrupt Catch

    # Shared memory handlers
    @staticmethod
    def _create_shared_memory_allocation(shape: tuple):
        nbytes = int(np.prod(shape) * np.dtype(np.float64).itemsize)
        shm = shared_memory.SharedMemory(create=True, size=nbytes)
        array = np.ndarray(shape, dtype=np.float64, buffer=shm.buf, order="F")
        array.fill(0.0)
        spec = {"name": shm.name, "shape": shape}
        return shm, spec, array

    def _reform_array_in_worker(self, shm_desc):
        object.__setattr__(self, "_shm_handles", {})
        for k, v in shm_desc.items():
            shm = shared_memory.SharedMemory(name=v["name"], create=False)
            object.__setattr__(
                self,
                k,
                np.ndarray(
                    shape=v["shape"], dtype=np.float64, buffer=shm.buf, order="F"
                ),
            )
            self._shm_handles[k] = shm
            del shm
        return

    def _worker_close_memory(self):
        for k, v in self._shm_desc.items():
            try:
                self._shm_handles[k].close()
            except Exception:
                pass
            finally:
                del self._shm_handles[k]
        return

    def _coordinator_unlink_memory(self):
        if bool(self._shm_desc):
            for k, v in self._shm_desc.items():
                try:
                    self.__setattr__(k, self.__getattr__(k).copy())
                    self._shm_handles[k].close()
                    self._shm_handles[k].unlink()
                except Exception:
                    print(f"Warning, unable to close and unlike {k}")
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

    def __init__(
        self,
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        starting_chemistry: np.ndarray = None,
        previous_model: AbstractModel = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["local", "managed", "external"] = "local",
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
            run_type,
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
            print(
                f"using "
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
        _, _, _, _, out_species_abundances_array, _, success_flag = wrap.cloud(
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
            heatarray=self.heat_array,
            abundancestart=self.starting_chemistry_array
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
            else None,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(
                out_species_abundances_array[:out_species_length]
            )
        return {
            "success_flag": success_flag,
            "out_species_abundances_array": out_species_abundances_array,
        }

    def _create_init_dict(self):
        return {
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug,
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

    def __init__(
        self,
        collapse: str = "BE1.1",
        physics_output: str = None,
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        starting_chemistry: np.ndarray = None,
        previous_model: AbstractModel = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["local", "managed", "external"] = "local",
    ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        if collapse not in ["filament", "ambipolar"] and physics_output is None:
            warnings.warn(
                "`physics_output` was None but `collapse` was `filament` or `ambipolar`. No output file will be created.",
                UserWarning,
            )
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file,
            run_type,
        )
        if read_file is None:
            collapse_dict = {"BE1.1": 1, "BE4": 2, "filament": 3, "ambipolar": 4}
            try:
                self.collapse = collapse_dict[collapse]
            except KeyError:
                raise ValueError(f"collapse must be in {collapse_dict.keys()}")
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
        _, _, _, _, out_species_abundances_array, _, success_flag = wrap.collapse(
            collapsein=self.collapse,
            collapsefilein=self.physics_output,
            writeout=self.write_physics,
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
            heatarray=self.heat_array,
            abundancestart=self.starting_chemistry_array
            if self.starting_chemistry_array is not None
            else None,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(
                out_species_abundances_array[:out_species_length]
            )
        return {
            "success_flag": success_flag,
            "out_species_abundances_array": out_species_abundances_array,
        }

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
            "_debug": self._debug,
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

    def __init__(
        self,
        temp_indx: int = 1,
        max_temperature: float = 300.0,
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        starting_chemistry: np.ndarray = None,
        previous_model: AbstractModel = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["local", "managed", "external"] = "local",
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
            run_type,
        )
        if read_file is None:
            if temp_indx is None or max_temperature is None:
                raise (
                    "temp_indx and max_temperature must be specified if not reading from file."
                )
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
        _, _, _, _, out_species_abundances_array, _, success_flag = wrap.hot_core(
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
            heatarray=self.heat_array,
            abundancestart=self.starting_chemistry_array
            if self.starting_chemistry_array is not None
            else None,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(
                out_species_abundances_array[:out_species_length]
            )
        return {
            "success_flag": success_flag,
            "out_species_abundances_array": out_species_abundances_array,
        }

    def _create_init_dict(self):
        return {
            "temp_indx": self.temp_indx,
            "max_temperature": self.max_temperature,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug,
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

    def __init__(
        self,
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
        run_type: Literal["local", "managed", "external"] = "local",
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
            run_type,
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
        _, _, _, _, out_species_abundances_array, dissipation_time, _, success_flag = (
            wrap.cshock(
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
        )
        if success_flag < 0:
            dissipation_time = None
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(
                out_species_abundances_array[:out_species_length]
            )
        return {
            "success_flag": success_flag,
            "dissipation_time": dissipation_time,
            "out_species_abundances_array": out_species_abundances_array,
        }

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
            "_debug": self._debug,
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

    def __init__(
        self,
        shock_vel: float = 10.0,
        param_dict: dict = None,
        out_species: list = ["H", "N", "C", "O"],
        starting_chemistry: np.ndarray = None,
        previous_model: AbstractModel = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["local", "managed", "external"] = "local",
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
            run_type,
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
        _, _, _, _, out_species_abundances_array, _, success_flag = wrap.jshock(
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
            heatarray=self.heat_array,
            abundancestart=self.starting_chemistry_array,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(
                out_species_abundances_array[:out_species_length]
            )
        return {
            "success_flag": success_flag,
            "out_species_abundances_array": out_species_abundances_array,
        }

    def _create_init_dict(self):
        return {
            "shock_vel": self.shock_vel,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug,
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

    def __init__(
        self,
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
        run_type: Literal["local", "managed", "external"] = "local",
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
            run_type,
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
                    assert len(array) == len(
                        time_array
                    ), "All arrays must be the same length"
                    # Ensure Fortran memory
                    array = np.asfortranarray(array, dtype=np.float64)
                    self.postprocess_arrays[key] = array
            self.time_array = time_array
            self.coldens_H_array = coldens_H_array
            self.usecoldens = self.coldens_H_array is not None
            if not self.give_start_abund:
                self.starting_chemistry_array = np.zeros(
                    shape=(self.gridPoints, n_species),
                    dtype=np.float64,
                    order="F",
                )
            if self.run_type != "external":
                self.run()
        elif time_array is None and read_file is None:
            raise ValueError(
                f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given."
            )
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, _, out_species_abundances_array, _, success_flag = wrap.postprocess(
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
            heatarray=self.heat_array,
            abundancestart=self.starting_chemistry_array,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(
                out_species_abundances_array[:out_species_length]
            )
        return {
            "success_flag": success_flag,
            "out_species_abundances_array": out_species_abundances_array,
        }

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

    def __init__(
        self,
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
        run_type: Literal["local", "managed", "external"] = "local",
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
            run_type,
        )
        if read_file is None and time_array is not None:
            self.time_array = time_array
            self.postprocess_arrays = dict(
                timegrid=time_array,
                densgrid=density_array,
                gastempgrid=gas_temperature_array,
                dusttempgrid=dust_temperature_array,
                radfieldgrid=radfield_array,
                zetagrid=zeta_array,
            )
            for key, array in self.postprocess_arrays.items():
                if array is not None:
                    # Convert single values into arrays that can be used
                    if isinstance(array, float):
                        array = np.ones(shape=time_array.shape) * array
                    # Assure lengths are correct
                    assert len(array) == len(
                        time_array
                    ), "All arrays must be the same length"
                    # Ensure Fortran memory
                    array = np.asfortranarray(array, dtype=np.float64)
                    self.postprocess_arrays[key] = array
            if not self.give_start_abund:
                self.starting_chemistry_array = np.zeros(
                    shape=(self.gridPoints, n_species),
                    dtype=np.float64,
                    order="F",
                )
            if self.run_type != "external":
                self.run()
        elif time_array is None and read_file is None:
            raise ValueError(
                f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given."
            )
        return

    def run_fortran(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.run(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
        _, _, _, _, out_species_abundances_array, _, success_flag = wrap.postprocess(
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
            heatarray=self.heat_array,
            abundancestart=self.starting_chemistry_array,
        )
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(
                out_species_abundances_array[:out_species_length]
            )
        return {
            "success_flag": success_flag,
            "out_species_abundances_array": out_species_abundances_array,
        }

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
    """The SequentialModel class allows for multiple models to be run back to back.

    By defining a specific dictionary to hold the information of each model class to run in sequence, SewuentialModel allows
    for the automatic running of multiple models as well as matching some physical parameters from one model to the next.

    Args:
        sequenced_model_parameters (Dict): The dictionary to pass to SequentialModel takes the format of
            {"<First Model Class>":{"param_dict":{<parameters>}, <other arguments>}, "<Second Model Class>:{"param_dict":{<parameters>}, <other arguments>}, ...}
        parameters_to_match (List, optional): The list provided to this argument decides which parameters should be matched from a previous model
            to the next model in the sequence. Currently, supports ["finalDens", "finalTemp"].
    """

    def __init__(
        self,
        sequenced_model_parameters: Dict,
        parameters_to_match: List = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        for model in list(sequenced_model_parameters.keys()):
            assert model != SequentialModel
        self.models = []
        self.sequenced_model_parameters = sequenced_model_parameters
        self.parameters_to_match = parameters_to_match
        self.run_type = run_type
        self.model_count = 0
        self._pickle_dict = {}
        if self.run_type == "managed":
            self.run()

    def run(self):
        previous_model = None
        for model_type, model_dict in self.sequenced_model_parameters.items():
            model_dict["param_dict"] = {
                k.lower(): v for k, v in model_dict["param_dict"].items()
            }
            if self.model_count > 0:
                model_dict["param_dict"] = {
                    **{k.lower(): v for k, v in previous_model._param_dict.items()},
                    **model_dict["param_dict"],
                }
                if self.parameters_to_match is not None:
                    for parameter in self.parameters_to_match:
                        if parameter == "finalDens":
                            model_dict["param_dict"]["initialdens"] = (
                                previous_model.physics_array[-1, 0, 1]
                            )
                            continue
                        elif parameter == "finalTemp":
                            model_dict["param_dict"]["initialtemp"] = (
                                previous_model.physics_arrayy[-1, 0, 2]
                            )
                        else:
                            print(
                                f"Parameter '{parameter}' has not been implemented for parameter matching"
                            )
                tmp_model = REGISTRY[model_type](
                    **model_dict, run_type=self.run_type, previous_model=previous_model
                )
                if self.run_type == "external":
                    tmp_model.run()
                self.models += [
                    {
                        "Model_Type": model_type,
                        "Model_Order": self.model_count,
                        "Model": tmp_model,
                    }
                ]
            else:
                tmp_model = REGISTRY[model_type](
                    **model_dict, run_type=self.run_type, previous_model=previous_model
                )
                if self.run_type == "external":
                    tmp_model.run()
                self.models += [
                    {
                        "Model_Type": model_type,
                        "Model_Order": self.model_count,
                        "Model": tmp_model,
                    }
                ]
                self.models[self.model_count]["Successful"] = (
                    True
                    if self.models[self.model_count]["Model"].success_flag == 0
                    else False
                )
                previous_model = self.models[self.model_count]["Model"]
            self.model_count += 1
        return

    def save_model(
        self,
        file: str,
        name: str = "default",
        engine: str = "h5netcdf",
        overwrite: bool = False,
    ):
        for model in self.models:
            model["Model"].save_model(
                file=file,
                name=f'{name}_{model["Model_Type"]}_{model["Model_Order"]}',
                engine=engine,
                overwrite=overwrite,
            )
        return

    def check_conservation(
        self, element_list: list = ["H", "N", "C", "O"], percent: bool = True
    ):
        """Utility method to check conservation of the chemical abundances

        Args:
            element_list (list, optional): List of elements to check conservation for. Defaults to self.out_species_lists.
            percent (bool, optional): Flag on if percentage values should be used. Defaults to True.
        """
        for model in self.models:
            conserve_dicts = []
            if model["Model"]._param_dict["points"] > 1:
                for i in range(model["Model"]._param_dict["points"]):
                    conserve_dicts += [
                        check_element_conservation(
                            model["Model"].get_dataframes(i), element_list, percent
                        )
                    ]
            else:
                conserve_dicts += [
                    check_element_conservation(
                        model["Model"].get_dataframes(0), element_list, percent
                    )
                ]
            conserved = True
            for i in conserve_dicts:
                conserved = (
                    True if all([float(x[:1]) < 1 for x in i.values()]) else False
                )
            model["elements_conserved"] = conserved

    def pickle(self):
        if not bool(self._pickle_dict):
            for model in self.models:
                model["Model"] = model["Model"].pickle()
                self._pickle_dict[f'{model["Model_Type"]}_{model["Model_Order"]}'] = (
                    model["Model"]._pickle_dict.copy()
                )
        return

    def un_pickle(self):
        if bool(self._pickle_dict):
            for model in self.models:
                model["Model"]._pickle_dict = self._pickle_dict[
                    f'{model["Model_Type"]}_{model["Model_Order"]}'
                ]
                model["Model"] = model["Model"].un_pickle()
        return

    def _coordinator_unlink_memory(self):
        for model in self.models:
            model["Model"]._coordinator_unlink_memory()


def _run_grid_model(model_id, model_type, pending_model):
    """
    Internal function to run a single model. This is used by the GridModels class
    """
    cls = REGISTRY.get(model_type)
    model_obj = cls(**pending_model, run_type="external")
    model_obj.run()
    model_obj._coordinator_unlink_memory()
    model_obj.pickle()
    return (model_id, model_obj)


class GridModels:
    """GridModels, like SequentialModel is not an actual uclchem model, instead it allows running multiple models on a grid of parameter space.

    Args
        model_type (str of model class to run):
        full_parameters (Dict): The dictionary passed to GridModels should nest into it, the param_dict argument that would
            be passed to any other model, with the addition of extra keys for the none param_dict variables of a model.
            Any variables that are turned into lists or arrays, will automatically be assumed to be used for the gridding.
        max_workers (int, optional): Maximum number of workers to use in parallel for the grid run. Defaults to 8.
        grid_file (str, optional): Name and path of the output file to which the models should be saved.
            Defaults to "./default_grid_out.h5".
        model_name_prefix (str, optional): Name prefix convention to use. The fifth model in the grid would have the name
            "<model_name_prefix>5>" assigned to it. Defaults to "" which would make the fifth model have the name "5".
        delay_run (bool, optional): Defaults to False.
    """

    def __init__(
        self,
        model_type: AnyStr,
        full_parameters: Dict,
        max_workers: int = 8,
        grid_file: str = "./default_grid_out.h5",
        model_name_prefix: str = "",
        delay_run: bool = False,
    ):
        assert model_type in REGISTRY
        self.model_type = model_type
        self.full_parameters = full_parameters
        self.max_workers = (
            max_workers - 1
            if (max_workers < int(os.cpu_count()) and int(os.cpu_count()) > 32)
            else int(os.cpu_count()) - 1
            if int(os.cpu_count()) > 32
            else int(os.cpu_count() / 2) - 1
        )
        self.grid_file = grid_file if ".h5" in grid_file else grid_file + ".h5"
        # TODO: Implement model appending to grid file
        # TODO: Implement option to append or overwrite grid file.
        # Initial placeholder statement to remove pre-existing grid files
        if os.path.isfile(self.grid_file):
            os.remove(self.grid_file)
        #
        self.model_name_prefix = model_name_prefix
        self._orig_sigint = signal.getsignal(signal.SIGINT)
        self.parameters_to_grid = {}

        if self.model_type == "SequentialModel":
            for model_type, model_full_params in self.full_parameters.items():
                if isinstance(model_full_params, dict):
                    for k, v in model_full_params.items():
                        if k == "param_dict":
                            for k_p, v_p in v.items():
                                self._grid_def(k_p, v_p, model_type)
                        else:
                            self._grid_def(k, v, model_type)
            grids = np.meshgrid(*self.parameters_to_grid.values(), indexing="xy")
        else:
            for k, v in self.full_parameters.items():
                if k == "param_dict":
                    for k_p, v_p in v.items():
                        self._grid_def(k_p, v_p)
                else:
                    self._grid_def(k, v)
            grids = np.meshgrid(*self.parameters_to_grid.values(), indexing="xy")

        self.flat_grids = np.reshape(
            grids,
            (
                len(self.parameters_to_grid),
                int(np.prod(np.shape(grids)) / len(self.parameters_to_grid)),
            ),
        )

        self.model_id_dict = {}
        self.models = []
        self.physics_values = None
        self.chemical_abun_values = None
        if not delay_run:
            self.run()
        return

    def _grid_def(self, key, value, model_type=None):
        if model_type is None:
            model_type = ""
        if isinstance(value, list):
            self.parameters_to_grid[model_type + key] = value
            self.parameters_to_grid[model_type + key] = np.array(value, dtype=object)
        elif isinstance(value, (np.ndarray, np.generic)):
            self.parameters_to_grid[model_type + key] = value.astype(dtype=object)

    def run(self):
        signal.signal(signal.SIGINT, self._handler)
        pending = self.grid_iter(
            self.full_parameters,
            list(self.parameters_to_grid.keys()),
            self.flat_grids,
            self.model_type,
        )
        with mp.Pool(self.max_workers) as pool:
            active = 0
            submitted = 0
            completed = 0

            def submit_next() -> bool:
                nonlocal active, submitted
                try:
                    pending_model = next(pending)
                except StopIteration:
                    return False

                active += 1
                submitted += 1
                model_id = pending_model.pop("id")
                try:
                    pool.apply_async(
                        _run_grid_model,
                        args=(model_id, self.model_type, pending_model),
                        callback=lambda result: on_result(result),
                        error_callback=lambda exc, _model_id=model_id: on_error(
                            exc, model_id
                        ),
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
                    save_name = f"{self.model_name_prefix}{model_id}"
                    model_object.un_pickle()
                    model_object.save_model(
                        file=self.grid_file,
                        name=save_name,
                        engine="h5netcdf",
                        overwrite=True,
                    )
                    self.model_id_dict[model_id] = save_name
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
                self.model_id_dict[_model_id]._coordinator_unlink_memory()
                print(f"error: {_exc}; for model: {_model_id}")
                if not submit_next() and active == 0:
                    print("No more models to submit")

            # Submit initial batch
            for _ in range(np.shape(self.flat_grids)[1]):
                if not submit_next():
                    break

            pool.close()
            pool.join()
        signal.signal(signal.SIGINT, self._orig_sigint)
        self.models = [
            {"Model": v}
            for _, v in sorted(self.model_id_dict.items(), key=lambda item: item[1])
        ]
        self._load_params()

    def load_phys(self, engine: str = "h5netcdf"):
        if self.model_type == "SequentialModel":
            warnings.warn("Sequantial Model physics loading not implemented")
            return
        for model in range(len(self.models)):
            loaded_data = self._load_model_data(model=self.models[model]["Model"])
            if self.physics_values is None:
                self.physics_values = json.loads(loaded_data["attributes_dict"].item())[
                    "physics_values"
                ]
            loaded_data = loaded_data.assign_coords(
                {"physics_values": self.physics_values}
            )
            self.models[model]["physics_array"] = loaded_data["physics_array"]
        return

    def load_chem(
        self, out_specie_list: list = ["H", "N", "C", "O"], engine: str = "h5netcdf"
    ):
        if self.model_type == "SequentialModel":
            warnings.warn("Sequantial Model chemistry loading not implemented")
            return
        for model in range(len(self.models)):
            loaded_data = self._load_model_data(model=self.models[model]["Model"])
            if self.chemical_abun_values is None:
                self.chemical_abun_values = json.loads(
                    loaded_data["attributes_dict"].item()
                )["chemical_abun_values"]
            loaded_data = loaded_data.assign_coords(
                {"chemical_abun_values": self.chemical_abun_values}
            )
            self.models[model]["out_species_abundances_array"] = loaded_data[
                "chemical_abun_array"
            ].sel(chemical_abun_values=out_specie_list)
        return

    def _load_params(self, engine: str = "h5netcdf"):
        if self.model_type == "SequentialModel":
            for model in range(len(self.models)):
                model_number = 0
                for mt_k, mt_v in self.full_parameters.items():
                    if isinstance(mt_v, dict):
                        tmp_model = self._load_model_data(
                            model=f'{self.models[model]["Model"]}_{mt_k}_{model_number}'
                        )

                        self.models[model][f"{mt_k}_{model_number}"] = {
                            **{
                                k.replace(mt_k, ""): tmp_model._param_dict[
                                    k.replace(mt_k, "").lower()
                                ]
                                for k in list(self.parameters_to_grid.keys())
                                if mt_k in k
                                and k.replace(mt_k, "").lower() in tmp_model._param_dict
                            },
                            **{
                                k.replace(mt_k, ""): tmp_model.__getattr__(
                                    k.replace(mt_k, "")
                                )
                                for k in list(self.parameters_to_grid.keys())
                                if mt_k in k
                                and k.replace(mt_k, "").lower()
                                in tmp_model._data.keys()
                            },
                        }
                        self.models[model][f"{mt_k}_{model_number}"]["Successful"] = (
                            True
                            if tmp_model.success_flag == 0
                            else tmp_model.success_flag
                        )
                        model_number += 1
        else:
            for model in range(len(self.models)):
                loaded_data = self._load_model_data(model=self.models[model]["Model"])
                loaded_dict = json.loads(loaded_data["_param_dict"].item())
                self.models[model] = {
                    **self.models[model],
                    **{
                        k: loaded_dict[k.lower()]
                        for k in list(self.parameters_to_grid.keys())
                    },
                }
                self.models[model]["Successful"] = (
                    True
                    if json.loads(loaded_data["attributes_dict"].item())["success_flag"]
                    == 0
                    else json.loads(loaded_data["attributes_dict"].item())[
                        "success_flag"
                    ]
                )

    def _load_model_data(self, model: str, engine: str = "h5netcdf"):
        tmp_model = load_model(file=self.grid_file, name=model, engine=engine)
        return tmp_model

    def check_conservation(
        self, element_list: list = ["H", "N", "C", "O"], percent: bool = True
    ):
        """Utility method to check conservation of the chemical abundances

        Args:
            element_list (list, optional): List of elements to check conservation for. Defaults to self.out_species_lists.
            percent (bool, optional): Flag on if percentage values should be used. Defaults to True.
        """
        for model in range(len(self.models)):
            tmp_model = load_model(
                file=self.grid_file, name=self.models[model]["Model"]
            )
            conserve_dicts = []
            if tmp_model._param_dict["points"] > 1:
                for i in range(tmp_model._param_dict["points"]):
                    conserve_dicts += [
                        check_element_conservation(
                            tmp_model.get_dataframes(i), element_list, percent
                        )
                    ]
            else:
                conserve_dicts += [
                    check_element_conservation(
                        tmp_model.get_dataframes(0), element_list, percent
                    )
                ]
            conserved = True
            for i in conserve_dicts:
                conserved = (
                    True if all([float(x[:1]) < 1 for x in i.values()]) else False
                )
            self.models[model]["elements_conserved"] = conserved
        return

    def _handler(self, signum, frame):
        try:
            self.on_interrupt()  # your “final steps”
        finally:
            # Restore default and re-raise KeyboardInterrupt to stop execution
            signal.signal(signal.SIGINT, self._orig_sigint)
            raise KeyboardInterrupt

    def on_interrupt(self):
        return

    @staticmethod
    def grid_iter(
        full_parameters: dict,
        param_keys: list,
        flattened_grids: np.ndarray,
        model_type: str,
    ) -> Iterator[Dict[str, Any]]:
        if model_type == "SequentialModel":
            for i in range(len(flattened_grids[0])):
                combo = ()
                for j in range(np.shape(flattened_grids)[0]):
                    combo += (flattened_grids[j][i],)
                yield_dict = {"id": i}
                run_dict = {}
                for model_type, model_full_parameters in full_parameters.items():
                    if isinstance(model_full_parameters, dict):
                        grid_param_dict = {
                            k.replace(model_type, ""): v
                            for k, v in zip(param_keys, combo)
                            if k.replace(model_type, "")
                            in model_full_parameters["param_dict"]
                        }
                        grid_dict = {
                            k.replace(model_type, ""): v
                            for k, v in zip(param_keys, combo)
                            if k.replace(model_type, "")
                            not in model_full_parameters["param_dict"]
                            and (
                                not any(mt in k for mt in list(full_parameters.keys()))
                                or model_type in k
                            )
                        }
                        run_dict[model_type] = {
                            **model_full_parameters,
                            "param_dict": {
                                **model_full_parameters["param_dict"],
                                **grid_param_dict,
                            },
                            **grid_dict,
                        }
                    else:
                        yield_dict[model_type] = model_full_parameters
                yield {**yield_dict, **{"sequenced_model_parameters": {**run_dict}}}
        else:
            for i in range(len(flattened_grids[0])):
                combo = ()
                for j in range(np.shape(flattened_grids)[0]):
                    combo += (flattened_grids[j][i],)
                grid_param_dict = {
                    k: v
                    for k, v in zip(param_keys, combo)
                    if k in full_parameters["param_dict"]
                }
                grid_dict = {
                    k: v
                    for k, v in zip(param_keys, combo)
                    if k not in full_parameters["param_dict"]
                }
                yield {
                    **full_parameters,
                    "param_dict": {**full_parameters["param_dict"], **grid_param_dict},
                    **grid_dict,
                    "id": i,
                }
