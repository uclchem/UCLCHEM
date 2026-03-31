"""UCLCHEM Model Module.

Core module for running gas-grain chemical models under different physical conditions.

This module provides both object-oriented and functional interfaces for running
UCLCHEM chemical models. The object-oriented API (model classes) is recommended
for most use cases as it provides better state management and built-in analysis tools.

**Model Classes (Object-Oriented API - Recommended):**

- :class:`Cloud` - Static or freefall collapsing cloud
- :class:`Collapse` - Collapsing cloud with various prescriptions (BE, filament, ambipolar)
- :class:`PrestellarCore` - Prestellar core with heating and chemistry
- :class:`CShock` - C-type shock model
- :class:`JShock` - J-type shock model
- :class:`Postprocess` - Custom physics from user-provided arrays
- :class:`SequentialRunner` - Chain multiple physical stages
- :class:`GridRunner` - Run parameter grids in parallel

**Legacy Functions (Functional API):**

Available in :mod:`uclchem.functional` for backward compatibility.
Returns arrays/DataFrames instead of model objects.

**Quick Example:**

.. code-block:: python

    import uclchem

    # Create a collapsing cloud model
    cloud = uclchem.model.Cloud(
        param_dict={
            \"initialDens\": 1e2,
            \"initialTemp\": 10.0,
            \"finalTime\": 1e6,
            \"freefall\": True
        },
        out_species=[\"CO\", \"H2O\", \"CH3OH\"]
    )

    # Check for errors and plot
    cloud.check_error()
    cloud.create_abundance_plot([\"CO\", \"$CO\"])

**Model Workflow:**

1. **Initialize**: Create model object with parameters
2. **Run**: Model runs automatically on initialization (or use `read_file` to load)
3. **Analyze**: Access results via attributes (`.final_abundances`, `.chemistry_dataframe`)
4. **Plot**: Use built-in plotting methods (`.create_abundance_plot()`)
5. **Chain**: Use as input to next stage (`.previous_model` parameter)

**Common Parameters:**

All models accept these key parameters in `param_dict`:

- ``initialDens`` (float): Initial density [cm⁻³]
- ``initialTemp`` (float): Initial temperature [K]
- ``finalTime`` (float): Simulation end time [years]
- ``freefall`` (bool): Enable freefall collapse (Cloud only)
- ``outputFile`` (str): Output file path (optional with OO API)

See the user guide for complete parameter list.

**Species Naming:**

- Gas phase: ``CO``, ``H2O``, ``CH3OH``
- Ice surface: ``$CO``, ``$H2O``, ``$CH3OH``
- Ice bulk: ``@CO``, ``@H2O``, ``@CH3OH``

**See Also:**

- :mod:`uclchem.analysis` - Analyze model outputs and reaction pathways
- :mod:`uclchem.advanced` - Advanced controls for heating and network state

**Note on Thread Safety:**

Model objects are **not thread-safe** when using advanced features that modify
Fortran module state. Use multiprocessing (not threading) for parallel runs.
"""

from __future__ import annotations

import contextlib
import json
import logging

# /UCLCHEM related imports
# Multiprocessing imports
import multiprocessing as mp
import os
import signal
import typing
import warnings
from abc import ABC, abstractmethod
from collections.abc import Iterator
from datetime import datetime
from multiprocessing import pool, shared_memory
from pathlib import Path
from typing import Any, AnyStr, Literal

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uclchemwrap
import xarray as xr

# UCLCHEM related imports
from uclchemwrap import uclchemwrap as wrap

from uclchem._coolant_utils import load_coolant_level_names
from uclchem._fortran_capture import capture_fortran_output
from uclchem.analysis import (
    check_element_conservation,
)
from uclchem.constants import (
    DVODE_STAT_NAMES,
    N_DVODE_STATS,
    N_PHYSICAL_PARAMETERS,
    N_SE_STATS_PER_COOLANT,
    N_TOTAL_LEVELS,
    NCOOLANTS,
    PHYSICAL_PARAMETERS,
    SE_STAT_NAMES,
    TIMEPOINTS,
    default_elements_to_check,
    default_param_dictionary,
    n_reactions,
    n_species,
)
from uclchem.plot import create_abundance_plot, plot_species
from uclchem.utils import UCLCHEM_ROOT_DIR

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


# Model registration is intended to prevent code injection during loading time.
REGISTRY: dict[str, type[AbstractModel]] = {}


def register_model(cls: type[AbstractModel]) -> type[AbstractModel]:
    """Register a new model in the model registry.

    Args:
        cls (type[AbstractModel]): class to register.

    Returns:
        cls (type[AbstractModel]): class

    Raises:
        ValueError: If a model with the same name as cls is already in the registry.

    """
    name = getattr(cls, "MODEL_NAME", cls.__name__)
    if name in REGISTRY and REGISTRY[name] is not cls:
        raise ValueError(f"Duplicate model registration for {name}")
    REGISTRY[name] = cls
    return cls


# /Global variables determining formats of write files


# Reaction and Species name retrieval classes to reduce file read repetition.
def reaction_line_formatter(line: list[str]) -> str:
    """Format a list of strings as a reaction, while filtering out "NAN"s.

    Args:
        line (list[str]): list of species involved in the reaction

    Returns:
        str: formatted reaction for printing.

    """
    reactants = list(filter(lambda x: not str(x).lower().endswith("nan"), line[0:3]))
    products = list(filter(lambda x: not str(x).lower().endswith("nan"), line[3:7]))
    return " + ".join(reactants) + " -> " + " + ".join(products)


class ReactionNamesStore:  # noqa: D101
    def __init__(self):
        self.reaction_names = None

    def __call__(self) -> list[str]:
        """Get a list of formatted reactions.

        Only load the reactions once, after that use the cached version

        Returns:
            list[str]: List of formatted reactions.

        """
        if self.reaction_names is None:
            reactions = pd.read_csv(UCLCHEM_ROOT_DIR / "reactions.csv")
            # format the reactions:
            self.reaction_names = [
                reaction_line_formatter(line) for idx, line in reactions.iterrows()
            ]
        return self.reaction_names


get_reaction_names = ReactionNamesStore()


class SpeciesNameStore:  # noqa: D101
    def __init__(self):
        self.species_names = None

    def __call__(self) -> list[str]:
        """Get the species names.

        Only loads the species once, after that use the cached version

        Returns:
            list[str]: List of species names
        """
        if self.species_names is None:
            species = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv")
            self.species_names = species["NAME"].tolist()
        return self.species_names


get_species_names = SpeciesNameStore()
# /Reaction and Species name retrieval classes to reduce file read repetition.


# Universal model loader
def load_model(
    *,
    file_obj: h5py.File | None = None,
    file: str | None = None,
    name: str = "default",
    debug: bool = False,
) -> AbstractModel:
    """Load a pre-existing model from a file. Bypasses `__init__`.

    Args:
        file_obj (h5py.File | None): open h5py file object.
        file (str | None): Path to a file that contains previously run and stored models.
        name (str): Name of the stored object, if none was provided `default` will have been used.
            Defaults to 'default'
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features

    Returns:
        obj (object): Loaded object that inherited from AbstractModel and has the class
            of to the model found in the loaded file.

    Raises:
        ValueError: If file_obj and file are both passed, or neither are passed.
        Exception: If the model with name `name` is not found in the file.

    """
    opened_file = False
    if (file_obj is None) == (file is None):
        raise ValueError("file_obj or file must be passed.")
    elif file_obj is None:
        file_obj = h5py.File(file, "r")
        opened_file = True

    if name not in file_obj:
        raise Exception(f"model {name} was not found in the save file that was passed.")
    model_group = file_obj[name]
    coords = {}
    if "_coords" in model_group:
        for name in model_group["_coords"]:
            coords[name] = _read_array(model_group["_coords"], name)
    data_vars = {}
    for name in model_group:
        if name == "_coords":
            continue
        data_vars[name] = _read_array(model_group, name)
    loaded_data = xr.Dataset(data_vars, coords=coords)

    if opened_file:
        file_obj.close()

    model_class = json.loads(loaded_data["attributes_dict"].item())["model_type"]
    cls = REGISTRY.get(model_class)
    if cls is None:
        raise ValueError(
            f"Unrecognized model type '{model_class}'. Not in trusted registry."
        )
    return cls.load_from_dataset(model_ds=loaded_data, debug=debug)


def _read_array(model_group: dict[str, xr.Dataset], name: str) -> xr.Variable:
    """Read an array from a model group.

    Args:
        model_group (dict[str, Dataset]): model group.
        name (str): key in model_group

    Returns:
        xr.Variable: xr array

    """
    ds = model_group[name]
    data = ds[()]
    if data.dtype.kind == "S":
        data = data.astype(str)
    dims = list(ds.attrs["_dims"])
    attrs = json.loads(ds.attrs.get("_attrs", "{}"))
    return xr.Variable(dims, data, attrs=attrs)


# /Universal model loader


# Worker entry for parallel jobs
def _worker_entry(
    model_class: str,
    init_kwargs: dict,
    shm_descs: dict,
    result_queue: mp.Queue,
    advanced_snapshot: dict = None,
):
    # Restore advanced settings captured in the coordinator process.
    if advanced_snapshot is not None:
        from uclchem.advanced.worker_state import restore_snapshot

        restore_snapshot(advanced_snapshot)

    cls = REGISTRY.get(model_class)
    if cls is None:
        raise ValueError(
            f"Unrecognized model type '{model_class}'. Not in trusted registry."
        )
    model = cls.worker_build(init_kwargs=init_kwargs, shm_desc=shm_descs)
    with capture_fortran_output(label="last_model", log_file="./last_model_fortran.log"):
        output = model.run_fortran()
    result_queue.put(output)
    return


# /Worker entry for parallel jobs


# Short compatibility helper for legacy parameter `endAtFinalDensity`
def _convert_legacy_stopping_param(param_dict: dict[str, Any]) -> dict:
    """Minimal conversion of legacy `endAtFinalDensity` to `parcelStoppingMode`.

    Args:
        param_dict (dict[str, Any]): parameter dictionary.

    Rules (short and strict):
      - If both keys are present: raise RuntimeError
      - If `endAtFinalDensity` is present and points>1: raise RuntimeError
      - If `endAtFinalDensity=True` and freefall=False: raise ValueError
      - Otherwise convert True->1, False->0 and remove the old key.

    Returns:
        dict[str, Any]: Converted dictionary

    Raises:
        RuntimeError: If `endAtFinaldensity` and `parcelStoppingMode` are both
            in `param_dict`.
        RuntimeError: If `endAtFinalDensity` is being used with a multi-point model.

    Note:
        This function assumes param_dict is already a copy and is case-normalized (lowercase keys).

    """
    if param_dict is None:
        return param_dict

    has_old = "endatfinaldensity" in param_dict
    has_new = "parcelstoppingmode" in param_dict
    if has_old and has_new:
        raise RuntimeError(
            "Cannot specify both 'endAtFinalDensity' and 'parcelStoppingMode'. Use 'parcelStoppingMode' only."
        )
    if has_old:
        points = param_dict.get("points", 1)
        if points > 1:
            raise RuntimeError(
                "endAtFinalDensity is no longer supported for multi-point models (points > 1). Use 'parcelStoppingMode' instead."
            )
        old_val = param_dict.pop("endatfinaldensity")
        new_val = 1 if old_val else 0
        param_dict["parcelstoppingmode"] = new_val

        # Check if parcelStoppingMode requires freefall validation
        # Defer full validation to AbstractModel.__init__() where model_type is known
        if new_val != 0 and not param_dict.get("freefall", False):
            param_dict["_needs_freefall_validation"] = True
    return param_dict


# TODO Add catch of ctrl+c or other aborts so that it saves model and a
# full output to files of year, month, day, time type.
class AbstractModel(ABC):
    """Base model class used for inheritance only.

    The AbstractModel class serves as an abstract class from which other model classes can be built.
    It is not intended to be used as a standalone class for running UCLCHEM.

    Args:
        param_dict (dict | None): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values if not provided.
        out_species_list (list[str] | None): List of species to focus on for outputs.
            If None, defaults to `uclchem.constants.default_elements_to_check`.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from AbstractModel,
            to use for the starting abundances of the new UCLCHEM model that will be run.
            Defaults to None.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
        run_type (Literal["managed", "external"]): Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as `model.run()`.

    """

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species_list: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: object | None = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        if out_species_list is None:
            out_species_list = default_elements_to_check
        self._data = xr.Dataset()
        self._pickle_dict = {}
        # Per-instance metadata containers (scalars and small values)
        object.__setattr__(self, "_meta", {})
        object.__setattr__(self, "_pickle_meta", {})
        # Set run_type into metadata
        self.run_type = run_type
        self.separate_worker_types = ["managed"]
        # Shared memory
        self._shm_desc = {}
        self._shm_handles = {}
        self._proc_handle = None
        # /Shared memory

        # Signal Interrupt
        self._was_interrupted = False
        self._orig_sigint = signal.getsignal(signal.SIGINT)
        # /Signal Interrupt

        self.model_type = str(self.__class__.__name__)
        self._param_dict = {}
        self.out_species_list = out_species_list
        self.out_species = ""
        self.full_array = None
        self._debug = debug
        self.success_flag = None
        # Note: specname is now accessed via get_species_names() global function
        # Note: PHYSICAL_PARAMETERS is now accessed via the global constant

        self.n_out = 0 if read_file is None else None
        self.timepoints = timepoints
        self.was_read = read_file is not None

        self._reform_inputs(param_dict, self.out_species_list)

        # Validate parcelStoppingMode usage after model_type is known
        if self._param_dict.get("_needs_freefall_validation", False):
            if self.model_type not in ["Collapse"]:  # Collapse models imply freefall
                raise ValueError(
                    "parcelStoppingMode != 0 can only be used with:\n"
                    "  - Cloud models with freefall=True\n"
                    "  - Collapse models (freefall is implied)\n"
                    f"Current model_type={self.model_type}, freefall={self._param_dict.get('freefall', False)}\n"
                    "Please either set freefall=True or set parcelStoppingMode=0"
                )
            self._param_dict.pop("_needs_freefall_validation", None)  # Clean up flag

        # If we were given a previously-written output file, populate the model
        # arrays and metadata from that file now so later initialization can rely on them.
        if read_file is not None:
            self.legacy_read_output_file(read_file)
        if "points" not in self._param_dict:
            self._param_dict["points"] = 1
        # Expose grid points as attribute for legacy code expecting `gridPoints`
        object.__setattr__(self, "gridPoints", self._param_dict["points"])

        self.outputFile = (
            self._param_dict.pop("outputfile")
            if "outputfile" in self._param_dict
            else None
        )
        self.abundSaveFile = (
            self._param_dict.pop("abundsavefile")
            if "abundsavefile" in self._param_dict
            else None
        )
        self.abundLoadFile = (
            self._param_dict.pop("abundloadfile")
            if "abundloadfile" in self._param_dict
            else None
        )

        self.starting_chemistry_array = None
        if previous_model is None and self.abundLoadFile is None:
            self._create_starting_array(starting_chemistry)
        elif self.abundLoadFile is not None:
            self.legacy_read_starting_chemistry()
        elif previous_model is not None and previous_model.has_attr(
            "next_starting_chemistry_array"
        ):
            # All models must use the same hard-coded PHYSICAL_PARAMETERS from constants.py
            # If previous_model was loaded from a legacy file with different parameters,
            # that's a compatibility issue that should have been caught during load.
            self._create_starting_array(previous_model.next_starting_chemistry_array)

        self.give_start_abund = self.starting_chemistry_array is not None
        assert not np.all(self.starting_chemistry_array == 0.0), (
            "Detected all zeros starting chemistry array."
        )

        # Only initialize next_starting_chemistry_array if we didn't load it from a file
        # (legacy_read_output_file sets it from the last timestep)
        if read_file is None:
            self.next_starting_chemistry_array = None

        # Only create new arrays if we didn't load them from a file
        if read_file is None:
            self.physics_array = None
            self.chemical_abun_array = None
            self._create_fortran_array()
            self.rate_constants_array = None
            self._create_rate_constants_array()
            self.heat_array = None
            self._create_heating_array()
            self.stats_array = None
            self._create_stats_array()
            self.level_populations_array = None
            self._create_level_populations_array()
            self.se_stats_array = None
            self._create_se_stats_array()
            self.out_species_abundances_array = None
        else:
            # When loading from file, arrays are already populated; just initialize
            # the arrays that weren't loaded
            self.rate_constants_array = None
            self.heat_array = None
            self.stats_array = None
            self.level_populations_array = None
            self.se_stats_array = None
            self.out_species_abundances_array = None
        return

    def __del__(self):
        if hasattr(self, "_shm_desc"):
            self._coordinator_unlink_memory()

    # Separate class building method(s)
    @classmethod
    def load_from_dataset(
        cls, model_ds: xr.Dataset, debug: bool = False
    ) -> AbstractModel:
        """Load an abstract model from an xr Dataset.

        Args:
            model_ds (xr.Dataset): Dataset to load
            debug (bool): Flag to set

        Returns:
            obj (AbstractModel): instantiated model
        """
        obj = cls.__new__(cls)
        obj._param_dict = json.loads(model_ds["_param_dict"].item())
        model_ds.__delitem__("_param_dict")
        obj._data = xr.Dataset()
        obj._data = model_ds.copy()
        model_ds.close()
        temp_attribute_dict = json.loads(obj._data["attributes_dict"].item())
        # Restore these values into the metadata dict rather than dataset variables
        object.__setattr__(obj, "_meta", temp_attribute_dict)
        obj._data.__delitem__("attributes_dict")
        obj.debug = debug
        obj._coord_assign()
        return obj

    @classmethod
    def worker_build(cls, init_kwargs, shm_desc):  # noqa: ANN001, D102
        obj = cls.__new__(cls)
        for k, v in init_kwargs.items():
            object.__setattr__(obj, k, v)
        obj._reform_array_in_worker(shm_desc)
        return obj

    @classmethod
    def from_file(
        cls,
        file: str,
        name: str = "default",
        debug: bool = False,
    ):
        """Load a model from a file.

        This is a convenience class method that wraps the module-level load_model function.

        Args:
            file (str): Path to a file that contains previously run and stored models.
            name (str): Name of the stored object. Defaults to 'default'.
            debug (bool): Flag if extra debug information should be printed to the terminal.
                Defaults to False. #TODO Add debug features

        Returns:
            Model object loaded from the file.

        """
        return load_model(file=file, name=name, debug=debug)

    # /Separate class building methods

    # Class utility methods
    def __getattr__(self, key: str) -> Any:
        # Internal attributes behave normally
        if key.startswith("_") and key != "_data":
            return super().__getattribute__(key)

        # First check instance metadata
        try:
            meta = super().__getattribute__("_meta")
        except Exception:
            meta = {}
        if key in meta:
            return meta[key]

        # Next check the xarray Dataset
        if key in super().__getattribute__("_data"):
            values = super().__getattribute__("_data")[key].values
            if np.shape(values) != ():
                if isinstance(values, tuple):
                    return values[1]
                else:
                    return values
            else:
                return self._data[key].item()

        raise AttributeError(
            f'{self.__class__.__name__} has no attribute of name: "{key}".'
        )

    def __setattr__(self, key: str, value: Any) -> None:
        # Underscored attributes are real attributes
        if key.startswith("_"):
            super().__setattr__(key, value)
            return

        # Determine ndim safely
        try:
            ndim = np.ndim(value)
        except Exception:
            ndim = None

        # If this looks like an array variable (name contains '_array' and ndim >= 1), store in _data
        if ndim is not None and "_array" in key and ndim >= 1:
            # Ensure value is a numpy array (convert lists, tuples, etc.)
            if not isinstance(value, np.ndarray):
                value = np.asarray(value)
                ndim = value.ndim  # Update ndim after conversion

            # Remove any existing conflicting metadata entry
            try:
                meta = super().__getattribute__("_meta")
                if key in meta:
                    del meta[key]
            except Exception:
                pass

            # Remove any existing dataset var of same name before inserting
            if key in self._data:
                with contextlib.suppress(Exception):
                    self._data = self._data.drop_vars(key)

            # Choose a time dimension name that avoids conflicts with existing dims
            time_dim = "time_step"
            existing_time = None
            try:
                # Use .sizes (mapping of dim name -> length) instead of .dims to avoid FutureWarning
                existing_time = self._data.sizes.get("time_step", None)
            except Exception:
                existing_time = None

            if existing_time is not None and existing_time != value.shape[0]:
                # create a per-variable time dim to avoid conflicts
                base_time_dim = f"time_step_{key}"
                time_dim = base_time_dim
                i = 1
                while time_dim in self._data.sizes:
                    time_dim = f"{base_time_dim}_{i}"
                    i += 1

            if ndim == 3:
                # Determine the 'values' dimension name for the variable
                values_dim = key.replace("array", "values")
                # If an existing coordinate with this name has a conflicting
                # length, create a per-variable values dimension to avoid
                # xarray merge conflicts (similar to per-variable time dims).
                try:
                    existing_values_len = self._data.sizes.get(values_dim, None)
                except Exception:
                    existing_values_len = None

                if (
                    existing_values_len is not None
                    and existing_values_len != value.shape[2]
                ):
                    base_values_dim = f"{values_dim}_{key}"
                    new_values_dim = base_values_dim
                    j = 1
                    while new_values_dim in self._data.sizes:
                        new_values_dim = f"{base_values_dim}_{j}"
                        j += 1
                    self._data[key] = ([time_dim, "point", new_values_dim], value)
                    # add numeric coords for the new dimension
                    self._data = self._data.assign_coords(
                        {new_values_dim: np.arange(value.shape[2])}
                    )
                else:
                    self._data[key] = ([time_dim, "point", values_dim], value)

                # Add coordinate for custom time dim if we created one
                if time_dim != "time_step":
                    self._data = self._data.assign_coords(
                        {time_dim: np.arange(value.shape[0])}
                    )
            elif ndim == 2:
                self._data[key] = (["point", key], value)
            elif ndim == 1:
                # For 1D arrays, use the array name as the dimension
                dim_name = key.replace("_array", "")
                self._data[key] = ([dim_name], value)
            else:
                self._data[key] = value
            return

        # Otherwise, store in metadata
        try:
            meta = super().__getattribute__("_meta")
        except Exception:
            object.__setattr__(self, "_meta", {})
            meta = self._meta

        # If a dataset variable exists with this name, drop it to avoid ambiguity
        if key in self._data:
            with contextlib.suppress(Exception):
                self._data = self._data.drop_vars(key)

        meta[key] = value
        return

    def has_attr(self, key: str) -> bool:
        """Method to check if the object has an attribute stored in self._meta or self._data.

        Args:
            key (str): name of attribute

        Returns:
            bool: whether the object has the attribute.

        """
        try:
            meta = super().__getattribute__("_meta")
        except Exception:
            meta = {}
        return (key in meta) or (key in self._data)

    # /Class utility method

    # UCLCHEM utility and analysis wrappers
    def check_conservation(
        self, element_list: list[str] | None = None, percent: bool = True
    ) -> None:
        """Utility method to check conservation of the chemical abundances.

        Args:
            element_list (list): List of elements to check conservation for.
                Defaults to `self.out_species_lists`.
            percent (bool): Flag on if percentage values should be used. Defaults to True.

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
                check_element_conservation(self.get_dataframes(0), element_list, percent)
            )

    def check_error(self, only_error: bool = False, raise_on_error: bool = True) -> None:
        """Checks the model error status and raises RuntimeError on failure.

        Args:
            only_error (bool): If True, only act when there was an error (skip success message).
            raise_on_error (bool): If True (default), raises RuntimeError on failure. If False, prints.

        """
        if self.success_flag is not None and self.success_flag < 0:
            from uclchem.utils import check_error as _check_error

            if raise_on_error:
                _check_error(self.success_flag, raise_on_error=True)
            else:
                msg = _check_error(self.success_flag, raise_on_error=False)
                print(msg)
        elif self.success_flag == 0 and not only_error:
            print("Model ran successfully.")
        elif self.success_flag is None:
            print("Model has not been run.")

    def create_abundance_plot(
        self,
        species: list[str] | None = None,
        figsize: tuple[2] = (16, 9),
        point: int = 0,
        plot_file: str | Path | None = None,
    ) -> tuple[plt.Figure, plt.Axes]:
        """`uclchem.plot.create_abundance_plot` wrapper method.

        Args:
            species (list[str] | None): List of species to plot. If None, uses self.out_species_list.
                Default = None.
            figsize (tuple[2]): The figure size to use for matplotlib Defaults to (16, 9).
            point (int): Integer referring to which point of the UCLCHEM model to use. Defaults to 0.
            plot_file (str | Path | None): if not None, save to a path. Default = None.

        Returns:
            tuple[fig, ax]: matplotlib figure and axis objects

        Raises:
            ValueError: If `point` is larger than the number of points in the model run.

        """
        if species is None:
            species = self.out_species_list

        if point > self._param_dict["points"]:
            raise ValueError("'point' must be less than number of modelled points.")
        return create_abundance_plot(
            self.get_dataframes(point), species, figsize, plot_file
        )

    def get_dataframes(
        self,
        point: int | None = None,
        joined: bool = True,
        with_rate_constants: bool = False,
        with_heating: bool = False,
        with_stats: bool = False,
        with_level_populations: bool = False,
        with_se_stats: bool = False,
    ) -> pd.DataFrame | tuple[pd.DataFrame, ...]:  # Returns joined DF or tuple of DFs
        """Converts the model physics and chemical_abun arrays from numpy to pandas arrays.

        Args:
            point (int | None): Integer referring to which point of the UCLCHEM model to return.
                If None, returns data for all points with a 'Point' column. Defaults to None.
            joined (bool): Flag on whether the returned pandas dataframe should be one, or if
                two dataframes should be
                returned. One physical, one chemical_abun dataframe. Defaults to True.
            with_rate_constants (bool): Flag on whether to include reaction rate constants
                in the dataframe, and/or as a separate dataframe depending on the value of `joined`.
                Defaults to False.
            with_heating (bool): Flag on whether to include heating/cooling rates in the dataframe,
                and/or as a separate dataframe depending on the value of `joined`. Defaults to False.
            with_stats (bool): Flag on whether to include DVODE solver statistics in the dataframe,
                and/or as a separate dataframe depending on the value of `joined`. Defaults to False.
            with_level_populations (bool): Flag on whether to include coolant level populations
                in the output. Defaults to False.
            with_se_stats (bool): Flag on whether to include SE solver statistics in the output.
                Defaults to False.

        Returns:
            return_df (pd.DataFrame): Dataframe of the joined arrays for point `point` if joined = True
            physics_df (pd.DataFrame): Dataframe of the physical parameters for point `point`
                if joined = False
            chemistry_df (pd.DataFrame): Dataframe of the chemical abundances  for point `point`
                if joined = False
            rate_constants_df (pd.DataFrame): Dataframe of the reaction rate constants for point `point`
                if joined = False and with_rate_constants = True
            heating_df (pd.DataFrame): Dataframe of the heating/cooling rates  for point `point`
                if joined = False and with_heating = True
            stats_df (pd.DataFrame): Dataframe of DVODE solver statistics for point `point`
                if joined = False and with_stats = True
            level_populations_df (pd.DataFrame): Dataframe of coolant level populations for point `point`
                if joined = False and with_level_populations = True
            se_stats_df (pd.DataFrame): Dataframe of SE solver statistics for point `point`
                if joined = False and with_se_stats = True

        """

        # Helper function to add Point column to a dataframe
        def add_point_column(df: pd.DataFrame, point_num: int) -> pd.DataFrame:
            if df is not None:
                df.insert(0, "Point", point_num)
            return df

        # Determine total number of points in model
        n_points = self._param_dict.get("points", 1)

        # Get dataframes for the requested points
        if point is None:
            # Collect dataframes for all points
            all_dfs = []
            for pt in range(n_points):
                dfs = self._get_single_point_dataframes(
                    pt,
                    with_rate_constants,
                    with_heating,
                    with_stats,
                    with_level_populations,
                    with_se_stats,
                )
                # Add Point columns to all dataframes (1-indexed)
                dfs = tuple(add_point_column(df, pt + 1) for df in dfs)
                all_dfs.append(dfs)

            # Transpose to group by dataframe type instead of by point
            # e.g., [[phys0, chem0, rates0], [phys1, chem1, rates1]] -> [[phys0, phys1], [chem0, chem1], [rates0, rates1]] # noqa: W505
            df_collections = list(zip(*all_dfs))

            # Concatenate each type vertically
            concatenated = [
                pd.concat([df for df in collection if df is not None], ignore_index=True)
                if any(df is not None for df in collection)
                else None
                for collection in df_collections
            ]
        else:
            # Single point mode
            concatenated = list(
                self._get_single_point_dataframes(
                    point,
                    with_rate_constants,
                    with_heating,
                    with_stats,
                    with_level_populations,
                    with_se_stats,
                )
            )
            # Add Point columns (1-indexed)
            concatenated = [add_point_column(df, point + 1) for df in concatenated]

        # Join horizontally if requested
        if joined:
            result_df = concatenated[0]  # Start with physics
            for df in concatenated[1:]:
                if df is not None:
                    # Drop duplicate Point column from subsequent dataframes
                    result_df = result_df.join(df.drop(columns=["Point"]))
            return result_df
        else:
            return tuple(concatenated)

    def _get_single_point_dataframes(
        self,
        point: int,
        with_rate_constants: bool,
        with_heating: bool,
        with_stats: bool,
        with_level_populations: bool,
        with_se_stats: bool,
    ) -> tuple[pd.DataFrame]:
        """Helper method to get dataframes for a single point without Point column.

        Args:
            point (int): Spatial point index (for multi-point models).
            with_rate_constants (bool): Flag on whether to include a reaction rate constant dataframe
                in the tuple.
            with_heating (bool): Flag on whether to include heating/cooling rates dataframe in the tuple.
            with_stats (bool): Flag on whether to include DVODE solver statistics dataframe in the tuple.
            with_level_populations (bool): Flag on whether to include coolant level
                populations in the tuple.
            with_se_stats (bool): Flag on whether to include SE solver statistics
                in the tuple

        Returns:
            tuple[pd.DataFrames]: a tuple of pd.DataFrame with physics_df, chemistry_df, and all
                additional information based off whether the flags were True.
        """
        # Create a physical parameter dataframe using global constants
        # Arrays are guaranteed to match these dimensions due to validation in legacy_read_output_file
        physics_df = pd.DataFrame(
            self.physics_array[:, point, :],
            index=None,
            columns=PHYSICAL_PARAMETERS,
        )
        # Create an abundances dataframe using global species names
        species_names = get_species_names()
        chemistry_df = pd.DataFrame(
            self.chemical_abun_array[:, point, :], index=None, columns=species_names
        )
        if self.rate_constants_array is not None and with_rate_constants:
            # Create a rate constants dataframe.
            rate_constants_df = pd.DataFrame(
                self.rate_constants_array[:, point, :],
                index=None,
                columns=get_reaction_names(),
            )
        else:
            rate_constants_df = None

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

        if self.stats_array is not None and with_stats:
            stats_df = pd.DataFrame(
                self.stats_array[:, point, :], index=None, columns=DVODE_STAT_NAMES
            )
        else:
            stats_df = None

        if self.level_populations_array is not None and with_level_populations:
            level_populations_df = self.get_level_populations_dataframe(point=point)
        else:
            level_populations_df = None

        if self.se_stats_array is not None and with_se_stats:
            se_stats_df = self.get_se_stats_dataframe(point=point)
        else:
            se_stats_df = None

        # Return tuple of all dataframes (some may be None)
        result = [physics_df, chemistry_df]
        if with_rate_constants:
            result.append(rate_constants_df)
        if with_heating:
            result.append(heating_df)
        if with_stats:
            result.append(stats_df)
        if with_level_populations:
            result.append(level_populations_df)
        if with_se_stats:
            result.append(se_stats_df)
        return tuple(result)

    def get_solver_stats_dataframe(self, point: int | None = None) -> pd.DataFrame | None:
        """Get all solver statistics including failed attempts.

        This method returns statistics for EVERY DVODE solver call,
        including failed attempts that were retried. This is different
        from the regular stats in get_dataframes() which only shows
        the final successful attempt per trajectory timestep.

        Args:
            point (int | None): Spatial point index (for multi-point models). If None, uses point 0.

        Returns:
            pd.DataFrame | None: DataFrame with columns from DVODE_STAT_NAMES,
                or None if stats not available.
                TRAJECTORY_INDEX column links solver attempts to trajectory timesteps.
                Rows where TRAJECTORY_INDEX=0 are filtered out (unused preallocated space).

        Example:
            >>> model = uclchem.model.Cloud(param_dict)
            >>> solver_stats = model.get_solver_stats_dataframe()
            >>> # Count failed attempts
            >>> failures = solver_stats[solver_stats['ISTATE'] < 0]
            >>> print(f"Failed attempts: {len(failures)}")

        """
        if self.stats_array is None:
            return None

        if point is None:
            point = 0

        # Extract stats for this point and filter out unused rows
        stats_data = self.stats_array[:, point, :]
        valid_mask = stats_data[:, 0] > 0  # TRAJECTORY_INDEX > 0
        stats_data = stats_data[valid_mask]

        if len(stats_data) == 0:
            return None

        df = pd.DataFrame(stats_data, columns=DVODE_STAT_NAMES)
        df["TRAJECTORY_INDEX"] = df["TRAJECTORY_INDEX"].astype(int)

        return df

    def get_failed_solver_attempts(self, point: int | None = None) -> pd.DataFrame | None:
        """Get only the failed solver attempts (ISTATE < 0).

        Returns a DataFrame containing only solver attempts that failed
        and required retry (ISTATE = -1, -2, -4, -5, etc.).

        Args:
            point (int | None): Spatial point index (for multi-point models). If None, uses point 0.

        Returns:
            pd.DataFrame | None: DataFrame of failed attempts,
                or None if no failures or stats unavailable.

        Example:
            >>> failures = model.get_failed_solver_attempts()
            >>> if failures is not None:
            >>>     print(f"Total retries needed: {len(failures)}")
            >>>     print(failures.groupby('ISTATE').size())

        """
        df = self.get_solver_stats_dataframe(point)
        if df is None:
            return None

        failed = df[df["ISTATE"] < 0]
        return failed if len(failed) > 0 else None

    def get_solver_efficiency_summary(
        self, point: int | None = None
    ) -> dict[str, int | float] | None:
        """Calculate solver efficiency metrics.

        Args:
            point (int | None): Spatial point index (for multi-point models). If None, uses point 0.

        Returns:
            dict[str, int | float] | None: dict with keys:
                - total_attempts: Total DVODE calls
                - successful_attempts: Calls that advanced the trajectory
                - failed_attempts: Calls that were retried
                - efficiency_ratio: successful / total (1.0 = no retries)
                - total_cpu_time: Sum of all CPU time
                - wasted_cpu_time: CPU time spent on failed attempts

        """
        df = self.get_solver_stats_dataframe(point)
        if df is None:
            return None

        total_attempts = len(df)
        failed_attempts = len(df[df["ISTATE"] < 0])
        successful_attempts = total_attempts - failed_attempts

        total_cpu = df["CPU_TIME"].sum()
        wasted_cpu = df[df["ISTATE"] < 0]["CPU_TIME"].sum()

        return {
            "total_attempts": total_attempts,
            "successful_attempts": successful_attempts,
            "failed_attempts": failed_attempts,
            "efficiency_ratio": successful_attempts / total_attempts
            if total_attempts > 0
            else 0.0,
            "total_cpu_time": total_cpu,
            "wasted_cpu_time": wasted_cpu,
            "wasted_fraction": wasted_cpu / total_cpu if total_cpu > 0 else 0.0,
        }

    def plot_species(
        self,
        ax: plt.Axes,
        species: list[str] | None = None,
        point: int = 0,
        legend: bool = True,
        **plot_kwargs,
    ) -> plt.Axes:
        """`uclchem.plot.plot(species)` wrapper method.

        Args:
            ax (pyplot.axis): An axis object to plot on
            df (pd.DataFrame): A dataframe created by `read_output_file`
            species (list[str]): A list of species names to be plotted.
                If species name starts with "$" instead of "#" or "@",
                plots the sum of surface and bulk abundances. If None, default to
                `self.out_species_list`.
            point (int): Grid point index. Default = 0.
            legend (bool): Whether to add a legend to the plot. Default = True.
            plot_kwargs (dict[str, Any]): keyword arguments passed to `ax.plot`.

        Returns:
            plt.Axes: Modified input axis is returned

        """
        if species is None:
            species = self.out_species_list
        return plot_species(
            ax, self.get_dataframes(point), species, legend, **plot_kwargs
        )

    # /UCLCHEM utility and analysis wrappers

    # Methods to start run of model
    def run(self) -> None:
        """Reset the Fortran arrays if the model was not read,
        allowing the arrays to be reused for new runs.

        Raises:
            RuntimeError: If the model was read.
            ValueError: If the model's run_type is invalid.

        """
        if self.was_read:
            raise RuntimeError("This model was read. It can not be run. ")
            self.physics_array = None
            self.chemical_abun_array = None
            self.rateConstantsArray = None
            self.heatArray = None
            self.statsArray = None
            self._create_fortran_array()
            self._create_rate_constants_array()
            self._create_heating_array()
            self._create_stats_array()

        def _handler(signum, frame):  # noqa: ARG001, ANN001
            try:
                self.on_interrupt()  # your “final steps”
            finally:
                # Restore default and re-raise KeyboardInterrupt to stop execution
                signal.signal(signal.SIGINT, self._orig_sigint)
                raise KeyboardInterrupt

        signal.signal(signal.SIGINT, _handler)
        if self.run_type not in self.separate_worker_types:
            output = self.run_fortran()
        elif self.run_type in self.separate_worker_types:
            from uclchem.advanced.worker_state import create_snapshot

            snapshot = create_snapshot()
            init_kwargs = self._create_init_dict()
            ctx = mp.get_context("spawn")
            result_queue = ctx.Queue()
            self._proc_handle = ctx.Process(
                target=_worker_entry,
                args=(
                    self.model_type,
                    init_kwargs,
                    self._shm_desc,
                    result_queue,
                    snapshot,
                ),
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

        signal.signal(signal.SIGINT, self._orig_sigint)

        if hasattr(self, "_shm_handles"):
            self._coordinator_unlink_memory()

        for k, v in output.items():
            self.__setattr__(k, v)

        self._array_clean()
        self.check_error(only_error=True)
        if self.outputFile is not None:
            logging.debug(f"Writing output file: {self.outputFile}")
            logging.debug(
                f"Physics array shape: {self.physics_array.shape if self.physics_array is not None else None}"
            )
            logging.debug(
                f"Chemical array shape: {self.chemical_abun_array.shape if self.chemical_abun_array is not None else None}"
            )
            try:
                self.legacy_write_full()
                logging.debug(f"Successfully wrote {self.outputFile}")
            except Exception as e:
                logging.error(f"Failed to write {self.outputFile}: {e}", exc_info=True)
                raise
        if self.abundSaveFile is not None:
            logging.debug(f"Writing abundance file: {self.abundSaveFile}")
            try:
                self.legacy_write_starting_chemistry()
                logging.debug(f"Successfully wrote {self.abundSaveFile}")
            except Exception as e:
                logging.error(f"Failed to write {self.abundSaveFile}: {e}", exc_info=True)
                raise
        return

    @abstractmethod
    def run_fortran(self) -> dict[str, int | list]:  # noqa: D102
        raise NotImplementedError

    # /Methods to start run of model

    # Model saving
    def save_model(
        self,
        *,
        file_obj: h5py.File | None = None,
        file: str | None = None,
        name: str = "default",
        overwrite: bool = False,
    ) -> None:
        """Save a model to file on disk. Multiple models can be saved into the same file
        if different names are used to store them.

        Args:
            file_obj (h5py.File | None): open file object
            file (str | None): Path to a file to store models.
            name (str): Name to use for the group of the object. Defaults to 'default'
            overwrite (bool): Boolean on whether to overwrite pre-existing models, or error out.
                Defaults to False

        Raises:
            ValueError: If file_obj and file are both passed, or neither are passed.
        """
        opened_file = False
        if file_obj is None and file is None:
            raise ValueError("file_obj or file must be passed.")
        elif file_obj is None:
            file_obj = h5py.File(file, "a")
            opened_file = True

        if name in file_obj:
            if not overwrite:
                warnings.warn(
                    f"Model with name: `{name}` already exists in save file but overwrite is set to False. Unable to save model."
                )
                return
            else:
                del file_obj[name]
        # TODO: Allow for toggling of saving float64 or float32 for the arrays
        temp_attribute_dict = {}
        with contextlib.suppress(Exception):
            temp_attribute_dict.update(super().__getattribute__("_meta"))

        # Work on a copy so save_model is non-destructive to self._data
        save_data = self._data.copy()
        # Collect remaining non-array dataset variables into attributes (same behaviour as before)
        for v in list(save_data.variables):
            if "_array" not in v and v not in ["_orig_sigint"]:
                if np.shape(save_data[v].values) != ():
                    if isinstance(save_data[v].values, tuple):
                        temp_attribute_dict[v] = save_data[v].values[1].tolist()
                        save_data = save_data.drop_vars(v)
                    else:
                        temp_attribute_dict[v] = save_data[v].values.tolist()
                        save_data = save_data.drop_vars(v)
                else:
                    temp_attribute_dict[v] = save_data[v].item()
                    save_data = save_data.drop_vars(v)
        #
        save_data["attributes_dict"] = xr.DataArray([json.dumps(temp_attribute_dict)])
        save_data["_param_dict"] = xr.DataArray([json.dumps(self._param_dict)])
        model_group = file_obj.create_group(name)
        coord_grp = model_group.create_group("_coords")
        for name, coord in save_data.coords.items():
            self._write_array(coord_grp, name, coord)
        for name, var in save_data.data_vars.items():
            self._write_array(model_group, name, var)
        if opened_file:
            file_obj.flush()
            file_obj.close()

    @staticmethod
    def _write_array(model_group: h5py.Group, name: str, xr_var: xr.DataArray) -> None:
        data = xr_var.values
        if data.dtype.kind == "U":
            data = data.astype(bytes)
        ds = model_group.create_dataset(name, data=data)
        ds.attrs["_dims"] = list(xr_var.dims)

    # /Model saving

    # Model Passing through Pickling
    def pickle(self) -> AbstractModel:
        """Pickle the model.

        Returns:
            AbstractModel

        """
        if self._data is not None and not bool(self._pickle_dict):
            for v in self._data.variables:
                if np.shape(self._data[v].values) != ():
                    if isinstance(self._data[v].values, tuple):
                        self._pickle_dict[v] = self._data[v].values[1].tolist()
                    else:
                        self._pickle_dict[v] = self._data[v].values.tolist()
                else:
                    self._pickle_dict[v] = self._data[v].item()
            # Save metadata separately for pickle roundtrip
            try:
                object.__setattr__(
                    self, "_pickle_meta", super().__getattribute__("_meta").copy()
                )
            except Exception:
                object.__setattr__(self, "_pickle_meta", {})
            self._data = None
            # Clear runtime metadata to reflect pickled state
            object.__setattr__(self, "_meta", {})
        return self

    def un_pickle(self) -> AbstractModel:
        """Un-pickle the model.

        Returns:
            AbstractModel

        """
        if self._data is None and bool(self._pickle_dict):
            self._data = xr.Dataset()
            for k, v in self._pickle_dict.items():
                if np.ndim(v) == 3 and "_array" in k:
                    # Avoid colliding with existing 'time_step' dim sizes by
                    # making a per-variable time dim if necessary
                    time_dim = "time_step"
                    existing_time = None
                    try:
                        existing_time = self._data.sizes.get("time_step", None)
                    except Exception:
                        existing_time = None
                    v_arr = np.asarray(v)
                    if existing_time is not None and existing_time != np.shape(v_arr)[0]:
                        base_time_dim = f"time_step_{k}"
                        time_dim = base_time_dim
                        i = 1
                        while time_dim in self._data.sizes:
                            time_dim = f"{base_time_dim}_{i}"
                            i += 1
                        self._data[k] = (
                            [time_dim, "point", k.replace("array", "values")],
                            v_arr,
                        )
                        self._data = self._data.assign_coords(
                            {time_dim: np.arange(np.shape(v_arr)[0])}
                        )
                    else:
                        self._data[k] = (
                            ["time_step", "point", k.replace("array", "values")],
                            v_arr,
                        )
                elif np.ndim(v) == 2 and "_array" in k:
                    self._data[k] = (["point", k], v)
                elif "_values" in k:
                    pass
                else:
                    self._data[k] = v
            self._pickle_dict = {}
            # Restore saved metadata if present
            try:
                if hasattr(self, "_pickle_meta") and self._pickle_meta:
                    object.__setattr__(self, "_meta", self._pickle_meta.copy())
            except Exception:
                pass
            finally:
                object.__setattr__(self, "_pickle_meta", {})
        self._coord_assign()
        return self

    # /Model Passing through Pickling

    # Legacy in & output support
    def legacy_read_output_file(
        self,
        read_file: str | Path,
        rate_constants_load_file: str | Path | None = None,  # noqa: ARG002
    ) -> None:
        """Perform classic output file reading.

        This reader constructs the internal :class:`xarray.Dataset` directly from the
        parsed header and data in the legacy ``.dat`` output format. The classic
        files place a `point` column between the physics and chemistry columns; we
        therefore use the location of `point` in the header to split the columns
        reliably and avoid using global constants as authoritative metadata.

        Args:
            read_file (str | Path): path to file
            rate_constants_load_file (str | Path | None): Not used

        Raises:
            ValueError: If there is any incompatibility error.

        """
        self.was_read = True
        # Read header and numeric data
        columns = np.char.strip(
            np.loadtxt(read_file, delimiter=",", max_rows=1, dtype=str, comments="%")
        )
        raw_array = np.loadtxt(read_file, delimiter=",", skiprows=1)

        # Determine where the `point` column is and how many points exist
        point_index = np.where(columns == "point")[0][0]
        self._param_dict["points"] = int(np.max(raw_array[:, point_index]))

        # Some legacy files include an additional metadata row; if more than one
        # point exists the legacy format contains that extra row which we must skip
        if self._param_dict["points"] > 1:
            array = np.loadtxt(read_file, delimiter=",", skiprows=2)
        else:
            array = raw_array

        row_count = int(np.shape(array)[0] / self._param_dict["points"])

        # Extract physics and species columns from legacy file header
        physics_cols_from_file = list(columns[:point_index].tolist())
        species_cols_from_file = list(columns[point_index + 1 :].tolist())

        # Validate compatibility with current UCLCHEM constants - these are non-negotiable
        # PHYSICAL_PARAMETERS are hard-coded and tied to the Fortran wrapper
        if physics_cols_from_file != list(PHYSICAL_PARAMETERS):
            # Special case: backward compatibility for missing 'dstep' parameter
            # Old UCLCHEM versions didn't have dstep; we can infer it if safe
            missing_params = set(PHYSICAL_PARAMETERS) - set(physics_cols_from_file)
            extra_params = set(physics_cols_from_file) - set(PHYSICAL_PARAMETERS)

            if missing_params == {"dstep"} and not extra_params:
                # Only dstep is missing - check if we can safely infer it
                # If there are no duplicate timesteps, we can assume dstep=1
                time_column_index = physics_cols_from_file.index("Time")
                time_values = array[:, time_column_index]

                if len(time_values) == len(np.unique(time_values)):
                    # No duplicate timesteps - safe to assume dstep=1
                    warnings.warn(
                        "LEGACY FILE COMPATIBILITY: The file is missing the 'dstep' parameter.\n"
                        "No duplicate timesteps detected, safely assuming dstep=1 for all rows.\n"
                        "This file was likely created with an older UCLCHEM version.\n"
                        "Consider regenerating the file with the current version for full compatibility.",
                        UserWarning,
                    )
                    # Add dstep=1 column to the array
                    dstep_column = np.ones((array.shape[0], 1))
                    array = np.hstack(
                        [array[:, :point_index], dstep_column, array[:, point_index:]]
                    )
                    physics_cols_from_file.append("dstep")
                    point_index += 1  # point column shifted by 1
                else:
                    raise ValueError(
                        f"INCOMPATIBLE LEGACY FILE: Cannot infer 'dstep' parameter.\n\n"
                        f"The file is missing 'dstep' and contains duplicate timesteps,\n"
                        f"making it impossible to safely infer the particle step values.\n\n"
                        f"  File contains:        {physics_cols_from_file}\n"
                        f"  Current UCLCHEM has:  {list(PHYSICAL_PARAMETERS)}\n\n"
                        f"To load this file, regenerate it with the current UCLCHEM version."
                    )
            else:
                # Other parameter mismatch - cannot fix automatically
                raise ValueError(
                    f"INCOMPATIBLE LEGACY FILE: Physical parameters mismatch.\n\n"
                    f"The file you are loading has different PHYSICAL_PARAMETERS than the currently installed UCLCHEM.\n"
                    f"This means the file was created with a different version of UCLCHEM and cannot be loaded.\n\n"
                    f"  File contains:        {physics_cols_from_file} ({len(physics_cols_from_file)} parameters)\n"
                    f"  Current UCLCHEM has:  {list(PHYSICAL_PARAMETERS)} ({len(PHYSICAL_PARAMETERS)} parameters)\n"
                    f"  Missing parameters:   {list(missing_params)}\n"
                    f"  Extra parameters:     {list(extra_params)}\n\n"
                    f"PHYSICAL_PARAMETERS are hard-coded in the Fortran wrapper and cannot be changed at runtime.\n"
                    f"To load this file, you must use the same UCLCHEM version that created it, or regenerate the file.\n"
                )

        if species_cols_from_file != list(get_species_names()):
            raise ValueError(
                f"INCOMPATIBLE LEGACY FILE: Species list mismatch.\n\n"
                f"The file you are loading has a different species network than the currently installed UCLCHEM.\n"
                f"This means the file was created with a different chemical network and cannot be loaded.\n\n"
                f"  File contains:        {len(species_cols_from_file)} species\n"
                f"  Current UCLCHEM has:  {len(get_species_names())} species\n\n"
                f"The species list is tied to the chemical network compiled into UCLCHEM.\n"
                f"To load this file, you must use the same UCLCHEM version/network that created it, or regenerate the file.\n"
            )

        # At this point, we've validated that dimensions match
        # Allocate arrays using the validated file dimensions (which equal current constants)
        self.physics_array = np.empty(
            (row_count, self._param_dict["points"], len(PHYSICAL_PARAMETERS))
        )
        self.chemical_abun_array = np.empty(
            (row_count, self._param_dict["points"], n_species)
        )

        for p in range(self._param_dict["points"]):
            sel = np.where(array[:, point_index] == p + 1)[0]
            self.physics_array[:, p, :] = array[sel, :point_index]
            self.chemical_abun_array[:, p, :] = array[sel, point_index + 1 :]

        # Construct Dataset using current UCLCHEM constants (validated to match file)
        import xarray as xr

        self._data = xr.Dataset(
            {
                "physics_array": (
                    ["time_step", "point", "physics_values"],
                    self.physics_array,
                ),
                "chemical_abun_array": (
                    ["time_step", "point", "chemical_abun_values"],
                    self.chemical_abun_array,
                ),
            },
            coords={
                "physics_values": PHYSICAL_PARAMETERS,
                "chemical_abun_values": get_species_names(),
                "time_step": np.arange(self.physics_array.shape[0]),
                "point": np.arange(self._param_dict["points"]),
            },
        )

        # Clean trailing zero-padded timesteps and finish
        self._array_clean()

        # Ensure timepoints is consistent with the dataset we just constructed
        # Fortran expects arrays of shape (timepoints+1, points, ...), so we
        # store `timepoints` as the number of simulated timesteps minus one.
        try:
            tp = int(self._data.sizes.get("time_step", 0))
            object.__setattr__(self, "timepoints", tp - 1 if tp > 0 else 0)
        except Exception:
            # Be defensive; if something goes wrong leave timepoints as-is
            pass

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

    def legacy_read_starting_chemistry(self) -> None:
        """Method to read the starting chemistry from the self.abundLoadFile provided in _param_dict."""
        self._create_starting_array(np.loadtxt(self.abundLoadFile, delimiter=","))
        return

    def legacy_write_full(self) -> None:
        """Perform classic output file writing to the file self.outputFile provided in _param_dict."""
        phys = self.physics_array.reshape(-1, self.physics_array.shape[-1])
        chem = self.chemical_abun_array.reshape(-1, self.chemical_abun_array.shape[-1])
        full_array = np.append(phys, chem, axis=1)
        # TODO Move away from the magic numbers seen here.
        species_names = get_species_names()
        string_fmt_string = f"{', '.join([PHYSICAL_PARAMETERS_HEADER_FORMAT] * (len(PHYSICAL_PARAMETERS)))}, {', '.join([SPECNAME_HEADER_FORMAT] * len(species_names))}"
        # Magic numbers here to match/improve the formatting of the classic version
        # TODO Move away from the magic numbers seen here.
        number_fmt_string = f"{PHYSICAL_PARAMETERS_VALUE_FORMAT}, {', '.join([SPECNAME_VALUE_FORMAT] * len(species_names))}"
        columns = np.array([PHYSICAL_PARAMETERS[:-1] + ["point"] + species_names])
        np.savetxt(self.outputFile, columns, fmt=string_fmt_string)
        with open(self.outputFile, "ab") as f:
            np.savetxt(f, full_array, fmt=number_fmt_string)
        return

    def legacy_write_starting_chemistry(self) -> None:
        """Perform classic starting abundance file writing to the file `self.abundSaveFile`
        provided in `_param_dict`.
        """
        last_timestep_index = self.chemical_abun_array[:, 0, 0].nonzero()[0][-1]
        # TODO Move away from the magic numbers seen here.
        species_names = get_species_names()
        number_fmt_string = f" {', '.join(['%9.5E'] * len(species_names))}"
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
        """Assign coordinate values derived from arrays and metadata.

        Be conservative: prefer existing dataset coordinates if they are consistent
        with array shapes. Only override coords from metadata when lengths match
        the corresponding array dimensions to avoid xarray dimension conflicts.
        """
        # physics coords
        try:
            phys_len = self.physics_array.shape[-1]
        except Exception:
            phys_len = None
        if phys_len is not None:
            if "physics_values" in self._data.coords:
                if len(self._data.coords["physics_values"]) != phys_len:
                    if len(PHYSICAL_PARAMETERS) == phys_len:
                        self._data = self._data.assign_coords(
                            {"physics_values": PHYSICAL_PARAMETERS}
                        )
                    else:
                        # fallback to numeric indices if dimensions don't match
                        self._data = self._data.assign_coords(
                            {"physics_values": np.arange(phys_len)}
                        )
            else:
                if len(PHYSICAL_PARAMETERS) == phys_len:
                    self._data = self._data.assign_coords(
                        {"physics_values": PHYSICAL_PARAMETERS}
                    )
                else:
                    self._data = self._data.assign_coords(
                        {"physics_values": np.arange(phys_len)}
                    )

        # chemical abundances coords
        try:
            chem_len = self.chemical_abun_array.shape[-1]
        except Exception:
            chem_len = None
        if chem_len is not None:
            species_names = get_species_names()
            if "chemical_abun_values" in self._data.coords:
                if len(self._data.coords["chemical_abun_values"]) != chem_len:
                    if len(species_names) == chem_len:
                        self._data = self._data.assign_coords(
                            {"chemical_abun_values": species_names}
                        )
                    else:
                        self._data = self._data.assign_coords(
                            {"chemical_abun_values": np.arange(chem_len)}
                        )
            else:
                if len(species_names) == chem_len:
                    self._data = self._data.assign_coords(
                        {"chemical_abun_values": species_names}
                    )
                else:
                    self._data = self._data.assign_coords(
                        {"chemical_abun_values": np.arange(chem_len)}
                    )
        return

    def _reform_inputs(self, param_dict: dict, out_species: list[str]) -> None:
        """Internal Method.
        Copies param_dict so as not to modify user's dictionary.
        Then reformats out_species from pythonic list
        to a string of space separated names for Fortran.

        Args:
            param_dict (dict): Parameter dictionary passed by the user to the model.
            out_species (list[str]): List of output species that are considered important for this model.

        Raises:
            ValueError: If an duplicate key is encountered in `param_dict`.
            ValueError: If an entry in `out_species` is not a valid species name.

        """
        if param_dict is None:
            # avoid mutating the shared default dictionary
            self._param_dict = default_param_dictionary.copy()
        else:
            # lower case (and conveniently copy so we don't edit) the user's dictionary
            # this is key to UCLCHEM's "case insensitivity"
            new_param_dict = {}
            for k, v in param_dict.items():
                if k.lower() in new_param_dict:
                    raise ValueError(
                        f"Duplcate lower case key {k} is already in the dict, stopping"
                    )
                if isinstance(v, Path):
                    v = str(v)
                new_param_dict[k.lower()] = v

            # Handle deprecated endAtFinalDensity parameter (after lowercasing)
            new_param_dict = _convert_legacy_stopping_param(new_param_dict)

            self._param_dict = {**default_param_dictionary, **new_param_dict.copy()}
            del new_param_dict
        # Remove keys with None values from the merged _param_dict
        # Check the merged dict, not the defaults, to preserve user-provided values
        keys_to_delete = [k for k, v in self._param_dict.items() if v is None]
        for k in keys_to_delete:
            del self._param_dict[k]
        if out_species is not None:
            # Validate out_species: list/tuple of strings and known species
            if not (
                isinstance(out_species, (list, tuple))
                and all(isinstance(s, str) for s in out_species)
                and all(s.strip() in get_species_names() for s in out_species)
            ):
                raise ValueError(
                    "out_species must be a list/tuple of valid species names; check available species via uclchem.model.get_species_names()"
                )
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
        return

    def _create_rate_constants_array(self):
        """Internal Method.
        Creates Fortran compliant np.array for rate constants that can
        be passed to the Fortran part of UCLCHEM.
        """
        # For shared memory:
        (
            self._shm_handles["rate_constants_array"],
            self._shm_desc["rate_constants_array"],
            self.rate_constants_array,
        ) = self._create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], n_reactions)
        )
        return

    def _create_heating_array(self):
        """Internal Method.
        Creates Fortran compliant np.array for heating/cooling rates that can
        be passed to the Fortran part of UCLCHEM.
        """
        try:
            heating_array_size = (
                2
                + uclchemwrap.heating.ncooling
                + uclchemwrap.heating.nheating
                + uclchemwrap.f2py_constants.ncoolants
            )
            # For shared memory:
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
        except AttributeError:
            # Heating module not available, likely compiled without heating support
            logging.debug("Heating module not available in uclchemwrap")
            self.heat_array = None
        return

    def _create_stats_array(self):
        """Internal Method.
        Creates Fortran compliant np.array for DVODE solver statistics.
        """
        (
            self._shm_handles["stats_array"],
            self._shm_desc["stats_array"],
            self.stats_array,
        ) = self._create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], N_DVODE_STATS)
        )
        return

    def _create_level_populations_array(self):
        """Internal Method.
        Creates Fortran compliant np.array for coolant level populations.
        Shape: (timepoints+1, gridpoints, total_levels).
        """
        (
            self._shm_handles["level_populations_array"],
            self._shm_desc["level_populations_array"],
            self.level_populations_array,
        ) = self._create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], N_TOTAL_LEVELS)
        )
        return

    def _create_se_stats_array(self):
        """Internal Method.
        Creates Fortran compliant np.array for SE solver statistics.
        Shape: (timepoints+1, gridpoints, NCOOLANTS*3).
        """
        n_stats = NCOOLANTS * N_SE_STATS_PER_COOLANT  # 35 * 3 = 105

        (
            self._shm_handles["se_stats_array"],
            self._shm_desc["se_stats_array"],
            self.se_stats_array,
        ) = self._create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], n_stats)
        )
        return

    def get_level_populations_dataframe(self, point: int = 0) -> pd.DataFrame | None:
        """Get level populations as a DataFrame for a specific grid point.

        Args:
            point: Grid point index (default 0)

        Returns:
            pd.DataFrame | None: DataFrame with columns for each level with meaningful
                coolant/level names

        """
        if (
            self.level_populations_array is None
            or self.level_populations_array.shape[0] < 3
        ):
            return None

        # Try to load meaningful level names from coolant files
        try:
            level_names_dict = load_coolant_level_names()
            # Build column names from parsed coolant files
            columns = []
            for i in range(NCOOLANTS):
                if i in level_names_dict:
                    columns.extend(level_names_dict[i])
        except (FileNotFoundError, ValueError, RuntimeError) as e:
            # Fallback to generic names if coolant files can't be loaded
            logging.warning(
                f"Could not load coolant level names: {e}. Using generic names."
            )
            columns = [f"LEVEL_{i}" for i in range(N_TOTAL_LEVELS)]

        return pd.DataFrame(self.level_populations_array[:, point, :], columns=columns)

    def get_se_stats_dataframe(self, point: int = 0) -> pd.DataFrame | None:
        """Get SE solver statistics as a DataFrame for a specific grid point.

        Args:
            point (int): Grid point index. Default = 0.

        Returns:
            pd.DataFrame | None: DataFrame with per-coolant SE solver statistics using
                actual coolant names

        """
        if self.se_stats_array is None or self.se_stats_array.shape[0] < 3:
            return None

        # Build meaningful column names using actual coolant names
        try:
            from uclchemwrap import f2py_constants

            coolant_names = [
                str(name.decode()).strip() for name in f2py_constants.coolantnames
            ]
            columns = []
            for coolant_name in coolant_names:
                columns.extend(
                    [
                        f"{coolant_name}_CONVERGED",
                        f"{coolant_name}_ITERATIONS",
                        f"{coolant_name}_MAX_REL_CHANGE",
                    ]
                )
        except Exception:
            # Fallback to generic names
            columns = SE_STAT_NAMES

        return pd.DataFrame(self.se_stats_array[:, point, :], columns=columns)

    def _create_starting_array(self, starting_chemistry: np.ndarray | None) -> None:
        if starting_chemistry is None:
            self.starting_chemistry_array = None
        else:
            if len(np.shape(starting_chemistry)) == 1:
                starting_chemistry = starting_chemistry[np.newaxis, :]
            # For shared memory:
            (
                self._shm_handles["starting_chemistry_array"],
                self._shm_desc["starting_chemistry_array"],
                self.starting_chemistry_array,
            ) = self._create_shared_memory_allocation(np.shape(starting_chemistry))
            np.copyto(self.starting_chemistry_array, starting_chemistry, casting="no")
        return

    # /Creation of arrays

    # Signal Interrupt Catch
    def on_interrupt(self, grid: bool = False, model_name: str | None = None) -> None:
        """Catch interruption. Save model to file.

        Args:
            grid (bool): whether the model was part of a grid
            model_name (str | None): the name of the model to save it under.
                If None, name is set to "interrupted".

        """
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

    def _reform_array_in_worker(self, shm_desc: dict[str, dict]) -> None:
        object.__setattr__(self, "_shm_handles", {})
        for k, v in shm_desc.items():
            shm = shared_memory.SharedMemory(name=v["name"], create=False)
            object.__setattr__(
                self,
                k,
                np.ndarray(shape=v["shape"], dtype=np.float64, buffer=shm.buf, order="F"),
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
        param_dict (dict): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in `defaultparameters.f90`.
        out_species (list | None): List of species whose abundances at the end of the model are
            returned. If None, defaults to `uclchem.constants.default_elements_to_check`.
            Default = None.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from AbstractModel,
            to use for the starting abundances of the new UCLCHEM model that will be run.
            Defaults to None.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.

    """

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.
        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            debug=debug,
            read_file=read_file,
            run_type=run_type,
        )
        if self.run_type != "external":
            self.run()

    def run_fortran(self) -> dict[str, int | list]:
        """Runs the UCLCHEM model, first by resetting the numpy arrays by using
        `AbstractModel.run()`, then running the model. `check_error` and `array_clean`
        are automatically called after the model run.

        Returns:
            dict[str, int | list]: Dictionary with two keys:
                "success_flag" with value the success flag
                "out_species_abundances_array" with value a list of the outspecies abundances.


        """
        # f2py returns all non-intent(in) values in Fortran signature order:
        # [0-6] physicsarray..sestatsarray (in,out, modified in-place),
        # [7] abundance_out, [8] specname_out, [9] successFlag
        result = wrap.cloud(
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrateconstants=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            rateconstantsarray=self.rate_constants_array,
            heatarray=self.heat_array,
            statsarray=self.stats_array,
            levelpopulationsarray=self.level_populations_array,
            sestatsarray=self.se_stats_array,
            abundancestart=self.starting_chemistry_array
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
            else None,
        )
        abundance_out, specname_out, success_flag = result[-3], result[-2], result[-1]
        if success_flag < 0:
            out_species_abundances_array = np.array([])
        else:
            out_species_length = (
                len(self.out_species_list) if self.out_species_list is not None else 0
            )
            out_species_abundances_array = list(abundance_out[:out_species_length])
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
        collapse (str): A string containing the collapse type.
            Options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'. Defaults to 'BE1.1'.
        physics_output (str): Filename to store physics output, only relevant for
            'filament' and 'ambipolar' collapses. If None, no physics output will be saved.
        param_dict (dict): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in `defaultparameters.f90`.
        out_species (list | None): List of species whose abundances at the end of the model are
            returned. If None, defaults to `uclchem.constants.default_elements_to_check`.
            Default = None.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from AbstractModel,
            to use for the starting abundances of the new UCLCHEM model that will be run.
            Defaults to None.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.

    """

    def __init__(
        self,
        collapse: Literal["BE1.1", "BE4", "filament", "ambipolar"] = "BE1.1",
        physics_output: str = None,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Raises:
            ValueError: If `collapse` is not one of `["BE1.1", "BE4", "filament", "ambipolar"]`.

        """
        if out_species is None:
            out_species = default_elements_to_check
        if collapse not in ["filament", "ambipolar"] and physics_output is None:
            warnings.warn(
                "`physics_output` was None but `collapse` was `filament` or `ambipolar`. No output file will be created.",
                UserWarning,
            )
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            debug=debug,
            read_file=read_file,
            run_type=run_type,
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

    def run_fortran(self) -> dict[str, int | list]:
        """Runs the UCLCHEM model, first by resetting the numpy arrays by using
        `AbstractModel.run()`, then running the model. `check_error` and `array_clean`
        are automatically called after the model run.

        Returns:
            dict[str, int | list]: Dictionary with two keys:
                "success_flag" with value the success flag
                "out_species_abundances_array" with value a list of the outspecies abundances.


        """
        result = wrap.collapse(
            collapsein=self.collapse,
            collapsefilein=self.physics_output,
            writeout=self.write_physics,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrateconstants=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            rateconstantsarray=self.rate_constants_array,
            heatarray=self.heat_array,
            statsarray=self.stats_array,
            levelpopulationsarray=self.level_populations_array,
            sestatsarray=self.se_stats_array,
            abundancestart=self.starting_chemistry_array
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
            else None,
        )
        out_species_abundances_array = result[-2]
        success_flag = result[-1]
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
    """PrestellarCore model class inheriting from AbstractModel.
    This model type was previously known as hot core.

    Args:
        temp_indx (int): Used to select the mass of the prestellar core from the following selection
            [1=1Msun, 2=5, 3=10, 4=15, 5=25,6=60]. Defaults to 1, which is 1 Msun
        max_temperature (float): Value at which gas temperature will stop increasing. Defaults to 300.0.
        param_dict (dict): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in `defaultparameters.f90`.
        out_species (list | None): List of species whose abundances at the end of the model are
            returned. If None, defaults to `uclchem.constants.default_elements_to_check`.
            Default = None.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from AbstractModel,
            to use for the starting abundances of the new UCLCHEM model that will be run.
            Defaults to None.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.

    """

    def __init__(
        self,
        temp_indx: int = 1,
        max_temperature: float = 300.0,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Raises:
            ValueError: If `read_file` is None, but `temp_idx` or `max_temperature` is also None.

        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            debug=debug,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None:
            if temp_indx is None or max_temperature is None:
                raise ValueError(
                    "temp_indx and max_temperature must be specified if not reading from file."
                )
            self.temp_indx = temp_indx
            self.max_temperature = max_temperature
            if self.run_type != "external":
                self.run()
        return

    def run_fortran(self) -> dict[str, int | list]:
        """Runs the UCLCHEM model, first by resetting the numpy arrays by using
        `AbstractModel.run()`, then running the model. `check_error` and `array_clean`
        are automatically called after the model run.

        Returns:
            dict[str, int | list]: Dictionary with two keys:
                "success_flag" with value the success flag
                "out_species_abundances_array" with value a list of the outspecies abundances.


        """
        _, _, _, _, _, _, _, out_species_abundances_array, _, success_flag = (
            wrap.hot_core(
                temp_indx=self.temp_indx,
                max_temp=self.max_temperature,
                dictionary=self._param_dict,
                outspeciesin=self.out_species,
                timepoints=self.timepoints,
                gridpoints=self._param_dict["points"],
                returnarray=True,
                returnrateconstants=True,
                givestartabund=self.give_start_abund,
                physicsarray=self.physics_array,
                chemicalabunarray=self.chemical_abun_array,
                rateconstantsarray=self.rate_constants_array,
                heatarray=self.heat_array,
                statsarray=self.stats_array,
                levelpopulationsarray=self.level_populations_array,
                sestatsarray=self.se_stats_array,
                abundancestart=self.starting_chemistry_array
                if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
                else None,
            )
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
        shock_vel (float): Velocity of the shock in km/s. Defaults to 10.0.
        timestep_factor (float): Whilst the time is less than 2 times the dissipation time of shock,
            timestep is timestep_factor*dissipation time. Essentially controls how well resolved the
            shock is in your model. Defaults to 0.01.
        minimum_temperature (float): Minimum post-shock temperature. Defaults to 0.0 (no minimum). The
            shocked gas typically cools to `initialTemp` if this is not set.
        param_dict (dict): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in `defaultparameters.f90`.
        out_species (list | None): List of species whose abundances at the end of the model are
            returned. If None, defaults to `uclchem.constants.default_elements_to_check`.
            Default = None.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from AbstractModel,
            to use for the starting abundances of the new UCLCHEM model that will be run.
            Defaults to None.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.

    """

    def __init__(
        self,
        shock_vel: float = 10.0,
        timestep_factor: float = 0.01,
        minimum_temperature: float = 0.0,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Raises:
            ValueError: If `read_file` is None, but `shock_vel` is also not set.

        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            debug=debug,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None:
            if shock_vel is None:
                raise ValueError("shock_vel must be specified if not reading from file.")
            self.shock_vel = shock_vel
            self.timestep_factor = timestep_factor
            self.minimum_temperature = minimum_temperature
            self.dissipation_time = -1
            if self.run_type != "external":
                self.run()
        return

    def run_fortran(self) -> dict[str, int | list]:
        """Runs the UCLCHEM model, first by resetting the numpy arrays by using
        `AbstractModel.run()`, then running the model. `check_error` and `array_clean`
        are automatically called after the model run.

        Returns:
            dict[str, int | list]: Dictionary with two keys:
                "success_flag" with value the success flag
                "out_species_abundances_array" with value a list of the outspecies abundances.


        """
        result = wrap.cshock(
            shock_vel=self.shock_vel,
            timestep_factor=self.timestep_factor,
            minimum_temperature=self.minimum_temperature,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrateconstants=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            rateconstantsarray=self.rate_constants_array,
            heatarray=self.heat_array,
            statsarray=self.stats_array,
            levelpopulationsarray=self.level_populations_array,
            sestatsarray=self.se_stats_array,
            abundancestart=self.starting_chemistry_array
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
            else None,
        )
        dissipation_time = result[-3]
        out_species_abundances_array = result[-2]
        success_flag = result[-1]
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
        param_dict (dict | None): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in `defaultparameters.f90`.
        out_species (list | None): List of species whose abundances at the end of the model are
            returned. If None, defaults to `uclchem.constants.default_elements_to_check`.
            Default = None.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from AbstractModel,
            to use for the starting abundances of the new UCLCHEM model that will be run.
            Defaults to None.
        timepoints (int): Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to `uclchem.constants.TIMEPOINTS`.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.

    """

    def __init__(
        self,
        shock_vel: float = 10.0,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        debug: bool = False,
        read_file: str = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Raises:
            ValueError: If `read_file` is None, but `shock_vel` is also not set.

        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            debug=debug,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None:
            if shock_vel is None:
                raise ValueError("shock_vel must be specified if not reading from file.")
            self.shock_vel = shock_vel
            if self.run_type != "external":
                self.run()
        return

    def run_fortran(self) -> dict[str, int | list]:
        """Runs the UCLCHEM model, first by resetting the numpy arrays by using
        `AbstractModel.run()`, then running the model. `check_error` and `array_clean`
        are automatically called after the model run.

        Returns:
            dict[str, int | list]: Dictionary with two keys:
                "success_flag" with value the success flag
                "out_species_abundances_array" with value a list of the outspecies abundances.


        """
        result = wrap.jshock(
            shock_vel=self.shock_vel,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrateconstants=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            rateconstantsarray=self.rate_constants_array,
            heatarray=self.heat_array,
            statsarray=self.stats_array,
            levelpopulationsarray=self.level_populations_array,
            sestatsarray=self.se_stats_array,
            abundancestart=self.starting_chemistry_array
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
            else None,
        )
        out_species_abundances_array = result[-2]
        success_flag = result[-1]
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

    Postprocess allows for additional controls of the time, density, gas temperature, radiation field,
    cosmic ray ionisation rate, atomic and molecular Hydrogen, CO and C column densities through the
    use of arrays. Using these arrays allows for experimental model crafting beyond the standard models
    in other model classes.

    Args:
        param_dict (dict | None): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in `defaultparameters.f90`.
        out_species (list | None): List of species whose abundances at the end of the model are
            returned. If None, defaults to `uclchem.constants.default_elements_to_check`.
            Default = None.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances
            to use for the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from AbstractModel,
            to use for the starting abundances of the new UCLCHEM model that will be run.
            Defaults to None.
        time_array (np.ndarray | None): Represents the time grid to be used for the model.
            This sets the target timesteps for which outputs will be stored.
        density_array (np.ndarray | None): Represents the value of the density at different
            timepoints found in time_array.
        gas_temperature_array (np.ndarray | None): Represents the value of the gas temperature
            at different timepoints found in time_array.
        dust_temperature_array (np.ndarray | None):Represents the value of the dust temperature
            at different timepoints found in time_array.
        zeta_array (np.ndarray | None): Represents the value of the cosmic ray ionisation rate
            at different timepoints found in time_array.
        radfield_array (np.ndarray | None): Represents the value of the UV radiation field
            at different timepoints found in time_array.
        coldens_H_array (np.ndarray | None): Represents the value of the column density of H
            at different timepoints found in time_array.
        coldens_H2_array (np.ndarray | None): Represents the value of the column density of H2
            at different timepoints found in time_array.
        coldens_CO_array (np.ndarray | None): Represents the value of the column density of CO
            at different timepoints found in time_array.
        coldens_C_array (np.ndarray | None): Represents the value of the column density of C
            at different timepoints found in time_array.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.

    """

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        time_array: np.ndarray | None = None,
        density_array: np.ndarray | None = None,
        gas_temperature_array: np.ndarray | None = None,
        dust_temperature_array: np.ndarray | None = None,
        zeta_array: np.ndarray | None = None,
        radfield_array: np.ndarray | None = None,
        visual_extinction_array: np.ndarray | None = None,
        coldens_H_array: np.ndarray | None = None,
        coldens_H2_array: np.ndarray | None = None,
        coldens_CO_array: np.ndarray | None = None,
        coldens_C_array: np.ndarray | None = None,
        debug: bool = False,
        read_file: str | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Raises:
            ValueError: If not all arrays have the same length.
            ValueError: If `read_file` is None, but `time_array` is not an array.

        """
        # Allocate 1.5x the input timesteps to give the DVODE solver
        # headroom for additional internal substeps.
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=int(1.5 * len(time_array)),
            debug=debug,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None and time_array is not None:
            n_input = len(time_array)
            object.__setattr__(
                self,
                "postprocess_arrays",
                {
                    "timegrid": time_array,
                    "densgrid": density_array,
                    "gastempgrid": gas_temperature_array,
                    "dusttempgrid": dust_temperature_array,
                    "radfieldgrid": radfield_array,
                    "zetagrid": zeta_array,
                    "avgrid": visual_extinction_array,
                    "nhgrid": coldens_H_array,
                    "nh2grid": coldens_H2_array,
                    "ncogrid": coldens_CO_array,
                    "ncgrid": coldens_C_array,
                },
            )
            for key, array in self.postprocess_arrays.items():
                if array is not None:
                    if isinstance(array, float):
                        array = np.ones(n_input) * array
                    if not len(array) == n_input:
                        raise ValueError("All arrays must be the same length")

                    # Pad to self.timepoints so Fortran gets
                    # correctly sized arrays.
                    padded = np.zeros(
                        self.timepoints,
                        dtype=np.float64,
                    )
                    padded[:n_input] = array
                    self.postprocess_arrays[key] = np.asfortranarray(padded)
            self.time_array = time_array
            # Column-density (coldens) and visual extinction (Av) arrays
            self.coldens_H_array = coldens_H_array
            self.visual_extinction_array = visual_extinction_array
            # Flags exposed for Fortran wrapper (mutually exclusive)
            self.usecoldens = self.coldens_H_array is not None
            self.useav = self.visual_extinction_array is not None
            assert not (self.usecoldens and self.useav), (
                "Cannot use both column density and visual extinction arrays simultaneously."
            )

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

    def run_fortran(self) -> dict[str, int | list]:
        """Runs the UCLCHEM model, first by resetting the numpy arrays by using
        `AbstractModel.run()`, then running the model. `check_error` and `array_clean`
        are automatically called after the model run.

        Returns:
            dict[str, int | list]: Dictionary with two keys:
                "success_flag" with value the success flag
                "out_species_abundances_array" with value a list of the outspecies abundances.


        """
        # Determine whether an Av grid was provided and set the flag expected by the Fortran wrapper
        # Only pass arrays that are present (not None) to the Fortran wrapper
        post_kwargs = {k: v for k, v in self.postprocess_arrays.items() if v is not None}
        result = wrap.postprocess(
            usecoldens=self.usecoldens,
            useav=self.useav,
            **post_kwargs,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrateconstants=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            rateconstantsarray=self.rate_constants_array,
            heatarray=self.heat_array,
            statsarray=self.stats_array,
            levelpopulationsarray=self.level_populations_array,
            sestatsarray=self.se_stats_array,
            abundancestart=self.starting_chemistry_array
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
            else None,
        )
        out_species_abundances_array = result[-2]
        success_flag = result[-1]
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
            "useav": self.useav,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
            "_debug": self._debug,
            "postprocess_arrays": self.postprocess_arrays,
            **self.postprocess_arrays,
        }


@register_model
class Model(AbstractModel):
    """Model, like Postprocess, represents a model class with additional controls.
    It inherits from AbstractModel.

    Model follows the same logic as Postprocess but without the coldens Arguments.
    It allows for additional controls of the time, density, gas temperature, radiation field,
    and cosmic ray ionisation rate through the use of arrays. Using these arrays allows for
    experimental model crafting beyond the standard models in other model classes.

    Args:
        param_dict (dict): Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in `defaultparameters.f90`.
        out_species (list | None): List of species whose abundances at the end of the model are
            returned. If None, defaults to `uclchem.constants.default_elements_to_check`.
            Default = None.
        starting_chemistry (np.ndarray | None): Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model (AbstractModel | None): Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model
            that will be run. Defaults to None.
        time_array (np.ndarray | None): Represents the time grid to be used for the model.
            This sets the target timesteps for which outputs will be stored.
        density_array (np.ndarray | None): Represents the value of the density at different
            timepoints found in time_array.
        gas_temperature_array (np.ndarray | None): Represents the value of the gas temperature
            at different timepoints found in time_array.
        dust_temperature_array (np.ndarray | None):Represents the value of the dust temperature
            at different timepoints found in time_array.
        zeta_array (np.ndarray | None): Represents the value of the cosmic ray ionisation rate
            at different timepoints found in time_array.
        radfield_array (np.ndarray | None): Represents the value of the UV radiation field at
            different timepoints found in time_array.
        debug (bool): Flag if extra debug information should be printed to the terminal.
            Defaults to False. #TODO Add debug features
        read_file (str | None): Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None.

    """

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        time_array: np.ndarray | None = None,
        density_array: np.ndarray | None = None,
        gas_temperature_array: np.ndarray | None = None,
        dust_temperature_array: np.ndarray | None = None,
        zeta_array: np.ndarray | None = None,
        radfield_array: np.ndarray | None = None,
        debug: bool = False,
        read_file: str | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Raises:
            ValueError: If not all arrays have the same length.
            ValueError: If `read_file` is None, but `time_array` is not an array.

        """
        # Allocate 1.5x the input timesteps to give the DVODE solver
        # headroom for additional internal substeps.
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=int(1.5 * len(time_array)),
            debug=debug,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None and time_array is not None:
            n_input = len(time_array)
            self.time_array = time_array
            self.postprocess_arrays = {
                "timegrid": time_array,
                "densgrid": density_array,
                "gastempgrid": gas_temperature_array,
                "dusttempgrid": dust_temperature_array,
                "radfieldgrid": radfield_array,
                "zetagrid": zeta_array,
            }
            for key, array in self.postprocess_arrays.items():
                if array is not None:
                    if isinstance(array, float):
                        array = np.ones(n_input) * array
                    if not len(array) == n_input:
                        raise ValueError("All arrays must be the same length")
                    # Pad to self.timepoints so Fortran gets
                    # correctly sized arrays.
                    padded = np.zeros(
                        self.timepoints,
                        dtype=np.float64,
                    )
                    padded[:n_input] = array
                    self.postprocess_arrays[key] = np.asfortranarray(padded)
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

    def run_fortran(self) -> dict[str, int | list]:
        """Runs the UCLCHEM model, first by resetting the numpy arrays by using
        `AbstractModel.run()`, then running the model. `check_error` and `array_clean`
        are automatically called after the model run.

        Returns:
            dict[str, int | list]: Dictionary with two keys:
                "success_flag" with value the success flag
                "out_species_abundances_array" with value a list of the outspecies abundances.


        """
        result = wrap.postprocess(
            usecoldens=False,
            **self.postprocess_arrays,
            dictionary=self._param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self._param_dict["points"],
            returnarray=True,
            returnrateconstants=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            rateconstantsarray=self.rate_constants_array,
            heatarray=self.heat_array,
            statsarray=self.stats_array,
            levelpopulationsarray=self.level_populations_array,
            sestatsarray=self.se_stats_array,
            abundancestart=self.starting_chemistry_array
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")
            else None,
        )
        out_species_abundances_array = result[-2]
        success_flag = result[-1]
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
class SequentialRunner:
    """The SequentialRunner class allows for multiple models to be run back to back.

    By defining a specific dictionary to hold the information of each model class to run in sequence, SewuentialModel allows
    for the automatic running of multiple models as well as matching some physical parameters from one model to the next.

    Args:
        sequenced_model_parameters (list[dict[str, Any]]): The List of dictionaries to pass to SequentialRunner takes the format of
            [{"<First Model Class>":{"param_dict":{<parameters>}, <other arguments>}}, {"<Second Model Class>:{"param_dict":{<parameters>}, <other arguments>}}, ...}]
        parameters_to_match (list[str]): The list provided to this argument decides which parameters
            should be matched from a previous model to the next model in the sequence.
            Currently, supports one or more of `["finalDens", "finalTemp"]`.
        run_type (Literal['managed', 'external']): Run type.
            Must be either "managed", or "external".


    Raises:
        NotImplementedError: If a parameter in `parameters_to_match` is not one of
            `["finalDens", "finalTemp"]`.

    """  # noqa: W505

    def __init__(
        self,
        sequenced_model_parameters: list,
        parameters_to_match: list = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        for model in sequenced_model_parameters:
            assert model[list(model.keys())[0]] != SequentialRunner
        self.models = []
        self.sequenced_model_parameters = sequenced_model_parameters
        self.parameters_to_match = parameters_to_match

        if self.parameters_to_match is not None:
            for parameter in self.parameters_to_match:
                if parameter not in ["finalTemp", "finalDens"]:
                    raise NotImplementedError(
                        f"Parameter '{parameter}' has not been implemented for parameter matching"
                    )

        self.run_type = run_type
        self.model_count = 0
        self._pickle_dict = {}
        self.success_flag = None
        if self.run_type == "managed":
            self.run()

    def run(self) -> None:
        """Run the sequential model.

        Raises:
            NotImplementedError: If a parameter in `parameters_to_match` is not one of
                `["finalDens", "finalTemp"]`.

        """
        previous_model = None
        for base_model_dict in self.sequenced_model_parameters:
            for model_type, model_dict in base_model_dict.items():
                model_dict["param_dict"] = {
                    k.lower(): v for k, v in model_dict["param_dict"].items()
                }
                if self.model_count > 0:
                    previous_param_dict = previous_model._param_dict.copy()
                    # Remove the converted stopping-mode key inherited
                    # from the previous stage before merging
                    previous_param_dict.pop("parcelstoppingmode", None)
                    model_dict["param_dict"] = {
                        **previous_param_dict,
                        **model_dict["param_dict"],
                    }
                    if self.parameters_to_match is not None:
                        for parameter in self.parameters_to_match:
                            if parameter == "finalDens":
                                model_dict["param_dict"]["initialdens"] = (
                                    previous_model.physics_array[-1, 0, 1].item()
                                )
                                continue
                            elif parameter == "finalTemp":
                                model_dict["param_dict"]["initialtemp"] = (
                                    previous_model.physics_array[-1, 0, 2].item()
                                )
                            else:
                                raise NotImplementedError(
                                    f"Parameter '{parameter}' has not been implemented for parameter matching"
                                )
                    tmp_model = REGISTRY[model_type](
                        **model_dict,
                        run_type=self.run_type,
                        previous_model=previous_model,
                    )
                else:
                    tmp_model = REGISTRY[model_type](
                        **model_dict,
                        run_type=self.run_type,
                        previous_model=previous_model,
                    )
                    if self.run_type == "external":
                        tmp_model.run()
                    self.models += [
                        {
                            "Model_Type": model_type,
                            "Model_Order": self.model_count,
                            "Model": tmp_model,
                            "Success": tmp_model.success_flag,
                        }
                    ]
                    self.models[self.model_count]["Successful"] = (
                        self.models[self.model_count]["Model"].success_flag == 0
                    )
                    previous_model = self.models[self.model_count]["Model"]
                self.model_count += 1
        self.success_flag = all(d["Successful"] for d in self.models)
        return

    def save_model(
        self,
        *,
        file_obj: h5py.File | None = None,
        file: str | None = None,
        name: str = "",
        overwrite: bool = False,
    ) -> None:
        """Save a model to an open file object or to a file.

        Args:
            file_obj (h5py.File | None): open file h5py file object. Default = None.
            file (str | None): file to write to. Default = None.
            name (str): name to save model under.
            overwrite (bool): Boolean on whether to overwrite pre-existing models, or error out.
                Defaults to False

        Raises:
            ValueError: If file_obj and file are both passed, or neither are passed.

        """
        opened_file = False
        if (file_obj is None) == (file is None):
            raise ValueError("file_obj or file must be passed.")
        elif file_obj is None:
            file_obj = h5py.File(file, "a")
            opened_file = True

        for model in self.models:
            model["Model"].save_model(
                file_obj=file_obj,
                name=f"{name}_{model['Model_Order']}_{model['Model_Type']}",
                overwrite=overwrite,
            )

        if opened_file:
            file_obj.close()

    def check_conservation(
        self, element_list: list[str] | None = None, percent: bool = True
    ) -> None:
        """Check conservation of the chemical abundances.

        Args:
            element_list (list[str] | None): List of elements to check conservation for.
                If None, use `uclchem.constants.default_elements_to_check`. Default = None.
            percent (bool): Flag on if percentage values should be used. Defaults to True.

        """
        if element_list is None:
            element_list = default_elements_to_check

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
                conserved = all(float(x[:1]) < 1 for x in i.values())
            model["elements_conserved"] = conserved

    def pickle(self) -> None:
        """Pickle the models."""
        if not bool(self._pickle_dict):
            for model in self.models:
                model["Model"] = model["Model"].pickle()
                self._pickle_dict[f"{model['Model_Order']}_{model['Model_Type']}"] = (
                    model["Model"]._pickle_dict.copy()
                )

    def un_pickle(self) -> None:
        """Un-pickle the models."""
        if bool(self._pickle_dict):
            for model in self.models:
                model["Model"]._pickle_dict = self._pickle_dict[
                    f"{model['Model_Order']}_{model['Model_Type']}"
                ]
                model["Model"] = model["Model"].un_pickle()

    def _coordinator_unlink_memory(self):
        for model in self.models:
            model["Model"]._coordinator_unlink_memory()


def _run_grid_model(
    model_id: int,
    model_type: str,
    pending_model: dict[str, Any],
    log_dir: str | Path | None = None,
) -> tuple[int, object]:
    """Run a single model. This is used by the GridRunner class.

    Args:
        model_id (int): id of model
        model_type (str): string representing the type of model
        pending_model (dict[str, Any]): dictionary with arguments necessary to initialize
            model.
        log_dir (str | Path | None): If not None, write logs to "model_{model_id}.log".
            If None, do not write logs. Default = None.

    Returns:
        model_id (int): model id of run model
        model_obj (object): pickled model object

    """
    log_file = None
    if log_dir is not None:
        log_file = os.path.join(log_dir, f"model_{model_id}.log")

    cls = REGISTRY.get(model_type)
    with capture_fortran_output(label=f"model_{model_id}", log_file=log_file):
        model_obj = cls(
            **pending_model,
            run_type=("managed" if model_type == "SequentialRunner" else "external"),
        )
        if model_type != "SequentialRunner":
            model_obj.run()
    model_obj._coordinator_unlink_memory()
    model_obj.pickle()
    return model_id, model_obj


# The following parameters from various chemical models, cannot be used as grid parameters.
NoGridParameters = [
    "out_species",
    "starting_chemistry",
    "time_array",
    "density_array",
    "gas_temperature_array",
    "dust_temperature_array",
    "zeta_array",
    "radfield_arrayvisual_extinction_array",
    "coldens_H_array",
    "coldens_H2_array",
    "coldens_CO_array",
    "coldens_C_array",
    "debug",
    "read_file",
    "run_type",
]


class _NoDaemonProcess(mp.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):  # noqa :ANN001
        pass


class NoDaemonPool(pool.Pool):  # noqa
    @staticmethod
    def Process(ctx, *args, **kwargs):  # noqa
        return _NoDaemonProcess(*args, **kwargs)


class GridRunner:
    """GridRunner, like SequentialRunner is not an actual uclchem model,
    instead it allows running multiple models on a grid of parameter space.

    Args:
        model_type (str of model class to run):
        full_parameters (Dict): The dictionary passed to GridRunner should nest into it,
            the param_dict argument that would be passed to any other model, with the addition
            of extra keys for the none param_dict variables of a model. Any variables that are
            turned into lists or arrays, will automatically be assumed to be used for the gridding.
        max_workers (int): Maximum number of workers to use in parallel for the grid run.
            Defaults to 8.
        grid_file (str): Name and path of the output file to which the models should be saved.
            Defaults to "./default_grid_out.h5".
        model_name_prefix (str): Name prefix convention to use. The fifth model in the grid
            would have the name "<model_name_prefix>5>" assigned to it. Defaults to "",
            which would make the fifth model have the name "5", for example.
        delay_run (bool): Whether to immediately start the models upon initialization,
            or delay until the user calls `self.run()`. Defaults to False (start immediately).
        log_dir (str | None): Where to write logs. If None, do not write logs. Default = None.
        model_ids (list | None): Optional subset of model indices (0-based column in flat_grids)
            to run. None means run all models in the grid. Default = None.

    """

    def __init__(
        self,
        model_type: AnyStr,
        full_parameters: dict | list,
        max_workers: int = 8,
        grid_file: str = "./default_grid_out.h5",
        model_name_prefix: str = "",
        overwrite_models: bool = False,
        delay_run: bool = False,
        log_dir: str | None = None,
        model_ids: list = None,
        create_grid: bool = True,
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
        self.overwrite_models = overwrite_models
        self.log_dir = log_dir
        if self.log_dir is not None:
            os.makedirs(self.log_dir, exist_ok=True)
            self._main_log = os.path.join(self.log_dir, "grid.log")
        else:
            self._main_log = None
        self._orig_sigint = signal.getsignal(signal.SIGINT)
        self.parameters_to_grid = {}

        if self.model_type == "SequentialRunner":
            if not isinstance(self.full_parameters, list):
                raise TypeError(
                    f"For SequentialRunner types, full_parameters must be a list. {type(self.full_parameters)} was passed."
                )
            for model_count in range(len(self.full_parameters)):
                for model_type, model_full_params in self.full_parameters[
                    model_count
                ].items():
                    if not isinstance(model_full_params, dict):
                        continue
                    for k, v in model_full_params.items():
                        if k == "param_dict":
                            for k_p, v_p in v.items():
                                self._grid_def(k_p, v_p, model_count)
                        else:
                            self._grid_def(k, v, model_count)

        #               grids = np.meshgrid(*self.parameters_to_grid.values(), indexing="xy")
        else:
            if not isinstance(self.full_parameters, dict):
                raise TypeError(
                    f"For none SequentialRunner types, full_parameters must be a dictionary. {type(self.full_parameters)} was passed."
                )
            for k, v in self.full_parameters.items():
                if k == "param_dict":
                    for k_p, v_p in v.items():
                        self._grid_def(k_p, v_p)
                else:
                    self._grid_def(k, v)
        if create_grid:
            grids = np.meshgrid(*self.parameters_to_grid.values(), indexing="xy")
            self.flat_grids = np.reshape(
                grids,
                (
                    len(self.parameters_to_grid),
                    int(np.prod(np.shape(grids)) / len(self.parameters_to_grid)),
                ),
            )
        else:
            assert len({len(v) for v in self.parameters_to_grid.values()}) == 1
            self.flat_grids = np.array(
                [
                    [p[i] for p in self.parameters_to_grid.values()]
                    for i in range(len(next(iter(self.parameters_to_grid.values()))))
                ]
            ).T

        # Optional subset of model indices (0-based column in flat_grids) to run.
        # None means run all. Accepts any iterable; stored as a frozenset for O(1) lookup.
        self.model_ids = frozenset(model_ids) if model_ids is not None else None

        self.model_id_dict = {}
        self.models = []
        self.physics_values = None
        self.chemical_abun_values = None
        if not delay_run:
            self.run()

    def _grid_def(self, key: str, value: Any, model_count: int | None = None) -> None:
        if model_count is None:
            model_count = ""
        else:
            model_count = f"{str(model_count)}_"
        if isinstance(value, list) and key not in NoGridParameters:
            self.parameters_to_grid[model_count + key] = value
            self.parameters_to_grid[model_count + key] = np.array(value, dtype=object)
        elif isinstance(value, (np.ndarray, np.generic)) and key not in NoGridParameters:
            self.parameters_to_grid[model_count + key] = value.astype(dtype=object)

    def _log_main(self, msg: str) -> None:
        """Append a timestamped line to the main grid log file."""
        if self._main_log is None:
            return
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(self._main_log, "a") as f:
            f.write(f"{ts} {msg}\n")

    def run(self) -> None:
        """Run the grid."""
        signal.signal(signal.SIGINT, self._handler)
        n_total = np.shape(self.flat_grids)[1]
        self._log_main(f"Grid started: {n_total} models, {self.max_workers} workers")

        # Capture advanced settings so spawned workers start with the same
        # Fortran module state as the coordinator process.
        from uclchem.advanced.worker_state import _pool_initializer, create_snapshot

        snapshot = create_snapshot()

        pending = self.grid_iter(
            self.full_parameters,
            list(self.parameters_to_grid.keys()),
            self.flat_grids,
            self.model_type,
        )
        file_obj = h5py.File(self.grid_file, "a")

        PoolClass = NoDaemonPool if self.model_type == "SequentialRunner" else mp.Pool

        with PoolClass(
            self.max_workers,
            initializer=_pool_initializer,
            initargs=(snapshot,),
            maxtasksperchild=1,
        ) as pool:
            completed = 0

            def on_result(result: tuple[int, object]) -> None:
                nonlocal completed
                completed += 1
                model_id, model_object = result
                try:
                    save_name = f"{self.model_name_prefix}{model_id}"
                    model_object.un_pickle()
                    failed = (
                        model_object.success_flag is not None
                        and model_object.success_flag < 0
                    )
                    if failed:
                        from uclchem.utils import check_error as _check_error

                        msg = _check_error(
                            model_object.success_flag, raise_on_error=False
                        )
                        self._log_main(
                            f"model_{model_id} failed ({completed}/{n_total}): {msg}"
                        )
                    else:
                        self._log_main(
                            f"model_{model_id} completed ({completed}/{n_total})"
                        )
                    model_object.save_model(
                        file_obj=file_obj,
                        name=save_name,
                        overwrite=self.overwrite_models,
                    )
                    self.model_id_dict[model_id] = save_name
                    file_obj.flush()
                except Exception as e:
                    print(f"Error saving model {model_id}: {e}")
                    import traceback

                    traceback.print_exc()

            def on_error(_exc: Any, _model_id: int) -> None:
                self._log_main(f"model_{_model_id} error: {_exc}")
                print(f"error: {_exc}; for model: {_model_id}")

            # Submit all tasks upfront; the pool limits actual concurrency to
            # max_workers. pool.close()/join() are called only after all tasks
            # have been submitted, so apply_async never races with a closed pool.
            for pending_model in pending:
                model_id = pending_model.pop("id")
                if self.model_ids is not None and model_id not in self.model_ids:
                    continue
                self._log_main(f"model_{model_id} started")
                pool.apply_async(
                    _run_grid_model,
                    args=(model_id, self.model_type, pending_model, self.log_dir),
                    callback=on_result,
                    error_callback=lambda exc, _mid=model_id: on_error(exc, _mid),
                )
            pool.close()
            pool.join()
            file_obj.close()

        signal.signal(signal.SIGINT, self._orig_sigint)
        self._log_main(f"Grid finished: {len(self.model_id_dict)} models completed")
        self.models = [
            {"Model": v}
            for _, v in sorted(self.model_id_dict.items(), key=lambda item: item[1])
        ]
        self._load_params()

    def load_phys(self) -> None:
        """Load the physics.

        Raises:
            NotImplementedError: If the model type is `SequentialRunner`.

        """
        if self.model_type == "SequentialRunner":
            raise NotImplementedError("Sequential Runner physics loading not implemented")
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

    def load_chem(self, out_species_list: list[str]) -> None:
        """Load the chemistry.

        Args:
            out_species_list (list[str]): list of species to load abundances for.

        Raises:
            NotImplementedError: If the model type is `SequentialRunner`.
        """
        if self.model_type == "SequentialRunner":
            raise NotImplementedError(
                "Sequential Runner chemistry loading not implemented"
            )
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
            ].sel(chemical_abun_values=out_species_list)

    def _load_params(self) -> None:
        """Loop through the models present in the grid models self.models
        attribute in order to load the changing physical parameters of the models.
        The method splits the loops into two cases. SequentialRunner, and other model cases.
        In both instances, the for loop loads the model data using the _load_model_data()
        method. Then, it matches the given parameters, with the changing parameters in
        order to populate the dictionary attribute self.models for users to view the
        models run for the given GridRunner object. The Differentiation of
        SequentialRunner stems from SequentialRunner instances being nested models,
        resulting in an additional required loop to take into account the
        additional nesting.

        """
        # The following loops perform the same actions, but for the different model types
        if self.model_type == "SequentialRunner":
            for model in range(len(self.models)):
                for model_count in range(len(self.full_parameters)):
                    for mt_k, mt_v in self.full_parameters[model_count].items():
                        if not isinstance(mt_v, dict):
                            continue
                        tmp_model = self._load_model_data(
                            model=f"{self.models[model]['Model']}_{model_count}_{mt_k}"
                        )

                        self.models[model][f"{model_count}_{mt_k}"] = {
                            **{
                                k.replace(f"{model_count}_", ""): tmp_model._param_dict[
                                    k.replace(f"{model_count}_", "").lower()
                                ]
                                for k in list(self.parameters_to_grid.keys())
                                if k[: len(str(model_count))] == str(model_count)
                                and k.replace(f"{model_count}_", "").lower()
                                in tmp_model._param_dict
                            },
                            **{
                                k.replace(f"{model_count}_", ""): tmp_model.__getattr__(
                                    k.replace(f"{model_count}_", "")
                                )
                                for k in list(self.parameters_to_grid.keys())
                                if mt_k in k
                                and k.replace(mt_k, "").lower() in tmp_model._data
                            },
                        }
                        self.models[model][f"{model_count}_{mt_k}"]["Successful"] = (
                            True
                            if tmp_model.success_flag == 0
                            else tmp_model.success_flag
                        )
        else:
            for model in range(len(self.models)):
                loaded_data = self._load_model_data(model=self.models[model]["Model"])
                loaded_dict = loaded_data._param_dict
                self.models[model] = {
                    **self.models[model],
                    **{
                        k: loaded_dict[k.lower()]
                        for k in list(self.parameters_to_grid.keys())
                    },
                }
                self.models[model]["Successful"] = (
                    True if loaded_data.success_flag == 0 else loaded_data.success_flag
                )

    def _load_model_data(self, model: str):
        tmp_model = load_model(file=self.grid_file, name=model)
        return tmp_model

    def check_conservation(
        self, element_list: list[str] | None = None, percent: bool = True
    ) -> None:
        """Check conservation of the chemical abundances.

        Args:
            element_list (list[str] | None): List of elements to check conservation for.
                If None, use `uclchem.constants.default_elements_to_check`. Default = None.
            percent (bool): Flag on if percentage values should be used.
                Defaults to True.

        """
        if element_list is None:
            element_list = default_elements_to_check
        for model in range(len(self.models)):
            tmp_model = load_model(file=self.grid_file, name=self.models[model]["Model"])
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
                conserved = all(float(x[:1]) < 1 for x in i.values())
            self.models[model]["elements_conserved"] = conserved

    def _handler(self, signum: Any, frame: Any) -> None:  # noqa: ARG002
        try:
            self.on_interrupt()  # your “final steps”
        finally:
            # Restore default and re-raise KeyboardInterrupt to stop execution
            signal.signal(signal.SIGINT, self._orig_sigint)
            raise KeyboardInterrupt

    @typing.overload
    def on_interrupt(self) -> None:
        return

    @staticmethod
    def grid_iter(
        full_parameters: dict | list,
        param_keys: list,
        flattened_grids: np.ndarray,
        model_type: str,
    ) -> Iterator[dict[str, Any]]:
        """Provide an iterable dictionary of parameters that can be used with the
        grid-based multiprocessing worker distribution.

        Args:
            full_parameters (dict | list): dictionary or list
                (if model_type == SequentialRunner)
                of the full parameters that will be used for the model
            param_keys (list[str]): list of parameters that are changing in
                this GridRunner object
            flattened_grids (np.ndarray): list of all the values for the changing parameters
                for this GridRunner object
            model_type (str): Type of model to use. 'SequentialRunner' results in
                alternative way of executing this phase as each SequentialRunner
                represents multiple models to be run in series.

        Yields:
            Next dictionary containing the parameter values for the model to run
                in the grid of models. Only offers one model per request.

        """
        if model_type == "SequentialRunner":
            # As SequentialRunner types contain multiple models,
            # sometimes of the same type, we split this type out to follow altered
            # logic to arrive at equivalently expected outputs.
            for i in range(len(flattened_grids[0])):
                combo = ()
                for j in range(np.shape(flattened_grids)[0]):
                    combo += (flattened_grids[j][i],)
                yield_dict = {"id": i}
                # run_list contains the dictionaries of each model to be run
                # in the sequence of models. This is the SequentialRunner
                # sequenced_model_parameters input.
                run_list = []
                for model_count in range(len(full_parameters)):
                    # run_dict contains all input parameters for a model,
                    # not just param_dict, for an individual model that is
                    # part of the SequentialRunner.
                    run_dict = {}
                    for model_type, model_full_parameters in full_parameters[
                        model_count
                    ].items():
                        if isinstance(model_full_parameters, dict):
                            # grid_param_dict is filled with the param_dict values of a model.
                            grid_param_dict = {
                                k.replace(f"{model_count}_", ""): v
                                for k, v in zip(param_keys, combo)
                                if k.replace(f"{model_count}_", "")
                                in model_full_parameters["param_dict"]
                                and k[: len(str(model_count))] == str(model_count)
                            }
                            # grid_dict is filled with the input parameters of a value,
                            # not part of param_dict
                            grid_dict = {
                                k.replace(f"{model_count}_", ""): v
                                for k, v in zip(param_keys, combo)
                                if k.replace(f"{model_count}_", "")
                                not in model_full_parameters["param_dict"]
                                and (
                                    k[: len(str(model_count))] == str(model_count)
                                    if k[: len(str(model_count))].isdigit()
                                    else False
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
                            run_list += [run_dict]
                        else:
                            yield_dict[model_type] = model_full_parameters
                yield {
                    "parameters_to_match": ["finalDens"],
                    **yield_dict,
                    **{"sequenced_model_parameters": run_list},
                }
        else:
            for i in range(len(flattened_grids[0])):
                combo = ()
                for j in range(np.shape(flattened_grids)[0]):
                    combo += (flattened_grids[j][i],)

                # grid_param_dict contains the param_dict values of the next model to run.
                grid_param_dict = {
                    k: v if not isinstance(v, float) else v.item()
                    for k, v in zip(param_keys, combo)
                    if k in full_parameters["param_dict"]
                }
                # grid_dict is filled with the input parameters of a value,
                # not part of param_dict, for the next model to run
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
