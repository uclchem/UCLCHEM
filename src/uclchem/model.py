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
    >>> import uclchem
    >>>
    >>> # Create a collapsing cloud model
    >>> cloud = uclchem.model.Cloud(
    ...     param_dict={
    ...         "initialDens": 1e2,
    ...         "initialTemp": 10.0,
    ...         "finalTime": 1e6,
    ...         "freefall": True
    ...     },
    ...     out_species=["CO", "H2O", "CH3OH"]
    ... )
    >>>
    >>> # Check for errors and plot
    >>> cloud.check_error()
    Model ran successfully
    >>> cloud.create_abundance_plot(["CO", "$CO"]) #doctest: +SKIP

**Model Workflow:**

1. **Initialize**: Create model object with parameters
2. **Run**: Model runs automatically on initialization (or use ``read_file`` to load)
3. **Analyze**: Access results via attributes (``.final_abundances``, ``.chemistry_dataframe``)
4. **Plot**: Use built-in plotting methods (``.create_abundance_plot()``)
5. **Chain**: Use as input to next stage (``.previous_model`` parameter)

**Common Parameters:**

All models accept these key parameters in ``param_dict``:

- ``initialDens`` (float): Initial density [cm⁻³]
- ``initialTemp`` (float): Initial temperature [K]
- ``finalTime`` (float): Simulation end time [years]
- ``freefall`` (bool): Enable freefall collapse (Cloud only)
- ``outputFile`` (str): Output file path (optional with OO API)

See the user guide for complete parameter list.

**Species Naming:**

- Gas phase: ``CO``, ``H2O``, ``CH3OH``
- Ice surface: ``#CO``, ``#H2O``, ``#CH3OH``
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
import traceback
import warnings
from abc import ABC, abstractmethod
from datetime import datetime
from multiprocessing import pool, shared_memory
from pathlib import Path
from time import perf_counter
from typing import TYPE_CHECKING, Any, Literal

import h5py
import numpy as np
import pandas as pd

# UCLCHEM related imports
import uclchemwrap
import xarray as xr
from uclchemwrap import uclchemwrap as wrap

from uclchem._coolant_utils import load_coolant_level_names
from uclchem._fortran_capture import capture_fortran_output
from uclchem.advanced.worker_state import (
    _pool_initializer,
    create_snapshot,
    restore_snapshot,
)
from uclchem.analysis import check_element_conservation
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
from uclchem.utils import (
    UCLCHEM_ROOT_DIR,
    ArrayLike,
    SuccessFlag,
    convert_keys_to_lowercase,
    get_dtype,
)

if TYPE_CHECKING:
    from collections.abc import Iterator

    import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# Global variables determining formats of write files
PHYSICAL_PARAMETERS_HEADER_FORMAT = "%10s"
# in the below variable, the outputs were chosen according to the spacing needed for
# "      Time,    Density,    gasTemp,   dustTemp,         Av,   radfield,       zeta,
#       point,    parcel_radius"
PHYSICAL_PARAMETERS_VALUE_FORMAT = (
    "%10.3E, %10.4E, %10.2f, %10.2f, %10.4E, %10.4E, %10.4E, %10i, %10.4E"
)
SPECNAME_HEADER_FORMAT = "%11s"
SPECNAME_VALUE_FORMAT = "%9.5E"


# Model registration is intended to prevent code injection during loading time.
REGISTRY: dict[str, type[AbstractModel]] = {}


def register_model(
    cls: type[AbstractModel],
) -> type[AbstractModel]:
    """Register a new model in the model registry.

    Parameters
    ----------
    cls : type[AbstractModel]
        class to register.

    Returns
    -------
    cls : type[AbstractModel]
        class

    Raises
    ------
    ValueError
        If a model with the same name as cls is already in the registry.

    """
    name = getattr(cls, "MODEL_NAME", cls.__name__)
    if name in REGISTRY and REGISTRY[name] is not cls:
        msg = f"Duplicate model registration for {name}"
        raise ValueError(msg)
    REGISTRY[name] = cls
    return cls


# /Global variables determining formats of write files


# Reaction and Species name retrieval classes to reduce file read repetition.
def reaction_line_formatter(line: list[str] | pd.Series) -> str:
    """Format a list of strings as a reaction, while filtering out "NAN"s.

    Parameters
    ----------
    line : list[str] | pd.Series
        list of species involved in the reaction

    Returns
    -------
    str
        formatted reaction for printing.

    Examples
    --------
    >>> print(reaction_line_formatter(["#OH", "#H", "LH", "#H2O", "NAN", "NAN", "NAN"]))
    #OH + #H + LH -> #H2O
    >>> print(reaction_line_formatter(["H2", "PHOTON", "NAN", "H", "H", "NAN", "NAN"]))
    H2 + PHOTON -> H + H

    """
    reactants = list(filter(lambda x: not str(x).lower().endswith("nan"), line[0:3]))
    products = list(filter(lambda x: not str(x).lower().endswith("nan"), line[3:7]))
    return " + ".join(reactants) + " -> " + " + ".join(products)


class ReactionNamesStore:
    """Way to only read reaction file once, and keep them loaded after."""

    def __init__(self):
        """Initialize the ReactionNamesStore."""
        self.reaction_names = None

    def __call__(self) -> list[str]:
        """Get a list of formatted reactions.

        Only load the reactions once, after that use the cached version.

        Returns
        -------
        list[str]
            List of formatted reactions.

        """
        if self.reaction_names is None:
            logger.debug(f"Reading reaction file {UCLCHEM_ROOT_DIR / 'reactions.csv'}")
            reactions = pd.read_csv(UCLCHEM_ROOT_DIR / "reactions.csv")
            # format the reactions:
            self.reaction_names = [
                reaction_line_formatter(line) for idx, line in reactions.iterrows()
            ]
        return self.reaction_names


get_reaction_names = ReactionNamesStore()


class SpeciesNamesStore:
    """Way to only read species file once, and keep them loaded after."""

    def __init__(self):
        """Initialize the SpeciesNamesStore."""
        self.species_names = None

    def __call__(self) -> list[str]:
        """Get the species names.

        Only loads the species once, after that use the cached version

        Returns
        -------
        list[str]
            List of species names

        """
        if self.species_names is None:
            logger.debug(f"Reading file {UCLCHEM_ROOT_DIR / 'species.csv'}")
            species = pd.read_csv(UCLCHEM_ROOT_DIR / "species.csv")
            self.species_names = species["NAME"].tolist()
        return self.species_names


get_species_names = SpeciesNamesStore()
# /Reaction and Species name retrieval classes to reduce file read repetition.


# Universal model loader
def load_model(
    file: h5py.File | str | Path,
    name: str = "default",
) -> AbstractModel:
    """Load a pre-existing model from a file. Bypasses ``__init__``.

    Parameters
    ----------
    file : h5py.File | str | Path
        open h5py file object, or Path to a file that contains
        previously run and stored models.
    name : str
        Name of the stored object. Defaults = "default".

    Returns
    -------
    obj : AbstractModel
        Loaded object that inherits from AbstractModel and has the class
        of to the model found in the loaded file.

    Raises
    ------
    Exception
        If the model with name ``name`` is not found in the file.
    TypeError
        if ``file`` is not a string, Path or ``h5py.File`` instance.
    ValueError
        If the model type is not in the model registry.

    """
    if isinstance(file, str | Path):
        opened_file = True
        file_obj = h5py.File(file, "r")
    elif isinstance(file, h5py.File):
        opened_file = False
        file_obj = file
    else:
        msg = f"Expected file to be type h5py.File, string or Path, but got type {type(file)}"
        raise TypeError(msg)

    if name not in file_obj:
        msg = f"model {name} was not found in the save file that was passed."
        raise Exception(msg)

    model_group = file_obj[name]
    coords = {}
    if "_coords" in model_group:
        for coord_name in model_group["_coords"]:
            coords[coord_name] = _read_array(model_group["_coords"], coord_name)
    data_vars = {}
    for coord_name in model_group:
        if coord_name == "_coords":
            continue
        data_vars[coord_name] = _read_array(model_group, coord_name)
    loaded_data = xr.Dataset(data_vars, coords=coords)

    if opened_file:
        file_obj.close()

    model_class = json.loads(loaded_data["attributes_dict"].item())["model_type"]
    cls = REGISTRY.get(model_class)
    if cls is None:
        msg = f"Unrecognized model type '{model_class}'. Not in trusted registry."
        raise ValueError(msg)
    return cls.load_from_dataset(model_ds=loaded_data)


def _read_array(model_group: dict[str, xr.Dataset], name: str) -> xr.Variable:
    """Read an array from a model group.

    Parameters
    ----------
    model_group : dict[str, xr.Dataset]
        model group.
    name : str
        key in model_group

    Returns
    -------
    xr.Variable
        xr array

    """
    ds = model_group[name]
    data = ds[()]
    if data.dtype.kind == "S":
        data = data.astype(str)
    dims = list(ds.attrs["_dims"])
    attrs = json.loads(ds.attrs.get("_attrs", "{}"))
    return xr.Variable(dims, data, attrs=attrs)


# /Universal model loader


def _write_array(
    model_group: h5py.Group,
    name: str,
    xr_var: xr.DataArray,
    array_dtype: np.typing.DTypeLike | np.dtype | str | None = None,
) -> None:
    data = xr_var.values
    if data.dtype.kind == "U":
        data = data.astype(bytes)

    dataset_kwargs = {}
    if array_dtype is not None and "_array" in name:
        dtype = get_dtype(array_dtype)
        dataset_kwargs["dtype"] = dtype

    ds = model_group.create_dataset(name, data=data, **dataset_kwargs)
    ds.attrs["_dims"] = list(xr_var.dims)


def _create_shared_memory_allocation(
    shape: tuple[int, ...],
) -> tuple[shared_memory.SharedMemory, dict, np.ndarray]:
    """Create a shared memory object.

    Parameters
    ----------
    shape : tuple[int, ...]
        shape of new array

    Returns
    -------
    shm : shared_memory.SharedMemory
        SharedMemory object.
    spec : dict
        Description of shared memory objects name and shape.
    array : np.ndarray
        Created numpy array (filled with zeros), with dtype ``np.float64``.

    """
    logger.debug(f"Creating shared memory allocation with shape {shape}")
    nbytes = int(np.prod(shape) * np.dtype(np.float64).itemsize)
    shm = shared_memory.SharedMemory(create=True, size=nbytes)
    array = np.ndarray(shape, dtype=np.float64, buffer=shm.buf, order="F")
    array.fill(0.0)
    spec = {"name": shm.name, "shape": shape}
    return shm, spec, array


# Worker entry for parallel jobs
def _worker_entry(
    model_class: str,
    init_kwargs: dict,
    shm_descs: dict,
    result_queue: mp.Queue,
    advanced_snapshot: dict | None = None,
) -> None:
    """Start the model.

    Parameters
    ----------
    model_class : str
        name of the model to run
    init_kwargs : dict
        keyword arguments to initialize the model
    shm_descs : dict
        dictionary with shared memory handles.
    result_queue : mp.Queue
        result queue where output of ``AbstractModel.run_fortran()``
        will be put.
    advanced_snapshot : dict | None
        snapshot to restore using func:`restore_snapshot`.
        If None, do not restore a snapshot. Default = None.

    Raises
    ------
    ValueError
        If ``model_class`` is not in the registry of models.

    """
    # Restore advanced settings captured in the coordinator process.
    if advanced_snapshot is not None:
        restore_snapshot(advanced_snapshot)

    cls = REGISTRY.get(model_class)
    if cls is None:
        msg = f"Unrecognized model type '{model_class}'. Not in trusted registry."
        raise ValueError(msg)
    model = cls.worker_build(init_kwargs=init_kwargs, shm_desc=shm_descs)
    with capture_fortran_output(label="last_model", log_file="./last_model_fortran.log"):
        output = model.run_fortran()
    result_queue.put(output)


# /Worker entry for parallel jobs


# Short compatibility helper for legacy parameter ``endAtFinalDensity``
def _convert_legacy_stopping_param(
    param_dict: dict[str, Any] | None,
) -> dict[str, Any] | None:
    """Minimal conversion of legacy ``endAtFinalDensity`` to ``parcelStoppingMode``.

    Parameters
    ----------
    param_dict : dict[str, Any] | None
        parameter dictionary.

    Returns
    -------
    dict[str, Any] | None
        Converted dictionary, or None if ``param_dict`` is None.

    Raises
    ------
    RuntimeError
        If ``endAtFinaldensity`` and ``parcelStoppingMode`` are both
        in ``param_dict``.
    RuntimeError
        If ``endAtFinalDensity`` is being used with a multi-point model.

    Notes
    -----
    This function assumes param_dict is already a copy and is case-normalized (lowercase keys).

    """
    if param_dict is None:
        return param_dict

    has_old = "endatfinaldensity" in param_dict
    has_new = "parcelstoppingmode" in param_dict
    if has_old and has_new:
        msg = "Cannot specify both 'endAtFinalDensity' and 'parcelStoppingMode'. Use 'parcelStoppingMode' only."
        raise RuntimeError(msg)
    if has_old:
        points = param_dict.get("points", 1)
        if points > 1:
            msg = "endAtFinalDensity is no longer supported for multi-point models (points > 1). Use 'parcelStoppingMode' instead."
            raise RuntimeError(msg)
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

    """

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species_list: list[str] | None = None,
        starting_chemistry: ArrayLike | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Create a new AbstractModel.

        Parameters
        ----------
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values if not provided.
        out_species_list : list[str] | None
            List of species to focus on for outputs.
            If None, defaults to ``uclchem.constants.default_elements_to_check``.
        starting_chemistry : ArrayLike | None
            Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model that will
            be run. Defaults to None.
        timepoints : int
            Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object
            prevents it from being run. Defaults to None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If ``run_type`` is not one of ``["managed", "external"]``.

        """
        if out_species_list is None:
            out_species_list = default_elements_to_check
        if run_type not in {"managed", "external"}:
            msg = (
                f"run_type should be one of ['managed', 'external'], but got '{run_type}'"
            )
            raise ValueError(msg)

        self._data = xr.Dataset()
        self._pickle_dict: dict[str, Any] = {}
        # Per-instance metadata containers (scalars and small values)
        object.__setattr__(self, "_meta", {})
        object.__setattr__(self, "_pickle_meta", {})
        # Set run_type into metadata
        self.run_type = run_type
        # Shared memory
        self._shm_desc: dict[str, dict] = {}
        self._shm_handles: dict[str, shared_memory.SharedMemory] = {}
        self._proc_handle: mp.context.SpawnProcess | None = None
        # /Shared memory

        # Signal Interrupt
        self._was_interrupted = False
        self._orig_sigint = signal.getsignal(signal.SIGINT)
        # /Signal Interrupt

        self.model_type = str(self.__class__.__name__)
        self._param_dict: dict[str, Any] = {}
        self.out_species_list = out_species_list
        self.out_species = ""
        self.success_flag: None | SuccessFlag = None
        # Note: specname is now accessed via get_species_names() global function
        # Note: PHYSICAL_PARAMETERS is now accessed via the global constant

        self.n_out = 0 if read_file is None else None
        self.timepoints = timepoints
        self.was_read = read_file is not None

        self._reform_inputs(param_dict, self.out_species_list)

        # Validate parcelStoppingMode usage after model_type is known
        if self._param_dict.get("_needs_freefall_validation", False):
            if self.model_type != "Collapse":  # Collapse models imply freefall
                msg = (
                    "parcelStoppingMode != 0 can only be used with:\n"
                    "  - Cloud models with freefall=True\n"
                    "  - Collapse models (freefall is implied)\n"
                    f"Current model_type={self.model_type}, freefall={self._param_dict.get('freefall', False)}\n"
                    "Please either set freefall=True or set parcelStoppingMode=0"
                )
                raise ValueError(msg)
            self._param_dict.pop("_needs_freefall_validation", None)  # Clean up flag

        # If we were given a previously-written output file, populate the model
        # arrays and metadata from that file now so later initialization can rely on them.
        if read_file is not None:
            self.legacy_read_output_file(read_file)
        if "points" not in self._param_dict:
            self._param_dict["points"] = 1
        # Expose grid points as attribute for legacy code expecting ``gridPoints``
        object.__setattr__(self, "gridPoints", self._param_dict["points"])

        self.outputFile: Path | None = (
            Path(self._param_dict.pop("outputfile"))
            if "outputfile" in self._param_dict
            else None
        )
        self.abundSaveFile: Path | None = (
            Path(self._param_dict.pop("abundsavefile"))
            if "abundsavefile" in self._param_dict
            else None
        )
        self.abundLoadFile: Path | None = (
            Path(self._param_dict.pop("abundloadfile"))
            if "abundloadfile" in self._param_dict
            else None
        )

        self.starting_chemistry_array: np.ndarray | None = None
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
        if np.all(self.starting_chemistry_array == 0.0):
            msg = "Detected all zeros in the starting chemistry array."
            raise ValueError(msg)

        # Only initialize next_starting_chemistry_array if we didn't load it from a file
        # (legacy_read_output_file sets it from the last timestep)
        self.next_starting_chemistry_array: np.ndarray | None
        if read_file is None:
            self.next_starting_chemistry_array = None

        # Only create new arrays if we didn't load them from a file
        if read_file is None:
            self.physics_array: np.ndarray | None = None
            self.chemical_abun_array: np.ndarray | None = None
            self._create_fortran_array()
            self.rate_constants_array: np.ndarray | None = None
            self._create_rate_constants_array()
            self.heat_array: np.ndarray | None = None
            self._create_heating_array()
            self.stats_array: np.ndarray | None = None
            self._create_stats_array()
            self.level_populations_array: np.ndarray | None = None
            self._create_level_populations_array()
            self.se_stats_array: np.ndarray | None = None
            self._create_se_stats_array()
            self.out_species_abundances_array: np.ndarray | None = None
        else:
            # When loading from file, arrays are already populated; just initialize
            # the arrays that weren't loaded
            self.rate_constants_array = None
            self.heat_array = None
            self.stats_array = None
            self.level_populations_array = None
            self.se_stats_array = None
            self.out_species_abundances_array = None

    def __del__(self):
        """Unlink all shared memory objects.

        If the AbstractModel object goes out of scope, ensure that no shared memory objects stick around.

        """
        if hasattr(self, "_shm_desc") and bool(self._shm_desc):
            self._coordinator_unlink_memory()

    # Separate class building method(s)
    @classmethod
    def load_from_dataset(cls, model_ds: xr.Dataset) -> AbstractModel:
        """Load an abstract model from an xr Dataset.

        Parameters
        ----------
        model_ds : xr.Dataset
            Dataset to load

        Returns
        -------
        obj : AbstractModel
            instantiated model

        """
        obj = cls.__new__(cls)
        obj._param_dict = json.loads(model_ds["_param_dict"].item())
        del model_ds["_param_dict"]
        obj._data = xr.Dataset()
        obj._data = model_ds.copy()
        model_ds.close()
        temp_attribute_dict = json.loads(obj._data["attributes_dict"].item())
        # Restore these values into the metadata dict rather than dataset variables
        object.__setattr__(obj, "_meta", temp_attribute_dict)  # noqa: PLC2801
        del obj._data["attributes_dict"]
        obj._coord_assign()
        if obj.success_flag is not None:
            obj.success_flag = SuccessFlag(obj.success_flag)
        return obj

    @classmethod
    def worker_build(cls, init_kwargs, shm_desc) -> AbstractModel:  # noqa: ANN001, D102
        obj = cls.__new__(cls)
        for k, v in init_kwargs.items():
            object.__setattr__(obj, k, v)  # noqa: PLC2801
        obj._reform_array_in_worker(shm_desc)
        return obj

    @classmethod
    def from_file(
        cls,
        file: str | Path | h5py.File,
        name: str = "default",
    ) -> AbstractModel:
        """Load a model from a file.

        This is a convenience class method that wraps the module-level load_model function.

        Parameters
        ----------
        file : str | Path | h5py.File
            Path to a file that contains previously run
            and stored models, or open ``h5py.File`` object.
        name : str
            Name of the stored object. Default = 'default'.

        Returns
        -------
        AbstractModel
            Model object loaded from the file.

        """
        return load_model(file, name=name)

    # /Separate class building methods

    # Class utility methods
    def __getattr__(self, key: str) -> Any:
        """Get attribute ``key``.

        Searches both ``self._meta`` and ``self._data``. If ``key`` starts with ``"_"``,
        just return the attribute.

        Parameters
        ----------
        key : str
            name of attribute

        Returns
        -------
        Any
            Attribute value

        Raises
        ------
        AttributeError
            If no attribute ``key`` can be found in ``self._meta``,
            ``self._data``.

        """
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

        msg = f'{self.__class__.__name__} has no attribute of name: "{key}".'
        raise AttributeError(msg)

    def __setattr__(self, key: str, value: Any) -> None:
        """Set attribute ``key`` to ``value``.

        Parameters
        ----------
        key : str
            attribute to set
        value : Any
            value to set attribute ``key`` to

        """
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
            with contextlib.suppress(Exception):
                meta = super().__getattribute__("_meta")
                if key in meta:
                    del meta[key]

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

    def has_attr(self, key: str) -> bool:
        """Check if the object has an attribute stored in self._meta or self._data.

        Parameters
        ----------
        key : str
            name of attribute

        Returns
        -------
        bool
            whether the object has the attribute.

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
        """Check the conservation of the elemental abundances.

        Parameters
        ----------
        element_list : list[str] | None
            List of elements to check conservation for.
            If None, uses ``self.out_species_lists``. Default = None.
        percent : bool
            Flag on if whether changes should be printed in percentages.
            Defaults to True.

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
                        self.get_joined_dataframes(i),
                        element_list=element_list,
                        percent=percent,
                    )
                )
        else:
            print("Element conservation report")
            print(
                check_element_conservation(
                    self.get_joined_dataframes(0),
                    element_list=element_list,
                    percent=percent,
                )
            )

    def check_error(self, only_error: bool = False, raise_on_error: bool = True) -> None:
        """Check the model error status and raises RuntimeError on failure.

        Parameters
        ----------
        only_error : bool
            If True, only act when there was an error (skip success message).
            Default = False.
        raise_on_error : bool
            If True (default), raises RuntimeError on failure. If False, prints.
            Default = True.

        """
        if self.success_flag is None:
            print("Model has not been run.")
            return
        msg = self.success_flag.check_error(
            only_error=only_error, raise_on_error=raise_on_error
        )
        if msg is not None:
            print(msg)

    def create_abundance_plot(
        self,
        species: list[str] | None = None,
        figsize: tuple[float, float] = (16, 9),
        point: int = 0,
        plot_file: str | Path | None = None,
    ) -> tuple[plt.Figure, plt.Axes]:
        """``uclchem.plot.create_abundance_plot`` wrapper method.

        Parameters
        ----------
        species : list[str] | None
            List of species to plot. If None, uses self.out_species_list.
            Default = None.
        figsize : tuple[float, float]
            The figure size to use for matplotlib.
            Default = (16, 9).
        point : int
            Integer referring to which point of the UCLCHEM model to use.
            Default = 0.
        plot_file : str | Path | None
            if not None, save to a path.
            Default = None.

        Returns
        -------
        tuple[plt.Figure, plt.Axes]
            matplotlib figure and axis objects

        Raises
        ------
        ValueError
            If ``point`` is larger than the number of points in the model run.

        """
        if species is None:
            species = self.out_species_list

        if point > self._param_dict["points"]:
            msg = "'point' must be less than number of modeled points."
            raise ValueError(msg)
        return create_abundance_plot(
            self.get_joined_dataframes(point),
            species,
            figsize,
            plot_file,
        )

    def get_dataframes(
        self,
        point: int | None = None,
        with_rate_constants: bool = False,
        with_heating: bool = False,
        with_stats: bool = False,
        with_level_populations: bool = False,
        with_se_stats: bool = False,
    ) -> tuple[pd.DataFrame, ...]:
        """Get the model physics and chemical abundances as multiple separate pandas DataFrames.

        Parameters
        ----------
        point : int | None
            Integer referring to which point of the UCLCHEM model to return.
            If None, returns data for all points with a 'Point' column. Defaults to None.
        with_rate_constants : bool
            Whether to include reaction rate constants
            as a separate dataframe. Default = False.
        with_heating : bool
            Flag on whether to include heating/cooling rates
            as a separate dataframe. Default = False.
        with_stats : bool
            Whether to include DVODE solver statistics
            as a separate dataframe. Default = False.
        with_level_populations : bool
            Whether to include coolant level populations
            in the output. Default = False.
        with_se_stats : bool
            Whether to include SE solver statistics in the output.
            Default = False.

        Returns
        -------
        physics_df : pd.DataFrame
            Dataframe of the physical parameters for point ``point``
        chemistry_df : pd.DataFrame
            Dataframe of the chemical abundances  for point ``point``
        rate_constants_df : pd.DataFrame
            Dataframe of the reaction rate constants for point
            ``point`` if with_rate_constants = True
        heating_df : pd.DataFrame
            Dataframe of the heating/cooling rates  for point ``point``
            if with_heating = True
        stats_df : pd.DataFrame
            Dataframe of DVODE solver statistics for point ``point``
            if with_stats = True
        level_populations_df : pd.DataFrame
            Dataframe of coolant level populations for point
            ``point`` if with_level_populations = True
        se_stats_df : pd.DataFrame
            Dataframe of SE solver statistics for point ``point``
            with_se_stats = True

        Notes
        -----
        Returns a tuple of at least 2 pandas DataFrames:
            (``physics_df`` and ``chemistry_df``), or 3 or more, depending on how many
            ``with_...`` are True.

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
                dfs = tuple(add_point_column(df, pt + 1) for df in dfs)  # type: ignore[assignment]
                all_dfs.append(dfs)

            # Transpose to group by dataframe type instead of by point
            # e.g., [[phys0, chem0, rates0], [phys1, chem1, rates1]] -> [[phys0, phys1], [chem0, chem1], [rates0, rates1]] # noqa: W505
            df_collections = list(zip(*all_dfs, strict=True))

            # Concatenate each type vertically
            concatenated: tuple[pd.DataFrame, ...] = tuple(
                pd.concat([df for df in collection if df is not None], ignore_index=True)
                for collection in df_collections
                if any(df is not None for df in collection)
            )
        else:
            # Single point mode
            concatenated = tuple(
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
            concatenated = tuple(add_point_column(df, point + 1) for df in concatenated)

        return concatenated

    def get_joined_dataframes(
        self,
        point: int | None = None,
        with_rate_constants: bool = False,
        with_heating: bool = False,
        with_stats: bool = False,
        with_level_populations: bool = False,
        with_se_stats: bool = False,
    ) -> pd.DataFrame:
        """Get the model physics and chemical abundances as one singular pandas DataFrame.

        Parameters
        ----------
        point : int | None
            Integer referring to which point of the UCLCHEM model to return.
            If None, returns data for all points with a 'Point' column. Default = None.
        with_rate_constants : bool
            Whether to include reaction rate constants
            in the dataframe. Default = False.
        with_heating : bool
            Whether to include heating/cooling rates in the dataframe.
            Default = False.
        with_stats : bool
            Whether to include DVODE solver statistics in the dataframe.
            Default = False.
        with_level_populations : bool
            Whether to include coolant level populations
            in the dataframe. Default = False.
        with_se_stats : bool
            Flag on whether to include SE solver statistics in the dataframe.
            Default = False.

        Returns
        -------
        result_df : pd.DataFrame
            Dataframe of the joined arrays for point ``point``

        """
        result_dfs = self.get_dataframes(
            point=point,
            with_rate_constants=with_rate_constants,
            with_heating=with_heating,
            with_stats=with_stats,
            with_level_populations=with_level_populations,
            with_se_stats=with_se_stats,
        )

        result_df = result_dfs[0]
        for df in result_dfs[1:]:
            if df is not None:
                # Drop duplicate Point column from subsequent dataframes
                result_df = result_df.join(df.drop(columns=["Point"]))
        return result_df

    def _get_single_point_dataframes(
        self,
        point: int,
        with_rate_constants: bool,
        with_heating: bool,
        with_stats: bool,
        with_level_populations: bool,
        with_se_stats: bool,
    ) -> tuple[pd.DataFrame, ...]:
        """Get the dataframes for a single point without Point column.

        Parameters
        ----------
        point : int
            Spatial point index (for multi-point models).
        with_rate_constants : bool
            Flag on whether to include a reaction rate constant dataframe
            in the tuple.
        with_heating : bool
            Flag on whether to include heating/cooling rates dataframe in the tuple.
        with_stats : bool
            Flag on whether to include DVODE solver statistics dataframe in the tuple.
        with_level_populations : bool
            Flag on whether to include coolant level
            populations in the tuple.
        with_se_stats : bool
            Flag on whether to include SE solver statistics
            in the tuple

        Returns
        -------
        tuple[pd.DataFrame, ...]
            a tuple of pd.DataFrame with physics_df, chemistry_df, and all
            additional information based off whether the flags were True.

        Raises
        ------
        ValueError
            If ``point`` is larger than or equal to the number of points
            that were run in the model.
        ValueError
            If ``self.physics_array`` or ``self.chemical_abun_array`` is None.

        """
        # Determine total number of points in model
        n_points = self._param_dict.get("points", 1)

        if point >= n_points:
            msg = f"point {point} was larger than the number of points in the model {n_points}"
            raise ValueError(msg)

        # Create a physical parameter dataframe using global constants
        # Arrays are guaranteed to match these dimensions due to validation in legacy_read_output_file
        if self.physics_array is None:
            msg = "physics_array is None, so cannot create the single point dataframes."
            raise ValueError(msg)
        physics_df = pd.DataFrame(
            self.physics_array[:, point, :],
            index=None,
            columns=PHYSICAL_PARAMETERS,
        )

        if self.chemical_abun_array is None:
            msg = "chemical_abun_array is None, so cannot create the single point dataframes."
            raise ValueError(msg)
        # Create an abundances dataframe using global species names
        species_names = get_species_names()
        chemistry_df = pd.DataFrame(
            self.chemical_abun_array[:, point, :],
            index=None,
            columns=species_names,
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

            heating_columns.extend(
                [f"{label} Line Cooling" for label in line_cooling_labels]
            )

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
        result: list[pd.DataFrame | None] = [physics_df, chemistry_df]
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
        return tuple(result)  # type: ignore[arg-type]

    def get_solver_stats_dataframe(self, point: int | None = None) -> pd.DataFrame | None:
        """Get all solver statistics including failed attempts.

        This method returns statistics for EVERY DVODE solver call,
        including failed attempts that were retried. This is different
        from the regular stats in get_dataframes() which only shows
        the final successful attempt per trajectory timestep.

        Parameters
        ----------
        point : int | None
            Spatial point index (for multi-point models).
            If None, uses point 0. Default = None.

        Returns
        -------
        pd.DataFrame | None
            DataFrame with columns from DVODE_STAT_NAMES,
            or None if stats not available.
            TRAJECTORY_INDEX column links solver attempts to trajectory timesteps.
            Rows where TRAJECTORY_INDEX=0 are filtered out (unused preallocated space).

        Examples
        --------
        >>> import uclchem
        >>> param_dict = {}
        >>> model = uclchem.model.Cloud(param_dict)
        >>> solver_stats = model.get_solver_stats_dataframe()
        >>> # Count failed attempts
        >>> failures = solver_stats[solver_stats['ISTATE'] < 0]
        >>> print(f"Failed attempts: {len(failures)}") # doctest: +ELLIPSIS
        Failed attempts: ...

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

        Parameters
        ----------
        point : int | None
            Spatial point index (for multi-point models).
            If None, uses point 0. Default = None.

        Returns
        -------
        pd.DataFrame | None
            DataFrame of failed attempts,
            or None if no failures or stats unavailable.

        Examples
        --------
        >>> import uclchem
        >>> param_dict = {}
        >>> model = uclchem.model.Cloud(param_dict)
        >>> failures = model.get_failed_solver_attempts()
        >>> if failures is not None:
        ...     print(f"Total retries needed: {len(failures)}")
        ...     print(failures.groupby('ISTATE').size())
        ... else:
        ...     print("No failures occurred.")
        ...
        No failures occurred.

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

        Parameters
        ----------
        point : int | None
            Spatial point index (for multi-point models).
            If None, uses point 0. Default = None.

        Returns
        -------
        dict[str, int | float] | None
            dict with keys:
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
        plot_kwargs: dict[str, Any] | None = None,
    ) -> plt.Axes:
        """``uclchem.plot.plot_species`` wrapper method.

        Parameters
        ----------
        ax : plt.Axes
            An axis object to plot on
        species : list[str] | None
            A list of species names to be plotted.
            If species name starts with "$" instead of "#" or "@",
            plots the sum of surface and bulk abundances. If None, default to
            ``self.out_species_list``. Default = None.
        point : int
            Grid point index. Default = 0.
        legend : bool
            Whether to add a legend to the plot. Default = True.
        plot_kwargs : dict[str, Any] | None
            keyword arguments passed to ``ax.plot``.
            Default = None.

        Returns
        -------
        plt.Axes
            Modified input axis is returned

        """
        if species is None:
            species = self.out_species_list
        return plot_species(
            ax,
            self.get_joined_dataframes(point),
            species,
            legend=legend,
            plot_kwargs=plot_kwargs,
        )

    # /UCLCHEM utility and analysis wrappers

    # Methods to start run of model
    def run(self) -> None:
        """Run the model.

        Reset the Fortran arrays if the model was not read, allowing the arrays to be reused
        for new runs.

        Raises
        ------
        RuntimeError
            If the model was read.
        RuntimeError
            If the dictionary returned by ``self.run_fortran()``
            does not contain a key ``"success_flag"``.
        TypeError
            If the dictionary returned by ``self.run_fortran()`` does not have
            a valid SuccessFlag type (integer)

        """
        if self.was_read:
            msg = "This model was read. It can not be run. "
            raise RuntimeError(msg)

        def _handler(signum: Any, frame: Any) -> None:  # noqa: ARG001
            """Handle a raised exception.

            Parameters
            ----------
            signum : Any
                Not used
            frame : Any
                Not used.

            Raises
            ------
            KeyboardInterrupt
                Always.

            """
            try:
                self.on_interrupt()  # your “final steps”
            finally:
                # Restore default and re-raise KeyboardInterrupt to stop execution
                signal.signal(signal.SIGINT, self._orig_sigint)
                raise KeyboardInterrupt

        logger.debug("Running model")
        start = perf_counter()
        signal.signal(signal.SIGINT, _handler)
        if self.run_type == "managed":
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
            output = self.run_fortran()

        logger.debug(f"Model finished. Took {perf_counter() - start:.2f} seconds.")

        signal.signal(signal.SIGINT, self._orig_sigint)

        if hasattr(self, "_shm_handles"):
            self._coordinator_unlink_memory()

        if "success_flag" not in output:
            msg = "Output dictionary from 'model.run_fortran()' does not contain 'success_flag' key."
            raise RuntimeError(msg)
        if output["success_flag"] is None or not isinstance(output["success_flag"], int):
            msg = "Output dictionary from 'model.run_fortran()' does not have valid 'success_flag' type (integer),"
            msg += f" but had type {type(output['success_flag'])}"
            raise TypeError(msg)

        output["success_flag"] = SuccessFlag(output["success_flag"])
        for k, v in output.items():
            self.__setattr__(k, v)  # noqa: PLC2801

        self._array_clean()
        self.check_error(only_error=True)
        if self.outputFile is not None:
            logger.debug(f"Writing output file: {self.outputFile}")
            logger.debug(
                f"Physics array shape: {self.physics_array.shape if self.physics_array is not None else None}"
            )
            logger.debug(
                f"Chemical array shape: {self.chemical_abun_array.shape if self.chemical_abun_array is not None else None}"
            )
            try:
                self.legacy_write_full()
                logger.debug(f"Successfully wrote {self.outputFile}")
            except Exception as e:
                logger.error(f"Failed to write {self.outputFile}: {e}", exc_info=True)
                raise
        if self.abundSaveFile is not None:
            logger.debug(f"Writing abundance file: {self.abundSaveFile}")
            try:
                self.legacy_write_starting_chemistry()
                logger.debug(f"Successfully wrote {self.abundSaveFile}")
            except Exception as e:
                logger.error(f"Failed to write {self.abundSaveFile}: {e}", exc_info=True)
                raise

    @abstractmethod
    def run_fortran(self) -> dict[str, int | list]:  # noqa: D102
        raise NotImplementedError

    # /Methods to start run of model

    # Model saving
    def save_model(
        self,
        file: h5py.File | str | Path,
        name: str = "default",
        overwrite: bool = False,
        array_dtype: np.typing.DTypeLike | np.dtype | str | None = None,
    ) -> None:
        """Save a model to file on disk.

        Multiple models can be saved into the same file if different names ``name``
        are used to store them.

        Parameters
        ----------
        file : h5py.File | str | Path
            Open h5py file object, or Path to a file that contains
            previously run and stored models.
        name : str
            Name to use for the group of the object. Default = "default".
        overwrite : bool
            Whether to overwrite pre-existing models, or just warn and
            not write if a dataset with ``name`` is already in the file. Default = False.
        array_dtype : np.typing.DTypeLike | np.dtype | str | None
            Precision to save arrays such
            as chemical abundances and physical conditions in. Can be used to save some storage.
            Default = None (infer dtype from arrays themselves).

        Raises
        ------
        TypeError
            if ``file`` is not a string, Path or ``h5py.File`` instance.

        Notes
        -----
        Saving arrays as ``np.float16`` is technically possible, but strongly discouraged.
            The maximum value that can be saved in a ``np.float16`` is less than $10^5$, so
            even densities of $10^6$ cm$^{-3}$ would cause it to overflow.

        """
        if isinstance(file, str | Path):
            opened_file = True
            file_obj = h5py.File(file, "a")
            logger.debug(f"Opened file {file}")
        elif isinstance(file, h5py.File):
            opened_file = False
            file_obj = file
        else:
            msg = f"Expected file to be type h5py.File, string or Path, but got type {type(file)}"
            raise TypeError(msg)

        if name in file_obj:
            if not overwrite:
                msg = f"Model with name '{name}' already exists in save file but overwrite is set to False. Unable to save model."
                warnings.warn(msg, stacklevel=2)
                return
            else:
                logger.debug(f"Deleting group {name} in file {file_obj.filename}")
                del file_obj[name]

        temp_attribute_dict = {}
        with contextlib.suppress(Exception):
            temp_attribute_dict.update(super().__getattribute__("_meta"))

        # Work on a copy so save_model is non-destructive to self._data
        save_data = self._data.copy()
        # Collect remaining non-array dataset variables into attributes (same behavior as before)
        v: str
        for v in list(save_data.variables):  # type: ignore[assignment]
            if "_array" not in v and v != "_orig_sigint":
                if np.shape(save_data[v].values) != ():
                    if isinstance(save_data[v].values, tuple):
                        temp_attribute_dict[v] = save_data[v].values[1].tolist()
                    else:
                        temp_attribute_dict[v] = save_data[v].values.tolist()
                else:
                    temp_attribute_dict[v] = save_data[v].item()
                save_data = save_data.drop_vars(v)
        for key, value in temp_attribute_dict.items():
            if isinstance(value, Path):
                temp_attribute_dict[key] = str(value)

        save_data["attributes_dict"] = xr.DataArray([json.dumps(temp_attribute_dict)])
        save_data["_param_dict"] = xr.DataArray([json.dumps(self._param_dict)])
        model_group = file_obj.create_group(name)
        coord_grp = model_group.create_group("_coords")

        save_name: str
        for save_name, coord in save_data.coords.items():  # type: ignore[assignment]
            _write_array(coord_grp, save_name, coord, array_dtype=array_dtype)
        for save_name, var in save_data.data_vars.items():  # type: ignore[assignment]
            _write_array(model_group, save_name, var, array_dtype=array_dtype)
        if opened_file:
            file_obj.flush()
            logger.debug(f"Closing file {file_obj.filename}")
            file_obj.close()

    # /Model saving

    # Model Passing through Pickling
    def pickle(self) -> AbstractModel:
        """Pickle the model.

        Returns
        -------
        AbstractModel
            Pickled model

        """
        v: str
        if self._data is not None and not bool(self._pickle_dict):
            for v in self._data.variables:  # type: ignore[assignment]
                if np.shape(self._data[v].values) != ():
                    if isinstance(self._data[v].values, tuple):
                        self._pickle_dict[v] = self._data[v].values[1].tolist()
                    else:
                        self._pickle_dict[v] = self._data[v].values.tolist()
                else:
                    self._pickle_dict[v] = self._data[v].item()
            # Save metadata separately for pickle roundtrip
            try:
                object.__setattr__(  # noqa: PLC2801
                    self, "_pickle_meta", super().__getattribute__("_meta").copy()
                )
            except Exception:
                object.__setattr__(self, "_pickle_meta", {})  # noqa: PLC2801
            self._data = xr.Dataset()
            # Clear runtime metadata to reflect pickled state
            object.__setattr__(self, "_meta", {})  # noqa: PLC2801
        return self

    def un_pickle(self) -> AbstractModel:
        """Un-pickle the model.

        Returns
        -------
        AbstractModel
            Unpickled model.

        """
        if (self._data is None or len(self._data.dims) == 0) and bool(self._pickle_dict):
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
                    object.__setattr__(self, "_meta", self._pickle_meta.copy())  # noqa: PLC2801
            except Exception as e:
                logger.exception(f"Exception occurred while restoring metadata: {e}")
            finally:
                object.__setattr__(self, "_pickle_meta", {})  # noqa: PLC2801
        else:
            warnings.warn("Un-pickling an object that was not pickled.", stacklevel=2)
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
        files place a ``point`` column between the physics and chemistry columns; we
        therefore use the location of ``point`` in the header to split the columns
        reliably and avoid using global constants as authoritative metadata.

        Parameters
        ----------
        read_file : str | Path
            path to file
        rate_constants_load_file : str | Path | None
            Not used. Default = None.

        Raises
        ------
        ValueError
            If there is any incompatibility error.

        """
        self.was_read = True
        # Read header and numeric data
        columns = np.char.strip(
            np.loadtxt(read_file, delimiter=",", max_rows=1, dtype=str, comments="%")
        )
        raw_array = np.loadtxt(read_file, delimiter=",", skiprows=1)

        # Determine where the ``point`` column is and how many points exist
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

            if missing_params <= {"dstep", "parcel_radius"} and not extra_params:
                # dstep and/or parcel_radius missing — check if we can safely infer
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
                        stacklevel=2,
                    )
                    if "dstep" in missing_params:
                        # Add dstep=1 column before point
                        dstep_column = np.ones((array.shape[0], 1))
                        array = np.hstack(
                            [array[:, :point_index], dstep_column, array[:, point_index:]]
                        )
                        physics_cols_from_file.append("dstep")
                        point_index += 1  # point column shifted by 1
                    if "parcel_radius" in missing_params:
                        # Add parcel_radius=0 column before point (collapse files only)
                        parcel_radius_column = np.zeros((array.shape[0], 1))
                        array = np.hstack(
                            [
                                array[:, :point_index],
                                parcel_radius_column,
                                array[:, point_index:],
                            ]
                        )
                        physics_cols_from_file.append("parcel_radius")
                        point_index += 1  # point column shifted by 1
                else:
                    msg = (
                        f"INCOMPATIBLE LEGACY FILE: Cannot infer 'dstep' parameter.\n\n"
                        f"The file is missing 'dstep' and contains duplicate timesteps,\n"
                        f"making it impossible to safely infer the particle step values.\n\n"
                        f"  File contains:        {physics_cols_from_file}\n"
                        f"  Current UCLCHEM has:  {list(PHYSICAL_PARAMETERS)}\n\n"
                        f"To load this file, regenerate it with the current UCLCHEM version."
                    )
                    raise ValueError(msg)
            else:
                # Other parameter mismatch - cannot fix automatically
                msg = (
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
                raise ValueError(msg)

        if species_cols_from_file != list(get_species_names()):
            msg = (
                f"INCOMPATIBLE LEGACY FILE: Species list mismatch.\n\n"
                f"The file you are loading has a different species network than the currently installed UCLCHEM.\n"
                f"This means the file was created with a different chemical network and cannot be loaded.\n\n"
                f"  File contains:        {len(species_cols_from_file)} species\n"
                f"  Current UCLCHEM has:  {len(get_species_names())} species\n\n"
                f"The species list is tied to the chemical network compiled into UCLCHEM.\n"
                f"To load this file, you must use the same UCLCHEM version/network that created it, or regenerate the file.\n"
            )
            raise ValueError(msg)

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
        # store ``timepoints`` as the number of simulated timesteps minus one.
        try:
            tp = int(self._data.sizes.get("time_step", 0))
            object.__setattr__(self, "timepoints", tp - 1 if tp > 0 else 0)  # noqa: PLC2801
        except Exception:
            # Be defensive; if something goes wrong leave timepoints as-is
            logger.debug("Could not set timepoints. Leave them unchanged.")

        nonzero_indices = self.physics_array[:, 0, 0].nonzero()[0]
        if len(nonzero_indices) > 0:
            last_timestep_index = nonzero_indices[-1]
            self.next_starting_chemistry_array = self.chemical_abun_array[
                last_timestep_index, :, :
            ]
        else:
            # Model failed immediately, no valid timesteps
            self.next_starting_chemistry_array = None

    def legacy_read_starting_chemistry(self) -> None:
        """Read the starting chemistry from the self.abundLoadFile provided in _param_dict.

        Raises
        ------
        ValueError
            If ``self.abundLoadFile`` is None.

        """
        if self.abundLoadFile is None:
            msg = "abundLoadFile was None, so cannot read the starting chemistry."
            raise ValueError(msg)
        logger.debug(f"Loading starting chemistry from file {self.abundLoadFile}")
        self._create_starting_array(np.loadtxt(self.abundLoadFile, delimiter=","))

    def legacy_write_full(self) -> None:
        """Perform classic output writing to file ``self.outputFile`` provided in ``_param_dict``.

        Raises
        ------
        RuntimeError
            If ``self.physics_array`` and ``self.chemical_abun_array``
            have not yet been initialized.
        ValueError
            If ``self.outputFile`` is None.

        """
        logger.debug(f"Writing output to {self.outputFile}")
        if self.physics_array is None or self.chemical_abun_array is None:
            msg = "Model arrays have not yet been initialized."
            raise RuntimeError(msg)

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

        if self.outputFile is None:
            msg = "outputFile was None, so cannot write the output."
            raise ValueError(msg)
        np.savetxt(self.outputFile, columns, fmt=string_fmt_string)
        with self.outputFile.open("ab") as f:
            np.savetxt(f, full_array, fmt=number_fmt_string)

    def legacy_write_starting_chemistry(self) -> None:
        """Perform classic starting abundance writing to file ``self.abundSaveFile``.

        Raises
        ------
        RuntimeError
            If ``self.chemical_abun_array`` is None.
        ValueError
            If ``self.abundSaveFile`` is None.

        """
        logger.debug(f"Writing starting chemistry to {self.abundSaveFile}")
        if self.chemical_abun_array is None:
            msg = "Chemical abundance array is None, so cannot write starting chemistry."
            raise RuntimeError(msg)
        last_timestep_index = self.chemical_abun_array[:, 0, 0].nonzero()[0][-1]
        # TODO Move away from the magic numbers seen here.
        species_names = get_species_names()
        number_fmt_string = f" {', '.join(['%9.5E'] * len(species_names))}"

        if self.abundSaveFile is None:
            msg = "abundSaveFile was None, so cannot write the output."
            raise ValueError(msg)
        with self.abundSaveFile.open("wb") as f:
            np.savetxt(
                f,
                self.chemical_abun_array[last_timestep_index, :, :],
                fmt=number_fmt_string,
            )

    # /Legacy in & output support

    # Cleaning of array & inptus
    def _array_clean(self):
        """Clean the arrays changed by UCLCHEM Fortran code."""
        logger.debug("Cleaning Fortran arrays")

        # Find the first element with all the zeros
        nonzero_indices = self.physics_array[:, 0, 0].nonzero()[0]
        if len(nonzero_indices) == 0:
            # Model failed immediately, keep only the first row to indicate failure
            logger.debug("Found 0 non-zero indices in physics_array")
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
            elif len(PHYSICAL_PARAMETERS) == phys_len:
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
            elif len(species_names) == chem_len:
                self._data = self._data.assign_coords(
                    {"chemical_abun_values": species_names}
                )
            else:
                self._data = self._data.assign_coords(
                    {"chemical_abun_values": np.arange(chem_len)}
                )

    def _reform_inputs(self, param_dict: dict | None, out_species: list[str]) -> None:
        """Reformat the input parameter dictionary.

        Copies param_dict so as not to modify user's dictionary.
        Then reformats out_species from pythonic list
        to a string of space separated names for Fortran.

        Parameters
        ----------
        param_dict : dict | None
            Parameter dictionary passed by the user to the model.
        out_species : list[str]
            List of output species that are considered important for this model.

        Raises
        ------
        ValueError
            If an duplicate key is encountered in ``param_dict``.
        ValueError
            If a key in ``param_dict`` is not in
            ``uclchem.constants.default_param_dictionary``.
        ValueError
            If an entry in ``out_species`` is not a valid species name.

        """
        if param_dict is None:
            # avoid mutating the shared default dictionary
            self._param_dict = default_param_dictionary.copy()
        else:
            # lower case (and conveniently copy so we don't edit) the user's dictionary
            # this is key to UCLCHEM's "case insensitivity"
            new_param_dict = convert_keys_to_lowercase(param_dict)
            for key, value in new_param_dict.items():
                if isinstance(value, Path):
                    new_param_dict[key] = str(value)

            # Handle deprecated endAtFinalDensity parameter (after lowercasing)
            new_param_dict = _convert_legacy_stopping_param(new_param_dict)  # type: ignore[assignment]

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
                msg = "out_species must be a list/tuple of valid species names; check available species via uclchem.model.get_species_names()"
                raise ValueError(msg)
            self.n_out = len(out_species)
            self._param_dict["outspecies"] = self.n_out
            self.out_species = " ".join(out_species)
        else:
            self.out_species = ""
            self.n_out = 0

    # /Cleaning of array & inptus

    # Creation of arrays
    def _create_fortran_array(self):
        """Create Fortran compliant np.arrays for physics and chemical abundances."""
        logger.debug("Creating fortran arrays for physics and abundances")
        # For shared memory:
        (
            self._shm_handles["physics_array"],
            self._shm_desc["physics_array"],
            self.physics_array,
        ) = _create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], N_PHYSICAL_PARAMETERS)
        )
        (
            self._shm_handles["chemical_abun_array"],
            self._shm_desc["chemical_abun_array"],
            self.chemical_abun_array,
        ) = _create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], n_species)
        )

    def _create_rate_constants_array(self):
        """Create Fortran compliant np.array for rate constants."""
        logger.debug("Creating fortran arrays for reaction rate constants")
        # For shared memory:
        (
            self._shm_handles["rate_constants_array"],
            self._shm_desc["rate_constants_array"],
            self.rate_constants_array,
        ) = _create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], n_reactions)
        )

    def _create_heating_array(self):
        """Create Fortran compliant np.array for heating/cooling rates."""
        logger.debug("Creating fortran arrays for heating")
        if not hasattr(uclchemwrap, "heating"):
            # Heating module not available, likely compiled without heating support
            logger.warning("Heating module not available in uclchemwrap")
            self.heat_array = None
            return

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
        ) = _create_shared_memory_allocation(
            (
                self.timepoints + 1,
                self._param_dict["points"],
                heating_array_size,
            )
        )

    def _create_stats_array(self):
        """Create Fortran compliant np.array for DVODE solver statistics."""
        logger.debug("Creating fortran arrays for DVODE solver statistics")
        (
            self._shm_handles["stats_array"],
            self._shm_desc["stats_array"],
            self.stats_array,
        ) = _create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], N_DVODE_STATS)
        )

    def _create_level_populations_array(self):
        """Create Fortran compliant np.array for coolant level populations.

        Shape: (timepoints+1, gridpoints, total_levels).

        """
        logger.debug("Creating fortran arrays for coolant level populations")
        (
            self._shm_handles["level_populations_array"],
            self._shm_desc["level_populations_array"],
            self.level_populations_array,
        ) = _create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], N_TOTAL_LEVELS)
        )

    def _create_se_stats_array(self):
        """Create Fortran compliant np.array for SE solver statistics.

        Shape: (timepoints+1, gridpoints, NCOOLANTS*3).

        """
        logger.debug("Creating fortran arrays for SE solver statistics")

        n_stats = NCOOLANTS * N_SE_STATS_PER_COOLANT  # 35 * 3 = 105

        (
            self._shm_handles["se_stats_array"],
            self._shm_desc["se_stats_array"],
            self.se_stats_array,
        ) = _create_shared_memory_allocation(
            (self.timepoints + 1, self._param_dict["points"], n_stats)
        )

    def get_level_populations_dataframe(self, point: int = 0) -> pd.DataFrame | None:
        """Get level populations as a DataFrame for a specific grid point.

        Parameters
        ----------
        point : int
            Grid point index. Default = 0.

        Returns
        -------
        pd.DataFrame | None
            DataFrame with columns for each level with meaningful
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
            logger.warning(
                f"Could not load coolant level names: {e}. Using generic names."
            )
            columns = [f"LEVEL_{i}" for i in range(N_TOTAL_LEVELS)]

        return pd.DataFrame(self.level_populations_array[:, point, :], columns=columns)

    def get_se_stats_dataframe(self, point: int = 0) -> pd.DataFrame | None:
        """Get SE solver statistics as a DataFrame for a specific grid point.

        Parameters
        ----------
        point : int
            Grid point index. Default = 0.

        Returns
        -------
        pd.DataFrame | None
            DataFrame with per-coolant SE solver statistics using
            actual coolant names

        """
        if self.se_stats_array is None or self.se_stats_array.shape[0] < 3:
            return None

        # Build meaningful column names using actual coolant names
        try:
            coolant_names = [
                str(name.decode()).strip()
                for name in uclchemwrap.f2py_constants.coolantnames
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

    def _create_starting_array(self, starting_chemistry: ArrayLike | None) -> None:
        if starting_chemistry is None:
            logger.debug(
                "starting_chemistry is None, so not creating starting_chemistry_array"
            )
            self.starting_chemistry_array = None
        else:
            logger.debug("Creating starting chemistry array")
            starting_chemistry = np.asarray(starting_chemistry)
            if len(np.shape(starting_chemistry)) == 1:
                starting_chemistry = starting_chemistry[np.newaxis, :]

            if np.shape(starting_chemistry)[1] != n_species:
                msg = f"Mismatch between number of species in the network ({n_species}) and number of entries in starting_chemistry ({np.shape(starting_chemistry)[1]})."
                msg += " The first model was possibly run with a different network."
                raise RuntimeError(msg)

            # For shared memory:
            (
                self._shm_handles["starting_chemistry_array"],
                self._shm_desc["starting_chemistry_array"],
                self.starting_chemistry_array,
            ) = _create_shared_memory_allocation(np.shape(starting_chemistry))
            np.copyto(self.starting_chemistry_array, starting_chemistry, casting="no")

    # /Creation of arrays

    # Signal Interrupt Catch
    def on_interrupt(self, grid: bool = False, model_name: str | None = None) -> None:
        """Catch interruption. Save model to file.

        Parameters
        ----------
        grid : bool
            whether the model was part of a grid. Default = False
        model_name : str | None
            the name of the model to save it under.
            If None, name is set to "interrupted". Default = None.

        """
        logger.info("Model was interrupted")

        if self._proc_handle is not None:
            logger.debug("Terminating process")
            self._proc_handle.terminate()
            self._proc_handle.join()
            self._proc_handle = None

        if bool(self._shm_desc):
            self._coordinator_unlink_memory()

        self._array_clean()

        error_time = datetime.now().strftime("%y_%m_%d_%H_%M")
        if not grid:
            if self.outputFile is None:
                self.outputFile = Path("./" + error_time + ".dat")
            elif "/" in str(self.outputFile):
                outputFile_str = str(self.outputFile)  # noqa: N806
                self.outputFile = Path(
                    outputFile_str[: outputFile_str.rfind("/") + 1]
                    + error_time
                    + outputFile_str[outputFile_str.rfind(".") :]
                )
            else:
                self.outputFile = Path("./" + error_time + ".dat")
            self.legacy_write_full()

        self._was_interrupted = True
        self.save_model(
            file="./" + error_time + ".h5"
            if not grid
            else "./grid_interrupted_models.h5",
            name=model_name if model_name is not None else "interrupted",
            overwrite=True,
        )

    # /Signal Interrupt Catch

    # Shared memory handlers

    def _reform_array_in_worker(self, shm_desc: dict[str, dict]) -> None:
        """Reform SharedMemory objects to allocate during the model.

        Parameters
        ----------
        shm_desc : dict[str, dict]
            Description of shared memory.

        """
        object.__setattr__(self, "_shm_handles", {})  # noqa: PLC2801
        for k, v in shm_desc.items():
            shm = shared_memory.SharedMemory(name=v["name"], create=False)
            object.__setattr__(  # noqa: PLC2801
                self,
                k,
                np.ndarray(shape=v["shape"], dtype=np.float64, buffer=shm.buf, order="F"),
            )
            self._shm_handles[k] = shm
            del shm

    def _worker_close_memory(self):
        for k in self._shm_desc:
            try:
                self._shm_handles[k].close()
            except Exception:  # noqa: S110
                pass
            finally:
                del self._shm_handles[k]

    def _coordinator_unlink_memory(self):
        logger.debug("Unlinking shared memory")
        if bool(self._shm_desc):
            for k in self._shm_desc:
                try:
                    self.__setattr__(k, self.__getattr__(k).copy())  # noqa: PLC2801
                    self._shm_handles[k].close()
                    self._shm_handles[k].unlink()
                except Exception:
                    print(f"Warning, unable to close and unlink {k}")
                    pass
                finally:
                    del self._shm_handles[k]
            del self._shm_desc
            self._shm_desc = {}
            del self._shm_handles
            self._shm_handles = {}

    # /Shared memory handlers

    @abstractmethod
    def _create_init_dict(self):
        raise NotImplementedError


@register_model
class Cloud(AbstractModel):
    """Cloud model class inheriting from AbstractModel."""

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: ArrayLike | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Create a new ``Cloud`` instance.

        Initiates the model first with ``AbstractModel.__init__()``,
        then with any additional commands needed for the model.

        Parameters
        ----------
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in ``defaultparameters.f90``. Default = None.
        out_species : list[str] | None
            List of species whose abundances at the end of the model are
            returned. If None, defaults to ``uclchem.constants.default_elements_to_check``.
            Default = None.
        starting_chemistry : ArrayLike | None
            Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model
            that will be run. Defaults to None.
        timepoints : int
            Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            read_file=read_file,
            run_type=run_type,
        )
        if self.run_type != "external" and not self.was_read:
            self.run()

    def run_fortran(self) -> dict[str, int | list]:
        """Run the fortran side of the UCLCHEM model.

        Returns
        -------
        dict[str, int | list]
            Dictionary with two keys:
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
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")  # noqa: PLC2801
            else None,
        )
        abundance_out, _specname_out, success_flag = result[-3], result[-2], result[-1]
        if success_flag < 0:
            out_species_abundances_array = []
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
        }


@register_model
class Collapse(AbstractModel):
    """Collapse model class inheriting from AbstractModel."""

    # Time (years) at which each collapse mode's density evolution ends and the fitting
    # functions become singular.
    _COLLAPSE_FINAL_TIMES = {
        "BE1.1": 1.173387e6,
        "BE4": 1.84265e5,
        "filament": 1.393761e6,
        "ambipolar": 1.6132984e7,
    }

    def __init__(
        self,
        collapse: Literal["BE1.1", "BE4", "filament", "ambipolar"] = "BE1.1",
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: ArrayLike | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Create a new ``Collapse`` model instance.

        Initiates the model first with ``AbstractModel.__init__()``,
        then with any additional commands needed for the model.

        Parameters
        ----------
        collapse : Literal['BE1.1', 'BE4', 'filament', 'ambipolar']
            A string containing
            the collapse type. Defaults to 'BE1.1'.
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in ``defaultparameters.f90``. Default = None.
        out_species : list[str] | None
            List of species whose abundances at the end of the model are
            returned. If None, defaults to ``uclchem.constants.default_elements_to_check``.
            Default = None.
        starting_chemistry : ArrayLike | None
            Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model
            that will be run. Defaults to None.
        timepoints : int
            Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If ``collapse`` is not one of `["BE1.1", "BE4", "filament", "ambipolar"]`.

        """
        collapse_dict = {"BE1.1": 1, "BE4": 2, "filament": 3, "ambipolar": 4}
        if collapse not in collapse_dict:
            msg = f"collapse must be in {collapse_dict.keys()}"
            raise ValueError(msg)

        collapse_final_time = self._COLLAPSE_FINAL_TIMES[collapse]

        if out_species is None:
            out_species = default_elements_to_check

        if param_dict is not None and "initialDens" in param_dict:
            warnings.warn(
                "initialDens is ignored for collapse models: the initial density is determined "
                "with fit functions, not the initialDens parameter.",
                UserWarning,
                stacklevel=2,
            )
        if (
            param_dict is not None
            and param_dict.get("points", 1) == 1
            and param_dict.get("r_in", 0.0) != 0.0
        ):
            msg = (
                "r_in has no effect when points=1: the single parcel is placed at rout. "
                "Either set points > 1 or remove r_in from param_dict."
            )
            raise ValueError(msg)

        # For collapse models, endAtFinalDensity controls finalTime behavior.
        # Reject parcelStoppingMode since collapse models don't use density-based stopping.
        _param = param_dict or {}
        if "parcelStoppingMode" in _param:
            msg = (
                "parcelStoppingMode is not supported for collapse models. "
                "Use endAtFinalDensity instead: True (default) to stop at collapse endpoint, "
                "False to extend chemistry until a custom finalTime with frozen density."
            )
            raise ValueError(msg)

        # endAtFinalDensity controls finalTime for collapse models:
        # - True (default): finalTime = collapseFinalTime (self-consistent density
        #   evolution until collapse)
        # - False: user must set finalTime > collapseFinalTime (density frozen,
        #   chemistry continues)
        user_final_time = _param.get("finalTime", None)
        end_at_final_density = _param.get(
            "endAtFinalDensity", True
        )  # Default True for collapse

        if end_at_final_density:
            # Scenario 1 (default): stop at collapse endpoint.
            if user_final_time is not None and user_final_time != collapse_final_time:
                msg = (
                    f"For {collapse!r} collapse with endAtFinalDensity=True, finalTime is fixed at "
                    f"collapseFinalTime={collapse_final_time:.3e} yr. "
                    f"To use a custom finalTime, set endAtFinalDensity=False."
                )
                raise ValueError(msg)
            param_dict = {**_param, "finalTime": collapse_final_time}
            # Remove endAtFinalDensity so _convert_legacy_stopping_param doesn't affect it.
            param_dict.pop("endAtFinalDensity", None)
        else:
            # Scenario 2: user extends chemistry past collapse endpoint with frozen density.
            if user_final_time is None:
                msg = (
                    f"endAtFinalDensity=False requires setting finalTime > collapseFinalTime "
                    f"({collapse_final_time:.3e} yr) to extend chemistry beyond the collapse endpoint."
                )
                raise ValueError(msg)
            if user_final_time < collapse_final_time:
                msg = (
                    f"endAtFinalDensity=False with finalTime={user_final_time:.3e} yr < "
                    f"collapseFinalTime={collapse_final_time:.3e} yr is invalid. "
                    f"Either set endAtFinalDensity=True (use default finalTime) or "
                    f"set finalTime > collapseFinalTime for extended chemistry."
                )
                raise ValueError(msg)
            # finalTime > collapseFinalTime is valid; density will freeze at collapseFinalTime.
            # Remove endAtFinalDensity so _convert_legacy_stopping_param doesn't affect it.
            param_dict = {**_param, "parcelStoppingMode": 0}
            param_dict.pop("endAtFinalDensity", None)

        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            read_file=read_file,
            run_type=run_type,
        )
        self.collapse_final_time = collapse_final_time
        if read_file is None:
            self.collapse = collapse_dict[collapse]
            if self.run_type != "external":
                self.run()

    def run_fortran(self) -> dict[str, int | list]:
        """Run the fortran side of the UCLCHEM model.

        Returns
        -------
        dict[str, int | list]
            Dictionary with two keys:
            "success_flag" with value the success flag
            "out_species_abundances_array" with value a list of the outspecies abundances.

        """
        result = wrap.collapse(
            collapsein=self.collapse,
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
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")  # noqa: PLC2801
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
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
        }


@register_model
class PrestellarCore(AbstractModel):
    """PrestellarCore model class inheriting from AbstractModel.

    This model type was previously known as hot core.

    """

    def __init__(
        self,
        temp_index: int = 1,
        max_temperature: float = 300.0,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: ArrayLike | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Create a new ``PrestellarCore`` model.

        Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Parameters
        ----------
        temp_index : int
            Used to select the mass of the prestellar core from the following selection
            [1=1Msun, 2=5, 3=10, 4=15, 5=25,6=60]. Defaults to 1, which is 1 Msun
        max_temperature : float
            Value at which gas temperature will stop increasing.
            Defaults to 300.0 K.
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in ``defaultparameters.f90``. Default = None.
        out_species : list[str] | None
            List of species whose abundances at the end of the model are
            returned. If None, defaults to ``uclchem.constants.default_elements_to_check``.
            Default = None.
        starting_chemistry : ArrayLike | None
            Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model that
            will be run. Defaults to None.
        timepoints : int
            Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If ``read_file`` is None, but ``temp_idx`` or ``max_temperature`` is also None.

        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None:
            if temp_index is None or max_temperature is None:
                msg = "temp_index and max_temperature must be specified if not reading from file."
                raise ValueError(msg)
            self.temp_index = temp_index
            self.max_temperature = max_temperature
            if self.run_type != "external":
                self.run()

    def run_fortran(self) -> dict[str, int | list]:
        """Run the fortran side of the UCLCHEM model.

        Returns
        -------
        dict[str, int | list]
            Dictionary with two keys:
            "success_flag" with value the success flag
            "out_species_abundances_array" with value a list of the outspecies abundances.

        """
        _, _, _, _, _, _, _, out_species_abundances_array, _, success_flag = (
            wrap.hot_core(
                temp_index=self.temp_index,
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
                if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")  # noqa: PLC2801
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
            "temp_index": self.temp_index,
            "max_temperature": self.max_temperature,
            "_param_dict": self._param_dict,
            "out_species_list": self.out_species_list,
            "out_species": self.out_species,
            "timepoints": self.timepoints,
            "give_start_abund": self.give_start_abund,
        }


@register_model
class CShock(AbstractModel):
    """C-Shock model class inheriting from AbstractModel."""

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
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Create a new ``CShock`` instance.

        Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Parameters
        ----------
        shock_vel : float
            Velocity of the shock in km/s. Defaults to 10.0.
        timestep_factor : float
            Whilst the time is less than 2 times the dissipation time of shock,
            timestep is timestep_factor*dissipation time. Essentially controls how well resolved the
            shock is in your model. Defaults to 0.01.
        minimum_temperature : float
            Minimum post-shock temperature. Defaults to 0.0 (no minimum).
            The shocked gas typically cools to ``initialTemp`` if this is not set.
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in ``defaultparameters.f90``. Default = None.
        out_species : list[str] | None
            List of species whose abundances at the end of the model are
            returned. If None, defaults to ``uclchem.constants.default_elements_to_check``.
            Default = None.
        starting_chemistry : np.ndarray | None
            Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model that
            will be run. Defaults to None.
        timepoints : int
            Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If ``read_file`` is None, but ``shock_vel`` is also not set.

        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None:
            if shock_vel is None:
                msg = "shock_vel must be specified if not reading from file."
                raise ValueError(msg)
            self.shock_vel = shock_vel
            self.timestep_factor = timestep_factor
            self.minimum_temperature = minimum_temperature
            self.dissipation_time = -1
            if self.run_type != "external":
                self.run()

    def run_fortran(self) -> dict[str, int | list]:
        """Run the fortran side of the UCLCHEM model.

        Returns
        -------
        dict[str, int | list]
            Dictionary with two keys:
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
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")  # noqa: PLC2801
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
        }


@register_model
class JShock(AbstractModel):
    """J-Shock model class inheriting from AbstractModel."""

    def __init__(
        self,
        shock_vel: float = 10.0,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        timepoints: int = TIMEPOINTS,
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Instantiate a new ``JShock`` model.

        Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Parameters
        ----------
        shock_vel : float
            Velocity of the shock. Defaults to 10.0.
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in ``defaultparameters.f90``. Default = None.
        out_species : list[str] | None
            List of species whose abundances at the end of the model are
            returned. If None, defaults to ``uclchem.constants.default_elements_to_check``.
            Default = None.
        starting_chemistry : np.ndarray | None
            Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model that
            will be run. Defaults to None.
        timepoints : int
            Integer value of how many timesteps should be calculated before
            aborting the UCLCHEM model. Defaults to ``uclchem.constants.TIMEPOINTS``.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If ``read_file`` is None, but ``shock_vel`` is also not set.

        """
        if out_species is None:
            out_species = default_elements_to_check
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=timepoints,
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is None:
            if shock_vel is None:
                msg = "shock_vel must be specified if not reading from file."
                raise ValueError(msg)
            self.shock_vel = shock_vel
            if self.run_type != "external":
                self.run()

    def run_fortran(self) -> dict[str, int | list]:
        """Run the fortran side of the UCLCHEM model.

        Returns
        -------
        dict[str, int | list]
            Dictionary with two keys:
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
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")  # noqa: PLC2801
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
        }


@register_model
class Postprocess(AbstractModel):
    """Postprocess represents a model class with additional controls. It inherits from AbstractModel.

    Postprocess allows for additional controls of the time, density, gas temperature, radiation field,
    cosmic ray ionization rate, atomic and molecular Hydrogen, CO and C column densities through the
    use of arrays. Using these arrays allows for experimental model crafting beyond the standard models
    in other model classes.

    """

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: ArrayLike | None = None,
        previous_model: AbstractModel | None = None,
        time_array: ArrayLike | None = None,
        density_array: ArrayLike | int | float | None = None,
        gas_temperature_array: ArrayLike | int | float | None = None,
        dust_temperature_array: ArrayLike | int | float | None = None,
        zeta_array: ArrayLike | int | float | None = None,
        radfield_array: ArrayLike | int | float | None = None,
        visual_extinction_array: ArrayLike | int | float | None = None,
        coldens_H_array: ArrayLike | int | float | None = None,  # noqa: N803
        coldens_H2_array: ArrayLike | int | float | None = None,  # noqa: N803
        coldens_CO_array: ArrayLike | int | float | None = None,  # noqa: N803
        coldens_C_array: ArrayLike | int | float | None = None,  # noqa: N803
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        r"""Initiate the postprocessing model.

        Initiates the model first with ``AbstractModel.__init__()``,
        then with any additional commands needed for the model.

        Any values passed as integers or floats are converted to arrays with the same
        length as ``time_array``. This means that if you for example pass
        ``density_array = 5e4``, the model runs with a constant density of
        $5 \times 10^{4}$ $cm$^{-3}$.

        Parameters
        ----------
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in ``defaultparameters.f90``. Default = None.
        out_species : list[str] | None
            List of species whose abundances at the end of the model are
            returned. If None, defaults to ``uclchem.constants.default_elements_to_check``.
            Default = None.
        starting_chemistry : ArrayLike | None
            Array containing the starting abundances
            to use for the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model
            that will be run. Defaults to None.
        time_array : ArrayLike | None
            Represents the time grid to be used for the model.
            This sets the target timesteps for which outputs will be stored. Default = None.
        density_array : ArrayLike | int | float | None
            Represents the value of the density at
            different timepoints found in time_array. Default = None.
        gas_temperature_array : ArrayLike | int | float | None
            Represents the value of the gas
            temperature at different timepoints found in time_array. Default = None.
        dust_temperature_array : ArrayLike | int | float | None
            Represents the value of the dust
            temperature at different timepoints found in time_array. Default = None.
        zeta_array : ArrayLike | int | float | None
            Represents the value of the cosmic ray
            ionization rate at different timepoints found in time_array. Default = None.
        radfield_array : ArrayLike | int | float | None
            Represents the value of the UV
            radiation field at different timepoints found in time_array. Default = None.
        visual_extinction_array : ArrayLike | int | float | None
            The value of the visual extinction
            Av at each of the timepoints in time_array. Default = None.
        coldens_H_array : ArrayLike | int | float | None
            Represents the value of the
            column density of H at different timepoints found in time_array. Default = None.
        coldens_H2_array : ArrayLike | int | float | None
            Represents the value of the
            column density of H2 at different timepoints found in time_array. Default = None.
        coldens_CO_array : ArrayLike | int | float | None
            Represents the value of the
            column density of CO at different timepoints found in time_array. Default = None.
        coldens_C_array : ArrayLike | int | float | None
            Represents the value of the
            column density of C at different timepoints found in time_array. Default = None.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None. Default = None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If not all arrays have the same length as ``time_array``.
        ValueError
            If ``read_file`` is None, but ``time_array`` is not an array.

        """
        # Allocate 1.5x the input timesteps to give the DVODE solver
        # headroom for additional internal substeps.
        if out_species is None:
            out_species = default_elements_to_check

        if read_file is None and time_array is None:
            msg = f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given."
            raise ValueError(msg)

        n_input = len(time_array) if time_array is not None else 1
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=int(n_input * 1.5),
            read_file=read_file,
            run_type=run_type,
        )
        if read_file is not None:
            return

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
                if isinstance(array, int | float):
                    array = np.ones(n_input) * array
                if len(array) != n_input:
                    msg = "All arrays must be the same length"
                    raise ValueError(msg)

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
        if self.usecoldens and self.useav:
            msg = "Cannot use both column density and visual extinction arrays simultaneously."
            raise ValueError(msg)

        if not self.give_start_abund:
            self.starting_chemistry_array = np.zeros(
                shape=(self.gridPoints, n_species),
                dtype=np.float64,
                order="F",
            )
        if self.run_type != "external":
            self.run()

    def run_fortran(self) -> dict[str, int | list]:
        """Run the fortran side of the UCLCHEM model.

        Returns
        -------
        dict[str, int | list]
            Dictionary with two keys:
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
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")  # noqa: PLC2801
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
            "postprocess_arrays": self.postprocess_arrays,
            **self.postprocess_arrays,
        }


@register_model
class Model(AbstractModel):
    """Model, like Postprocess, represents a model class with additional controls.

    Model follows the same logic as Postprocess but without the coldens Arguments.
    It allows for additional controls of the time, density, gas temperature, radiation field,
    and cosmic ray ionization rate through the use of arrays. Using these arrays allows for
    experimental model crafting beyond the standard models in other model classes.

    """

    def __init__(
        self,
        param_dict: dict | None = None,
        out_species: list[str] | None = None,
        starting_chemistry: np.ndarray | None = None,
        previous_model: AbstractModel | None = None,
        time_array: ArrayLike | None = None,
        density_array: ArrayLike | int | float | None = None,
        gas_temperature_array: ArrayLike | int | float | None = None,
        dust_temperature_array: ArrayLike | int | float | None = None,
        zeta_array: ArrayLike | int | float | None = None,
        radfield_array: ArrayLike | int | float | None = None,
        read_file: str | Path | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Create a new ``Model`` instance.

        Initiates the model first with AbstractModel.__init__(),
        then with any additional commands needed for the model.

        Parameters
        ----------
        param_dict : dict | None
            Dictionary containing the parameters to use for the UCLCHEM model.
            Uses UCLCHEM default values found in ``defaultparameters.f90``. Default = None.
        out_species : list[str] | None
            List of species whose abundances at the end of the model are
            returned. If None, defaults to ``uclchem.constants.default_elements_to_check``.
            Default = None.
        starting_chemistry : np.ndarray | None
            Array containing the starting abundances to use for
            the UCLCHEM model. Defaults to None.
        previous_model : AbstractModel | None
            Model object, a class that inherited from
            AbstractModel, to use for the starting abundances of the new UCLCHEM model
            that will be run. Defaults to None.
        time_array : ArrayLike | None
            Represents the time grid to be used for the model.
            This sets the target timesteps for which outputs will be stored. Default = None.
        density_array : ArrayLike | int | float | None
            Represents the value of the density at
            different timepoints found in time_array. Default = None.
        gas_temperature_array : ArrayLike | int | float | None
            Represents the value of the gas
            temperature at different timepoints found in time_array. Default = None.
        dust_temperature_array : ArrayLike | int | float | None
            Represents the value of the dust
            temperature at different timepoints found in time_array. Default = None.
        zeta_array : ArrayLike | int | float | None
            Represents the value of the cosmic ray
            ionization rate at different timepoints found in time_array. Default = None.
        radfield_array : ArrayLike | int | float | None
            Represents the value of the UV radiation
            field at different timepoints found in time_array. Default = None.
        read_file : str | Path | None
            Path to the file to be read. Reading a file to a model object,
            prevents it from being run. Defaults to None. Default = None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If not all arrays have the same length.
        ValueError
            If ``read_file`` is None, but ``time_array`` is not an array.

        """
        if out_species is None:
            out_species = default_elements_to_check

        if time_array is None and read_file is None:
            msg = f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given."
            raise ValueError(msg)

        # Allocate 1.5x the input timesteps to give the DVODE solver
        # headroom for additional internal substeps.
        n_input = len(time_array) if time_array is not None else 1
        super().__init__(
            param_dict=param_dict,
            out_species_list=out_species,
            starting_chemistry=starting_chemistry,
            previous_model=previous_model,
            timepoints=int(n_input * 1.5),
            read_file=read_file,
            run_type=run_type,
        )

        if read_file is not None:
            return

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
                if isinstance(array, int | float):
                    array = np.ones(n_input) * array
                if len(array) != n_input:
                    msg = "All arrays must be the same length"
                    raise ValueError(msg)

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

    def run_fortran(self) -> dict[str, int | list]:
        """Run the fortran side of the UCLCHEM model.

        Returns
        -------
        dict[str, int | list]
            Dictionary with two keys:
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
            if "starting_chemistry_array" in object.__getattribute__(self, "__dict__")  # noqa: PLC2801
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
            **self.postprocess_arrays,
        }


class SequentialRunner:
    """The SequentialRunner class allows for multiple models to be run back to back.

    By defining a specific dictionary to hold the information of each model class to run in sequence,
    SequentialModel allows for the automatic running of multiple models as well as matching some
    physical parameters from one model to the next.

    """

    def __init__(
        self,
        sequenced_model_parameters: list[dict[str, Any]],
        parameters_to_match: list[str] | None = None,
        run_type: Literal["managed", "external"] = "managed",
    ):
        """Initialize the SequentialRunner.

        Parameters
        ----------
        sequenced_model_parameters : list[dict[str, Any]]
            List of dictionaries with format
            [
                {"<First Model Class>":{"param_dict":{<parameters>}, <other arguments>}},
                {"<Second Model Class>":{"param_dict":{<parameters>}, <other arguments>}},
                {...},
            ]
        parameters_to_match : list[str] | None
            The list provided to this argument decides which
            parameters should be matched from a previous model to the next model in the sequence.
            Currently, supports one or more of ``["finalDens", "finalTemp"]``. Default = None.
        run_type : Literal['managed', 'external']
            Run type. "external" means that the model is not
            run directly after instantiation, but can instead be run as ``model.run()``.
            Default = "managed".

        Raises
        ------
        ValueError
            If any of the model names in ``sequenced_model_parameters``
            are ``"SequentialRunner"``.
        NotImplementedError
            If a parameter in ``parameters_to_match`` is not one of
            ``["finalDens", "finalTemp"]``.

        """
        for model in sequenced_model_parameters:
            if model[list(model.keys())[0]] == "SequentialRunner":
                msg = "Cannot run a SequentialRunner within a SequentialRunner"
                raise ValueError(msg)

        self.models: list[dict[str, Any]] = []
        self.sequenced_model_parameters = sequenced_model_parameters
        self.parameters_to_match = parameters_to_match

        if self.parameters_to_match is not None:
            for parameter in self.parameters_to_match:
                if parameter not in {"finalTemp", "finalDens"}:
                    msg = f"Parameter '{parameter}' has not been implemented for parameter matching"
                    raise NotImplementedError(msg)

        self.run_type = run_type
        self.model_count = 0
        self._pickle_dict: dict[str, dict] = {}
        self.success_flag: bool | None = None
        if self.run_type == "managed":
            self.run()

    def run(self) -> None:
        """Run the sequential model.

        Raises
        ------
        RuntimeError
            If ``previous_model`` is somehow None after the first model has been run.
        NotImplementedError
            If a parameter in ``parameters_to_match`` is not one of
            ``["finalDens", "finalTemp"]``.

        """
        previous_model: AbstractModel | None = None
        for base_model_dict in self.sequenced_model_parameters:
            for model_type, model_dict in base_model_dict.items():
                model_dict["param_dict"] = convert_keys_to_lowercase(
                    model_dict["param_dict"]
                )
                if self.model_count > 0:
                    if previous_model is None:
                        msg = "previous_model is still None, but expected it to be an AbstractModel"
                        raise RuntimeError(msg)
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
                                    previous_model.physics_array[-1, 0, 1].item()  # type: ignore[index]
                                )
                                continue
                            elif parameter == "finalTemp":
                                model_dict["param_dict"]["initialtemp"] = (
                                    previous_model.physics_array[-1, 0, 2].item()  # type: ignore[index]
                                )
                            else:
                                msg = f"Parameter '{parameter}' has not been implemented for parameter matching"
                                raise NotImplementedError(msg)

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

    def save_model(
        self,
        file: h5py.File | str | Path,
        name: str = "",
        overwrite: bool = False,
        array_dtype: np.typing.DTypeLike | np.dtype | str | None = None,
    ) -> None:
        """Save a model to an open file object or to a file.

        Parameters
        ----------
        file : h5py.File | str | Path
            open h5py file object, or Path to a file that contains
            previously run and stored models.
        name : str
            name to save model under. Default = "".
        overwrite : bool
            Boolean on whether to overwrite pre-existing models, or error out.
            Default = False
        array_dtype : np.typing.DTypeLike | np.dtype | str | None
            Precision to save arrays such
            as chemical abundances and physical conditions in. Can be used to save some storage.
            Default = None (infer dtype from arrays themselves).

        Raises
        ------
        TypeError
            if ``file`` is not a string, Path or ``h5py.File`` instance.

        """
        if isinstance(file, str | Path):
            opened_file = True
            file_obj = h5py.File(file, "a")
        elif isinstance(file, h5py.File):
            opened_file = False
            file_obj = file
        else:
            msg = f"Expected file to be type h5py.File, string or Path, but got type {type(file)}"
            raise TypeError(msg)

        for model in self.models:
            model["Model"].save_model(
                file=file_obj,
                name=f"{name}_{model['Model_Order']}_{model['Model_Type']}",
                overwrite=overwrite,
                array_dtype=array_dtype,
            )

        if opened_file:
            file_obj.close()

    @classmethod
    def load_from_dataset(  # noqa: D102
        cls, model_ds: xr.Dataset
    ) -> SequentialRunner:
        msg = "Loading a SequentialRunner from a dataset has not been implemented yet."
        raise NotImplementedError(msg)

    @classmethod
    def worker_build(cls, init_kwargs, shm_desc):  # noqa: ANN001, D102
        msg = "SequentialRunner.worker_build has not been implemented yet."
        raise NotImplementedError(msg)

    def check_conservation(self, element_list: list[str] | None = None) -> None:
        """Check conservation of the chemical abundances.

        Adds an entry in the self.models[model] dictionary with key
            ``"elements_conserved"``, which indicates whether all the points in model with
            index ``model_idx`` had changes in elemental abundances below 1%.

        Parameters
        ----------
        element_list : list[str] | None
            List of elements to check conservation for.
            If None, use ``uclchem.constants.default_elements_to_check``. Default = None.

        """
        if element_list is None:
            element_list = default_elements_to_check

        for model in self.models:
            conserve_dicts: list[dict[str, str]] = []
            if model["Model"]._param_dict["points"] > 1:
                for point in range(model["Model"]._param_dict["points"]):
                    conserve_dicts += [
                        check_element_conservation(
                            model["Model"].get_dataframes(point),
                            element_list=element_list,
                            percent=True,
                        )
                    ]
            else:
                conserve_dicts += [
                    check_element_conservation(
                        model["Model"].get_dataframes(0),
                        element_list=element_list,
                        percent=True,
                    )
                ]
            conserved = True
            for conserve_dict in conserve_dicts:
                conserved = conserved and all(
                    float(x[:1]) < 1 for x in conserve_dict.values()
                )
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

    Parameters
    ----------
    model_id : int
        id of model
    model_type : str
        string representing the type of model
    pending_model : dict[str, Any]
        dictionary with arguments necessary to initialize
        model.
    log_dir : str | Path | None
        If not None, write logs to "model_{model_id}.log".
        If None, do not write logs. Default = None.

    Returns
    -------
    model_id : int
        model id of run model
    model_obj : object
        pickled model object

    Raises
    ------
    ValueError
        If no model called ``model_type`` is in the model registry.

    """
    log_file = None
    if log_dir is not None:
        log_file = Path(log_dir) / f"model_{model_id}.log"

    cls = REGISTRY.get(model_type)
    if cls is None:
        msg = f"No model type {model_type} in registry."
        raise ValueError(msg)

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
    """Allows running multiple models on a grid of parameter space."""

    def __init__(
        self,
        model_type: str,
        full_parameters: dict | list,
        max_workers: int = 8,
        grid_file: str | Path = "./default_grid_out.h5",
        model_name_prefix: str = "",
        overwrite_models: bool = False,
        array_dtype: np.typing.DTypeLike | np.dtype | str | None = None,
        delay_run: bool = False,
        log_dir: str | Path | None = None,
        model_ids: list[int] | None = None,
        create_grid: bool = True,
    ):
        """Create a new ``GridRunner`` object.

        Parameters
        ----------
        model_type : str
            Type of model class to run
        full_parameters : dict | list
            The dictionary passed to GridRunner should nest into it,
            the param_dict argument that would be passed to any other model, with the addition
            of extra keys for the none param_dict variables of a model. Any variables that are
            turned into lists or arrays, will automatically be assumed to be used for the gridding.
        max_workers : int
            Maximum number of workers to use in parallel for the grid run.
            Defaults to 8.
        grid_file : str | Path
            Name and path of the output file to which the models should be saved.
            Defaults to "./default_grid_out.h5".
        model_name_prefix : str
            Name prefix convention to use. The fifth model in the grid
            would have the name "<model_name_prefix>5>" assigned to it. Defaults to "",
            which would make the fifth model have the name "5", for example.
        overwrite_models : bool
            Whether to overwrite the models in ``grid_file`` if they have
            the same name. Default = False.
        array_dtype : np.typing.DTypeLike | np.dtype | str | None
            Precision to save arrays such
            as chemical abundances and physical conditions in. Can be used to save some storage.
            Default = None (infer dtype from arrays themselves).
        delay_run : bool
            Whether to immediately start the models upon initialization,
            or delay until the user calls ``self.run()``. Defaults to False (start immediately).
        log_dir : str | Path | None
            Where to write logs. If None, do not write logs. Default = None.
        model_ids : list[int] | None
            Optional subset of model indices (0-based column in flat_grids)
            to run. None means run all models in the grid. Default = None.
        create_grid : bool
            whether to create the grid from the arrays in ``full_parameters``
            (using ``np.meshgrid``) or not. Default = True.

        Raises
        ------
        ValueError
            If ``model_type`` is not in the model registry.
        TypeError
            If ``model_type`` is ``"SequentialRunner"``, but ``full_parameters``
            is not a list of dictionaries.
        TypeError
            If ``model_type`` is not ``"SequentialRunner"``, but ``full_parameters``
            is not a dictionary.

        """
        if model_type not in REGISTRY:
            msg = f"Model type {model_type} not in model registry. Available model types: {REGISTRY.keys()}"
            raise ValueError(msg)
        self.model_type = model_type
        self.full_parameters = full_parameters

        os_cpu_count = os.cpu_count()
        if os_cpu_count is None:
            warnings.warn(
                f"Could not determine number of physical CPU cores. Cannot determine whether the number of workers `max_workers` ({max_workers}) is larger than that or not. Using it directly",
                stacklevel=2,
            )
            self.max_workers = max_workers - 1
        else:
            self.max_workers = (
                max_workers - 1
                if (max_workers < os_cpu_count and os_cpu_count > 32)
                else os_cpu_count - 1
                if os_cpu_count > 32
                else int(os_cpu_count / 2) - 1
            )
        self.grid_file = (
            Path(grid_file) if ".h5" in str(grid_file) else Path(str(grid_file) + ".h5")
        )
        # TODO: Implement model appending to grid file
        # TODO: Implement option to append or overwrite grid file.
        # Initial placeholder statement to remove pre-existing grid files
        if self.grid_file.is_file():
            self.grid_file.unlink()
        self.model_name_prefix = model_name_prefix
        self.overwrite_models = overwrite_models
        self.array_dtype = array_dtype
        if log_dir is not None:
            self.log_dir: Path | None = Path(log_dir)
            self.log_dir.mkdir(exist_ok=True, parents=True)
            self._main_log: Path | None = self.log_dir / "grid.log"
        else:
            self.log_dir = None
            self._main_log = None
        self._orig_sigint = signal.getsignal(signal.SIGINT)
        self.parameters_to_grid: dict[str, np.ndarray] = {}

        if self.model_type == "SequentialRunner":
            if not isinstance(self.full_parameters, list):
                msg = f"For SequentialRunner types, full_parameters must be a list. {type(self.full_parameters)} was passed."
                raise TypeError(msg)
            for model_count in range(len(self.full_parameters)):
                for model_full_params in self.full_parameters[model_count].values():
                    if not isinstance(model_full_params, dict):
                        continue
                    for k, v in model_full_params.items():
                        if k == "param_dict":
                            for k_p, v_p in v.items():
                                self._grid_def(k_p, v_p, model_count)
                        else:
                            self._grid_def(k, v, model_count)
        else:
            if not isinstance(self.full_parameters, dict):
                msg = f"For non-SequentialRunner types, full_parameters must be a dictionary. {type(self.full_parameters)} was passed."
                raise TypeError(msg)
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
            if len({len(v) for v in self.parameters_to_grid.values()}) != 1:
                msg = "Not all passed arrays had the same length"
                raise ValueError(msg)
            self.flat_grids = np.array(
                [
                    [p[i] for p in self.parameters_to_grid.values()]
                    for i in range(len(next(iter(self.parameters_to_grid.values()))))
                ]
            ).T
        logger.debug(f"Created grid with {len(self.flat_grids)} gridpoints")

        # Optional subset of model indices (0-based column in flat_grids) to run.
        # None means run all. Accepts any iterable; stored as a frozenset for O(1) lookup.
        self.model_ids = frozenset(model_ids) if model_ids is not None else None

        self.model_id_dict: dict[int, str] = {}
        self.models: list[dict[str, Any]] = []
        self.physics_values = None
        self.chemical_abun_values = None
        if not delay_run:
            self.run()

    def _grid_def(self, key: str, value: Any, model_count: int | None = None) -> None:
        """Define a grid.

        Parameters
        ----------
        key : str
            name of parameter.
        value : Any
            value of parameter.
        model_count : int | None
            count of model. If not None,
            prepend the key ``key`` in ``self.parameters_to_grid`` with ``str(model_count)_``.
            Default = None.

        """
        if model_count is None:
            model_count_string = ""
        else:
            model_count_string = f"{str(model_count)}_"
        if isinstance(value, list) and key not in NoGridParameters:
            self.parameters_to_grid[model_count_string + key] = np.array(
                value, dtype=object
            )
        elif isinstance(value, (np.ndarray, np.generic)) and key not in NoGridParameters:
            self.parameters_to_grid[model_count_string + key] = value.astype(dtype=object)

    def _log_main(self, msg: str) -> None:
        """Append a timestamped line to the main grid log file.

        Parameters
        ----------
        msg : str
            message to log

        """
        if self._main_log is None:
            return
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with self._main_log.open("a") as f:
            f.write(f"{ts} {msg}\n")

    def run(self) -> None:
        """Run the grid."""

        def _handler(signum: Any, frame: Any) -> None:  # noqa: ARG001
            """Handle a raised exception.

            Parameters
            ----------
            signum : Any
                Not used
            frame : Any
                Not used.

            Raises
            ------
            KeyboardInterrupt
                Always.

            """
            try:
                self.on_interrupt()  # your “final steps”
            finally:
                # Restore default and re-raise KeyboardInterrupt to stop execution
                signal.signal(signal.SIGINT, self._orig_sigint)
                raise KeyboardInterrupt

        signal.signal(signal.SIGINT, _handler)
        n_total = np.shape(self.flat_grids)[1]
        self._log_main(f"Grid started: {n_total} models, {self.max_workers} workers")

        # Capture advanced settings so spawned workers start with the same
        # Fortran module state as the coordinator process.

        snapshot = create_snapshot()

        pending = self.grid_iter(
            self.full_parameters,
            list(self.parameters_to_grid.keys()),
            self.flat_grids,
            self.model_type,
        )
        file_obj = h5py.File(self.grid_file, "a")

        PoolClass = NoDaemonPool if self.model_type == "SequentialRunner" else mp.Pool  # noqa: N806

        with PoolClass(
            self.max_workers,
            initializer=_pool_initializer,
            initargs=(snapshot,),
            maxtasksperchild=1,
        ) as pool:
            completed = 0

            def on_result(result: tuple[int, AbstractModel]) -> None:
                """Log a successful model, and save it.

                Parameters
                ----------
                result : tuple[int, AbstractModel]
                    tuple of model index
                    and the AbstractModel object.

                """
                nonlocal completed
                completed += 1
                model_id, model_object = result
                try:
                    save_name = f"{self.model_name_prefix}{model_id}"
                    model_object.un_pickle()
                    if model_object.success_flag is not None:
                        if model_object.success_flag != SuccessFlag.SUCCESS:
                            msg = model_object.success_flag.check_error(
                                raise_on_error=False
                            )
                            self._log_main(
                                f"model_{model_id} failed ({completed}/{n_total}): {msg}"
                            )
                        else:
                            self._log_main(
                                f"model_{model_id} completed ({completed}/{n_total})"
                            )
                    else:
                        self._log_main(
                            f"model_{model_id} did not complete ({completed}/{n_total})"
                        )
                    model_object.save_model(
                        file=file_obj,
                        name=save_name,
                        overwrite=self.overwrite_models,
                        array_dtype=self.array_dtype,
                    )
                    self.model_id_dict[model_id] = save_name
                    file_obj.flush()
                except Exception as e:
                    print(f"Error saving model {model_id}: {e}")

                    traceback.print_exc()

            def on_error(_exc: Any, _model_id: int) -> None:
                """Log when an error occurs.

                Parameters
                ----------
                _exc : Any
                    Error
                _model_id : int
                    Model index

                """
                self._log_main(f"model_{_model_id} error: {_exc}")
                print(f"error: {_exc}; for model: {_model_id}")

            # Submit all tasks upfront; the pool limits actual concurrency to
            # max_workers. pool.close()/join() are called only after all tasks
            # have been submitted, so apply_async never races with a closed pool.
            for pending_model in pending:
                model_id = pending_model.pop("id")
                if self.model_ids is not None and model_id not in self.model_ids:
                    logger.debug(f"Skipping model id {model_id}, not in self.model_ids")
                    continue
                self._log_main(f"model_{model_id} started")
                pool.apply_async(
                    _run_grid_model,
                    args=(model_id, self.model_type, pending_model, self.log_dir),
                    callback=on_result,
                    error_callback=lambda exc, _mid=model_id: on_error(exc, _mid),  # type: ignore[misc]
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

        Raises
        ------
        NotImplementedError
            If the model type is ``SequentialRunner``.
        RuntimeError
            If no models were run yet.

        """
        if self.model_type == "SequentialRunner":
            msg = "SequentialRunner physics loading not implemented"
            raise NotImplementedError(msg)
        if not self.models:
            msg = "No models were run yet, so cannot load their physics arrays."
            raise RuntimeError(msg)

        for model_idx in range(len(self.models)):
            logger.debug(f"Loading physics of model with index {model_idx}")
            loaded_data = self._load_model_data(model=self.models[model_idx]["Model"])
            if self.physics_values is None:
                self.physics_values = json.loads(loaded_data["attributes_dict"].item())[
                    "physics_values"
                ]
            loaded_data = loaded_data.assign_coords(
                {"physics_values": self.physics_values}
            )
            self.models[model_idx]["physics_array"] = loaded_data["physics_array"]

    def load_chemistry(self, out_species_list: list[str]) -> None:
        """Load the chemistry.

        Parameters
        ----------
        out_species_list : list[str]
            list of species to load abundances for.

        Raises
        ------
        NotImplementedError
            If the model type is ``SequentialRunner``.
        RuntimeError
            If no models were run yet.

        """
        if self.model_type == "SequentialRunner":
            msg = "SequentialRunner chemistry loading not implemented"
            raise NotImplementedError(msg)
        if not self.models:
            msg = "No models were run yet, so cannot load their chemistry arrays."
            raise RuntimeError(msg)

        for model_idx in range(len(self.models)):
            logger.debug(f"Loading chemistry of model with index {model_idx}")
            loaded_data = self._load_model_data(model=self.models[model_idx]["Model"])
            if self.chemical_abun_values is None:
                self.chemical_abun_values = json.loads(
                    loaded_data["attributes_dict"].item()
                )["chemical_abun_values"]
            loaded_data = loaded_data.assign_coords(
                {"chemical_abun_values": self.chemical_abun_values}
            )
            self.models[model_idx]["out_species_abundances_array"] = loaded_data[
                "chemical_abun_array"
            ].sel(chemical_abun_values=out_species_list)

    def _load_params(self) -> None:
        """Loop through the models present in ``self.models`` to load the changing physical parameters.

        The method splits the loops into two cases: ``SequentialRunner``, and other model cases.
        In both instances, the for loop loads the model data using the ``_load_model``
        method. Then, it matches the given parameters, with the changing parameters in
        order to populate the dictionary attribute ``self.models`` for users to view the
        models that were run for the given ``GridRunner`` object.

        The differentiation of ``SequentialRunner`` stems from ``SequentialRunner``
        instances being nested models, resulting in an additional required loop to take into
        account the additional nesting.

        """
        # The following loops perform the same actions, but for the different model types
        if self.model_type == "SequentialRunner":
            for model_idx in range(len(self.models)):
                for model_count in range(len(self.full_parameters)):
                    for mt_k, mt_v in self.full_parameters[model_count].items():
                        if not isinstance(mt_v, dict):
                            continue
                        tmp_model = self._load_model(
                            model=f"{self.models[model_idx]['Model']}_{model_count}_{mt_k}"
                        )

                        self.models[model_idx][f"{model_count}_{mt_k}"] = {
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
                                k.replace(f"{model_count}_", ""): tmp_model.__getattr__(  # noqa: PLC2801
                                    k.replace(f"{model_count}_", "")
                                )
                                for k in list(self.parameters_to_grid.keys())
                                if mt_k in k
                                and k.replace(mt_k, "").lower() in tmp_model._data
                            },
                        }
                        self.models[model_idx][f"{model_count}_{mt_k}"]["Successful"] = (
                            True
                            if tmp_model.success_flag == 0
                            else tmp_model.success_flag
                        )
        else:
            for model_idx in range(len(self.models)):
                loaded_data = self._load_model(model=self.models[model_idx]["Model"])
                loaded_dict = loaded_data._param_dict
                self.models[model_idx] = {
                    **self.models[model_idx],
                    **{
                        k: loaded_dict[k.lower()]
                        for k in list(self.parameters_to_grid.keys())
                    },
                }
                self.models[model_idx]["Successful"] = (
                    True if loaded_data.success_flag == 0 else loaded_data.success_flag
                )

    def _load_model(self, model: str) -> AbstractModel:
        """Load an entire model from a file within the grid file.

        Parameters
        ----------
        model : str
            Name of the model in ``self.grid_file``.

        Returns
        -------
        model_object : AbstractModel
            Loaded model.

        """
        logger.debug(f"Loading model in file {self.grid_file} with name {model}")
        model_object = load_model(file=self.grid_file, name=model)
        return model_object

    def _load_model_data(self, model: str) -> xr.Dataset:
        """Load the ``_data`` attribute of a model, and get a copy.

        Parameters
        ----------
        model : str
            Name of the model in ``self.grid_file``.

        Returns
        -------
        xr.Dataset
            Copy of loaded model's ``_data`` attribute.

        """
        return self._load_model(model)._data.copy()

    def check_conservation(self, element_list: list[str] | None = None) -> None:
        """Check conservation of the chemical abundances.

        Adds an entry in the self.models[model_idx] dictionary with key
            ``"elements_conserved"``, which indicates whether all the points in model with
            index ``model_idx`` had changes in elemental abundances below 1%.

        Parameters
        ----------
        element_list : list[str] | None
            List of elements to check conservation for.
            If None, use ``uclchem.constants.default_elements_to_check``. Default = None.

        """
        if element_list is None:
            element_list = default_elements_to_check
        for model_idx in range(len(self.models)):
            tmp_model = self._load_model(self.models[model_idx]["Model"])
            conserve_dicts = []
            if tmp_model._param_dict["points"] > 1:
                for i in range(tmp_model._param_dict["points"]):
                    conserve_dicts += [
                        check_element_conservation(
                            tmp_model.get_dataframes(i)[1],
                            element_list=element_list,
                            percent=True,
                        )
                    ]
            else:
                conserve_dicts += [
                    check_element_conservation(
                        tmp_model.get_dataframes(0)[1],
                        element_list=element_list,
                        percent=True,
                    )
                ]
            conserved = True
            for conserve_dict in conserve_dicts:
                conserved = all(float(x[:1]) < 1 for x in conserve_dict.values())
            self.models[model_idx]["elements_conserved"] = conserved

    def on_interrupt(self) -> None:  # noqa: PLR6301
        """Catch interruption. Does nothing."""
        return

    @staticmethod
    def grid_iter(
        full_parameters: dict | list[dict],
        param_keys: list[str],
        flattened_grids: np.ndarray,
        model_type: str,
    ) -> Iterator[dict[str, Any]]:
        """Iterate over dictionaries of parameters for each model in the grid.

        Provide an iterable dictionary of parameters that can be used with the
        grid-based multiprocessing worker distribution.

        Parameters
        ----------
        full_parameters : dict | list[dict]
            dictionary or list
            (if model_type == SequentialRunner)
            of the full parameters that will be used for the model
        param_keys : list[str]
            list of parameters that are changing in
            this GridRunner object
        flattened_grids : np.ndarray
            list of all the values for the changing parameters
            for this GridRunner object.
        model_type : str
            Type of model to use. 'SequentialRunner' results in
            alternative way of executing this phase as each SequentialRunner
            represents multiple models to be run in series.

        Yields
        ------
        dict[str, Any]
            Next dictionary containing the parameter values for the model to run
            in the grid of models. Only offers one model per request.

        Raises
        ------
        TypeError
            If ``model_type`` is ``"SequentialRunner"``, but ``full_parameters``
            is not a list of dictionaries.

        """
        if model_type == "SequentialRunner":
            if not isinstance(full_parameters, list):
                msg = "If model_type == 'SequentialRunner', full_parameters needs to be a list of dictionaries,"
                msg += f" but got type {type(full_parameters)}"
                raise TypeError(msg)
            # As SequentialRunner types contain multiple models,
            # sometimes of the same type, we split this type out to follow altered
            # logic to arrive at equivalently expected outputs.
            for i in range(len(flattened_grids[0])):
                combo: tuple[Any, ...] = ()  # type: ignore[assignment]
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
                    for current_model_type, model_full_parameters in full_parameters[
                        model_count
                    ].items():
                        if isinstance(model_full_parameters, dict):
                            # grid_param_dict is filled with the param_dict values of a model.
                            grid_param_dict = {
                                k.replace(f"{model_count}_", ""): v
                                for k, v in zip(param_keys, combo, strict=True)
                                if k.replace(f"{model_count}_", "")
                                in model_full_parameters["param_dict"]
                                and k[: len(str(model_count))] == str(model_count)
                            }
                            # grid_dict is filled with the input parameters of a value,
                            # not part of param_dict
                            grid_dict = {
                                k.replace(f"{model_count}_", ""): v
                                for k, v in zip(param_keys, combo, strict=True)
                                if k.replace(f"{model_count}_", "")
                                not in model_full_parameters["param_dict"]
                                and (
                                    k[: len(str(model_count))] == str(model_count)
                                    if k[: len(str(model_count))].isdigit()
                                    else False
                                )
                            }
                            run_dict[current_model_type] = {
                                **model_full_parameters,
                                "param_dict": {
                                    **model_full_parameters["param_dict"],
                                    **grid_param_dict,
                                },
                                **grid_dict,
                            }
                            run_list += [run_dict]
                        else:
                            yield_dict[current_model_type] = model_full_parameters
                yield {
                    "parameters_to_match": ["finalDens"],
                    **yield_dict,
                    **{"sequenced_model_parameters": run_list},
                }
        else:
            if not isinstance(full_parameters, dict):
                msg = "If model_type != 'SequentialRunner', full_parameters needs to be a dictionary,"
                msg += f" but got type {type(full_parameters)}"
                raise TypeError(msg)
            for i in range(len(flattened_grids[0])):
                combo = ()
                for j in range(np.shape(flattened_grids)[0]):
                    combo += (flattened_grids[j][i],)

                # grid_param_dict contains the param_dict values of the next model to run.
                grid_param_dict = {
                    k: v
                    if not isinstance(v, float)
                    else (v.item() if hasattr(v, "item") else v)
                    for k, v in zip(param_keys, combo, strict=True)
                    if k in full_parameters["param_dict"]
                }
                # grid_dict is filled with the input parameters of a value,
                # not part of param_dict, for the next model to run
                grid_dict = {
                    k: v
                    for k, v in zip(param_keys, combo, strict=True)
                    if k not in full_parameters["param_dict"]
                }
                yield {
                    **full_parameters,
                    "param_dict": {**full_parameters["param_dict"], **grid_param_dict},
                    **grid_dict,
                    "id": i,
                }
