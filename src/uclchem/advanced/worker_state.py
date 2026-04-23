"""Snapshot and restore of advanced settings for multiprocessing propagation.

When UCLCHEM runs grid models or managed-mode models, worker processes are
spawned via ``mp.Pool`` or ``mp.Process`` with the ``spawn`` context.  These
fresh processes import ``uclchemwrap`` from scratch, losing any runtime
modifications made through ``GeneralSettings``, ``HeatingSettings``, or
``NetworkState``.

This module provides :func:`create_snapshot` / :func:`restore_snapshot` to
capture the current Fortran module state into a picklable dict and re-apply
it in a worker process before the model runs.
"""

import contextlib
import logging
from typing import Any

import numpy as np
import uclchemwrap
from uclchemwrap import f2py_constants as f2py_constants_module
from uclchemwrap import heating as heating_module
from uclchemwrap import network as network_module

from uclchem.advanced.constants import (
    FILE_PATH_PARAMETERS,
    FORTRAN_PARAMETERS,
    INTERNAL_PARAMETERS,
)

logger = logging.getLogger(__name__)

# Module names mirroring GeneralSettings._discover_modules()
_MODULE_NAMES = [
    "defaultparameters",
    "network",
    "heating",
    "physicscore",
    "constants",
    "cloud_mod",
    "collapse_mod",
    "cshock_mod",
    "jshock_mod",
    "hotcore",
    "chemistry",
    "rates",
    "photoreactions",
    "surfacereactions",
    "io",
    "f2py_constants",
    "postprocess_mod",
    "sputtering",
]

# Modules where *all* 0-d array attributes are Fortran PARAMETERs (compile-time
# constants compiled into read-only pages).  Their runtime state is either
# handled by dedicated snapshot sections (e.g. f2py_constants → heating section)
# or never needs propagation to workers.  Attempting setattr on these in a
# freshly-spawned worker causes SIGBUS on macOS.
_MODULES_SKIP_0D = frozenset(
    {
        "constants",  # physical constants (c, k_boltz, …) – all PARAMETERs
        "f2py_constants",  # build-time counts (nspec, nReac, …) – all PARAMETERs
        "surfacereactions",  # grain/surface constants – all PARAMETERs
    }
)


def create_snapshot() -> dict[str, Any]:
    """Capture the current Fortran module state into a picklable dict.

    Reads directly from Fortran memory (not cached Python values) to ensure
    accuracy even when settings were modified outside the wrapper classes.

    The returned dict has three sections:

    * ``"general"`` – scalar settings from all uclchemwrap sub-modules
      (excluding PARAMETERs, INTERNAL, FILE_PATH, and arrays).
    * ``"heating"`` – heating/cooling boolean arrays, scalars, and coolant
      configuration.
    * ``"network"`` – reaction-rate and binding-energy arrays.

    Returns:
        dict[str, Any]: Fully picklable dict suitable for passing to :func:`restore_snapshot`.

    """
    logger.debug("Creating snapshot")

    snapshot: dict[str, Any] = {}

    # --- General settings (scalars only) ---
    general: dict[str, dict[str, Any]] = {}
    for mod_name in _MODULE_NAMES:
        if not hasattr(uclchemwrap, mod_name):
            continue
        mod = getattr(uclchemwrap, mod_name)
        mod_snapshot: dict[str, Any] = {}
        for attr in dir(mod):
            if attr.startswith("_"):
                continue
            try:
                value = getattr(mod, attr)
            except Exception as e:
                logger.exception(
                    f"Exception occurred when trying to get attribute {attr} from module {mod_name}:\n",
                    e,
                )
                continue
            if callable(value):
                continue
            # Skip arrays (handled by heating/network sections, or immutable),
            # but keep 0-d arrays (f2py scalar wrappers like cloud_mod scalars).
            if isinstance(value, np.ndarray) and value.ndim > 0:
                continue
            # Skip 0-d arrays for modules that are entirely Fortran PARAMETERs —
            # these are read-only in spawned workers and cause SIGBUS on macOS.
            if (
                isinstance(value, np.ndarray)
                and value.ndim == 0
                and mod_name in _MODULES_SKIP_0D
            ):
                continue
            # Skip parameters that cannot or should not be set
            attr_lower = attr.lower()
            if attr_lower in FORTRAN_PARAMETERS:
                continue
            if attr_lower in INTERNAL_PARAMETERS:
                continue
            if attr_lower in FILE_PATH_PARAMETERS:
                continue
            # Convert 0-d numpy arrays to Python scalars for clean pickling
            if isinstance(value, np.ndarray) and value.ndim == 0:
                value = value.item()
            mod_snapshot[attr] = value
        if mod_snapshot:
            general[mod_name] = mod_snapshot
    snapshot["general"] = general

    # --- Heating / cooling settings ---
    heating: dict[str, Any] = {
        "heating_modules": np.copy(heating_module.heating_modules),
        "cooling_modules": np.copy(heating_module.cooling_modules),
        "dust_gas_coupling_method": int(heating_module.dust_gas_coupling_method),
        "line_solver_attempts": int(heating_module.line_solver_attempts),
        "pahabund": float(heating_module.pahabund),
        "coolantdatadir": np.copy(f2py_constants_module.coolantdatadir),
        "coolant_active": np.copy(f2py_constants_module.coolant_active),
    }

    # Coolant restart mode – accessor pattern varies between builds
    if hasattr(uclchemwrap, "get_coolant_restart_mode_wrap"):
        heating["coolant_restart_mode"] = int(uclchemwrap.get_coolant_restart_mode_wrap())
    elif hasattr(
        getattr(uclchemwrap, "uclchemwrap", None), "get_coolant_restart_mode_wrap"
    ):
        heating["coolant_restart_mode"] = int(
            uclchemwrap.uclchemwrap.get_coolant_restart_mode_wrap()
        )
    snapshot["heating"] = heating

    # --- Network state (rate parameters + binding energies) ---
    snapshot["network"] = {
        "alpha": np.copy(network_module.alpha),
        "beta": np.copy(network_module.beta),
        "gama": np.copy(network_module.gama),
        "bindingenergy": np.copy(network_module.bindingenergy),
    }

    return snapshot


def restore_snapshot(snapshot: dict[str, Any]) -> None:
    """Apply a previously captured snapshot to the current process.

    Must be called **before** running any model in the worker process.

    Args:
        snapshot (dict[str, Any]): Dict produced by :func:`create_snapshot`.

    """
    logger.debug("Regenerating snapshot")
    # --- General settings ---
    # If uclchem hangs here, the last debug line printed shows which Fortran
    # PARAMETER is blocking. Add it to src/uclchem/advanced/fortran_metadata.yaml
    # and regenerate with: uclchem-generate-metadata
    for mod_name, settings_dict in snapshot.get("general", {}).items():
        if not hasattr(uclchemwrap, mod_name):
            continue
        mod = getattr(uclchemwrap, mod_name)
        for attr, value in settings_dict.items():
            # Uncomment next line to debug hangs (last printed line is the blocker):
            # print(f"[DEBUG] setattr({mod_name}, {attr}, {value!r})", flush=True, file=sys.stderr) # noqa: ERA001
            with contextlib.suppress(AttributeError, TypeError):
                # read-only or incompatible – skip silently
                setattr(mod, attr, value)

    # --- Heating / cooling settings ---
    heating = snapshot.get("heating", {})
    if "heating_modules" in heating:
        heating_module.heating_modules[:] = heating["heating_modules"]
    if "cooling_modules" in heating:
        heating_module.cooling_modules[:] = heating["cooling_modules"]
    if "dust_gas_coupling_method" in heating:
        heating_module.dust_gas_coupling_method = heating["dust_gas_coupling_method"]
    if "line_solver_attempts" in heating:
        heating_module.line_solver_attempts = heating["line_solver_attempts"]
    if "pahabund" in heating:
        heating_module.pahabund = heating["pahabund"]
    if "coolantdatadir" in heating:
        f2py_constants_module.coolantdatadir = heating["coolantdatadir"]
    if "coolant_active" in heating:
        f2py_constants_module.coolant_active[:] = heating["coolant_active"]
    if "coolant_restart_mode" in heating:
        mode = heating["coolant_restart_mode"]
        if hasattr(uclchemwrap, "set_coolant_restart_mode_wrap"):
            uclchemwrap.set_coolant_restart_mode_wrap(mode)
        elif hasattr(
            getattr(uclchemwrap, "uclchemwrap", None), "set_coolant_restart_mode_wrap"
        ):
            uclchemwrap.uclchemwrap.set_coolant_restart_mode_wrap(mode)

    # --- Network state ---
    net = snapshot.get("network", {})
    if "alpha" in net:
        np.copyto(network_module.alpha, net["alpha"])
    if "beta" in net:
        np.copyto(network_module.beta, net["beta"])
    if "gama" in net:
        np.copyto(network_module.gama, net["gama"])
    if "bindingenergy" in net:
        np.copyto(network_module.bindingenergy, net["bindingenergy"])


def _pool_initializer(snapshot: dict[str, Any]) -> None:
    """``mp.Pool`` initializer that restores advanced settings in each worker.

    Args:
        snapshot (dict[str, Any]): Snapshot created by func:`create_snapshot`.

    Example:
        >>> import multiprocessing as mp
        >>>
        >>> # Take a snapshot of the current Fortran module
        >>> snapshot = create_snapshot()
        >>>
        >>> # A pool can then be initialized as
        >>> n_workers = 2
        >>> mp.Pool(n_workers, initializer=_pool_initializer, initargs=(snapshot,))

    """
    restore_snapshot(snapshot)
