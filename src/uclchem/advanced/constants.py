"""
Constants and metadata for UCLCHEM advanced module.

This module loads Fortran parameter classifications from YAML and provides
them as Python sets for use in GeneralSettings.
"""

from typing import Set

import yaml

from uclchem.utils import UCLCHEM_ROOT_DIR


def _load_fortran_metadata() -> tuple[Set[str], Set[str], Set[str]]:
    """Load Fortran parameter classifications from YAML file.

    Returns:
        Tuple of (fortran_parameters, internal_parameters, file_path_parameters) as lowercase sets
    """
    yaml_path = UCLCHEM_ROOT_DIR / "advanced" / "fortran_metadata.yaml"

    with open(yaml_path, "r") as f:
        metadata = yaml.safe_load(f)

    # Flatten nested structure and convert to lowercase sets
    fortran_params = set()
    for module_params in metadata.get("fortran_parameters", {}).values():
        if isinstance(module_params, list):
            fortran_params.update(p.lower() for p in module_params if isinstance(p, str))

    internal_params = set()
    for module_params in metadata.get("internal_parameters", {}).values():
        if isinstance(module_params, list):
            internal_params.update(p.lower() for p in module_params if isinstance(p, str))

    file_path_params = set()
    for module_params in metadata.get("file_path_parameters", {}).values():
        if isinstance(module_params, list):
            file_path_params.update(
                p.lower() for p in module_params if isinstance(p, str)
            )

    return fortran_params, internal_params, file_path_params


# Load metadata on module import
FORTRAN_PARAMETERS, INTERNAL_PARAMETERS, FILE_PATH_PARAMETERS = _load_fortran_metadata()

__all__ = ["FORTRAN_PARAMETERS", "INTERNAL_PARAMETERS", "FILE_PATH_PARAMETERS"]
