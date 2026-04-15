"""Collection of tests for UCLCHEM.

Deprecated
"""

import numpy as np

import uclchem
from uclchem.constants import default_elements_to_check
from uclchem.utils import UCLCHEM_ROOT_DIR


def test_ode_conservation(
    element_list: list[str] | None = None,
) -> dict[str, str]:
    """Test whether the ODEs conserve elements. Useful to run each time you change network.
    Integrator errors may still cause elements not to be conserved but they cannot be conserved
    if the ODEs are not correct.

    Args:
        element_list (list[str]): A list of elements for which to check the conservation.
            If None, use ``uclchem.constants.default_elements_to_check``. Default = None.

    Returns:
        result (dict[str, str]): A dictionary of the elements in element list with values
            representing the total rate of change of each element.

    """
    if element_list is None:
        element_list = default_elements_to_check

    species_array = np.loadtxt(
        UCLCHEM_ROOT_DIR / "species.csv",
        usecols=[0],
        dtype=str,
        skiprows=1,
        unpack=True,
        delimiter=",",
        comments="%",
    )
    species_list = list(species_array)
    param_dict = {
        "endatfinaldensity": False,
        "freefall": True,
        "initialdens": 1e4,
        "initialtemp": 10.0,
        "finaldens": 1e5,
        "finaltime": 1.0e3,
        "outspecies": len(species_list),
    }
    model = uclchem.model.Cloud(param_dict, out_species=species_list)
    model.check_error()

    physics_df, abundances_df = model.get_dataframes()
    result = uclchem.analysis.check_element_conservation(
        abundances_df,  # type: ignore[arg-type]
        element_list=element_list,
    )
    return result
