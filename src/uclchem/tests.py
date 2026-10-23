"""Collection of tests for UCLCHEM.

Deprecated

"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


import uclchem
from uclchem.constants import ELEMENTS_TO_CHECK


def test_ode_conservation(
    element_list: Iterable[str] = ELEMENTS_TO_CHECK,  # ruff: ignore[pytest-parameter-with-default-argument]
) -> dict[str, float]:
    """Test whether the ODEs conserve elements.

    Useful to run each time you change network. Integrator errors may still cause
    elements not to be conserved but they cannot be conserved if the ODEs are not correct.

    Parameters
    ----------
    element_list : Iterable[str]
        A list of elements for which to check the conservation.
        Default = ``uclchem.constants.ELEMENTS_TO_CHECK``.

    Returns
    -------
    result : dict[str, float]
        A dictionary of the elements in element list with values
        representing the total change of each element.

    """
    param_dict = {
        "endatfinaldensity": False,
        "freefall": True,
        "initialdens": 1e4,
        "initialtemp": 10.0,
        "finaldens": 1e5,
        "finaltime": 1.0e3,
    }
    model = uclchem.model.Cloud(param_dict)
    model.check_error()

    _physics_df, abundances_df = model.get_split_dataframes()
    result = uclchem.analysis.check_element_conservation(
        abundances_df,
        element_list=element_list,
        percent=False,
    )
    result_floats = {
        element: float(discrepancy) for element, discrepancy in result.items()
    }
    return result_floats
