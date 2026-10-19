"""Guard against requesting starting abundances without supplying them.

Before PR #215, chained models set givestartabund=True but never passed the
abundance array to Fortran, so every species silently started at the
abundance floor. These tests make sure that situation now fails loudly with
MISSING_START_ABUNDANCE_ERROR instead of integrating an empty gas.
"""

import pytest

import uclchem
from uclchem.utils import SuccessFlag

PARAMS = {"finalTime": 1e2}
TIMEPOINTS = 200


@pytest.fixture(scope="module")
def first_stage():
    model = uclchem.model.Cloud(param_dict=PARAMS, timepoints=TIMEPOINTS)
    model.check_error()
    return model


def _second_stage(**kwargs):
    return uclchem.model.Cloud(
        param_dict=PARAMS, timepoints=TIMEPOINTS, run_type="external", **kwargs
    )


def test_chained_model_with_starting_abundances_succeeds(first_stage):
    model = _second_stage(previous_model=first_stage)
    model.run()
    assert model.success_flag == SuccessFlag.SUCCESS


def test_start_abundances_requested_but_not_passed_fails():
    model = _second_stage()
    assert model.starting_chemistry_array is None
    model.give_start_abund = True  # the pre-#215 state: flag set, no array passed
    with pytest.raises(RuntimeError, match="MISSING_START_ABUNDANCE_ERROR"):
        model.run()


def test_all_zero_starting_abundances_fail(first_stage):
    model = _second_stage(previous_model=first_stage)
    model.starting_chemistry_array[0, :] = 0.0
    with pytest.raises(RuntimeError, match="MISSING_START_ABUNDANCE_ERROR"):
        model.run()
