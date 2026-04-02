import pytest

from uclchem.advanced.advanced_settings import GeneralSettings
from uclchem.utils import SuccessFlag


def test_match_fortran_python_success_flags_match():
    settings = GeneralSettings()
    error_flags_dict = settings.search(
        "_ERROR", include_internal=True, include_parameters=False
    )
    error_flags_dict = {
        k.upper().split(".")[1]: int(val.get())
        for k, val in error_flags_dict.items()
        if "constants" in k
    }
    error_flags_dict["SUCCESS"] = 0

    assert {e.name: e.value for e in SuccessFlag} == error_flags_dict


def test_check_error_success_flags():
    success_flag_dict = {e.name: e.value for e in SuccessFlag}
    success_flags = [SuccessFlag(val) for key, val in success_flag_dict.items()]

    for success_flag in success_flags:
        if success_flag.value < 0:
            with pytest.raises(RuntimeError):
                success_flag.check_error(raise_on_error=True)
            success_flag.check_error(raise_on_error=False)
        else:
            success_flag.check_error()
