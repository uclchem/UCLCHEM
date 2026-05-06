"""Test that tempa and tempb are correctly set for both hotcore model variants.

tempa/tempb live in the Fortran hotcore module and are set during initializePhysics.
Using run_type="external" prevents PrestellarCore from auto-running in a managed
subprocess; calling model.run() manually then executes run_fortran() in-process,
making the Fortran state readable directly after the call.

Setup: 1e3 year cloud phase followed by a 1e2 year hotcore phase.
"""

import importlib.util

import numpy as np
import pytest

uclchem_imported = importlib.util.find_spec("uclchem") is not None
pytestmark = pytest.mark.skipif(not uclchem_imported, reason="uclchem not installed")

# Expected coefficient arrays from hotcore.f90 initializePhysics
TEMPA_0D = np.array([1.927e-1, 4.8560e-2, 7.8470e-3, 9.6966e-4, 1.706e-4, 4.74e-7])
TEMPB_0D = np.array([0.5339, 0.6255, 0.8395, 1.085, 1.289, 1.98])

TEMPA_1D = np.array([3.1417e-2, 3.5495e-2, 4.9653e-4, 9.5928e-4, 1.4158e-3, 2.817e-3])
TEMPB_1D = np.array([0.5329, 0.5324, 0.9, 0.9, 0.9, 0.9])

BASE_PARAMS = {
    "initialDens": 1e4,
    "initialTemp": 10.0,
    "finalDens": 1e4,
    "writeStep": 1,
}


def _read_coefficients(extra_params: dict) -> tuple[np.ndarray, np.ndarray]:
    """Run a cloud + hotcore in-process and return (tempa, tempb).

    PrestellarCore is constructed with run_type="external" so it does not
    auto-spawn a managed subprocess. Calling model.run() then executes
    run_fortran() in-process, leaving the Fortran hotcore module state
    readable in the current process.
    """
    from uclchemwrap import hotcore as hc

    import uclchem

    params = {**BASE_PARAMS, **extra_params}

    cloud = uclchem.model.Cloud(
        param_dict={**params, "finalTime": 1e3}, out_species=["CO"]
    )
    cloud.check_error()

    core = uclchem.model.PrestellarCore(
        temp_indx=1,
        max_temperature=300.0,
        param_dict={**params, "finalTime": 1e2},
        out_species=["CO"],
        previous_model=cloud,
        run_type="external",
    )
    core.run()

    return np.array(hc.tempa), np.array(hc.tempb)


class TestTempCoefficients0D:
    """0D hotcore (enable_radiative_transfer=False, points=1)."""

    @pytest.fixture(scope="class")
    def coefficients(self):
        return _read_coefficients({"enable_radiative_transfer": False, "points": 1})

    def test_tempa_0d(self, coefficients):
        tempa, _ = coefficients
        np.testing.assert_allclose(
            tempa, TEMPA_0D, rtol=1e-4, err_msg="tempa does not match expected 0D values"
        )

    def test_tempb_0d(self, coefficients):
        _, tempb = coefficients
        np.testing.assert_allclose(
            tempb, TEMPB_0D, rtol=1e-4, err_msg="tempb does not match expected 0D values"
        )


class TestTempCoefficients1D:
    """1D hotcore (enable_radiative_transfer=True, points=2)."""

    @pytest.fixture(scope="class")
    def coefficients(self):
        return _read_coefficients(
            {
                "enable_radiative_transfer": True,
                "points": 2,
                "rout": 0.05,
                "density_scale_radius": 0.05,
                "density_power_index": 2.4,
                "lum_star": 1e6,
                "temp_star": 4.5e4,
            }
        )

    def test_tempa_1d(self, coefficients):
        tempa, _ = coefficients
        np.testing.assert_allclose(
            tempa, TEMPA_1D, rtol=1e-4, err_msg="tempa does not match expected 1D values"
        )

    def test_tempb_1d(self, coefficients):
        _, tempb = coefficients
        np.testing.assert_allclose(
            tempb, TEMPB_1D, rtol=1e-4, err_msg="tempb does not match expected 1D values"
        )
