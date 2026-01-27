"""
Test multi-stage model runs with IN-MEMORY return modes only.

This test uses return_array and return_dataframe with starting_chemistry
to ensure all model stages work with Python in-memory arrays.
"""

import shutil
import tempfile
from pathlib import Path

import pytest

try:
    import uclchem

    uclchem_imported = True
except ImportError:
    uclchem_imported = False


def test_import_uclchem():
    if not uclchem_imported:
        pytest.fail(
            "uclchem module could not be imported, "
            "make sure your environment is loaded and UCLCHEM is installed."
        )


@pytest.fixture(scope="function")
def test_output_directory(request):
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


def test_static_model_return_array(test_output_directory):
    """Test static cloud model with return_array"""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
    }
    physics, chemistry, rates, heating, abundances_start, return_code = (
        uclchem.functional.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_array=True,
            return_rates=True,
        )
    )
    assert (
        return_code == 0
    ), f"Static model returned with nonzero exit code {return_code}"

    # Verify finalTime is respected (within 10% tolerance)
    max_time = physics[:, 0].max()  # Time is first column
    assert (
        max_time <= 1.1 * params["finalTime"]
    ), f"Model exceeded finalTime tolerance: {max_time:.2e} > {1.1 * params['finalTime']:.2e}"
    assert (
        max_time >= 0.9 * params["finalTime"]
    ), f"Model stopped too early: {max_time:.2e} < {0.9 * params['finalTime']:.2e}"


def test_static_model_return_dataframe(test_output_directory):
    """Test static cloud model with return_dataframe"""
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
    }
    physics, chemistry, rates, heating, abundances_start, return_code = (
        uclchem.functional.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_dataframe=True,
        )
    )
    assert (
        return_code == 0
    ), f"Static model returned with nonzero exit code {return_code}"

    # Verify finalTime is respected (within 10% tolerance)
    max_time = physics["Time"].max()
    assert (
        max_time <= 1.1 * params["finalTime"]
    ), f"Model exceeded finalTime tolerance: {max_time:.2e} > {1.1 * params['finalTime']:.2e}"
    assert (
        max_time >= 0.9 * params["finalTime"]
    ), f"Model stopped too early: {max_time:.2e} < {0.9 * params['finalTime']:.2e}"


def test_collapse_hotcore_return_array(test_output_directory):
    """Test collapse -> hot core chained models with return_array"""
    # Stage 1: Collapse with endAtFinalDensity=True
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "finalDens": 1e6,
        "finalTime": 1e5,
    }
    # return_array with return_rates=True returns 6 values: physics, chemistry, rates, heating(None), abundances, flag
    physics, chemistry, rates, heating, abundances_start, return_code = (
        uclchem.functional.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_array=True,
            return_rates=True,
        )
    )
    assert return_code == 0, f"Stage 1 returned with nonzero exit code {return_code}"

    # Verify endAtFinalDensity=True behavior: stops at finalTime OR finalDens
    max_time = physics[:, 0].max()
    # Check if physics array has density column (some models may not output it in all modes)
    if physics.shape[1] > 1:
        max_density = physics[:, 1].max()
        assert (
            max_time <= 1.1 * params["finalTime"]
            or max_density >= 0.9 * params["finalDens"]
        ), f"Collapse should stop at finalTime OR finalDens: time={max_time:.2e}, density={max_density:.2e}"
    else:
        # If no density column, just check time
        assert (
            max_time <= 1.1 * params["finalTime"]
        ), f"Collapse should stop at or before finalTime: time={max_time:.2e}, finalTime={params['finalTime']:.2e}"

    # Stage 2: Hot core using starting_chemistry
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    # prestellar_core with return_array also returns 6 values
    physics, chemistry, rates, heating, abundances_start, return_code = (
        uclchem.functional.prestellar_core(
            3,
            300.0,
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_array=True,
            starting_chemistry=abundances_start,
        )
    )
    assert return_code == 0, f"Stage 2 returned with nonzero exit code {return_code}"

    # Verify finalTime is respected (within 10% tolerance)
    max_time = physics[:, 0].max()
    assert (
        max_time <= 1.1 * params["finalTime"]
    ), f"Hot core exceeded finalTime tolerance: {max_time:.2e} > {1.1 * params['finalTime']:.2e}"
    assert (
        max_time >= 0.9 * params["finalTime"]
    ), f"Hot core stopped too early: {max_time:.2e} < {0.9 * params['finalTime']:.2e}"


def test_collapse_hotcore_return_dataframe(test_output_directory):
    """Test collapse -> hot core chained models with return_dataframe"""
    # Stage 1: Collapse with endAtFinalDensity=True
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "finalDens": 1e6,
        "finalTime": 1e5,
    }
    # return_dataframe returns 6 values: physics, chemistry, rates(None), heating(None), abundances, flag
    physics, chemistry, rates, heating, abundances_start, return_code = (
        uclchem.functional.cloud(
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_dataframe=True,
        )
    )
    assert return_code == 0, f"Stage 1 returned with nonzero exit code {return_code}"

    # Verify endAtFinalDensity=True behavior: stops at finalTime OR finalDens
    max_time = physics["Time"].max()
    max_density = physics["Density"].max()
    assert (
        max_time <= 1.1 * params["finalTime"]
        or max_density >= 0.9 * params["finalDens"]
    ), f"Collapse should stop at finalTime OR finalDens: time={max_time:.2e}, density={max_density:.2e}"

    # Stage 2: Hot core using starting_chemistry
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
    }
    # prestellar_core with return_dataframe also returns 6 values
    physics, chemistry, rates, heating, abundances_start, return_code = (
        uclchem.functional.prestellar_core(
            3,
            300.0,
            param_dict=params,
            out_species=["OH", "OCS", "CO", "CS", "CH3OH"],
            return_dataframe=True,
            starting_chemistry=abundances_start,
        )
    )
    assert return_code == 0, f"Stage 2 returned with nonzero exit code {return_code}"

    # Verify finalTime is respected (within 10% tolerance)
    max_time = physics["Time"].max()
    assert (
        max_time <= 1.1 * params["finalTime"]
    ), f"Hot core exceeded finalTime tolerance: {max_time:.2e} > {1.1 * params['finalTime']:.2e}"
    assert (
        max_time >= 0.9 * params["finalTime"]
    ), f"Hot core stopped too early: {max_time:.2e} < {0.9 * params['finalTime']:.2e}"


def test_cshock_return_dataframe(test_output_directory):
    """Test C-shock model with return_dataframe and starting_chemistry"""
    # Pre-shock cloud
    param_dict = {
        "endAtFinalDensity": False,
        "freefall": True,
        "initialDens": 1e2,
        "finalDens": 1e4,
        "initialTemp": 10.0,
        "finalTime": 6.0e6,
        "rout": 0.1,
        "baseAv": 1.0,
    }
    # return_dataframe returns 6 values: physics, chemistry, rates(None), heating(None), abundances, flag
    (
        df_stage1_physics,
        df_stage1_chemistry,
        rates,
        heating,
        final_abundances,
        return_code,
    ) = uclchem.functional.cloud(param_dict=param_dict, return_dataframe=True)
    assert (
        return_code == 0
    ), f"Pre-cshock cloud returned with nonzero exit code {return_code}"

    # Verify finalTime is respected (within 10% tolerance)
    max_time_stage1 = df_stage1_physics["Time"].max()
    assert (
        max_time_stage1 <= 1.1 * param_dict["finalTime"]
    ), f"Pre-shock cloud exceeded finalTime: {max_time_stage1:.2e} > {1.1 * param_dict['finalTime']:.2e}"

    # C-shock with starting_chemistry
    param_dict["initialDens"] = 1e4
    param_dict["finalTime"] = 1e6
    # CShock returns 7 values: physics, chemistry, rates(None), heating(None), dissipation_time, abundances, flag
    (
        df_stage2_physics,
        df_stage2_chemistry,
        rates,
        heating,
        dissipation_time,
        final_abundances2,
        return_code,
    ) = uclchem.functional.cshock(
        shock_vel=40,
        param_dict=param_dict,
        return_dataframe=True,
        starting_chemistry=final_abundances,
    )
    assert return_code == 0, f"C-shock returned with nonzero exit code {return_code}"

    # Verify finalTime is respected (shock models may stop early)
    max_time_stage2 = df_stage2_physics["Time"].max()
    assert (
        max_time_stage2 <= 1.1 * param_dict["finalTime"]
    ), f"C-shock exceeded finalTime: {max_time_stage2:.2e} > {1.1 * param_dict['finalTime']:.2e}"


def test_endAtFinalDensity_with_collapse(test_output_directory):
    """Test that endAtFinalDensity=True works correctly with Collapse models"""
    # Test Collapse model which should allow endAtFinalDensity without freefall
    # Create a unique output file for this test to avoid file conflicts
    output_file = str(test_output_directory / "collapse_test.dat")
    params = {
        "initialDens": 1e4,
        "finalDens": 1e6,
        "finalTime": 1e6,  # Longer finalTime for collapse models
        "initialTemp": 10.0,
        "endAtFinalDensity": True,
    }

    physics, chemistry, rates, heating, abundances, return_code = (
        uclchem.functional.collapse(
            collapse="BE4",  # BE collapse mode
            physics_output=output_file,
            param_dict=params,
            out_species=["OH", "CO"],
            return_array=True,
        )
    )
    assert (
        return_code == 0
    ), f"Collapse model returned with nonzero exit code {return_code}"

    # Should stop at finalTime OR finalDens (whichever comes first)
    max_time = physics[:, 0].max()
    # Only check density if physics array has multiple columns (when physics_output is provided)
    if physics.shape[1] > 1:
        max_density = physics[:, 1].max()
        assert (
            max_time <= 1.1 * params["finalTime"]
            or max_density >= 0.9 * params["finalDens"]
        ), f"Collapse should stop at finalTime OR finalDens: time={max_time:.2e}, density={max_density:.2e}"
    else:
        # When only time is output, just verify model ran successfully
        assert (
            max_time > 0
        ), f"Collapse model did not produce valid output: time={max_time:.2e}"


def test_endAtFinalDensity_validation(test_output_directory):
    """Test that endAtFinalDensity=True raises error without freefall for Cloud"""
    params = {
        "endAtFinalDensity": True,  # Invalid without freefall
        "freefall": False,
        "initialDens": 1e4,
        "finalTime": 1e5,
    }

    # Should raise ValueError for Cloud without freefall
    with pytest.raises(ValueError, match="endAtFinalDensity=True can only be used"):
        physics, chemistry, rates, heating, abundances, return_code = (
            uclchem.functional.cloud(
                param_dict=params,
                out_species=["OH", "CO"],
                return_array=True,
            )
        )


def main():
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
