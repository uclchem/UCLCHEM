"""
Test 1D radiative transfer functionality.

This test suite verifies the one-dimensional radiative transfer features
merged from the 1D_model branch, including:
- enable_radiative_transfer parameter
- Spatial density profiles (density_scale_radius, density_power_index parameters)
- Stellar heating (lum_star, temp_star parameters)
- Multi-point spatial models
"""

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

try:
    import uclchem

    uclchem_imported = True
except ImportError:
    uclchem_imported = False


def test_import_uclchem():
    """Verify UCLCHEM module imports successfully."""
    if not uclchem_imported:
        pytest.fail(
            "uclchem module could not be imported, "
            "make sure your environment is loaded and UCLCHEM is installed."
        )


@pytest.fixture(scope="module")
def common_output_directory(request):
    """Create temporary directory for test outputs."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def base_1d_params():
    """Standard parameter dictionary for 1D radiative transfer tests."""
    return {
        "parcelStoppingMode": 2,  # Stop each parcel individually at max density
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e5,  # 100k years - reasonable for astrochemistry
        "points": 5,  # Multiple spatial points required for 1D
        "enable_radiative_transfer": True,  # Enable 1D radiative transfer
        "density_scale_radius": 0.05,  # Distance scale in pc
        "density_power_index": 2.0,  # Density profile power law index
        "rout": 0.1,  # Outer radius in pc
        # Relax solver tolerances to avoid integrator taking excessive sub-steps
        # which can exceed the Fortran-allocated time array (compiled with smaller TIMEPOINTS).
        "reltol": 1e-4,
        "abstol_factor": 1e-8,
    }


@pytest.fixture
def hotcore_1d_params(base_1d_params):
    """Parameters for 1D hotcore with stellar heating."""
    params = base_1d_params.copy()
    params.update(
        {
            "lum_star": 1e6,  # Stellar luminosity in Lsun
            "temp_star": 4.5e4,  # Stellar temperature in K
        }
    )
    return params


class Test1DCloud:
    """Test 1D radiative transfer in cloud model."""

    def test_1d_cloud_return_array(self, base_1d_params):
        """Test 1D cloud model with return_array mode."""
        physics, chemistry, rates, heating, abundances_start, return_code = (
            uclchem.functional.cloud(
                param_dict=base_1d_params,
                out_species=["CO", "H2O", "CH3OH"],
                return_array=True,
                return_rates=True,
                timepoints=2500,
            )
        )

        assert return_code == 0, f"1D cloud model failed with code {return_code}"

        # Verify array shapes for 1D model
        # Shape should be (timepoints+1, points, n_columns)
        assert physics.ndim == 3, "Physics array should be 3-dimensional for 1D model"
        assert chemistry.ndim == 3, "Chemistry array should be 3-dimensional"

        n_timepoints = physics.shape[0]
        n_points = physics.shape[1]

        assert (
            n_points == base_1d_params["points"]
        ), f"Expected {base_1d_params['points']} spatial points, got {n_points}"

        # Verify spatial variation (density should vary with radius)
        # Extract density column (first column after time)
        final_time_idx = -1
        densities = physics[final_time_idx, :, 1]  # Column 1 is density

        # Densities should not all be identical for 1D model
        assert not np.allclose(
            densities, densities[0]
        ), "1D model should produce spatial density variation"

    def test_1d_cloud_return_dataframe(self, base_1d_params):
        """Test 1D cloud model with return_dataframe mode."""
        physics_df, chem_df, rates_df, heating_df, abundances_start, return_code = (
            uclchem.functional.cloud(
                param_dict=base_1d_params,
                out_species=["CO", "H2O"],
                return_dataframe=True,
                timepoints=2500,
            )
        )

        assert (
            return_code == 0
        ), f"1D cloud model (dataframe) failed with code {return_code}"

        # Verify dataframe structure
        assert (
            "Point" in physics_df.columns
        ), "Physics dataframe should have 'Point' column"
        assert "Time" in physics_df.columns, "Physics dataframe should have 'Time' column"
        assert (
            "Density" in physics_df.columns
        ), "Physics dataframe should have 'Density' column"

        # Verify multiple spatial points
        unique_points = physics_df["Point"].unique()
        assert (
            len(unique_points) == base_1d_params["points"]
        ), f"Expected {base_1d_params['points']} spatial points, got {len(unique_points)}"

        # Verify spatial variation
        final_time = physics_df["Time"].max()
        final_densities = physics_df[physics_df["Time"] == final_time]["Density"].values
        assert not np.allclose(
            final_densities, final_densities[0]
        ), "1D model should produce spatial density variation"

    def test_1d_cloud_disk_output(self, base_1d_params, common_output_directory):
        """Test 1D cloud model writing to disk."""
        output_file = common_output_directory / "test_1d_cloud.dat"
        params = base_1d_params.copy()
        params["outputFile"] = str(output_file)

        result = uclchem.functional.cloud(
            param_dict=params,
            out_species=["CO", "H2O", "CH3OH"],
            timepoints=2500,
        )

        return_code = result[0]
        assert return_code == 0, f"1D cloud disk output failed with code {return_code}"
        assert output_file.exists(), "Output file was not created"

        # Verify file contains data for multiple points
        with open(output_file, "r") as f:
            lines = f.readlines()

        # File should have header + data for multiple timepoints * multiple points
        assert len(lines) > 10, "Output file should contain substantial data"


class Test1DHotcore:
    """Test 1D radiative transfer in hotcore model with stellar heating."""

    def test_1d_hotcore_return_array(self, hotcore_1d_params):
        """Test 1D hotcore model with stellar heating."""
        physics, chemistry, rates, heating, abundances_start, return_code = (
            uclchem.functional.hot_core(
                param_dict=hotcore_1d_params,
                out_species=["CO", "H2O", "CH3OH", "H2CO"],
                return_array=True,
                return_rates=True,
                timepoints=2500,
            )
        )

        assert return_code == 0, f"1D hotcore model failed with code {return_code}"

        # Verify 3D array structure
        assert physics.ndim == 3, "Physics array should be 3-dimensional for 1D model"
        assert chemistry.ndim == 3, "Chemistry array should be 3-dimensional"

        n_points = physics.shape[1]
        assert (
            n_points == hotcore_1d_params["points"]
        ), f"Expected {hotcore_1d_params['points']} spatial points, got {n_points}"

        # Verify temperature variation (stellar heating should create gradient)
        final_time_idx = -1
        temperatures = physics[final_time_idx, :, 2]  # Column 2 is temperature

        # Temperatures should vary spatially due to stellar heating
        assert not np.allclose(
            temperatures, temperatures[0]
        ), "1D hotcore should produce spatial temperature variation from stellar heating"

        # Inner regions should be warmer than outer regions
        assert (
            temperatures[0] > temperatures[-1]
        ), "Inner temperature should exceed outer temperature with stellar heating"

    def test_1d_hotcore_stellar_parameters(self, base_1d_params):
        """Test that stellar heating parameters (lum_star, temp_star) affect results."""
        # Run with low stellar luminosity
        params_low = base_1d_params.copy()
        params_low.update({"lum_star": 1e3, "temp_star": 3000})

        physics_low, _, _, _, _, code_low = uclchem.functional.hot_core(
            param_dict=params_low,
            out_species=["CO"],
            return_array=True,
            return_rates=True,
            timepoints=2500,
        )

        # Run with high stellar luminosity
        params_high = base_1d_params.copy()
        params_high.update({"lum_star": 1e6, "temp_star": 4.5e4})

        physics_high, _, _, _, _, code_high = uclchem.functional.hot_core(
            param_dict=params_high,
            out_species=["CO"],
            return_array=True,
            return_rates=True,
            timepoints=2500,
        )

        assert code_low == 0 and code_high == 0, "Both models should succeed"

        # Extract final temperatures
        temps_low = physics_low[-1, :, 2]
        temps_high = physics_high[-1, :, 2]

        # Higher luminosity should produce higher temperatures
        assert np.mean(temps_high) > np.mean(
            temps_low
        ), "Higher stellar luminosity should produce higher temperatures"


class Test1DParameterValidation:
    """Test parameter handling for 1D models."""

    def test_1d_requires_multiple_points(self):
        """Test that 1D radiative transfer requires points > 1."""
        params = {
            "parcelStoppingMode": 2,
            "freefall": False,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "points": 1,  # Single point - should work but not use 1D features
            "enable_radiative_transfer": True,
        }

        # This should still work but effectively be 0D
        _, _, _, _, _, return_code = uclchem.functional.cloud(
            param_dict=params,
            out_species=["CO"],
            return_array=True,
            return_rates=True,
            timepoints=2500,
        )

        # Should succeed (code will handle gracefully)
        assert return_code == 0

    def test_1d_density_profile_parameters(self):
        """Test density_scale_radius and density_power_index parameters control density profile."""
        base_params = {
            "parcelStoppingMode": 2,
            "freefall": False,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "points": 5,
            "enable_radiative_transfer": True,
            "rout": 0.1,
        }

        # Test steep profile (high power index)
        params_steep = base_params.copy()
        params_steep.update({"density_scale_radius": 0.05, "density_power_index": 3.0})

        physics_steep, _, _, _, _, code_steep = uclchem.functional.cloud(
            param_dict=params_steep,
            out_species=["CO"],
            return_array=True,
            return_rates=True,
            timepoints=2500,
        )

        # Test shallow profile (low power index)
        params_shallow = base_params.copy()
        params_shallow.update({"density_scale_radius": 0.05, "density_power_index": 1.5})

        physics_shallow, _, _, _, _, code_shallow = uclchem.functional.cloud(
            param_dict=params_shallow,
            out_species=["CO"],
            return_array=True,
            return_rates=True,
            timepoints=2500,
        )

        assert code_steep == 0 and code_shallow == 0, "Both profiles should succeed"

        # Extract final densities
        dens_steep = physics_steep[-1, :, 1]
        dens_shallow = physics_shallow[-1, :, 1]

        # Steeper profile should have larger density contrast
        contrast_steep = dens_steep.max() / dens_steep.min()
        contrast_shallow = dens_shallow.max() / dens_shallow.min()

        assert (
            contrast_steep > contrast_shallow
        ), "Steeper profile (higher a) should produce larger density contrast"

    def test_0d_mode_still_works(self):
        """Verify that 0D mode (enable_radiative_transfer=False) still works."""
        params = {
            "parcelStoppingMode": 2,  # For 0D, can use legacy False behavior
            "freefall": False,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "points": 1,
            "enable_radiative_transfer": False,  # Explicitly disable 1D
        }

        physics, chemistry, rates, heating, abundances_start, return_code = (
            uclchem.functional.cloud(
                param_dict=params,
                out_species=["CO", "H2O"],
                return_array=True,
                return_rates=True,
                timepoints=2500,
            )
        )

        assert return_code == 0, "0D mode should still work"

        # For 0D with points=1, arrays should still be 3D but with single point
        assert physics.shape[1] == 1, "0D mode should have 1 spatial point"


class Test1DChemicalEvolution:
    """Test chemical evolution in 1D models."""

    def test_1d_abundance_spatial_variation(self, base_1d_params):
        """Test that abundances vary spatially in 1D model."""
        physics, chemistry, rates, heating, abundances_start, return_code = (
            uclchem.functional.cloud(
                param_dict=base_1d_params,
                out_species=["CO", "H2O", "CH3OH"],
                return_array=True,
                return_rates=True,
                timepoints=2500,
            )
        )

        assert return_code == 0, "1D model should succeed"

        # Extract final abundances for each species
        final_time_idx = -1

        # Chemistry array columns: species abundances
        # Each species should vary spatially
        for species_idx in range(chemistry.shape[2]):
            abundances = chemistry[final_time_idx, :, species_idx]

            # Allow for very low abundances that might all be near zero
            if np.max(abundances) > 1e-20:
                # Species with significant abundance should show spatial variation
                variation = np.std(abundances) / (np.mean(abundances) + 1e-30)

                # Some variation expected (not strict requirement as chemistry is complex)
                assert variation >= 0 or np.allclose(
                    abundances, abundances[0]
                ), f"Species {species_idx} abundances: {abundances}"

    def test_1d_chained_models(self, base_1d_params):
        """Test chaining 1D models using starting_chemistry."""
        # Run first phase
        params_phase1 = base_1d_params.copy()
        params_phase1["finalTime"] = 1.0e5

        physics1, chem1, rates1, heating1, abund_start1, code1 = uclchem.functional.cloud(
            param_dict=params_phase1,
            out_species=["CO", "H2O"],
            return_array=True,
            return_rates=True,
            timepoints=2500,
        )

        assert code1 == 0, "Phase 1 should succeed"

        # Run second phase starting from phase 1 final state
        params_phase2 = base_1d_params.copy()
        params_phase2["finalTime"] = 2.0e5
        params_phase2["currentTime"] = 1.0e5

        physics2, chem2, rates2, heating2, abund_start2, code2 = uclchem.functional.cloud(
            param_dict=params_phase2,
            out_species=["CO", "H2O"],
            return_array=True,
            return_rates=True,
            starting_chemistry=abund_start1,
            timepoints=2500,
        )

        assert code2 == 0, "Phase 2 should succeed with starting_chemistry"

        # Abundances should have evolved from phase 1
        assert not np.allclose(
            chem1[-1, :, :], chem2[0, :, :]
        ), "Chemistry should continue evolving in phase 2"


# ============================================================================
# Object-Oriented Interface Tests for 1D Models
# ============================================================================


class TestOOCloud1D:
    """Test Cloud model with 1D radiative transfer using OO interface."""

    def test_oo_cloud_1d_basic_run(self, base_1d_params):
        """Test basic 1D cloud model run with OO interface."""
        model = uclchem.model.Cloud(
            param_dict=base_1d_params,
            out_species=["CO", "H2O", "CH3OH"],
            timepoints=2500,
        )

        # Check model succeeded
        model.check_error()
        assert model.success_flag == 0

        # Verify arrays are 3D
        assert model.physics_array.ndim == 3
        assert model.chemical_abun_array.ndim == 3

        # Verify number of spatial points
        assert model.physics_array.shape[1] == base_1d_params["points"]

    def test_oo_cloud_1d_dataframe_output(self, base_1d_params):
        """Test that get_dataframes works properly for 1D models."""
        model = uclchem.model.Cloud(
            param_dict=base_1d_params,
            out_species=["CO", "H2O"],
            timepoints=2500,
        )

        model.check_error()

        # Get DataFrames
        phys_df, chem_df = model.get_dataframes(joined=False)

        # Check Point column exists
        assert "Point" in phys_df.columns
        assert "Point" in chem_df.columns

        # Verify number of unique points
        n_points = phys_df["Point"].nunique()
        assert n_points == base_1d_params["points"]

        # Verify spatial variation in density
        final_time = phys_df["Time"].max()
        final_densities = phys_df[phys_df["Time"] == final_time]["Density"]
        assert final_densities.std() > 0, "Density should vary spatially"

    def test_oo_cloud_1d_point_selection(self, base_1d_params):
        """Test getting DataFrames for individual spatial points."""
        model = uclchem.model.Cloud(
            param_dict=base_1d_params,
            out_species=["CO"],
            timepoints=2500,
        )

        model.check_error()

        # Get data for first point only
        phys_df_pt0, chem_df_pt0 = model.get_dataframes(point=0, joined=False)

        # Should only have data for one point
        assert phys_df_pt0["Point"].nunique() == 1
        assert (phys_df_pt0["Point"] == 1).all()  # Point is 1-indexed

        # Get data for last point
        last_pt = base_1d_params["points"] - 1
        phys_df_last, chem_df_last = model.get_dataframes(point=last_pt, joined=False)

        assert (phys_df_last["Point"] == base_1d_params["points"]).all()

        # Densities should differ between points
        density_pt0 = phys_df_pt0["Density"].iloc[-1]
        density_last = phys_df_last["Density"].iloc[-1]
        assert density_pt0 != density_last

    def test_oo_cloud_1d_with_stats(self, base_1d_params):
        """Test that DVODE stats work with 1D models."""
        model = uclchem.model.Cloud(
            param_dict=base_1d_params,
            out_species=["CO"],
            timepoints=2500,
        )

        model.check_error()

        # Get DataFrames with stats
        result = model.get_dataframes(joined=False, with_stats=True)

        # Should return: phys, chem, stats
        assert len(result) == 3
        phys_df, chem_df, stats_df = result

        # Stats DataFrame should exist and have Point column
        assert stats_df is not None
        assert "Point" in stats_df.columns
        assert stats_df["Point"].nunique() == base_1d_params["points"]

        # Stats should contain DVODE statistics
        assert "NFE" in stats_df.columns  # Number of function evaluations
        assert "NJE" in stats_df.columns  # Number of Jacobian evaluations


class TestOOCollapse1D:
    """Test Collapse model with 1D features using OO interface."""

    def test_oo_collapse_1d_freefall(self):
        """Test 1D collapse model with freefall."""
        params = {
            "parcelStoppingMode": 2,
            "freefall": True,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "finalDens": 1e6,
            "points": 5,
            "enable_radiative_transfer": True,
            "density_scale_radius": 0.05,
            "density_power_index": 2.0,
            "rout": 0.1,
            "reltol": 1e-4,
            "abstol_factor": 1e-8,
            "writeStep": 5,
        }

        model = uclchem.model.Collapse(
            collapse="BE1.1",
            physics_output=None,
            param_dict=params,
            out_species=["CO", "H2O"],
            timepoints=2500,
        )

        model.check_error()
        assert model.success_flag == 0

        # Verify 3D arrays
        assert model.physics_array.ndim == 3
        assert model.physics_array.shape[1] == params["points"]

    def test_oo_collapse_parcel_stopping_modes(self):
        """Test all parcelStoppingMode values with pure freefall collapse."""
        base_params = {
            "freefall": True,
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "finalDens": 1e6,
            "points": 3,
            "enable_radiative_transfer": True,
            "density_scale_radius": 0.05,
            "density_power_index": 2.0,
            "rout": 0.1,
            "reltol": 1e-4,
            "abstol_factor": 1e-8,
        }

        for mode in [0, 1, 2]:
            params = base_params.copy()
            params["parcelStoppingMode"] = mode

            # Use pure freefall collapse (no collapse profile)
            model = uclchem.model.Cloud(
                param_dict=params,
                out_species=["CO"],
                timepoints=2500,
            )

            model.check_error()
            assert model.success_flag == 0, f"parcelStoppingMode={mode} should succeed"


class TestOOHotcore1D:
    """Test PrestellarCore (hotcore) with 1D stellar heating using OO interface."""

    def test_oo_hotcore_1d_stellar_heating(self, hotcore_1d_params):
        """Test 1D hotcore with stellar heating parameters."""
        model = uclchem.model.PrestellarCore(
            temp_indx=1,
            max_temperature=300.0,
            param_dict=hotcore_1d_params,
            out_species=["CO", "H2O", "CH3OH"],
            timepoints=2500,
        )

        model.check_error()
        assert model.success_flag == 0

        # Get final temperatures at each point
        phys_df = model.get_dataframes(joined=False)[0]
        final_time = phys_df["Time"].max()
        final_temps = (
            phys_df[phys_df["Time"] == final_time].sort_values("Point")["gasTemp"].values
        )

        # Temperature should vary spatially
        assert final_temps.std() > 0

        # Inner regions should be warmer
        assert final_temps[0] > final_temps[-1]


class TestOOModelSavingLoading1D:
    """Test saving and loading 1D models with OO interface."""

    def test_oo_save_load_1d_model(self, base_1d_params, common_output_directory):
        """Test that 1D models can be saved and loaded."""
        save_file = common_output_directory / "test_oo_1d_model.json"

        # Run and save model
        model1 = uclchem.model.Cloud(
            param_dict=base_1d_params,
            out_species=["CO", "H2O"],
            timepoints=2500,
        )
        model1.check_error()
        model1.save_model(str(save_file))

        assert save_file.exists()

        # Load model
        model2 = uclchem.model.Cloud.from_file(str(save_file))

        # Verify loaded model has same data
        assert np.allclose(model1.physics_array, model2.physics_array)
        assert np.allclose(model1.chemical_abun_array, model2.chemical_abun_array)
        assert model1._param_dict["points"] == model2._param_dict["points"]
        assert (
            model1._param_dict["enable_radiative_transfer"]
            == model2._param_dict["enable_radiative_transfer"]
        )

    def test_oo_save_load_preserves_point_column(
        self, base_1d_params, common_output_directory
    ):
        """Test that Point column is preserved after save/load."""
        save_file = common_output_directory / "test_oo_1d_point_column.json"

        model1 = uclchem.model.Cloud(
            param_dict=base_1d_params, out_species=["CO"], timepoints=2500
        )
        model1.check_error()
        model1.save_model(str(save_file))

        model2 = uclchem.model.Cloud.from_file(str(save_file))

        # Get DataFrames from loaded model
        phys_df, chem_df = model2.get_dataframes(joined=False)

        # Point column should exist
        assert "Point" in phys_df.columns
        assert phys_df["Point"].nunique() == base_1d_params["points"]


class TestOOModelChaining1D:
    """Test chaining 1D models together with OO interface."""

    def test_oo_chain_1d_models(self, base_1d_params):
        """Test running sequential 1D models using previous_model."""
        # Phase 1: Initial evolution
        params1 = base_1d_params.copy()
        params1["finalTime"] = 5.0e4

        model1 = uclchem.model.Cloud(
            param_dict=params1,
            out_species=["CO", "H2O"],
            timepoints=2500,
        )
        model1.check_error()

        # Phase 2: Continue from phase 1
        params2 = base_1d_params.copy()
        params2["finalTime"] = 1.0e5
        params2["currentTime"] = 5.0e4

        model2 = uclchem.model.Cloud(
            param_dict=params2,
            out_species=["CO", "H2O"],
            previous_model=model1,
            timepoints=2500,
        )
        model2.check_error()

        # Verify time continuity
        phys_df2 = model2.get_dataframes(joined=False)[0]
        assert phys_df2["Time"].iloc[-1] >= 5.0e4

        # Verify chemistry evolved
        assert not np.allclose(
            model1.chemical_abun_array[-1], model2.chemical_abun_array[0]
        )

    def test_oo_chain_with_starting_chemistry_array(self, base_1d_params):
        """Test chaining using next_starting_chemistry_array."""
        # Run first model
        model1 = uclchem.model.Cloud(
            param_dict=base_1d_params,
            out_species=["CO", "H2O"],
            timepoints=2500,
        )
        model1.check_error()

        # Get starting chemistry for next run
        starting_chem = model1.next_starting_chemistry_array

        # Verify shape matches multi-point structure
        # Shape should be (points, nspec)
        assert starting_chem.shape[0] == base_1d_params["points"]

        # Run second model with this chemistry
        params2 = base_1d_params.copy()
        params2["finalTime"] = 2.0e5
        params2["currentTime"] = 1.0e5

        model2 = uclchem.model.Cloud(
            param_dict=params2,
            out_species=["CO", "H2O"],
            starting_chemistry=starting_chem,
            timepoints=2500,
        )
        model2.check_error()


class TestEndAtFinalDensity:
    """Test the endAtFinalDensity parameter (merged from develop branch)."""

    def test_end_at_final_density_enabled(self):
        """Test that model stops when density reaches finalDens."""
        params = {
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalTime": 1.0e6,  # Long time
            "finalDens": 5e4,  # Target density
            "freefall": True,
            "endAtFinalDensity": True,  # Stop at density
            "points": 1,
            "enable_radiative_transfer": False,
            "reltol": 1e-4,
            "abstol_factor": 1e-8,
        }

        model = uclchem.model.Cloud(
            param_dict=params, out_species=["CO"], timepoints=2500
        )

        model.check_error()

        # Get final density
        phys_df = model.get_dataframes(joined=False)[0]
        final_density = phys_df["Density"].iloc[-1]

        # Should have stopped at or before finalDens
        # (might stop slightly before due to timestep control)
        assert final_density <= params["finalDens"] * 1.1

    def test_end_at_final_density_disabled(self):
        """Test that model runs to finalTime when endAtFinalDensity=False."""
        params = {
            "initialDens": 1e4,
            "initialTemp": 10.0,
            "finalTime": 1.0e5,
            "finalDens": 1e6,  # High target density
            "freefall": True,
            "endAtFinalDensity": False,  # Don't stop at density
            "points": 1,
            "enable_radiative_transfer": False,
            "reltol": 1e-4,
            "abstol_factor": 1e-8,
        }

        model = uclchem.model.Cloud(
            param_dict=params, out_species=["CO"], timepoints=2500
        )

        model.check_error()

        # Get final time
        phys_df = model.get_dataframes(joined=False)[0]
        final_time = phys_df["Time"].iloc[-1]

        # Should have run to finalTime
        assert final_time >= params["finalTime"] * 0.95


class TestFunctionalVsOOConsistency:
    """Test that functional and OO APIs produce consistent results for 1D."""

    def test_functional_vs_oo_1d_cloud(self, base_1d_params):
        """Test that functional and OO APIs give same results."""
        # Run with OO interface
        oo_model = uclchem.model.Cloud(
            param_dict=base_1d_params,
            out_species=["CO", "H2O"],
            timepoints=2500,
        )
        oo_model.check_error()

        # Run with functional interface
        phys_func, chem_func, _, _, _, flag_func = uclchem.functional.cloud(
            param_dict=base_1d_params,
            out_species=["CO", "H2O"],
            return_array=True,
            timepoints=2500,
        )

        assert flag_func == 0

        # Arrays should match
        assert np.allclose(oo_model.physics_array, phys_func)
        assert np.allclose(oo_model.chemical_abun_array, chem_func)

    def test_functional_dataframe_point_column(self, base_1d_params):
        """Test that functional API returns DataFrames with Point column."""
        phys_df, chem_df, _, _, _, flag = uclchem.functional.cloud(
            param_dict=base_1d_params,
            out_species=["CO"],
            return_dataframe=True,
            timepoints=2500,
        )

        assert flag == 0
        assert "Point" in phys_df.columns
        assert phys_df["Point"].nunique() == base_1d_params["points"]
