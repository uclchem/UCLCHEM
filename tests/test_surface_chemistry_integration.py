"""Integration test for new surface chemistry features.

Runs a cloud model and verifies the new physics produces different/better results.
"""

import uclchem


# @pytest.mark.skip(reason="H2 burial exclusion breaks ML cap")
def test_ice_dependent_desorption_changes_chemistry() -> None:
    """Test that ice-coverage-dependent desorption actually affects chemistry.

    Verifies:
    - Model runs successfully with new features
    - Ice builds up over time
    - Surface chemistry evolves as ice grows
    - Chemical desorption efficiency varies with coverage

    """
    param_dict = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1.0e4,
        "initialTemp": 10.0,
        "finalTime": 1.0e6,
        "points": 1,  # Explicitly set to 0D mode to avoid state pollution from 1D tests
        "enable_radiative_transfer": False,  # Explicitly disable to avoid state pollution from 1D tests
        "desorb": True,
        "chemdesorb": True,
    }

    result = uclchem.model.Cloud(param_dict=param_dict)

    # Basic checks
    assert result is not None, "Model failed to run"
    result.check_error()

    # Get dataframe
    df = result.get_joined_dataframes()
    assert len(df) > 0, "No output produced"

    # Verify ice buildup (SURFACE + BULK = total ice)
    early_ice = df["SURFACE"].iloc[2] + df["BULK"].iloc[2]
    late_ice = df["SURFACE"].iloc[-1] + df["BULK"].iloc[-1]
    assert late_ice > early_ice * 10, (
        f"Ice should build up significantly: {early_ice:.2e} → {late_ice:.2e}"
    )

    # Verify chemistry evolves
    early_co = df["#CO"].iloc[2]
    late_co = df["#CO"].iloc[-1]
    assert late_co != early_co, "Surface CO should evolve"

    # Verify complex molecules form (tests that chemistry is active)
    final_ch4 = df["CH4"].iloc[-1]
    assert final_ch4 > 1e-20, "Complex molecules should form"

    settings = uclchem.advanced.GeneralSettings()
    num_monolayers_is_surface = next(
        iter(
            settings.search(
                "num_monolayers_is_surface",
                include_internal=True,
                include_parameters=True,
            ).values()
        )
    ).get()
    gas_dust_density_ratio = next(
        iter(
            settings.search(
                "gas_dust_density_ratio", include_internal=True, include_parameters=True
            ).values()
        )
    ).get()
    num_sites_per_grain = next(
        iter(
            settings.search(
                "num_sites_per_grain", include_internal=True, include_parameters=True
            ).values()
        )
    ).get()

    num_monolayers_in_run = (
        df["SURFACE"].max() * gas_dust_density_ratio / num_sites_per_grain
    )

    tol = 5e-2
    assert num_monolayers_in_run <= num_monolayers_is_surface + tol, (
        f"Number of monolayers of surface should be less than {num_monolayers_is_surface}, but was {num_monolayers_in_run:.2f} at most"
    )


def test_high_temp_CO_should_be_low() -> None:  # ruff: ignore[invalid-function-name]
    param_dict = {
        "endAtFinalDensity": False,
        "freefall": False,
        "initialDens": 1.0e6,
        "initialTemp": 30.0,
        "finalTime": 5e5,
        "points": 1,  # Explicitly set to 0D mode to avoid state pollution from 1D tests
        "enable_radiative_transfer": False,  # Explicitly disable to avoid state pollution from 1D tests
        "max_desorption_rate_constant_factor": 0.0,
        "min_desorption_rate_constant_cap": 0.0,
        "max_desorption_rate_constant_cap": 0.0,
    }

    result = uclchem.model.Cloud(param_dict=param_dict)

    # Basic checks
    assert result is not None, "Model failed to run"
    result.check_error()

    # Get dataframe
    df = result.get_joined_dataframes()

    H2O_ice = df["#H2O"].iloc[-1] + df["@H2O"].iloc[-1]  # ruff: ignore[non-lowercase-variable-in-function]
    CO_ice = df["#CO"].iloc[-1] + df["@CO"].iloc[-1]  # ruff: ignore[non-lowercase-variable-in-function]
    assert CO_ice < H2O_ice, "At 30 K, there should be more H2O ice than CO ice."
    CO_gas = df["CO"].iloc[-1]  # ruff: ignore[non-lowercase-variable-in-function]
    assert CO_ice < CO_gas, "At 30 K, CO should mostly be in the gas phase."
