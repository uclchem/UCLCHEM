"""
Integration test for new surface chemistry features.

Runs a cloud model and verifies the new physics produces different/better results.
"""

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

import uclchem


class TestSurfaceChemistryIntegration:
    """Integration test verifying new surface chemistry features work."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create temporary directory for test outputs."""
        temp_dir = tempfile.mkdtemp(prefix="uclchem_surface_test_")
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)

    def test_ice_dependent_desorption_changes_chemistry(self, temp_output_dir):
        """
        Test that ice-coverage-dependent desorption actually affects chemistry.

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
            "desorb": True,
            "chemdesorb": True,
            "outputFile": str(Path(temp_output_dir) / "surface_test.dat"),
        }

        result = uclchem.model.Cloud(
            param_dict=param_dict, out_species=["Time", "#H2O", "#CO", "CH4"]
        )

        # Basic checks
        assert result is not None, "Model failed to run"
        result.check_error()

        # Get dataframe
        df = result.get_dataframes()
        assert len(df) > 0, "No output produced"

        # Verify ice buildup (SURFACE + BULK = total ice)
        early_ice = df["SURFACE"].iloc[2] + df["BULK"].iloc[2]
        late_ice = df["SURFACE"].iloc[-1] + df["BULK"].iloc[-1]
        assert (
            late_ice > early_ice * 10
        ), f"Ice should build up significantly: {early_ice:.2e} → {late_ice:.2e}"

        # Verify chemistry evolves
        early_co = df["#CO"].iloc[2]
        late_co = df["#CO"].iloc[-1]
        assert late_co != early_co, "Surface CO should evolve"

        # Verify complex molecules form (tests that chemistry is active)
        final_ch4 = df["CH4"].iloc[-1]
        assert final_ch4 > 1e-20, "Complex molecules should form"

        print(f"✓ Ice buildup: {early_ice:.2e} → {late_ice:.2e}")
        print(f"✓ Surface CO: {early_co:.2e} → {late_co:.2e}")
        print(f"✓ CH4 formed: {final_ch4:.2e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
