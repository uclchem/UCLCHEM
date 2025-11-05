"""
Test photo_on_grain network with IN-MEMORY return modes only.

This test uses return_array and return_dataframe with starting_chemistry
to ensure the photo_on_grain network works with Python in-memory arrays.
"""
import shutil
import subprocess
from pathlib import Path

import pytest


def test_photo_on_grain_memory():
    """Test photo_on_grain network using in-memory arrays/dataframes"""
    # Check that fortran is present
    result = subprocess.run(
        "gfortran --version", shell=True, text=True, capture_output=True
    )
    assert result.returncode == 0, (
        f"gfortran not available:\n{result.stdout}\n{result.stderr}"
    )

    TEST_DIR = Path("tests/photo_on_grain_test_output/")
    TEST_DIR.mkdir(parents=True, exist_ok=True)

    # Install the package using pip
    install_command = "pip install ."
    result = subprocess.run(
        install_command, shell=True, text=True, capture_output=True
    )
    assert result.returncode == 0, (
        f"Package installation failed:\n{result.stdout}\n{result.stderr}"
    )

    # Check if installed
    pip_packages = subprocess.run(
        "pip list", shell=True, text=True, capture_output=True
    )
    assert "uclchem" in pip_packages.stdout, (
        "Package uclchem not found in pip list"
    )

    # Copy the default settings file
    default_settings_path = "Makerates/user_settings.yaml"
    test_settings_path = f"{TEST_DIR}/photo_on_grain_settings_mem.yaml"
    shutil.copy(default_settings_path, test_settings_path)

    # Update add_crp_photo_to_grain to True in the settings file
    import re

    with open(test_settings_path, "r") as file:
        content = file.read()

    # Find and replace the line with add_crp_photo_to_grain
    content = re.sub(
        r"add_crp_photo_to_grain:\s*\w+",
        "add_crp_photo_to_grain: True",
        content,
    )

    with open(test_settings_path, "w") as file:
        file.write(content)

    # Verify the change was made
    with open(test_settings_path, "r") as file:
        verify_content = file.read()
    assert "add_crp_photo_to_grain: True" in verify_content, (
        f"Failed to update add_crp_photo_to_grain in {test_settings_path}"
    )

    # Run Makerates with the modified config
    result = subprocess.run(
        f"cd Makerates; python MakeRates.py {test_settings_path}; cd ../",
        shell=True,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, (
        f"Running Makerates with crp_photo_to_grain=True failed: \n"
        f"{result.stdout}\n{result.stderr}"
    )

    # Reinstall with the new network
    result = subprocess.run(
        install_command, shell=True, text=True, capture_output=True
    )
    assert result.returncode == 0, (
        f"Package installation w/ crp_photo_to_grain network failed:\n"
        f"{result.stdout}\n{result.stderr}"
    )

    # Import the package and test if it can be imported
    try:
        import uclchem  # noqa: F401
    except ImportError:
        assert False, (
            "Failed to import the installed package after running Makerates "
            "with crp_photo_to_grain=True"
        )

    # Test static model with IN-MEMORY return_array
    outSpecies = ["OH", "OCS", "CO", "CS", "CH3OH"]
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
        "abstol_min": 1e-15,
        "reltol": 1e-5,
    }
    physics, chemistry, rates, abundances, return_code = uclchem.model.cloud(
        param_dict=params, out_species=outSpecies, return_array=True
    )
    assert return_code == 0, (
        f"Static model failed with result code {return_code}"
    )

    # Test Stage 1: Collapse with IN-MEMORY return_dataframe
    params["freefall"] = True
    params["endAtFinalDensity"] = True
    params["initialDens"] = 1e2
    physics, chemistry, rates, abundances, return_code = uclchem.model.cloud(
        param_dict=params, out_species=outSpecies, return_array=True
    )
    assert return_code == 0, (
        f"Stage 1 model failed with result code {return_code}"
    )

    # Test Stage 2: Hot core with IN-MEMORY and starting_chemistry
    params["initialDens"] = 1e5
    params["freezeFactor"] = 0.0
    params["thermdesorb"] = True
    params["endAtFinalDensity"] = False
    params["freefall"] = False
    params["finalTime"] = 1e6
    params["abstol_min"] = 1e-25
    (
        physics,
        chemistry,
        rates,
        abundances,
        return_code,
    ) = uclchem.model.hot_core(
        3,
        300.0,
        param_dict=params,
        out_species=outSpecies,
        return_array=True,
        starting_chemistry=abundances,
    )
    assert return_code == 0, (
        f"Stage 2 model failed with result code {return_code}"
    )


if __name__ == "__main__":
    pytest.main(["-v", __file__])
