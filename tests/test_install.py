import subprocess
from pathlib import Path

import pytest

from uclchem.makerates import run_makerates


@pytest.mark.timeout(600)
def test_package_installation():
    # Check that fortran is present:
    result = subprocess.run(
        "gfortran --version", shell=True, text=True, capture_output=True
    )

    assert result.returncode == 0, (
        f"Package installation failed:\n{result.stdout}\n{result.stderr}"
    )

    # Install the package using pip
    install_command = "pip install ."
    result = subprocess.run(install_command, shell=True, text=True, capture_output=True)

    # Check if the installation was successful
    assert result.returncode == 0, (
        f"Package installation failed:\n{result.stdout}\n{result.stderr}"
    )

    # Import the package and test if it can be imported
    try:
        import uclchem
    except ImportError as e:
        assert False, f"Failed to import the installed package, with ImportError: {e}"

    # Test generating a network with the small_chemistry configuration
    try:
        settings_path = (
            Path(__file__).parent / "networks" / "small_chemistry" / "user_settings.yaml"
        )
        # run_makerates will automatically find the project src/ directory
        run_makerates(str(settings_path), write_files=True)
    except Exception as e:
        assert False, f"Installing an alternative network failed: {e}"

    result = subprocess.run(install_command, shell=True, text=True, capture_output=True)

    # Check if the installation was successful
    assert result.returncode == 0, (
        f"Package installation w/ small network failed:\n{result.stdout}\n{result.stderr}"
    )

    # Import the package and test if it can be imported
    try:
        import uclchem  # noqa: F401
    except ImportError:
        assert False, "Failed to import the installed package"


if __name__ == "__main__":
    import pytest

    pytest.main()
