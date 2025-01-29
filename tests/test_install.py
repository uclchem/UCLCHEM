import subprocess


def test_package_installation():
    # Check that fortran is present:
    result = subprocess.run(
        "gfortran --version", shell=True, text=True, capture_output=True
    )

    assert (
        result.returncode == 0
    ), f"Package installation failed:\n{result.stdout}\n{result.stderr}"

    # Install the package using pip
    install_command = "pip install -e ."
    result = subprocess.run(install_command, shell=True, text=True, capture_output=True)

    # Check if the installation was successful
    assert (
        result.returncode == 0
    ), f"Package installation failed:\n{result.stdout}\n{result.stderr}"

    # Import the package and test if it can be imported
    try:
        import uclchem  # noqa: F401
    except ImportError:
        assert False, "Failed to import the installed package"

    result = subprocess.run(
        "cd Makerates; python Makerates.py data/small_chemistry/user_settings.yaml; cd ../ ",
        shell=True,
        text=True,
        capture_output=True,
    )

    assert (
        result.returncode == 0
    ), f"Installing an alternative network failed: \n{result.stdout}\n{result.stderr}"

    result = subprocess.run(install_command, shell=True, text=True, capture_output=True)

    # Check if the installation was successful
    assert (
        result.returncode == 0
    ), f"Package installation w/ small network failed:\n{result.stdout}\n{result.stderr}"

    # Import the package and test if it can be imported
    try:
        import uclchem  # noqa: F401
    except ImportError:
        assert False, "Failed to import the installed package"


if __name__ == "__main__":
    import pytest

    pytest.main()
