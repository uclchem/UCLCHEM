import pytest
import subprocess
import tempfile


def test_compilation():
    # Check that fortran is present:
    result = subprocess.run(
        "gfortran --version", shell=True, text=True, capture_output=True
    )

    assert result.returncode == 0, (
        f"gfortran installation failed:\n{result.stdout}\n{result.stderr}"
    )

    result = subprocess.run(
        "make", shell=True, text=True, capture_output=True, cwd="src/fortran_src"
    )
    assert result.returncode == 0, f"'make' failed:\n{result.stdout}\n{result.stderr}"

    param_dict = {"initialTemp": 10, "initialDens": 10000.0, "finalTime": 1e4}
    with tempfile.NamedTemporaryFile(mode="w") as param_dict_file:
        param_dict_file.write(str(param_dict).lower())
        param_dict_file.flush()

        result = subprocess.run(
            f"./uclchem CLOUD {param_dict_file.name}",
            shell=True,
            text=True,
            capture_output=True,
            cwd=".",
        )
    assert result.returncode == 0, (
        f"'./uclchem CLOUD <tempfile>' failed:\n{result.stdout}\n{result.stderr}"
    )
