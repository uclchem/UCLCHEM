"""Collection of useful functions to run a set of notebooks in a directory.

Example:
    `python3 run_notebooks.py directory`

    will run all the notebooks in `"directory"`.

"""

import argparse
import os
import shutil
import subprocess
import sys

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor


def get_python_executable() -> str:
    """Get a valid Python executable path.

    Returns:
        str: Path to python executable

    Raises:
        RuntimeError: If no valid python executable can be found.

    """
    # Try sys.executable first
    if os.path.exists(sys.executable):
        return sys.executable

    # Fall back to finding python3 in PATH
    python3_path = shutil.which("python3")
    if python3_path and os.path.exists(python3_path):
        return python3_path

    # Fall back to finding python in PATH
    python_path = shutil.which("python")
    if python_path and os.path.exists(python_path):
        return python_path

    msg = "Could not find a valid python executable"
    raise RuntimeError(msg)


def ensure_ipykernel_installed(python_exec: str) -> bool:
    """Ensure ipykernel is installed in the current environment.

    If it is not installed, tries to install it using `pip`.

    Args:
        python_exec (str): Path to python executable

    Returns:
        bool: Whether ipykernel is installed

    """
    try:
        # Check if ipykernel is already installed
        result = subprocess.run(
            [python_exec, "-c", "import ipykernel"], capture_output=True, text=True
        )
        if result.returncode == 0:
            print("ipykernel is already installed")
            return True
    except Exception:
        print(
            f"Could not determine whether ipykernel is installed for python executable {python_exec}"
        )

    # Install ipykernel if not present
    print("Installing ipykernel...")
    try:
        result = subprocess.run(
            [python_exec, "-m", "pip", "install", "ipykernel"],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print("ipykernel installed successfully")
            return True
        else:
            print(f"Failed to install ipykernel: {result.stderr}")
            return False
    except Exception as e:
        print(f"Error installing ipykernel: {e}")
        return False


def find_existing_kernel_spec(python_exec: str) -> None:
    """Check if a kernel spec already exists for this Python executable.

    Args:
        python_exec (str): Path to python executable

    """
    try:
        # List all available kernel specs
        result = subprocess.run(
            ["jupyter", "kernelspec", "list", "--json"], capture_output=True, text=True
        )

        if result.returncode == 0:
            import json

            kernel_specs = json.loads(result.stdout)

            # Check each kernel spec to see if it uses our Python executable
            specs = kernel_specs.get("kernelspecs", {}).items()
            for kernel_name, spec_info in specs:
                spec_file = os.path.join(spec_info["resource_dir"], "kernel.json")
                if os.path.exists(spec_file):
                    try:
                        with open(spec_file) as f:
                            kernel_config = json.load(f)
                            # Check if argv[0] (Python executable) matches ours
                            argv = kernel_config.get("argv", [])
                            if argv and os.path.samefile(argv[0], python_exec):
                                print(f"Found existing kernel spec: {kernel_name}")
                                return kernel_name
                    except (json.JSONDecodeError, OSError, FileNotFoundError):
                        continue
    except Exception as e:
        print(f"Error checking existing kernel specs: {e}")
    return None


def create_kernel_spec(python_exec: str) -> str | None:
    """Create a kernel spec for the current Python executable.

    Args:
        python_exec (str): Path to python executable

    Returns:
        str | None: name of the kernel, or None if it failed to create kernel

    """
    # Get the current working directory name for a human-readable kernel name
    cwd = os.getcwd()
    dir_name = os.path.basename(cwd)

    # If the directory is just "uclchem" or "UCLCHEM", make it more specific
    if dir_name.lower() == "uclchem":
        dir_name = "uclchem_notebooks"

    kernel_name = f"{dir_name}_kernel"

    try:
        # Install kernel spec
        result = subprocess.run(
            [
                python_exec,
                "-m",
                "ipykernel",
                "install",
                "--user",
                "--name",
                kernel_name,
                "--display-name",
                f"Python ({dir_name})",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            print(f"Created kernel spec: {kernel_name}")
            return kernel_name
        else:
            print(f"Failed to create kernel spec: {result.stderr}")
            return None
    except Exception as e:
        print(f"Error creating kernel spec: {e}")
        return None


def cleanup_kernel_spec(kernel_name: str) -> None:
    """Remove the temporary kernel spec.

    Args:
        kernel_name (str): name of the kernel to clean up

    """
    if kernel_name:
        try:
            subprocess.run(
                ["jupyter", "kernelspec", "remove", kernel_name, "-f"],
                capture_output=True,
            )
            print(f"Cleaned up kernel spec: {kernel_name}")
        except Exception as e:
            print(
                f"Exception occurred while cleaning up kernel spec: {kernel_name}. Exception: {e}"
            )


def run_all_notebooks(notebooks_dir: str) -> None:
    """Run all jupyter notebooks in a directory.

    Args:
        notebooks_dir (str): directory with notebooks.

    Raises:
        RuntimeError: If the ipykernel cannot be installed.
        KeyboardInterrupt: If the user presses Ctrl+c to interrupt the notebooks.

    """
    python_exec = get_python_executable()
    print(f"Using Python executable: {python_exec}")

    # Generate notebooks from .py files if notebooks.sh exists
    notebooks_script = os.path.join(notebooks_dir, "notebooks.sh")
    if os.path.exists(notebooks_script):
        print("Generating .ipynb from .py files...")
        try:
            subprocess.run(
                ["bash", notebooks_script, "--generate"],
                cwd=notebooks_dir,
                check=True,
                capture_output=True,
                text=True,
            )
            print("Notebooks generated successfully")
        except subprocess.CalledProcessError as e:
            print(f"Warning: Failed to generate notebooks: {e}")
            print(f"Stderr: {e.stderr}")

    # Ensure ipykernel is installed
    if not ensure_ipykernel_installed(python_exec):
        msg = "Failed to install ipykernel"
        raise RuntimeError(msg)

    # Check if a kernel spec already exists for this environment
    existing_kernel = find_existing_kernel_spec(python_exec)
    if existing_kernel:
        kernel_name = existing_kernel
        print(f"Using existing kernel: {kernel_name}")
    else:
        # Create a kernel spec for this environment
        kernel_name = create_kernel_spec(python_exec)
        if not kernel_name:
            msg = "Failed to create kernel spec"
            raise RuntimeError(msg)
        print(f"Created new kernel: {kernel_name}")

    try:
        for root, _, files in os.walk(notebooks_dir):
            for file in files:
                if file.endswith(".ipynb"):
                    notebook_path = os.path.join(root, file)
                    print(f"Running {notebook_path}...")

                    with open(notebook_path) as f:
                        nb = nbformat.read(f, as_version=4)

                        ep = ExecutePreprocessor(timeout=3600, kernel_name=kernel_name)

                        ep.preprocess(nb, {"metadata": {"path": root}})

                    with open(notebook_path, "w") as f:
                        nbformat.write(nb, f)
    except KeyboardInterrupt:
        print("\nNotebook execution interrupted by user")
        raise
    except Exception as e:
        print(f"Error running notebooks: {e}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run all Jupyter notebooks in a directory."
    )
    parser.add_argument(
        "notebooks_dir",
        type=str,
        help="Path to the directory containing Jupyter notebooks.",
    )
    args = parser.parse_args()

    run_all_notebooks(args.notebooks_dir)
