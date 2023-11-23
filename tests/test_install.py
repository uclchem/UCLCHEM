import subprocess

def test_package_installation():
    
    # Check that fortran is present:
    result = subprocess.run("gfortran --version", shell=True, text=True, capture_output=True)
    
    assert result.returncode == 0, f"Package installation failed:\n{result.stdout}\n{result.stderr}"

    # Install the package using pip
    install_command = f'pip install -e .'
    result = subprocess.run(install_command, shell=True, text=True, capture_output=True)
    
    # Check if the installation was successful
    assert result.returncode == 0, f"Package installation failed:\n{result.stdout}\n{result.stderr}"
    
    # Import the package and test if it can be imported
    try:
        import uclchem
    except ImportError:
        assert False, "Failed to import the installed package"

if __name__ == '__main__':
    import pytest
    pytest.main()
