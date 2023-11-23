from time import perf_counter
import os
import pytest
from pathlib import Path
import tempfile
import shutil
import uclchem
try:
    import uclchem
    uclchem_imported = True
except ImportError:
    uclchem_imported = False

# Define a test to check if uclchem was imported successfully
def test_import_uclchem():
    if not uclchem_imported:
        pytest.fail("uclchem module could not be imported, make sure your environment is loaded and UCLCHEM is installed correctly.")

@pytest.fixture(scope="module")
def common_output_directory(request):
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    # Clean up the temporary directory after all tests in the module
    shutil.rmtree(temp_dir, ignore_errors=True)

# Test for ODE element conservation
def test_element_conservation():
    result = uclchem.tests.test_ode_conservation()
    for key, value in result.items():
        assert abs(value) < 1e-12, f"{key} not conserved with total rate of change {value:.2e}"

# Test for running the static model (Stage 1)
def test_static_model(common_output_directory):
    params = {
        "endAtFinalDensity": False,
        "freefall": False,
        "writeStep": 1,
        "initialDens": 1e4,
        "initialTemp": 10.0,
        "finalDens": 1e5,
        "finalTime": 5.0e6,
        "outputFile": common_output_directory / "static-full.dat",
        "abundSaveFile": common_output_directory / "startstatic.dat",
    }
    start = perf_counter()
    print("Running the model")
    return_code = uclchem.model.cloud(param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"])
    print("Finished running!")
    print(return_code)
    stop = perf_counter()
    elapsed_time = stop - start
    assert return_code[0] == 1, f"Static model returned with nonzero exit code {return_code[0]}"

# Test for running Stage 1
def test_stage_1(common_output_directory):
    params = {
        "freefall": True,
        "endAtFinalDensity": True,
        "initialDens": 1e2,
        "abundSaveFile": common_output_directory / "startstage1.dat",
        "outputFile": common_output_directory / "stage1-full.dat",
        "columnFile": common_output_directory / "stage1-column.dat",
    }

    start = perf_counter()
    return_code = uclchem.model.cloud(param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"])
    stop = perf_counter()
    elapsed_time = stop - start
    assert return_code[0] == 1, f"Stage 1 returned with nonzero exit code {return_code[0]}"

# Test for running Stage 2
# Test for running Stage 2
def test_stage_2(common_output_directory):
    params = {
        "initialDens": 1e5,
        "freezeFactor": 0.0,
        "thermdesorb": True,
        "endAtFinalDensity": False,
        "freefall": False,
        "finalTime": 1e6,
        "outputFile": common_output_directory / "stage2-full.dat",
        "abundLoadFile": common_output_directory / "startstage1.dat",
    }
    start = perf_counter()
    return_code = uclchem.model.hot_core(3, 300.0, param_dict=params, out_species=["OH", "OCS", "CO", "CS", "CH3OH"])
    stop = perf_counter()
    elapsed_time = stop - start
    assert return_code[0] == 1, f"Stage 2 returned with nonzero exit code {return_code[0]}"

def main():
    import uclchem
    # Run the tests using pytest
    pytest.main(["-v", __file__])

if __name__ == "__main__":
    main()


# if __name__=="__main__":
#     import uclchem
#     test_element_conservation()
#     test_static_model()
#     test_stage_1()
#    test_stage_2()
