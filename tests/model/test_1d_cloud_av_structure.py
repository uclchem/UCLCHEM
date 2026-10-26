import numpy as np
import pytest

import uclchem


@pytest.mark.parametrize("run_type", ["managed", "external"])
def test_1d_cloud_correct_av(run_type):
    param_dict = {
        "freefall": False,
        "initialDens": 1e4,  # cm^-3
        "initialTemp": 10.0,  # K
        "finalTime": 1e1,  # yr
        "rout": 0.1,  # pc
        "baseAv": 1.0,
        "points": 5,
        "parcelStoppingMode": 0,
    }

    cloud = uclchem.model.Cloud(param_dict=param_dict, run_type=run_type)
    if run_type == "external":
        cloud.run()

    # time step 0, for all points, and show Av
    av = cloud.physics_array[0, :, 4]

    # Should be sorted in reverse order, because 1 is core, N is edge.
    gradient = np.diff(av)
    assert np.all(gradient < 0)  # Monotonically decreasing
    assert np.allclose(gradient, gradient[0])  # With constant spacing between each point
