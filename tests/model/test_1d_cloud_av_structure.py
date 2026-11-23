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

    # Should be sorted in reverse order, because 1 is core, N is edge.
    for timepoint in range(np.shape(cloud.physics_array)[0]):
        time = cloud.physics_array[timepoint, 0, 0]
        av = cloud.physics_array[timepoint, :, 4]
        gradient = np.diff(av)
        assert np.all(gradient < 0), (
            f"av should be monotonically decreasing, but at {time=} years (idx {timepoint}) got {gradient=}"
        )
        assert np.allclose(gradient, gradient[0]), (
            f"av should have consistent spacing, but at {time=} years (idx {timepoint}) got {gradient=}"
        )
