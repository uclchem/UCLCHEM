from itertools import product

import pytest

import uclchem

temps = [10, 20, 40, 50, 80, 100, 120, 150, 200]
densities = [1e3, 1e4, 1e5, 1e6, 1e7]

param_dicts = [
    {"initialTemp": temp, "initialDens": dens, "finalTime": 1e6}
    for temp, dens in product(temps, densities)
]


@pytest.mark.parametrize("param_dict", param_dicts)
def test_cloud_no_freefall(param_dict):
    model = uclchem.model.Cloud(param_dict=param_dict)
    model.check_error()

    df = model.get_joined_dataframes()

    assert df["Time"].iloc[-1] == param_dict["finalTime"]
    assert all(df["Density"] == param_dict["initialDens"])
    assert all(df["gasTemp"] == param_dict["initialTemp"])
    assert all(df["dustTemp"] == param_dict["initialTemp"])
