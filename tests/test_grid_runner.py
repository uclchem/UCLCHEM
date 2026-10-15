import numpy as np
import pytest

import uclchem

prefixes = [
    "",
    "prefix",
]


@pytest.mark.parametrize("prefix", prefixes)
def test_grid_matches_individual(prefix, tmp_path):
    dictionaries = {
        "param_dict": {
            "initialTemp": [10, 20, 30],
            "initialDens": [1e4, 1e3, 1e2],
            "finalTime": [1e2, 1e2, 1e2],
        }
    }
    grid = uclchem.model.GridRunner(
        "Cloud",
        dictionaries,
        create_grid=False,
        max_workers=3,
        grid_file=tmp_path,
        model_name_prefix=prefix,
    )

    for idx, grid_model in zip(
        range(len(dictionaries["param_dict"]["initialTemp"])),
        grid.get_models(),
        strict=True,
    ):
        param_dict = {
            key: values[idx] for key, values in dictionaries["param_dict"].items()
        }
        model = uclchem.model.Cloud(param_dict)

        assert np.allclose(model.physics_array, grid_model.physics_array)
        assert np.allclose(model.chemical_abun_array, grid_model.chemical_abun_array)


grid_data = [
    ([10], [1e3]),
    ([10, 20], [1e3, 1e4]),
    ([10, 20, 30, 40], [1e3, 1e4]),
    ([10, 20, 30, 40], [1e3]),
]


@pytest.mark.parametrize(("temperatures", "densities"), grid_data)
def test_grid_runner_create_grid(temperatures, densities):
    grid = uclchem.model.GridRunner(
        "Cloud",
        {"param_dict": {"initialTemp": temperatures, "initialDensity": densities}},
        create_grid=True,
        delay_run=True,
    )
    pending = list(
        uclchem.model.GridRunner.grid_iter(
            grid.full_parameters,
            list(grid.parameters_to_grid.keys()),
            grid.flat_grids,
            grid.model_type,
        )
    )
    assert len(pending) == len(temperatures) * len(densities)
    for param_dict in pending:
        assert param_dict["param_dict"]["initialTemp"] in temperatures
        assert param_dict["param_dict"]["initialDensity"] in densities


def test_invalid_model_name_raises():
    dictionaries = {
        "param_dict": {
            "initialTemp": [10, 20, 30],
            "initialDens": [1e4, 1e3, 1e2],
            "finalTime": [1e2, 1e2, 1e2],
        }
    }
    with pytest.raises(uclchem.model.UnknownModelError):
        _ = uclchem.model.GridRunner(
            "UNKNOWN_MODEL",
            dictionaries,
        )
