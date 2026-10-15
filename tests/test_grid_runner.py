import numpy as np

import uclchem


def test_grid_matches_individual(tmp_path):
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
