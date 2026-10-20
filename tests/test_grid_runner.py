from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import uclchem
from uclchem.constants import SECONDS_PER_YEAR

grid_runner_test_data = list(product(("", "prefix"), (False, True)))


@pytest.mark.parametrize(("prefix", "in_memory"), grid_runner_test_data)
def test_grid_matches_individual(prefix, in_memory, tmp_path):
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
        grid_file=None if in_memory else tmp_path,
        model_name_prefix=prefix,
    )

    assert grid.num_models == 3

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


create_grid_data = [
    ([10], [1e3]),
    ([10, 20], [1e3, 1e4]),
    ([10, 20, 30, 40], [1e3, 1e4]),
    ([10, 20, 30, 40], [1e3]),
]


@pytest.mark.parametrize(("temperatures", "densities"), create_grid_data)
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
    assert len(pending) == grid.num_models == len(temperatures) * len(densities)
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
            delay_run=True,
        )


NEATH_COLUMNS = [
    "time",
    "x",
    "y",
    "z",
    "density",
    "vx",
    "vy",
    "vz",
    "Tgas",
    "xH2",
    "xCO",
    "N_H",
    "N_H2",
    "N_CO",
]


def test_grid_runner_postprocessing(tmp_path):
    uclchem_root = Path(__file__).parent.parent
    df = pd.read_csv(
        uclchem_root / "examples/fortran_cli/neath_small_data.out",
        sep=r"\s+",
        header=None,
        dtype=np.float64,
    )
    df.columns = NEATH_COLUMNS
    df["time"] /= SECONDS_PER_YEAR  # Convert time to years
    df = df.where(df["time"] < 1e4)  # Only take short time
    df["particle_id"] = (df["time"] == 0).astype(int).cumsum()
    particle_ids = df["particle_id"].unique()

    models = []
    for particle_id in particle_ids:
        particle_df = df.query(f"particle_id == {particle_id}")

        model = uclchem.model.Postprocess(
            param_dict={},
            time_array=particle_df["time"].to_numpy(),
            density_array=particle_df["density"].to_numpy(),
            gas_temperature_array=particle_df["Tgas"].to_numpy(),
            dust_temperature_array=particle_df["Tgas"].to_numpy(),
            zeta_array=1.0,
            radfield_array=1.0,
            coldens_H_array=particle_df["N_H"].to_numpy(),
            coldens_H2_array=particle_df["N_H2"].to_numpy(),
            coldens_CO_array=particle_df["N_CO"].to_numpy(),
            coldens_C_array=0.0,
        )
        model.check_error(only_error=False)
        models.append(model)

    def _create_parameter_dictionary(
        dataframes: list[pd.DataFrame],
    ) -> dict[str, list[np.ndarray]]:
        param_dict = defaultdict(list)
        for df in dataframes:
            ones = np.full(df["time"].to_numpy().shape, 1.0)
            zeros = np.full(df["time"].to_numpy().shape, 0.0)

            param_dict["time_array"].append(df["time"].to_numpy())
            param_dict["density_array"].append(df["density"].to_numpy())
            param_dict["gas_temperature_array"].append(df["Tgas"].to_numpy())
            param_dict["dust_temperature_array"].append(df["Tgas"].to_numpy())
            param_dict["zeta_array"].append(ones)
            param_dict["radfield_array"].append(ones)
            param_dict["coldens_H_array"].append(df["N_H"].to_numpy())
            param_dict["coldens_H2_array"].append(df["N_H2"].to_numpy())
            param_dict["coldens_CO_array"].append(df["N_CO"].to_numpy())
            param_dict["coldens_C_array"].append(zeros)
        return param_dict

    param_dict = _create_parameter_dictionary(
        [df.query(f"particle_id == {particle_id}") for particle_id in particle_ids]
    )

    grid = uclchem.model.GridRunner(
        model_type="Postprocess",
        full_parameters={"param_dict": {}, **param_dict},
        create_grid=False,
        max_workers=5,
        delay_run=True,
        grid_file=tmp_path / "grid.h5",
    )
    assert grid.num_models == len(particle_ids)

    grid.run()

    grid_models = list(grid.get_models())
    assert len(grid_models) == grid.num_models
    for grid_model, model in zip(grid_models, models, strict=True):
        assert np.allclose(grid_model.physics_array, model.physics_array)
        assert np.allclose(grid_model.chemical_abun_array, grid_model.chemical_abun_array)
