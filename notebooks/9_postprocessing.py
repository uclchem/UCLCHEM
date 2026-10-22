# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: UCLCHEM_312
#     language: python
#     name: python3
# ---

# %%
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import uclchem
from uclchem.constants import SECONDS_PER_YEAR


def plot_postprocess_models(models: Iterable[uclchem.model.Postprocess]) -> list[plt.Axes]:
    max_time = 0.0
    _fig, axs = plt.subplots(3, 1, sharex=True, sharey=False)
    for model in models:
        model_df = model.get_joined_dataframes()
        axs[0].plot(model_df['Time'], model_df['CO'])
        axs[1].plot(model_df['Time'], model_df['Density'])
        axs[2].plot(model_df['Time'], model_df['gasTemp'])

        max_time = max(max_time, model_df['Time'].max())
    axs[0].set_ylabel("CO abund (w.r.t. H)")
    axs[1].set_ylabel(r"$n_{\mathrm{H}}$ (cm$^{-3}$)")
    axs[2].set_ylabel(r"$T_{\mathrm{gas}}$ (K)")

    axs[0].set_xlim([1e4, max_time])
    axs[0].set_ylim([1e-14, 1e-2])
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[1].set_yscale('log')
    axs[2].set_xlabel("Time (years)")

    return axs


# %%
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
df = pd.read_csv(
    "../examples/fortran_cli/neath_small_data.out",
    sep=r"\s+",
    header=None,
    dtype=np.float64,
)
df.columns = NEATH_COLUMNS
df['time'] /= SECONDS_PER_YEAR  # Convert time to years
df["particle_id"] = (df["time"] == 0).astype(int).cumsum()
particle_ids = df["particle_id"].unique()
df.head()


# %% [markdown]
# We can now run some models where we use the physical structure from the hydrodynamics, and evolve the chemistry using `uclchem.model.Postprocess`.

# %%
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

# %%
_ = plot_postprocess_models(models)
plt.tight_layout()
plt.show()


# %% [markdown]
# We can also make use of `uclchem.model.GridRunner` to instead run the models in parallel, if we do just a bit more work.

# %%
def create_parameter_dictionary(dataframes: list[pd.DataFrame]) -> dict[str, list[np.ndarray]]:
    param_dict = defaultdict(list)
    for df in dataframes:
        ones = np.full(df['time'].to_numpy().shape, 1.0)
        zeros = np.full(df['time'].to_numpy().shape, 0.0)

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


param_dict = create_parameter_dictionary(
    [df.query(f"particle_id == {particle_id}") for particle_id in particle_ids]
)

output_dir = Path('output_9')
if not output_dir.is_dir():
    output_dir.mkdir()

grid = uclchem.model.GridRunner(
    model_type="Postprocess",
    full_parameters={"param_dict": {}, **param_dict},
    create_grid=False,
    max_workers=5,
    grid_file=output_dir / "grid.h5",
)


# %%
_ = plot_postprocess_models(grid.get_models())
plt.tight_layout()
plt.show()
