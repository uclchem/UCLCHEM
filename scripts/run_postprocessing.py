"""Test some postprocessing models."""

import numpy as np
import pandas as pd

import uclchem

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
    # "Tdust",  # ruff: ignore[commented-out-code], Not present in the sample file
]

if __name__ == "__main__":
    df = pd.read_csv(
        "examples/fortran_cli/neath_small_data.out",
        sep=r"\s+",
        header=None,
        dtype=np.float64,
    )
    df.columns = NEATH_COLUMNS
    df["particle_id"] = (df["time"] == 0).astype(int).cumsum()
    for particle_id in df["particle_id"].unique():
        particle_df = df.query(f"particle_id == {particle_id}")
        model_nocoldens = uclchem.model.Postprocess(
            param_dict={},
            time_array=particle_df["time"].to_numpy(),
            density_array=particle_df["density"].to_numpy(),
            gas_temperature_array=particle_df["Tgas"].to_numpy(),
            dust_temperature_array=particle_df["Tgas"].to_numpy(),
            zeta_array=1.0,
            radfield_array=1.0,
            coldens_H_array=None,
            coldens_H2_array=None,
            coldens_CO_array=None,
            coldens_C_array=None,
        )
        model_nocoldens.check_error(only_error=True)

        _dfs_nocoldens = model_nocoldens.get_dataframes(joined=False)
        if not isinstance(_dfs_nocoldens, tuple):
            msg = "Expected tuple from get_dataframes(joined=False)"
            raise TypeError(msg)
        physics_df_nocoldens, abundances_df_nocoldens = (
            pd.DataFrame(_dfs_nocoldens[0]),
            pd.DataFrame(_dfs_nocoldens[1]),
        )
        physics_df_nocoldens.to_csv(
            "physics_nocoldens.csv",
            index=False,
        )
        abundances_df_nocoldens.to_csv("abunds_nocoldens.csv", index=False)

        model_coldens = uclchem.model.Postprocess(
            param_dict={
                #     outputfile="postprocess.dat", # ruff: ignore[commented-out-code]
            },
            out_species=["H2"],
            time_array=particle_df["time"].to_numpy(),
            density_array=particle_df["density"].to_numpy(),
            gas_temperature_array=particle_df["Tgas"].to_numpy(),
            dust_temperature_array=particle_df["Tgas"].to_numpy(),
            zeta_array=np.array([1.0]),
            radfield_array=np.array([1.0]),
            coldens_H_array=particle_df["N_H"].to_numpy(),
            coldens_H2_array=particle_df["N_H2"].to_numpy(),
            coldens_CO_array=particle_df["N_CO"].to_numpy(),
            coldens_C_array=np.array([0.0]),
        )
        model_coldens.check_error(only_error=False)
        _dfs_coldens = model_coldens.get_dataframes(joined=False)
        if not isinstance(_dfs_coldens, tuple):
            msg = "Expected tuple from get_dataframes(joined=False)"
            raise TypeError(msg)
        physics_df_coldens, abundances_df_coldens = (
            pd.DataFrame(_dfs_coldens[0]),
            pd.DataFrame(_dfs_coldens[1]),
        )
        physics_df_coldens.to_csv(
            "physics_coldens.csv",
            index=False,
        )

        abundances_df_coldens.to_csv("abunds_coldens.csv", index=False)
        break
