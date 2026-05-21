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
    # "Tdust",  # noqa: ERA001, Not present in the sample file
]

if __name__ == "__main__":
    df = pd.read_csv(
        "examples/fortran_cli/neath_small_data.out",
        sep=r"\s+",
        header=None,
        dtype=np.float64,
    )
    df.columns = NEATH_COLUMNS
    df["particle_id"] = (df["time"] == 0.0).astype(int).cumsum()
    for particle_id in df["particle_id"].unique():
        particle_df = df.query(f"particle_id == {particle_id}")
        model_nocoldens = uclchem.model.Postprocess(
            param_dict={},
            out_species=["H2"],
            time_array=particle_df["time"],
            density_array=particle_df["density"],
            gas_temperature_array=particle_df["Tgas"],
            dust_temperature_array=particle_df["Tgas"],
            zeta_array=1.0,
            radfield_array=1.0,
            coldens_H_array=None,
            coldens_H2_array=None,
            coldens_CO_array=None,
            coldens_C_array=None,
        )
        model_nocoldens.check_error(only_error=True)

        physics_df_nocoldens, abundances_df_nocoldens = model_nocoldens.get_dataframes()
        physics_df_nocoldens.to_csv(
            "physics_nocoldens.csv",
            index=False,
        )
        abundances_df_nocoldens.to_csv("abunds_nocoldens.csv", index=False)

        model_coldens = uclchem.model.Postprocess(
            param_dict={
                #     outputfile="postprocess.dat", # noqa: ERA001
            },
            out_species=["H2"],
            time_array=particle_df["time"],
            density_array=particle_df["density"],
            gas_temperature_array=particle_df["Tgas"],
            dust_temperature_array=particle_df["Tgas"],
            zeta_array=1.0,
            radfield_array=1.0,
            coldens_H_array=particle_df["N_H"],
            coldens_H2_array=particle_df["N_H2"],
            coldens_CO_array=particle_df["N_CO"],
            coldens_C_array=0.0,
        )
        model_coldens.check_error(only_error=False)
        physics_df_coldens, abundances_df_coldens = model_coldens.get_dataframes()
        physics_df_coldens.to_csv(
            "physics_coldens.csv",
            index=False,
        )

        abundances_df_coldens.to_csv("abunds_coldens.csv", index=False)
        break
