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
    # "Tdust",  # Not present in the sample file
]

if __name__ == "__main__":
    df = pd.read_csv(
        "examples/fortran_cli/neath_small_data.out",
        sep="\s+",
        header=None,
        dtype=np.float64,
    )
    df.columns = NEATH_COLUMNS
    df["particle_id"] = (df["time"] == 0.0).astype(int).cumsum()
    for particle_id in df["particle_id"].unique():
        particle_df = df.query(f"particle_id == {particle_id}")
        (
            physicsArray,
            chemicalAbunArray,
            abundanceStart,
            success_flag,
        ) = uclchem.model.postprocess(
            param_dict={},
            out_species=["H2"],
            return_array=True,
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
        pd.DataFrame(physicsArray[:, 0, :]).to_csv("physics.csv", index=False)
        pd.DataFrame(chemicalAbunArray[:, 0, :]).to_csv(
            "chemical_abundances.csv", index=False
        )
        break
