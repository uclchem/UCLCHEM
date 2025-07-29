import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from uclchemwrap import uclchemwrap as wrap
from uclchem.analysis import check_element_conservation, create_abundance_plot, plot_species
from uclchem.constants import (
    N_PHYSICAL_PARAMETERS,
    PHYSICAL_PARAMETERS,
    TIMEPOINTS,
    n_reactions,
    n_species,
)


def reaction_line_formatter(line):
    reactants = list(filter(lambda x: not str(x).lower().endswith("nan"), line[0:3]))
    products = list(filter(lambda x: not str(x).lower().endswith("nan"), line[3:7]))
    return " + ".join(reactants) + " -> " + " + ".join(products)


class ReactionNamesStore:
    def __init__(self):
        self.reaction_names = None

    def __call__(self):
        # Only load the reactions once, after that use the cached version
        if self.reaction_names is None:
            reactions = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "reactions.csv"))
            # format the reactions:
            self.reaction_names = [reaction_line_formatter(line) for idx, line in reactions.iterrows()]
        return self.reaction_names


get_reaction_names = ReactionNamesStore()


class AbstractModel:
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: object = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        self.param_dict = {}
        self.physics_array = None
        self.chemical_abun_array = None
        self.ratesArray = None
        self.full_array = None
        self.debug = debug
        self.success_flag = None
        specname = wrap.get_specname()
        self.specname = np.array([x.strip() for x in specname.astype(str) if x != ""])

        if read_file is not None:
            self.n_out = None
            self.out_species = None
            self.timepoints = None
            self.was_read = True
            self.PHYSICAL_PARAMETERS = None

            self.give_start_abund = False
            self.next_starting_chemistry = None
            self.read_output_file(read_file)
        else:
            self.n_out = 0
            self.out_species = []
            self.timepoints = timepoints
            self.was_read = False
            self.PHYSICAL_PARAMETERS = PHYSICAL_PARAMETERS
            self.outputFile = param_dict.pop("outputFile") if "outputFile" in param_dict else None
            self.abundSaveFile = param_dict.pop("abundSaveFile") if "abundSaveFile" in param_dict else None
            self.abundLoadFile = param_dict.pop("abundLoadFile") if "abundLoadFile" in param_dict else None
            self._reform_inputs(param_dict, out_species)
            if "points" not in self.param_dict:
                self.param_dict["points"] = 1
            if previous_model is None and self.abundLoadFile is None:
                self.starting_chemistry = starting_chemistry
            elif 'next_starting_chemistry' in [attr for attr in dir(previous_model) if
                                               not callable(getattr(previous_model, attr)) and
                                               not attr.startswith("__") and
                                               attr is not None]:
                self.starting_chemistry = previous_model.next_starting_chemistry
            elif self.abundLoadFile is not None:
                self.read_starting_chemistry_output_file()
            self.give_start_abund = self.starting_chemistry is not None
            self.next_starting_chemistry = None
        return

    def __run__(self):
        if not self.was_read:
            self.physics_array = None
            self.chemical_abun_array = None
            self.ratesArray = None
            self._create_fortran_array()
            self._create_rates_array()
        else:
            raise ("This model was read. It can not be run. ")
        return

    def check_conservation(self,
                           element_list: list = ["H", "N", "C", "O"],
                           percent: bool = True
                           ):
        if self.param_dict["points"] > 1:
            for i in range(self.param_dict["points"]):
                print(f"Element conservation report for point {i + 1} of {self.param_dict['points']}")
                print(check_element_conservation(self.get_dataframes(i), element_list, percent))
        else:
            print(f"Element conservation report")
            print(check_element_conservation(self.get_dataframes(0), element_list, percent))

    def check_error(self):
        """Prints the error message of the model"""
        if self.success_flag != 0 and self.success_flag is not None:
            errors = {
                -1: "Parameter read failed. Likely due to a mispelled parameter name, compare your dictionary to the parameters docs.",
                -2: "Physics intiialization failed. Often due to user chosing unacceptable parameters such as hot core masses or collapse modes that don't exist. Check the docs for your model function.",
                -3: "Chemistry initialization failed",  # this doesn't exist yet
                -4: "Unrecoverable integrator error, DVODE failed to integrate the ODEs in a way that UCLCHEM could not fix. Run UCLCHEM tests to check your network works at all then try to see if bad parameter combination is at play.",
                -5: "Too many integrator fails. DVODE failed to integrate the ODE and UCLCHEM repeatedly altered settings to try to make it pass but tried too many times without success so code aborted to stop infinite loop.",
                -6: "The model was stopped because there are not enough time points allocated in the time array. Increase the number of time points in the time array in constants.py and try again.",
            }
            try:
                print(f'{errors[self.success_flag]}')
            except KeyError:
                raise ValueError(f"Unknown error code: {self.success_flag}")
        elif self.success_flag == 0:
            print(f'Model ran successfully.')
        elif self.success_flag is None:
            print(f'Model has not been run.')

    def create_abundance_plot(self,
                              species: list = ["H", "N", "C", "O"],
                              figsize: tuple[2] = (10, 7),
                              point: int = 0,
                              plot_file=None):
        if point > self.param_dict["points"]:
            raise Exception("'point' must be less than number of modelled points.")
        return create_abundance_plot(self.get_dataframes(point), species, figsize, plot_file)

    def get_dataframes(self, point: int = 0, joined: bool = True, with_rates: bool = False):
        """Convert the output array of "point" to a pandas dataframe
        Returns:
            _type_: _description_
        """
        # Create a physical parameter dataframe
        physics_df = pd.DataFrame(
            self.physics_array[:, point, : len(self.PHYSICAL_PARAMETERS)],
            index=None,
            columns=self.PHYSICAL_PARAMETERS,
        )
        # Create an abundances dataframe.
        chemistry_df = pd.DataFrame(
            self.chemical_abun_array[:, point, :], index=None, columns=self.specname
        )
        if self.ratesArray is not None:
            # Create a rates dataframe.
            rates_df = pd.DataFrame(
                self.ratesArray[:, point, :], index=None, columns=get_reaction_names()
            )
        else:
            rates_df = None
        if joined:
            if with_rates:
                return_df = physics_df.join(chemistry_df.join(rates_df))
            else:
                return_df = physics_df.join(chemistry_df)
            return return_df
        else:
            if with_rates:
                return physics_df, chemistry_df, rates_df
            else:
                return physics_df, chemistry_df

    def plot_species(self, ax: plt.axes, species: list[str], point: int = 0, legend: bool = True, **plot_kwargs):
        return plot_species(ax, self.get_dataframes(point), species, legend, **plot_kwargs)

    def read_output_file(self, read_file: str, rates_load_file: str = None):
        '''Perform classic output file reading.
        Args:
            abund_load_file (str): path to file containing a full UCLCHEM output
        '''
        self.was_read = True
        columns = np.char.strip(np.loadtxt(read_file, delimiter=",", max_rows=1, dtype=str, comments='%'))
        array = np.loadtxt(read_file, delimiter=",", skiprows=1)
        point_index = np.where(columns == 'point')[0][0]
        self.param_dict["points"] = int(np.max(array[:, point_index]))
        if self.param_dict["points"] > 1:
            array = np.loadtxt(read_file, delimiter=",", skiprows=2)
        row_count = int(np.shape(array)[0] / self.param_dict["points"])

        self.PHYSICAL_PARAMETERS = [p for p in PHYSICAL_PARAMETERS if p in columns]
        specname = wrap.get_specname()
        self.specname = [c for c in np.array([x.strip() for x in specname.astype(str) if x != ""]) if c in columns]

        self.physics_array = np.empty((row_count, self.param_dict["points"], len(self.PHYSICAL_PARAMETERS) + 1))
        self.chemical_abun_array = np.empty((row_count, self.param_dict["points"], len(self.specname)))
        for p in range(self.param_dict["points"]):
            self.physics_array[:, p, :] = \
            array[np.where(array[:, point_index] == p + 1), :len(self.PHYSICAL_PARAMETERS) + 1][0]
            self.chemical_abun_array[:, p, :] = \
            array[np.where(array[:, point_index] == p + 1), (len(self.PHYSICAL_PARAMETERS) + 1):][0]
        self._array_clean()
        last_timestep_index = self.physics_array[:, 0, 0].nonzero()[0][-1]
        self.next_starting_chemistry = self.chemical_abun_array[last_timestep_index, :, :]
        return

    def read_starting_chemistry_output_file(self):
        self.starting_chemistry = np.loadtxt(self.abundLoadFile, delimiter=",")

    def write_full_output_file(self):
        '''Perform classic output file writing.

        '''
        phys = self.physics_array.reshape(-1, self.physics_array.shape[-1])
        chem = self.chemical_abun_array.reshape(-1, self.chemical_abun_array.shape[-1])
        full_array = np.append(phys, chem, axis=1)
        string_fmt_string = f'{", ".join(["%10s"] * (len(self.PHYSICAL_PARAMETERS)))}, {", ".join(["%11s"] * len(self.specname.tolist()))}'
        # Magic numbers here to match/improve the formatting of the classic version
        number_fmt_string = f'%10.3E, %10.4E, %10.2f, %10.2f, %10.4E, %10.4E, %10.4E, %10i, {", ".join(["%9.5E"] * len(self.specname.tolist()))}'
        columns = np.array([self.PHYSICAL_PARAMETERS[:-1] + ["point"] + self.specname.tolist()])
        np.savetxt(self.outputFile, columns, fmt=string_fmt_string)
        with open(self.outputFile, "ab") as f:
            np.savetxt(f, full_array, fmt=number_fmt_string)
        return

    def write_starting_chemistry_output_file(self):
        last_timestep_index = self.chemical_abun_array[:, 0, 0].nonzero()[0][-1]
        number_fmt_string = f' {", ".join(["%9.5E"] * len(self.specname.tolist()))}'
        with open(self.abundSaveFile, "wb") as f:
            np.savetxt(f, self.chemical_abun_array[last_timestep_index, :, :], fmt=number_fmt_string)
        return

    def _array_clean(self):
        """
        Clean the array
        """
        # Find the first element with all the zeros
        last_timestep_index = self.physics_array[:, 0, 0].nonzero()[0][-1]
        # Get the arrays for only the simulated timesteps (not the zero padded ones)
        self.physics_array = self.physics_array[: last_timestep_index + 1, :, :]
        self.chemical_abun_array = self.chemical_abun_array[: last_timestep_index + 1, :, :]
        if self.ratesArray is not None:
            self.ratesArray = self.ratesArray[: last_timestep_index + 1, :, :]
        # Get the last arrays simulated, easy for starting another model.
        self.next_starting_chemistry = self.chemical_abun_array[last_timestep_index, :, :]
        return

    def _create_fortran_array(self):
        # fencepost problem, need to add 1 to timepoints to account for the 0th timestep
        self.physics_array = np.zeros(
            shape=(self.timepoints + 1, self.param_dict["points"], N_PHYSICAL_PARAMETERS),
            dtype=np.float64,
            order="F",
        )
        self.chemical_abun_array = np.zeros(
            shape=(self.timepoints + 1, self.param_dict["points"], n_species),
            dtype=np.float64,
            order="F",
        )
        return

    def _create_rates_array(self):
        self.ratesArray = np.zeros(shape=(self.timepoints + 1, self.param_dict["points"], n_reactions),
                                   dtype=np.float64, order="F")
        return

    def _reform_inputs(self,
                       param_dict: dict,
                       out_species: list):
        """Copies param_dict so as not to modify user's dictionary. Then reformats out_species from pythonic list
        to a string of space separated names for Fortran.
        """
        if param_dict is None:
            self.param_dict = {}
        else:
            # lower case (and conveniently copy so we don't edit) the user's dictionary
            # this is key to UCLCHEM's "case insensitivity"
            new_param_dict = {}
            for k, v in param_dict.items():
                assert (
                        k.lower() not in new_param_dict
                ), f"Lower case key {k} is already in the dict, stopping"
                if isinstance(v, Path):
                    v = str(v)
                new_param_dict[k.lower()] = v
            self.param_dict = new_param_dict.copy()
            del new_param_dict
        if out_species is not None:
            self.n_out = len(out_species)
            self.param_dict["outspecies"] = self.n_out
            self.out_species = " ".join(out_species)
        else:
            self.out_species = ""
            self.n_out = 0
        return

    def _xarray_conversion(self):
        return


class Cloud(AbstractModel):
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file
        )
        if read_file is None:
            self.run()

    def run(self):
        super().__run__()
        _, _, _, abunds, _, self.success_flag = wrap.cloud(
            dictionary=self.param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self.param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            abundancestart=self.starting_chemistry,
        )
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class Collapse(AbstractModel):
    def __init__(self,
                 collapse: str,
                 physics_output: str,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file
        )
        if read_file is None:
            collapse_dict = {"BE1.1": 1, "BE4": 2, "filament": 3, "ambipolar": 4}
            try:
                self.collapse = collapse_dict[collapse]
            except KeyError:
                raise ValueError(
                    f"collapse must be in {collapse_dict.keys()}"
                )
            self.physics_output = physics_output
            self.write_physics = self.physics_output is not None
            if not self.write_physics:
                self.physics_output = ""
            self.run()

    def run(self):
        super().__run__()
        _, _, _, abunds, _, self.success_flag = wrap.collapse(
            collapseIn=self.collapse,
            collapseFileIn=self.physics_output,
            writeOut=self.write_physics,
            dictionary=self.param_dict,
            outspeciesin=self.out_species,
            returnarray=True,
            returnrates=True,
            givesstartabund=self.give_start_abund,
            timepoints=self.timepoints,
            gridpoints=self.param_dict["points"],
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            abundanceStart=self.starting_chemistry,
        )
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class PrestellarCore(AbstractModel):
    def __init__(self,
                 temp_indx,
                 max_temperature,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file
        )
        if read_file is None:
            if temp_indx is None or max_temperature is None:
                raise ("temp_indx and max_temperature must be specified if not reading from file.")
            self.temp_indx = temp_indx
            self.max_temperature = max_temperature
            self.run()

    def run(self):
        super().__run__()
        _, _, _, abunds, _, self.success_flag = wrap.hot_core(
            temp_indx=self.temp_indx,
            max_temp=self.max_temperature,
            dictionary=self.param_dict,
            outspeciesin=self.out_species,
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            timepoints=self.timepoints,
            gridpoints=self.param_dict["points"],
            physicsarray=self.physics_array,
            ratesarray=self.ratesArray,
            chemicalabunarray=self.chemical_abun_array,
            abundancestart=self.starting_chemistry,
        )
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class CShock(AbstractModel):
    def __init__(self,
                 shock_vel: float = None,
                 timestep_factor: float = 0.01,
                 minimum_temperature: float = 0.0,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file
        )
        if read_file is None:
            if shock_vel is None:
                raise ("shock_vel must be specified if not reading from file.")
            self.shock_vel = shock_vel
            self.timestep_factor = timestep_factor
            self.minimum_temperature = minimum_temperature
            self.dissipation_time = -1
            self.run()

    def run(self):
        super().__run__()
        _, _, _, abunds, self.dissipation_time, _, self.success_flag = wrap.cshock(
            shock_vel=self.shock_vel,
            timestep_factor=self.timestep_factor,
            minimum_temperature=self.minimum_temperature,
            dictionary=self.param_dict,
            outspeciesin=self.out_species,
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            timepoints=self.timepoints,
            gridpoints=self.param_dict["points"],
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            abundancestart=self.starting_chemistry,
        )
        if self.success_flag < 0:
            self.dissipation_time = None
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class JShock(AbstractModel):
    def __init__(self,
                 shock_vel: float = None,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            timepoints,
            debug,
            read_file
        )
        if read_file is None:
            if shock_vel is None:
                raise ("shock_vel must be specified if not reading from file.")
            self.shock_vel = shock_vel
            self.run()

    def run(self):
        super().__run__()
        _, _, _, abunds, _, self.success_flag = wrap.jshock(
            shock_vel=self.shock_vel,
            dictionary=self.param_dict,
            outspeciesin=self.out_species,
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            timepoints=self.timepoints,
            gridpoints=self.param_dict["points"],
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            abundancestart=self.starting_chemistry,
        )
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class Model(AbstractModel):
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 time_array: np.array = None,
                 density_array: np.array = None,
                 gas_temperature_array: np.array = None,
                 dust_temperature_array: np.array = None,
                 zeta_array: np.array = None,
                 radfield_array: np.array = None,
                 debug: bool = False,
                 read_file: str = None
                 ):
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            len(time_array),
            debug,
            read_file
        )
        if read_file is None:
            self.postprocess_arrays = dict(
                timegrid=time_array,
                densgrid=density_array,
                gastempgrid=gas_temperature_array,
                dusttempgrid=dust_temperature_array,
                radfieldgrid=radfield_array,
                zetagrid=zeta_array
            )
            for key, array in self.postprocess_arrays.items():
                if array is not None:
                    # Convert single values into arrays that can be used
                    if isinstance(array, float):
                        array = np.ones(shape=time_array.shape) * array
                    # Assure lengths are correct
                    assert len(array) == len(time_array), "All arrays must be the same length"
                    # Ensure Fortran memory
                    array = np.asfortranarray(array, dtype=np.float64)
                    self.postprocess_arrays[key] = array
            self.time_array = time_array
            if not self.give_start_abund:
                self.starting_chemistry = np.zeros(
                    shape=n_species,
                    dtype=np.float64,
                    order="F",
                )
            self.run()

    def run(self):
        super().__run__()
        _, _, _, abunds, _, self.success_flag = wrap.postprocess(
            dictionary=self.param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self.param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            abundancestart=self.starting_chemistry,
            usecoldens=False,
            **self.postprocess_arrays,
        )
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class Postprocess(AbstractModel):
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = None,
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 time_array: np.array = None,
                 density_array: np.array = None,
                 gas_temperature_array: np.array = None,
                 dust_temperature_array: np.array = None,
                 zeta_array: np.array = None,
                 radfield_array: np.array = None,
                 coldens_H_array: np.array = None,
                 coldens_H2_array: np.array = None,
                 coldens_CO_array: np.array = None,
                 coldens_C_array: np.array = None,
                 debug: bool = False,
                 read_file: str = None
                 ):
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            len(time_array),
            debug,
            read_file
        )
        if read_file is None:
            self.postprocess_arrays = dict(
                timegrid=time_array,
                densgrid=density_array,
                gastempgrid=gas_temperature_array,
                dusttempgrid=dust_temperature_array,
                radfieldgrid=radfield_array,
                zetagrid=zeta_array,
                nhgrid=coldens_H_array,
                nh2grid=coldens_H2_array,
                ncogrid=coldens_CO_array,
                ncgrid=coldens_C_array,
            )
            for key, array in self.postprocess_arrays.items():
                if array is not None:
                    # Convert single values into arrays that can be used
                    if isinstance(array, float):
                        array = np.ones(shape=time_array.shape) * array
                    # Assure lengths are correct
                    assert len(array) == len(time_array), "All arrays must be the same length"
                    # Ensure Fortran memory
                    array = np.asfortranarray(array, dtype=np.float64)
                    self.postprocess_arrays[key] = array
            self.time_array = time_array
            self.coldens_H_array = coldens_H_array
            if not self.give_start_abund:
                self.starting_chemistry = np.zeros(
                    shape=n_species,
                    dtype=np.float64,
                    order="F",
                )
            self.run()

    def run(self):
        super().__run__()
        _, _, _, abunds, _, self.success_flag = wrap.postprocess(
            dictionary=self.param_dict,
            outspeciesin=self.out_species,
            timepoints=self.timepoints,
            gridpoints=self.param_dict["points"],
            returnarray=True,
            returnrates=True,
            givestartabund=self.give_start_abund,
            physicsarray=self.physics_array,
            chemicalabunarray=self.chemical_abun_array,
            abundancestart=self.starting_chemistry,
            usecoldens=self.coldens_H_array is not None,
            **self.postprocess_arrays,
        )
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()

