import warnings
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
        if read_file is not None:
            self.n_out = None
            self.param_dict = {}
            self.out_species = None
            self.timepoints = None
            self.physics_array = None
            self.chemical_abun_array = None
            self.ratesArray = None
            self.debug = debug
            self.success_flag = None
            self.full_array = None
            self.was_read = True
            self.PHYSICAL_PARAMETERS = None

            specname = wrap.get_specname()
            self.specname = np.array([x.strip() for x in specname.astype(str) if x != ""])

            self.give_start_abund = False
            self.next_starting_chemistry = None
            self.read_output_file(read_file)
        else:
            self.n_out = 0
            self.param_dict = {}
            self.out_species = []
            self.timepoints = timepoints
            self.physics_array = None
            self.chemical_abun_array = None
            self.ratesArray = None
            self.debug = debug
            self.success_flag = None
            self.was_read = False
            self.PHYSICAL_PARAMETERS = PHYSICAL_PARAMETERS

            specname = wrap.get_specname()
            self.specname = np.array([x.strip() for x in specname.astype(str) if x != ""])

            self._reform_inputs(param_dict, out_species)
            if "points" not in self.param_dict:
                self.param_dict["points"] = 1
            if previous_model is None:
                self.starting_chemistry = starting_chemistry
            elif 'next_starting_chemistry' in [attr for attr in dir(previous_model) if
                                               not callable(getattr(previous_model, attr)) and
                                               not attr.startswith("__") and
                                               attr is not None]:
                self.starting_chemistry = previous_model.next_starting_chemistry
            self.give_start_abund = self.starting_chemistry is not None
            self.next_starting_chemistry = None
        return

    def run(self,):
        if not self.was_read:
            self.physics_array = None
            self.chemical_abun_array = None
            self.ratesArray = None
            self._create_fortran_array()
            self._create_rates_array()
        else:
            raise("This model was read. It can not be run. ")
        return

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


    def check_conservation(self,
                            element_list: list = ["H", "N", "C", "O"],
                            percent: bool = True
                            ):
        if self.param_dict["points"] > 1:
            for i in range(self.param_dict["points"]):
                print(f"Element conservation report for point {i+1} of {self.param_dict['points']}")
                print(check_element_conservation(self.get_dataframes(i), element_list, percent))
        else:
            print(f"Element conservation report")
            print(check_element_conservation(self.get_dataframes(0), element_list, percent))

    def create_abundance_plot(self,
                       species: list = ["H", "N", "C", "O"],
                       figsize: tuple[2] = (10, 7),
                       point: int = 0,
                       plot_file = None):
        if point > self.param_dict["points"]:
            raise Exception("'point' must be less than number of modelled points.")
        return create_abundance_plot(self.get_dataframes(point), species, figsize, plot_file)

    def plot_species(self, ax: plt.axes, species: list[str], point: int = 0, legend: bool=True, **plot_kwargs):
        return plot_species(ax, self.get_dataframes(point), species, legend, **plot_kwargs)

    def write_output_file(self, output_file: str):
        '''Perform classic output file writing.

        '''

        phys = self.physics_array.reshape(-1, self.physics_array.shape[-1])
        chem = self.chemical_abun_array.reshape(-1, self.chemical_abun_array.shape[-1])

        appended = np.append(phys, chem, axis=1)
        df = pd.DataFrame(appended, columns=[self.PHYSICAL_PARAMETERS + ["point"] + self.specname])
        np.savetxt(output_file, np.array([a[0] for a in df.columns.to_list()])[None, :], delimiter=',\t', fmt='%s')
        with open(output_file, "ab") as f:
            np.savetxt(f, df, delimiter=',\t', fmt='%.5g')
        return

    def read_output_file(self, abund_load_file: str, rates_load_file: str = None):
        '''Perform classic output file reading.
        Args:
            abund_load_file (str): path to file containing a full UCLCHEM output
        '''
        self.was_read = True
        loaded = pd.read_csv(abund_load_file)
        loaded.columns = loaded.columns.str.strip()
        self.param_dict["points"] = loaded["point"].max()
        row_count = int(len(loaded) / self.param_dict["points"]) + 1
        self.PHYSICAL_PARAMETERS = [p for p in PHYSICAL_PARAMETERS if p in loaded.columns]
        self.specname = [c for c in self.specname if c in loaded.columns]
        phys_output_array = np.empty((row_count, self.param_dict["points"], len(self.PHYSICAL_PARAMETERS)+1))
        chem_output_array = np.empty((row_count,  self.param_dict["points"], len(self.specname)))
        for p in range(self.param_dict["points"]):
            if p == 0:
                phys_output_array[:, p, :] = loaded[loaded['point'] == (p + 1)][self.PHYSICAL_PARAMETERS + ["point"]].to_numpy()
                chem_output_array[:, p, :] = loaded[loaded['point'] == (p + 1)][self.specname].to_numpy()
            else:
                phys_output_array[:, p, :] = np.append(np.zeros((1, len(self.PHYSICAL_PARAMETERS + ["point"]))),
                                                       loaded[loaded['point'] == (p + 1)][self.PHYSICAL_PARAMETERS + ["point"]].to_numpy(), axis=0)
                chem_output_array[:, p, :] = np.append(np.zeros((1, len(self.specname))),
                                                       loaded[loaded['point'] == (p + 1)][self.specname].to_numpy(), axis=0)
        self.physics_array = phys_output_array
        self.chemical_abun_array = chem_output_array

        self._array_clean()
        last_timestep_index = self.physics_array[:, 0, 0].nonzero()[0][-1]
        self.next_starting_chemistry = self.chemical_abun_array[last_timestep_index, 0, :]
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
        self.ratesArray = np.zeros(shape=(self.timepoints + 1, self.param_dict["points"], n_reactions), dtype=np.float64, order="F")
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
        self.next_starting_chemistry = self.chemical_abun_array[last_timestep_index, 0, :]
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
        super().run()
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
        super().run()
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
                raise("temp_indx and max_temperature must be specified if not reading from file.")
            self.temp_indx = temp_indx
            self.max_temperature = max_temperature

            self.run()

    def run(self):
        super().run()
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
            chemicalabunarray=self.chemical_abun_array,
            abundancestart=self.starting_chemistry,
        )
        self._array_clean()


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
        super().run()
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
        super().run()
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
        super().run()
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
        super().run()
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


'''
Beyond here lie deprecated functions for model.py
'''



def _reform_inputs(param_dict, out_species):
    """Copies param_dict so as not to modify user's dictionary. Then reformats out_species from pythonic list
    to a string of space separated names for Fortran.
    """
    if param_dict is None:
        param_dict = {}
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
        param_dict = new_param_dict.copy()
        del new_param_dict
    if out_species is not None:
        n_out = len(out_species)
        param_dict["outspecies"] = n_out
        out_species = " ".join(out_species)
    else:
        out_species = ""
        n_out = 0
    return n_out, param_dict, out_species


def _format_output(n_out, abunds, success_flag):
    if success_flag < 0 or n_out == 0:
        abunds = []
    else:
        abunds = list(abunds[:n_out])
    return [success_flag] + abunds


def _create_fortranarray(param_dict, nPhysParam, timepoints=TIMEPOINTS):
    # fencepost problem, need to add 1 to timepoints to account for the 0th timestep
    physicsArray = np.zeros(
        shape=(timepoints + 1, param_dict["points"], nPhysParam),
        dtype=np.float64,
        order="F",
    )
    chemicalAbunArray = np.zeros(
        shape=(timepoints + 1, param_dict["points"], n_species),
        dtype=np.float64,
        order="F",
    )
    return physicsArray, chemicalAbunArray


def _create_ratesarray(points, nReacs, timepoints=TIMEPOINTS):
    return np.zeros(shape=(timepoints + 1, points, nReacs), dtype=np.float64, order="F")


def _return_array_checks(params):
    if any([key.endswith("File") for key in list(params.keys())]):
        raise RuntimeError(
            "return_array or return_dataframe cannot be used if any output of input file is specified."
        )


def _array_clean(
    physicalParameterArray, chemicalAbundanceArray, specname, ratesArray = None
    
):
    """Clean the array

    Args:
        physicalParameterArray (np.ndarray): Array with the UCLCHEM physical parameters
        chemicalAbundanceArray (np.ndarray): Array with the output chemical abundances
        specname (np.ndarray): Numpy array with the names of all the species
        nPhysParam (int): The number of physical parameters you are interested in.

    Returns:
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray:
        The physical parameters, the abundances, overtime species names, the last
        nonzero physical parameters, the last abundances respectively.
    """
    specname_new = specname.astype(str)
    specname_new = np.array([x.strip() for x in specname_new if x != ""])

    # Find the first element with all the zeros
    last_timestep_index = physicalParameterArray[:, 0, 0].nonzero()[0][-1]
    # Get the arrays for only the simulated timesteps (not the zero padded ones)
    physicsArray = physicalParameterArray[: last_timestep_index + 1, :, :]
    chemArray = chemicalAbundanceArray[: last_timestep_index + 1, :, :]
    # Also clean the rates array, only if we have it:
    if ratesArray is not None:
        ratesArray = ratesArray[: last_timestep_index + 1, :, :]
    # Get the last arrays simulated, easy for starting another model.
    abundanceStart = chemicalAbundanceArray[last_timestep_index, 0, :]
    return physicsArray, chemArray, specname_new, ratesArray, abundanceStart


def outputArrays_to_DataFrame(
    physicalParameterArray, chemicalAbundanceArray, specname, ratesArray, 
):
    """Convert the output arrays to a pandas dataframe

    Args:
        physicalParameterArray (np.array): Array with the output physical parameters
        chemicalAbundanceArray (np.array): Array with the output chemical abundances
        specname (list): List with the names of all the species
        physParameter (list): Array with all the physical parameter names

    Returns:
        _type_: _description_
    """
    # Create a physical parameter dataframe
    physics_df = pd.DataFrame(
        physicalParameterArray[:, 0, : N_PHYSICAL_PARAMETERS],
        index=None,
        columns=PHYSICAL_PARAMETERS,
    )
    # Create a abundances dataframe.
    chemistry_df = pd.DataFrame(
        chemicalAbundanceArray[:, 0, :], index=None, columns=specname
    )
    if ratesArray is not None:
        # Create a rates dataframe.
        rates_df = pd.DataFrame(
            ratesArray[:, 0, :], index=None, columns=get_reaction_names()
        )
    else: 
        rates_df = None
    return physics_df, chemistry_df, rates_df


def cloud(
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    warnings.warn("The cloud function is deprecated and may not be up to date, please use the Cloud class.", DeprecationWarning)
    give_start_abund = starting_chemistry is not None
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    # Check to make sure no output files are specified, if so, halt the execution.
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=timepoints
    )
    ratesArray = _create_ratesarray(param_dict["points"], n_reactions, timepoints=timepoints)
    _, _, _, abunds, specname, success_flag = wrap.cloud(
        dictionary=param_dict,
        outspeciesin=out_species,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, ratesArray, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
        ) + (abundanceStart, success_flag,)
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def collapse(
    collapse,
    physics_output,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run collapse model from UCLCHEM based on Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)

    Args:
        collapse (str): A string containing the collapse type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'
        physics_output(str): Filename to store physics output, only relevant for 'filament' and 'ambipolar' collapses. If None, no physics output will be saved.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    warnings.warn("The collapse function is deprecated and may not be up to date, please use the Collapse class.", DeprecationWarning)
    collapse_dict = {"BE1.1": 1, "BE4": 2, "filament": 3, "ambipolar": 4}
    try:
        collapse = collapse_dict[collapse]
    except KeyError:
        raise ValueError(
            "collapse must be one of 'BE1.1', 'BE4', 'filament', or 'ambipolar'"
        )
    write_physics = physics_output is not None
    if not write_physics:
        physics_output = ""
    give_start_abund = starting_chemistry is not None

    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, len(PHYSICAL_PARAMETERS), timepoints=timepoints)
    ratesArray = _create_ratesarray(param_dict["points"], n_reactions, timepoints=timepoints)
    _, _, _, abunds, specname, success_flag = wrap.collapse(
        collapseIn=collapse,
        collapseFileIn=physics_output,
        writeOut=write_physics,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givesstartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        abundanceStart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, ratesArray, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
        ) + (abundanceStart, success_flag,)
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def hot_core(
    temp_indx,
    max_temperature,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run hot core model from UCLCHEM, based on Viti et al. 2004 and Collings et al. 2004

    Args:
        temp_indx (int): Used to select the mass of hot core. 1=1Msun,2=5, 3=10, 4=15, 5=25,6=60]
        max_temperature (float): Value at which gas temperature will stop increasing.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    warnings.warn("The hot_core function is deprecated and may not be up to date, please use the PrestellarCore class.", DeprecationWarning)
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    if return_array or return_dataframe:
        # Check to make sure no output files are specified, if so, halt the execution.
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, len(PHYSICAL_PARAMETERS), timepoints=timepoints
    )
    ratesArray = _create_ratesarray(param_dict["points"], n_reactions, timepoints=timepoints)
    give_start_abund = starting_chemistry is not None
    _, _, _, abunds, specname, success_flag = wrap.hot_core(
        temp_indx=temp_indx,
        max_temp=max_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        ratesarray=ratesArray,
        chemicalabunarray=chemicalAbunArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, ratesArray, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
        ) + (abundanceStart, success_flag,)
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def cshock(
    shock_vel,
    timestep_factor=0.01,
    minimum_temperature=0.0,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run C-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock in km/s
        timestep_factor (float, optional): Whilst the time is less than 2 times the dissipation time of shock, timestep is timestep_factor*dissipation time. Essentially controls
        how well resolved the shock is in your model. Defaults to 0.01.
        minimum_temperature (float, optional): Minimum post-shock temperature. Defaults to 0.0 (no minimum). The shocked gas typically cools to `initialTemp` if this is not set.
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, the second element is the dissipation time and further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - disspation_time (float): dissipation time in years
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - disspation_time (float): dissipation time in years
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    warnings.warn("The cshock function is deprecated and may not be up to date, please use the CShock class.", DeprecationWarning)
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=timepoints
    )
    ratesArray=_create_ratesarray(param_dict["points"], n_reactions, timepoints=timepoints)
    give_start_abund = starting_chemistry is not None
    _, _, _, abunds, disspation_time, specname, success_flag = wrap.cshock(
        shock_vel=shock_vel,
        timestep_factor=timestep_factor,
        minimum_temperature=minimum_temperature,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        ratesarray=ratesArray,
        chemicalabunarray=chemicalAbunArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates:
        ratesArray = None
    if success_flag < 0:
        disspation_time = None
        abunds = []
    else:
        abunds = list(abunds[:n_out])
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, ratesArray, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
        ) + (disspation_time, abundanceStart, success_flag,)
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            disspation_time,
            abundanceStart,
            success_flag,
        )
    else:
        if success_flag < 0:
            disspation_time = None
            abunds = []
        else:
            abunds = list(abunds[:n_out])
        result = [success_flag, disspation_time] + abunds
        return result


def jshock(
    shock_vel,
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    starting_chemistry=None,
    timepoints=TIMEPOINTS,
):
    """Run J-type shock model from UCLCHEM

    Args:
        shock_vel (float): Velocity of the shock
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem

    Returns:if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the model succeeded, the second element is the dissipation time and further elements are the abundances of all species in `out_species`.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.

    """
    warnings.warn("The jshock function is deprecated and may not be up to date, please use the JShock class.", DeprecationWarning)
    if "points" not in param_dict:
        param_dict["points"] = 1
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=timepoints
    )
    give_start_abund = starting_chemistry is not None
    ratesArray = _create_ratesarray(param_dict["points"], n_reactions, timepoints=timepoints)

    _, _, _, abunds, specname, success_flag = wrap.jshock(
        shock_vel=shock_vel,
        dictionary=param_dict,
        outspeciesin=out_species,
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        timepoints=timepoints,
        gridpoints=param_dict["points"],
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        abundancestart=starting_chemistry,
    )
    # Overwrite the ratesArray with None if its not used:
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, ratesArray, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, ratesArray,
        )
    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
        ) + (abundanceStart, success_flag,)
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)


def postprocess(
    param_dict=None,
    out_species=None,
    return_array=False,
    return_dataframe=False,
    return_rates=False,
    starting_chemistry=None,
    time_array=None,
    density_array=None,
    gas_temperature_array=None,
    dust_temperature_array=None,
    zeta_array=None,
    radfield_array=None,
    coldens_H_array=None,
    coldens_H2_array=None,
    coldens_CO_array=None,
    coldens_C_array=None,
):
    """Run cloud model from UCLCHEM

    Args:
        param_dict (dict,optional): A dictionary of parameters where keys are any of the variables in defaultparameters.f90 and values are value for current run.
        out_species (list, optional): A list of species for which final abundance will be returned. If None, no abundances will be returned.. Defaults to None.
        return_array (bool, optional): A boolean on whether a np.array should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        return_dataframe (bool, optional): A boolean on whether a pandas.DataFrame should be returned to a user, if both return_array and return_dataframe are false, this function will default to writing outputs to a file
        starting_chemistry (array, optional): np.array containing the starting chemical abundances needed by uclchem
    Returns:
        if return_array and return_dataframe are False:
            - A list where the first element is always an integer which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details. If the `out_species` parametere is provided, the remaining elements of this list will be the final abundances of the species in out_species.
        if return_array is True:
            - physicsArray (array): array containing the physical outputs for each written timestep
            - chemicalAbunArray (array): array containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
        if return_dataframe is True:
            - physicsDF (pandas.DataFrame): DataFrame containing the physical outputs for each written timestep
            - chemicalDF (pandas.DataFrame): DataFrame containing the chemical abundances for each written timestep
            - abundanceStart (array): array containing the chemical abundances of the last timestep in the format uclchem needs in order to perform an additional run after the initial model
            - success_flag (integer): which is negative if the model failed to run and can be sent to `uclchem.utils.check_error()` to see more details.
    """
    warnings.warn("The postprocess function is deprecated and may not be up to date, please use the Postprocess class.", DeprecationWarning)
    # Assure that every array is cast to a fortran array
    postprocess_arrays = dict(
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
    for key, array in postprocess_arrays.items():
        if array is not None:
            # Convert single values into arrays that can be used
            if isinstance(array, float):
                array = np.ones(shape=time_array.shape) * array
            # Assure lengths are correct
            assert len(array) == len(time_array), "All arrays must be the same length"
            # Ensure Fortran memory
            array = np.asfortranarray(array, dtype=np.float64)
            postprocess_arrays[key] = array
    give_start_abund = starting_chemistry is not None
    if not give_start_abund:
        starting_chemistry = np.zeros(
            shape=(n_species),
            dtype=np.float64,
            order="F",
        )
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    # Check to make sure no output files are specified, if so, halt the execution.
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=len(time_array)
    )
    ratesArray = _create_ratesarray(param_dict["points"], n_reactions, timepoints=len(time_array))
    _, _, _, abunds, specname, success_flag = wrap.postprocess(
        dictionary=param_dict,
        outspeciesin=out_species,
        timepoints=len(time_array),
        gridpoints=param_dict["points"],
        returnarray=return_array or return_dataframe,
        returnrates=return_rates,
        givestartabund=give_start_abund,
        physicsarray=physicsArray,
        chemicalabunarray=chemicalAbunArray,
        ratesarray=ratesArray,
        abundancestart=starting_chemistry,
        usecoldens=coldens_H_array is not None,
        **postprocess_arrays,
    )
    if not return_rates or not (return_array or return_dataframe):
        ratesArray = None
    if return_array or return_dataframe:
        physicsArray, chemicalAbunArray, specname, abundanceStart = _array_clean(
            physicsArray, chemicalAbunArray, specname, N_PHYSICAL_PARAMETERS
        )

    if return_dataframe:
        return outputArrays_to_DataFrame(
            physicsArray,
            chemicalAbunArray,
            specname,
            ratesArray,
        ) + (abundanceStart, success_flag,)
    elif return_array:
        return (
            physicsArray,
            chemicalAbunArray,
            ratesArray,
            abundanceStart,
            success_flag,
        )
    else:
        return _format_output(n_out, abunds, success_flag)
        return _format_output(n_out, abunds, success_flag)
