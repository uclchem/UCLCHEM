import warnings
import os
from pathlib import Path
from abc import ABC

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.f2py.auxfuncs import throw_error

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


class AbstractModel(ABC):
    """Base model class used for inheritance only

    The AbstractModel class serves as an abstract class from which other model classes can be built. It is not intended
    to be used as a standalone class for running UCLCHEM.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values if not provided.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_specie_list: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: object = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        """Initiates the model with all common factors found within models"""
        self.param_dict = {}
        self.out_species_list = out_specie_list
        self.out_species = ""
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
            self.timepoints = None
            self.was_read = True
            self.PHYSICAL_PARAMETERS = None
            self.give_start_abund = False
            self.next_starting_chemistry = None
            self._reform_inputs(param_dict, self.out_species_list)
            self.read_output_file(read_file)
        else:
            self.n_out = 0
            self.timepoints = timepoints
            self.was_read = False
            self.PHYSICAL_PARAMETERS = PHYSICAL_PARAMETERS
            self.outputFile = param_dict.pop("outputFile") if "outputFile" in param_dict else None
            self.abundSaveFile = param_dict.pop("abundSaveFile") if "abundSaveFile" in param_dict else None
            self.abundLoadFile = param_dict.pop("abundLoadFile") if "abundLoadFile" in param_dict else None
            self._reform_inputs(param_dict, self.out_species_list)
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
        """__run__ resets the Fortran arrays if the model was not read, allowing the arrays to be reused for new runs."""
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
                           element_list: list = None,
                           percent: bool = True
                           ):
        """Utility method to check conservation of the chemical abundances

        Args:
            element_list (list, optional): List of elements to check conservation for. Defaults to self.out_species_lists.
            percent (bool, optional): Flag on if percentage values should be used. Defaults to True.
        """
        if element_list is None:
            element_list = self.out_species_list

        if self.param_dict["points"] > 1:
            for i in range(self.param_dict["points"]):
                print(f"Element conservation report for point {i + 1} of {self.param_dict['points']}")
                print(check_element_conservation(self.get_dataframes(i), element_list, percent))
        else:
            print(f"Element conservation report")
            print(check_element_conservation(self.get_dataframes(0), element_list, percent))

    def check_error(self, only_error: bool = False):
        """
        Prints the error message of the model based on self.success_flag, this method was originally an uclchem.utils function.

        Args:
            only_error (bool, optional): Flag to inform check_error to only print a message if self.success_flag was
                not 0
        """
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
        elif self.success_flag == 0 and not only_error:
            print(f'Model ran successfully.')
        elif self.success_flag is None:
            print(f'Model has not been run.')

    def create_abundance_plot(self,
                              element_list: list = None,
                              figsize: tuple[2] = (16, 9),
                              point: int = 0,
                              plot_file=None):
        """uclchem.analysis.create_abundance_plot wrapper method
        Args:
            element_list (list, optional): List of elements to check conservation for. Defaults to  self.out_species_list.
            figsize (tuple[2], optional): The figure size to use for matplotlib Defaults to (16, 9).
            point (int, optional): Integer referring to which point of the UCLCHEM model to use. Defaults to 0.
        Returns:
            fig,ax: matplotlib figure and axis objects
        """
        if element_list is None:
            element_list = self.out_species_list

        if point > self.param_dict["points"]:
            raise Exception("'point' must be less than number of modelled points.")
        return create_abundance_plot(self.get_dataframes(point), element_list, figsize, plot_file)

    def get_dataframes(self, point: int = 0, joined: bool = True, with_rates: bool = False):
        """Converts the model physics and chemical_abun arrays from numpy to pandas arrays.
        Args:
            point (int, optional): Integer referring to which point of the UCLCHEM model to return as pandas does not support higher
                than 2D data structures. Defaults to 0.
            joined (bool, optional): Flag on whether the returned pandas dataframe should be one, or if two dataframes should be
                returned. One physical, one chemical_abun dataframe. Defaults to True.
            with_rates (bool, optional): Flag on whether to include reaction rates in the dataframe, and/or as a separate
                dataframe depending on the value of `joined`. Defaults to False. #TODO Add the code to read rates.
        Returns:
            return_df (pandas.DataFrame): Dataframe of the joined arrays for point 'point' if joined = True
            physics_df (pandas.DataFrame): Dataframe of the physical parameters for point 'point' if joined = False
            chemistry_df (pandas.DataFrame): Dataframe of the chemical abundances  for point 'point' if joined = False
            rates_df (pandas.DataFrame): Dataframe of the reaction rates  for point 'point' if joined = False and with_rates = True
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

    def plot_species(self, ax: plt.axes, species: list[str] = None, point: int = 0, legend: bool = True, **plot_kwargs):
        """uclchem.analysis.plot(species) wrapper method
        Args:
            ax (plt.axes):
            species (list, optional):
            point (int, optional):
            legend (bool, optional):
            plot_kwargs (dict, optional):
        """
        if species is None:
            species = self.out_species_list
        return plot_species(ax, self.get_dataframes(point), species, legend, **plot_kwargs)

    def read_output_file(self, read_file: str, rates_load_file: str = None):
        '''Perform classic output file reading.
        Args:
            read_file (str): path to file containing a full UCLCHEM output
            rates_load_file (str, optional): path to file containing the reaction rates output from UCLCHEM. Defaults
                to None. #TODO Add the code to read the rates files.
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
        """Method to read the starting chemistry from the self.abundLoadFile provided in param_dict"""
        self.starting_chemistry = np.loadtxt(self.abundLoadFile, delimiter=",")

    def write_full_output_file(self):
        """Perform classic output file writing to the file self.outputFile provided in param_dict"""
        phys = self.physics_array.reshape(-1, self.physics_array.shape[-1])
        chem = self.chemical_abun_array.reshape(-1, self.chemical_abun_array.shape[-1])
        full_array = np.append(phys, chem, axis=1)
        #TODO Move away from the magic numbers seen here.
        string_fmt_string = f'{", ".join(["%10s"] * (len(self.PHYSICAL_PARAMETERS)))}, {", ".join(["%11s"] * len(self.specname.tolist()))}'
        # Magic numbers here to match/improve the formatting of the classic version
        #TODO Move away from the magic numbers seen here.
        number_fmt_string = f'%10.3E, %10.4E, %10.2f, %10.2f, %10.4E, %10.4E, %10.4E, %10i, {", ".join(["%9.5E"] * len(self.specname.tolist()))}'
        columns = np.array([self.PHYSICAL_PARAMETERS[:-1] + ["point"] + self.specname.tolist()])
        np.savetxt(self.outputFile, columns, fmt=string_fmt_string)
        with open(self.outputFile, "ab") as f:
            np.savetxt(f, full_array, fmt=number_fmt_string)
        return

    def write_starting_chemistry_output_file(self):
        """Perform classic starting abundance file writing to the file self.abundSaveFile provided in param_dict"""
        last_timestep_index = self.chemical_abun_array[:, 0, 0].nonzero()[0][-1]
        #TODO Move away from the magic numbers seen here.
        number_fmt_string = f' {", ".join(["%9.5E"] * len(self.specname.tolist()))}'
        with open(self.abundSaveFile, "wb") as f:
            np.savetxt(f, self.chemical_abun_array[last_timestep_index, :, :], fmt=number_fmt_string)
        return

    def _array_clean(self):
        """Internal Method.
        Clean the arrays changed by UCLCHEM Fortran code.
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
        """Internal Method.
        Creates Fortran compliant np.arrays that can be passed to the Fortran part of UCLCHEM.
        """
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
        """Internal Method.
        Creates Fortran compliant np.array for rates that can be passed to the Fortran part of UCLCHEM.
        """
        self.ratesArray = np.zeros(shape=(self.timepoints + 1, self.param_dict["points"], n_reactions),
                                   dtype=np.float64, order="F")
        return

    def _reform_inputs(self,
                       param_dict: dict,
                       out_species: list):
        """Internal Method.
        Copies param_dict so as not to modify user's dictionary. Then reformats out_species from pythonic list
        to a string of space separated names for Fortran.

        Args:
            param_dict (dict): Parameter dictionary passed by the user to the model.
            out_species (list): List of output species that are considered important for this model.
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
        """Internal Method. #TODO for Issue #94 add xarray support and conversion.
                Creates Fortran compliant np.array for rates that can be passed to the Fortran part of UCLCHEM.
        """
        return


class Cloud(AbstractModel):
    """Cloud model class inheriting from AbstractModel.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
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
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.__run__(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
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
        self.check_error(only_error=True)
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class Collapse(AbstractModel):
    """Collapse model class inheriting from AbstractModel.

    Args:
        collapse (str, optional):A string containing the collapse type, options are 'BE1.1', 'BE4', 'filament', or 'ambipolar'.
            Defaults to 'BE1.1'.
        physics_output (str, optional): Filename to store physics output, only relevant for 'filament' and 'ambipolar' collapses.
            If None, no physics output will be saved.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 collapse: str = "BE1.1",
                 physics_output: str = None,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        if collapse not in ["filament", "ambipolar"] and physics_output is None:
            warnings.warn("`physics_output` was None but `collapse` was `filament` or `ambipolar`. No output file will be created.", UserWarning)
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
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.__run__(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
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
        self.check_error(only_error=True)
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class PrestellarCore(AbstractModel):
    """PrestellarCore model class inheriting from AbstractModel. This model type was previously known as hot core.

    Args:
        temp_indx (int, optional): Used to select the mass of the prestellar core from the following selection
            [1=1Msun, 2=5, 3=10, 4=15, 5=25,6=60]. Defaults to 1, which is 1 Msun
        max_temperature (float, optional): Value at which gas temperature will stop increasing. Defaults to 300.0.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 temp_indx: int = 1,
                 max_temperature: float = 300.0,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
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
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.__run__(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
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
        self.check_error(only_error=True)
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class CShock(AbstractModel):
    """C-Shock model class inheriting from AbstractModel.

    Args:
        shock_vel (float, optional): Velocity of the shock in km/s. Defaults to 10.0.
        timestep_factor (float, optional): Whilst the time is less than 2 times the dissipation time of shock,
            timestep is timestep_factor*dissipation time. Essentially controls how well resolved the shock is
            in your model. Defaults to 0.01.
        minimum_temperature (float, optional): Minimum post-shock temperature. Defaults to 0.0 (no minimum). The
            shocked gas typically cools to `initialTemp` if this is not set.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 shock_vel: float = 10.0,
                 timestep_factor: float = 0.01,
                 minimum_temperature: float = 0.0,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
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
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.__run__(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
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
        self.check_error(only_error=True)
        if self.success_flag < 0:
            self.dissipation_time = None
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class JShock(AbstractModel):
    """J-Shock model class inheriting from AbstractModel.

    Args:
        shock_vel (float): Velocity of the shock. Defaults to 10.0.
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        timepoints (int, optional): Integer value of how many timesteps should be calculated before aborting the UCLCHEM model.
            Defaults to uclchem.constants.TIMEPOINTS
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 shock_vel: float = 10.0,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
                 starting_chemistry: np.ndarray = None,
                 previous_model: AbstractModel = None,
                 timepoints: int = TIMEPOINTS,
                 debug: bool = False,
                 read_file: str = None
                 ):
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
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
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.__run__(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
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
        self.check_error(only_error=True)
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class Postprocess(AbstractModel):
    """Postprocess represents a model class with additional controls. It inherits from AbstractModel.

    Postprocess allows for additional controls of the time, density, gas temperature, radiation field, cosmic ray
    ionisation rate, atomic and molecular Hydrogen, CO and C column densities through the use of arrays. Using these
    arrays allows for experimental model crafting beyond the standard models in other model classes.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        time_array (np.array, optional): Represents the time grid to be used for the model. This sets the target timesteps for
            which outputs will be stored.
        density_array (np.array, optional): Represents the value of the density at different timepoints found in time_array.
        gas_temperature_array (np.array, optional): Represents the value of the gas temperature at different timepoints found in time_array.
        dust_temperature_array (np.array, optional):Represents the value of the dust temperature at different timepoints found in time_array.
        zeta_array (np.array, optional): Represents the value of the cosmic ray ionisation rate at different timepoints found in time_array.
        radfield_array (np.array, optional): Represents the value of the UV radiation field at different timepoints found in time_array.
        coldens_H_array (np.array, optional): Represents the value of the column density of H at different timepoints found in time_array.
        coldens_H2_array (np.array, optional): Represents the value of the column density of H2 at different timepoints found in time_array.
        coldens_CO_array (np.array, optional): Represents the value of the column density of CO at different timepoints found in time_array.
        coldens_C_array (np.array, optional): Represents the value of the column density of C at different timepoints found in time_array.
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
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
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            len(time_array),
            debug,
            read_file
        )
        if read_file is None and time_array is not None:
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
        elif time_array is None and read_file is None:
            throw_error(f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given.")

    def run(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.__run__(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
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
        self.check_error(only_error=True)
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


class Model(AbstractModel):
    """Model, like Postprocess, represents a model class with additional controls. It inherits from AbstractModel.

    Model follows the same logic as Postprocess but without the coldens Arguments. It allows for additional controls of
    the time, density, gas temperature, radiation field, and cosmic ray ionisation rate through the use of arrays.
    Using these arrays allows for experimental model crafting beyond the standard models in other model classes.

    Args:
        param_dict (dict, optional): Dictionary containing the parameters to use for the UCLCHEM model. Uses UCLCHEM
            default values found in defaultparameters.f90.
        out_species (list, optional): List of chemicals to focus on for outputs such as conservation check, if no other values are
            provided. Defaults to ["H", "N", "C", "O"].
        starting_chemistry (np.ndarray, optional): Numpy ndarray containing the starting abundances to use for the UCLCHEM model.
            Defaults to None.
        previous_model (object, optional): Model object, a class that inherited from AbstractModel, to use for the starting abundances
            of the new UCLCHEM model that will be run. Defaults to None
        time_array (np.array, optional): Represents the time grid to be used for the model. This sets the target timesteps for
            which outputs will be stored. While listed as optional, this is only done to allow
        density_array (np.array, optional): Represents the value of the density at different timepoints found in time_array.
        gas_temperature_array (np.array, optional): Represents the value of the gas temperature at different timepoints found in time_array.
        dust_temperature_array (np.array, optional):Represents the value of the dust temperature at different timepoints found in time_array.
        zeta_array (np.array, optional): Represents the value of the cosmic ray ionisation rate at different timepoints found in time_array.
        radfield_array (np.array, optional): Represents the value of the UV radiation field at different timepoints found in time_array.
        debug (bool, optional): Flag if extra debug information should be printed to the terminal. Defaults to False. #TODO Add debug features
        read_file (str, optional): Path to the file to be read. Reading a file to a model object, prevents it from
            being run. Defaults to None.
    """
    def __init__(self,
                 param_dict: dict = None,
                 out_species: list = ["H", "N", "C", "O"],
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
        """Initiates the model first with AbstractModel.__init__(), then with any additional commands needed for the model."""
        super().__init__(
            param_dict,
            out_species,
            starting_chemistry,
            previous_model,
            len(time_array),
            debug,
            read_file
        )
        if read_file is None and time_array is not None:
            self.time_array = time_array
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
            if not self.give_start_abund:
                self.starting_chemistry = np.zeros(
                    shape=n_species,
                    dtype=np.float64,
                    order="F",
                )
            self.run()
        elif time_array is None and read_file is None:
            throw_error(f"time_array must be an array if read_file is None. A value of {time_array} with type {type(time_array)} was given.")

    def run(self):
        """
        Runs the UCLCHEM model, first by resetting the np.arrays by using AbstractModel.__run__(), then running the model.
        check_error, and array_clean are automatically called post model run.
        """
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
        self.check_error(only_error=True)
        self._array_clean()
        if self.outputFile is not None:
            self.write_full_output_file()
        if self.abundSaveFile is not None:
            self.write_starting_chemistry_output_file()


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
        physicalParameterArray, chemicalAbundanceArray, specname, ratesArray=None

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
    warnings.warn("The cloud function is deprecated and may not be up to date, please use the Cloud class.",
                  DeprecationWarning)
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
    warnings.warn("The collapse function is deprecated and may not be up to date, please use the Collapse class.",
                  DeprecationWarning)
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
    warnings.warn("The hot_core function is deprecated and may not be up to date, please use the PrestellarCore class.",
                  DeprecationWarning)
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
    warnings.warn("The cshock function is deprecated and may not be up to date, please use the CShock class.",
                  DeprecationWarning)
    n_out, param_dict, out_species = _reform_inputs(param_dict, out_species)
    if "points" not in param_dict:
        param_dict["points"] = 1
    if return_array or return_dataframe:
        _return_array_checks(param_dict)
    physicsArray, chemicalAbunArray = _create_fortranarray(
        param_dict, N_PHYSICAL_PARAMETERS, timepoints=timepoints
    )
    ratesArray = _create_ratesarray(param_dict["points"], n_reactions, timepoints=timepoints)
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
    warnings.warn("The jshock function is deprecated and may not be up to date, please use the JShock class.",
                  DeprecationWarning)
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
    warnings.warn("The postprocess function is deprecated and may not be up to date, please use the Postprocess class.",
                  DeprecationWarning)
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
