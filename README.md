# UCLCHEM v1.3
UCLCHEM is a gas-grain chemical code written in Fortran 95. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. Included in the repository is MakeRates, a python script to combine a species list, UMIST reaction file and user-defined reaction file into a consistent network with all files required by UCLCHEM.

**************************************************************
Usage Instructions
**************************************************************

Full documentation is available from the website: uclchem.github.io

To build UCLCHEM, edit the Makefile in uclchem/src to use a compiler available on your machine. Then use "make" to create the executable.
- The Makefile also contains the choice of physics module.
- Building requires odes.f90 and network.f90 which are outputs of Makerates.
- uclchem/Makerates/ contains the Makerates python script to produce a network from the files in uclchem/Makerates/inputFiles

To run UCLCHEM, create an input file with the desired parameters. Any variable in defaultparameters.f90 can be set, any that are not will take the value given in defaultparameters.f90.
A full explanation of each parameter is given in defaultparameters.f90 and an example input file is given in example.inp
Call uclchem with the filename as an argument: "./uclchem example.inp"

**************************************************************
Python
**************************************************************
"Make python" builds a python module from the source code and wrap.f90. This can be imported into a python script and any subroutine in wrap.f90 is a function in the module. Currently, this is just uclchem.general() which takes a dictionary of any parameters in defaultparameters.f90 and runs the code.
An example script, grid.py in the scripts folder runs a grid of models by repeatedly calling that function. This demonstrates the basic use of the wrapper and how to use python Pool objects to parallize runing a grid.

Currently, the wrapper must be recompiled for different physics modules.

**************************************************************
Change Log
**************************************************************
**New python wrapper**

**Updated grain treatment to make all grain parameters self-consistent**

**Continuous thermal desorption added.** Makerates has a new parameter therm_flag which defaults to False. If set to True, thermal desorption reactions are added to the network and material will continually desorb from the grains. UCLCHEM's physics modules are the preferred method for controlling thermal desorption but the functionality is available nonetheless.

**General Code Improvements**
- Small changes to the code have been made to improve readability. uvy is now uv_yield, grain_area is grain_area_per_H to better represent its physical meaning and constants have been moved to constants.f90
- An error message now informs users that the common error of writing a file to a folder that doesn't exist has occurred. Improving on the basic fortran error message which says the *file* doesn't exist.
- Double precision variables declared in correct fashion for modern fortran

*************************************************************
Contributing
*************************************************************
This is an open source science code for the community and we are happy to accept pull requests. We are also happy to work with you to produce a physics module if none of the ones available in the repository suit the modelling work you wish to do. If you are contributing ,please try to work with our current code style. We have the following general guidelines:

- camelCase variable and subroutines names that are self-explanatory where possible 

- CAPITALIZED fortran built in functions to make code structure apparent.