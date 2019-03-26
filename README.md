# UCLCHEM v1.3
UCLCHEM is a gas-grain chemical code written in Fortran 95. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. Included in the repository is MakeRates, a python script to combine a species list, UMIST reaction file and user-defined reaction file into a consistent network with all files required by UCLCHEM.

**************************************************************
Usage Instructions
**************************************************************

A full manual is available from the website: uclchem.github.io

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
Support for python wrapping is limited. "Make python" builds a python library from the source code and wrap.f90. grid.py in the scripts folder runs a grid of models by repeatedly calling that python library. It is difficult to write a python wrapper for all possible use cases so wrap.f90 and grid.py are examples only.

 By changing the inputs and outputs of the subroutine in wrap.f90 or writing new subroutines, the user should be able to run anything they need. Set all parameters that are the same in every model in defaultparameters.f90 and then make variables into inputs to the subroutine.


**************************************************************
Change Log
**************************************************************
Various changes have been made with a view to improving readability of the code and improving sputtering in cshock.f90:
- Network moved to module allowing physics module to access the properties but NOT the abundances of species/reactions in the network.
- Sublimation subroutine added to all physics modules with a call in main.f90. 
- Thermal evaporation treatment moved to cloud.f90 sublimation subroutine
- cshock.f90 now has sputtering treatment from Jimenez-Serra et al. 2008
- Variable names changed to fit standard
- initializePhysics: sublimation flags etc reset so wrap.f90 can run in grids

*************************************************************
Contributing
*************************************************************
This is an open source science code for the community and we are happy to accept pull requests. We are also happy to work with you to produce a physics module if none of the ones available in the repository suit the modelling work you wish to do. If you are contributing ,please try to work with our current code style. We have the following general guidelines:

- camelCase variable and subroutines names that are self-explanatory where possible 

- CAPITALIZED fortran built in functions to make code structure apparent.


**************************************************************
To Do / General Inprovements
**************************************************************
- Create a keyword based python wrap in a similar to readparameters.f90 so that is more or less general purpose and usable by all.
- Trial 3 phase chemistry (gas, grain surface, bulk ice)