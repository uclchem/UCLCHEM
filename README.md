# UCLCHEM v1.1
UCLCHEM is a gas-grain chemical code written in Fortran 95. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. Included in the repository is MakeRates, a python script to combine a species list, UMIST reaction file and user-define reaction file into a consistent network with all files required by UCLCHEM.

**************************************************************
Usage Instructions
**************************************************************

A full manual is available from the website: uclchem.github.io

To build UCLCHEM, edit the Makefile in uclchem/src to use a compiler available on your machine. Then use "make" to create the executable.
- The Makefile also contains the choice of physics module.
- Building requires odes.f90 and network.f90 which are outputs of Makerates.
- uclchem/Makerates/ contains the Makerates python script to produce a network from the files in uclchem/Makerates/inputFiles

**************************************************************
Change Log
**************************************************************
Various changes have been made with the goal of making UCLCHEM more efficient and starting development of a python module. These changes include:
- Moving to F95 version of DVODE
- Removing all file inputs, the code needed to be compiled with odes.f90 anyway so reactions.csv and species.csv offered no benefits
- On going code clean up and bug fixes

*************************************************************
Contributing
*************************************************************
This is an open source science code for the community and we are happy to accept pull requests. We are also happy to work with you to produce a physics module if none of the ones available in the repository suit the modelling work you wish to do. If you are contributing ,please try to work with our current code style. We have the following general guidelines:

- camelCase variable and subroutines names that are self-explanatory where possible 

- CAPITALIZED fortran built in functions to make code structure apparent.

*************************************************************
General issues/ to do list
*************************************************************
- f2py wrapper