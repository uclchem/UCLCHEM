# UCLCHEM
UCLCHEM is a gas-grain chemical code written in Fortran 95. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. Included in the repository is MakeRates, a python script to combine a species list, various reaction files into a consistent network with all files required by UCLCHEM.

**************************************************************
Usage Instructions
**************************************************************

A manual is available from the website: uclchem.github.io

*************************************************************
Contributing
*************************************************************
This is an open source science code for the community and we are happy to accept pull requests. We are also happy to work with you to produce a physics module if none of the ones available in the repository suit the modelling work you wish to do. If you are contributing ,please try to work with our current code style. We have the following general guidelines:

-camelCase variable and subroutines names that are self-explanatory where possible 
-CAPITALIZED fortran built in functions to make code structure apparent.

*************************************************************
General issues/ to do list
*************************************************************
-f2py wrapper
-move to DVODE module
-General code clean up