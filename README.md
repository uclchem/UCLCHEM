# UCL_CHEM
UCLCHEM is a gas-grain chemical code written in Fortran 95. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. Included in the repository is MakeRates, a python script to combine a species list, UMIST reaction file and user-define reaction file into a consistent network with all files required by UCLCHEM.

**************************************************************
Usage Instructions
**************************************************************

A manual is available from the website: uclchem.github.io


*************************************************************
General issues/ to do list
*************************************************************
Zeros
	-rate(j)=1.0d-30 when it is really zero. This was to avoid computing issues. However, in many cases a fortran 95 code on a modern pc will run without errors. This could be changed.
