# UCL_CHEM
UCLCHEM is a gas-grain chemical code written in Fortran 95. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. Included in the repository is MakeRates, a python script to combine a species list, UMIST reaction file and user-define reaction file into a consistent network with all files required by UCLCHEM.

**************************************************************
Usage Instructions
**************************************************************

UCL_CHEM only works with my version of makerates. This writes nspec and nreac into the output files for array allocation as well as generating list of evaporation types.
	-Open Makerates folder and check Makerates.py to see if input files are correct
	-Run "python Makerates.py" and copy odes.f90,evaplists.csv,reactions.csv and species.csv into main folder
	-check freeze out alphas in reactions.csv, they default to 1
	-Nothing else needs to be changed, everything is updated automatically.


*************************************************************
General issues/ to do list
*************************************************************
Zeros
	rate(j)=1.0d-30 when it is really zero. This was to avoid computing issues. However, in many cases a fortran 95 code on a modern pc will run without errors. Investigate this.

