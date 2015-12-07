# UCL_CHEM
Fortran 95 version of uclchem for eventual public release


**************************************************************
Usage Instructions
**************************************************************
Eventually there will be a proper readme but for now here's the need to know stuff.

UCL_CHEM only works with my version of makerates. This writes nspec and nreac into the output files for array allocation as well as generating list of evaporation types.
	-Open Makerates folder and check Makerates.py to see if input files are correct
	-Run "python Makerates.py" and copy odes.f90,evaplists.csv,final_rates.csv and final_species.csv into main folder
	-check freeze out alphas in final_rates.csv, they default to 1
	-Nothing else needs to be changed, everything is updated automatically.

Whenever you change anything, including parameters type "make" into terminal. Currently, there is a parameter file called parameters.f90. As there is a makefile, recompiling after changing parameters takes a negligible amount of time. For that reason, having a parameter file that is just read in is fairly low priority.


*************************************************************
General issues/ to do list
*************************************************************
Freeze out
	Beta=1 for ions. Alpha is set as an 'efficiency' basically setting the fraction of a species
	that freezes down a particular route. Beta and gamma are then set to zero as they are useless.
	However, beta is set to 1 for ions so that ucl_chem can pick out the ions and freeze them out 
	differently. Either need to alert users or find a way to improve on this.

Zeros
	rate(j)=1.0d-30 when it is really zero. This was to avoid computing issues. However, in many cases a fortran 95 code on a modern pc will run without errors. Investigate this.



