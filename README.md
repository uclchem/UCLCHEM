# UCL_CHEM
Fortran 95 version of uclchem for eventual public release

Eventually there will be a proper readme but for now here's the need to know stuff.

newnew will not run with the wrong number of species etc. Nspec and Nreac are given by Makerates. Ngas has to be worked out. It's either the sum of the gas species or Nspec-Ngrain-1. The latter is easier as you can grep for #s in the species file to count the number of grain species.

Currently, there is a parameter file called parameters.f90. As there is a makefile, recompiling after changing parameters takes a negligible amount of time. For that reason, having a parameter file that is just read in is fairly low priority. It will be in the final release.

Not everything is named exactly like old ucl_chem but most things are close. Let me know if anything isn't clear so I can add comments.


