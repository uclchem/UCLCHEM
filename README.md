# UCL_CHEM
Fortran 95 version of uclchem for eventual public release

Eventually there will be a proper readme but for now here's the need to know stuff.

UCL_CHEM only works with my version of makerates. This writes nspec and nreac into the output files for array allocation as well as generating list of evaporation types.

Currently, there is a parameter file called parameters.f90. As there is a makefile, recompiling after changing parameters takes a negligible amount of time. For that reason, having a parameter file that is just read in is fairly low priority. It will be in the final release.

Not everything is named exactly like old ucl_chem but most things are close. Let me know if anything isn't clear so I can add comments.


