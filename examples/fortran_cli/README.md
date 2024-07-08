## Contents
Warning: this part of the code is not actively maintained as of now, it will be refactored at some point in the future and
only functions to keep use cases around until we can improve on them. DO NOT USE THIS INTERFACE FOR ANY SCIENCE.

`.inp files` - It is possible to run UCLCHEM as a binary. The inp files are how we set parameters when using the binary. 
    -   `static.inp` A static standalone model, it uses CLOUD
    -   `phase1.inp` A first stage model, it uses CLOUD
    -   `phase2.inp` It needs first stage to be run succesfully, it uses HOTCORE
    -   `shock.inp`  A shock model, it is set up also run after the first stage
    -   `postprocess.inp` This interfaces with the postprocessing code as proposed in [1]

## How to Use This
Go into the `fortran_src` directory, compile UCLCHEM as a binary using `make main`, go back into
the root UCLCHEM folder, and now you can run any of the codes using the following commands:
```bash
./uclchem CLOUD examples/fortran_cli/static.inp
./uclchem CLOUD examples/fortran_cli/phase1.inp
./uclchem HOTCORE examples/fortran_cli/phase2.inp 3 300.0
./uclchem CSHOCK examples/fortran_cli/shock.inp 20.0 0.01 10.0
./uclchem POSTPROCESS examples/fortran_cli/postprocess.inp
```


[1] Non-Equilibrium Abundances Treated Holistically (NEATH): the molecular composition of star-forming clouds 
     https://ui.adsabs.harvard.edu/link_gateway/2023MNRAS.524.5971P/doi:10.1093/mnras/stad2278
     https://github.com/fpriestley/neath/
