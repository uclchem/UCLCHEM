# Examples
This directory contains a few simple example runs of UCLCHEM. You can use it to test out the code, it's particularly helpful when you change network or modify the code because we focus on several major species so unless you've changed something absolutely fundamental, the abundances shouldn't differ from the examples.

## Contents
`example-output` - contains the output of UCLCHEM for a set of models. These can be taken as "correct" in the sense that they represent the output of a reasonable network from UCLCHEM working as described in the docs.

`test-output` - is a folder where your output files will go, we have various ways to generate equivalent files to those in `example-output`.

`fortran_cli` - Are files that document the way that one can interact with the LEGACY Fortran CLI. 


## How to Use This
In the `scripts/` directory, there is are two files: `scripts/run_uclchem_tests.py` and `scripts/plot_uclchem_tests.py`, run these both from the main directory after installing the code. This will first check your ODE conservation, then run three test cases of UCLCHEM. The plot_uclchem_tests.py script will plot the results of the test cases and compare them to the example outputs, creating a file called `examples/example-comparison.png` where you will be able to compare results between our examples and your tests. If they match well, you're all set! If they don't, you probably want to find out why rather than continuing.