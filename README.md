# UCLCHEM v3.0
`UCLCHEM` is a gas-grain chemical code that propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. We provide several physical models to enable the modelling of different astrophysical environments and a utility script `MakeRates` to help the user produce a chemical network from simple lists of reactions and species.


**************************************************************
## Installation Instructions
**************************************************************

Full documentation is available from the website: [uclchem.github.io](https://uclchem.github.io)

UCLCHEM is intended to be used as a python module. Once you've run MakeRates, you can simply run
```
pip install .
```
in the main directory to install Uclchem. You can then `import uclchem` in any python script.

To see the contents of this python module, check our [Python API docs](https://uclchem.github.io/docs/pythonapi). To see some example notebooks, check the tutorial section of the docs or the notebooks in `Tutorials/`.


If you want to build an executable from the Fortran source, head to `src/fortran_src` and run `make`.

### Prerequisites
To build UCLCHEM, you'll need gfortran and Cmake.

To run the python module, you'll need the python modules listed in `requirements.txt`


**************************************************************
## Change Log
**************************************************************
See change.log! We've made a large number of improvements for v3.0. The code has been restructured to be Python first in its intended use, different physical models can be accessed without recompilation, and MakeRates is more helpful than ever before.

*************************************************************
## Contributing
*************************************************************
This is an open source science code for the community and we are happy to accept pull requests. We are also happy to work with you to produce a physics module if none of the ones available in the repository suit the modelling work you wish to do. If you are contributing ,please try to work with our current code style. We have the following general guidelines:

- camelCase variable and subroutines names that are self-explanatory where possible 

- CAPITALIZED fortran built in functions to make code structure apparent.
