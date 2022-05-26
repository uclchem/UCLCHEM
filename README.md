# UCLCHEM v3.0
`UCLCHEM` is a gas-grain chemical code that propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. We provide several physical models to enable the modelling of different astrophysical environments and a utility script `MakeRates` to help the user produce a chemical network from simple lists of reactions and species.


**************************************************************
## Installation Instructions
**************************************************************

Full documentation is available from the website: [uclchem.github.io](https://uclchem.github.io)

UCLCHEM is intended to be used as a python module but must be installed from source rather than an online index such as Pypi. This is because users are expected to modify the source code, at least by creating their own networks. To obtain and install the code simply run:

```bash
git clone https://github.com/uclchem/UCLCHEM.git
cd UCLCHEM
pip install -r requirements.txt
pip install .
```

You can then `import uclchem` in any python script. You need to `pip install .` whenever you change your network. 

To see the contents of this python module, check our [Python API docs](https://uclchem.github.io/docs/pythonapi). To see some example notebooks, check the tutorial section of the docs or the notebooks in `Tutorials/`.


If you want to build an executable from the Fortran source, head to `src/fortran_src` and run `make`. You can then run the executable with `./uclchem CLOUD input_file.inp` where there examples of input files in the `examples/` directory. We do not suggest users use the code this way unt

### Prerequisites
To build UCLCHEM, you'll need gfortran, make and python 3.6+.

To run the python module, you'll need the python modules listed in `requirements.txt`


**************************************************************
## Change Log
**************************************************************
See change.log! We've made a large number of improvements for v3.0. The code has been restructured to be Python first in its intended use, different physical models can be accessed without recompilation, and MakeRates is more helpful than ever before.

*************************************************************
## Contributing
*************************************************************
This is an open source science code for the community and are open to pull requests. We are also happy to work with you to produce a physics module if none of the models available in the python module `uclchem.model` suit the modelling work you wish to do. If you are contributing, please try to work with our current code style. We have the following general guidelines:

### Python
- Use [Black](https://github.com/psf/black) to format your code.
- snake_case variables and functions with self-explanatory names
- Docstrings for all functions, they're used to produce the online docs!

### Fortran
- camelCase variable and subroutines names that are self-explanatory where possible 
- CAPITALIZED fortran built in functions to make code structure apparent.
- Modularization, related subroutines should be added as modules. Small tweaks should be inserted into relevant module

