# UCLCHEM
`UCLCHEM` is a gas-grain chemical code that propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. We provide several physical models to enable the modeling of different astrophysical environments and a utility script `MakeRates` to help the user produce a chemical network from simple lists of reactions and species.


> [!TIP]
>  THe UCLCHEM 4.0 preprint is available [online](https://doi.org/10.48550/arXiv.2606.20265)

## License
We use the [MIT License](https://github.com/uclchem/UCLCHEM/blob/main/LICENSE.txt), allowing the user liberties with how they use UCLCHEM. To foster our efforts and open science, we however kindly request the user to:
1. Cite the paper if you use our code.
2. Make public any flavours of UCLCHEM you derive and use for articles, in order of preference: as a pull request on this reposistory, as a fork, code uploaded somewhere.

## Installation Instructions

Full documentation is available from the website: [uclchem.github.io](https://uclchem.github.io)

UCLCHEM is intended to be used as a python module but must be installed from source rather than an online index such as Pypi. This is because users are expected to modify the source code, at least by creating their own networks. To obtain and install the code simply run:

```bash
git clone https://github.com/uclchem/UCLCHEM.git
cd UCLCHEM
pip install .
```

You can then `import uclchem` in any python script. You need to `pip install .` whenever you change your network with Makerates. 

To see the contents of this python module, check our [Python API docs](https://uclchem.github.io/develop/api/index.html). To see some example notebooks, check the tutorial section of the docs or the notebooks in `Tutorials/`.

### Prerequisites
To build UCLCHEM, you'll need gfortran, make and python 3.12+. On MacOS, make sure to have xcode installed. See more detailed installation intructions [here](https://uclchem.github.io/develop/getting-started/installation.html).



## Changes and releases
You can check a broad description of changes in each of the [releases](https://github.com/uclchem/UCLCHEM/releases).


## Contributing
This is an open source science code for the community and are open to pull requests. We are also happy to work with you to produce a physics module if none of the models available in the python module `uclchem.model` suit the modeling work you wish to do. If you are contributing, please try to work with our current code style. Feel free to checkout the latest developments with `git fetch; git checkout develop` We have the following general guidelines:

### Editable development Setup
After cloning the repository, install the development dependencies and set up pre-commit hooks:

```bash
pip install meson-python meson ninja # dependencies 
pip install -e ".[dev]" --no-build-isolation # uclchem
pre-commit install # linting, type checking and formatting
```
This last command will automatically run linting and formatting checks before each commit and save you time with failed CI/CD.
This might introduce some weird behaviour with fortran .mod files lingering if you quit compilation at any point, so beware!


### Github
- Work in a personal branch or fork to your own Github to develop features.
- Make sure you base your new work on the develop branch.
- Pull requests should be opened with the `develop` branch as target.
- In principle, squash and merge is preferred over keeping the entire git commit history when merging into develop.
- Releases will be coordinated by the core developers, this will constitute a push to main and creating a release.

### Python
- Use [Ruff](https://docs.astral.sh/ruff/) to format your code.
- snake_case variables and functions with self-explanatory names
- Docstrings for all functions (verified using [pymend](https://github.com/JanEricNitschke/pymend), they're used to produce the online docs!

### Fortran
- camelCase variable and subroutines names that are self-explanatory where possible 
- CAPITALIZED fortran built in functions to make code structure apparent.
- Modularization, related subroutines should be added as modules. Small tweaks should be inserted into relevant module

## Citing UCLCHEM

If you use UCLCHEM for your research, please cite the following [paper](https://ui.adsabs.harvard.edu/abs/2026arXiv260620265V/abstract) for the framework
and any important papers mentioned therin that contributed to specific modules you used. 

# Developers
## Developer tools

After adding or renaming Fortran `PARAMETER` declarations, regenerate `src/uclchem/advanced/fortran_metadata.yaml` by running `uclchem-generate-metadata` (or `--dry-run` to preview changes). The CI workflow `check-fortran-metadata.yml` will fail on pull requests if the file is out of date.
