# UCLCHEM

![GitHub Release](https://img.shields.io/github/v/release/uclchem/UCLCHEM?include_prereleases&style=for-the-badge)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fuclchem%2FUCLCHEM%2Fmain%2Fpyproject.toml&style=for-the-badge&logo=python&logoSize=auto)

<img alt="Read the Docs" src="https://img.shields.io/readthedocs/uclchem?style=for-the-badge&link=!%5BRead%20the%20Docs%5D(https%3A%2F%2Fimg.shields.io%2Freadthedocs%2Fuclchem%3Fstyle%3Dfor-the-badge%26link%3Duclchem.github.io)">

![Read the Docs](https://img.shields.io/readthedocs/uclchem?style=for-the-badge&logo=readthedocs&link=https%3A%2F%2Fuclchem.github.io)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/uclchem/UCLCHEM/ci-lint.yml?style=for-the-badge&logo=precommit&label=pre-commit)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/uclchem/UCLCHEM/.github%2Fworkflows%2Fci-test-suite.yml?style=for-the-badge&logo=pytest&logoSize=auto)
![GitHub License](https://img.shields.io/github/license/uclchem/UCLCHEM?style=for-the-badge&link=https%3A%2F%2Fgithub.com%2Fuclchem%2FUCLCHEM%2FLICENSE)

A gas-grain chemical code that propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas. We provide several physical models to enable the modeling of different astrophysical environments and a utility script `MakeRates` to help the user produce a chemical network from simple lists of reactions and species.

> [!TIP]
>  The UCLCHEM 4.0 preprint is available [online](https://doi.org/10.48550/arXiv.2606.20265)

## Installation Instructions

Full documentation is available from the website: [uclchem.github.io](https://uclchem.github.io)

UCLCHEM is intended to be used as a python module but must be installed from source rather than an online index such as Pypi. This is because users are expected to modify the source code, at least by creating their own networks. To obtain and install the code simply run:

```bash
git clone https://github.com/uclchem/UCLCHEM.git
cd UCLCHEM
pip install .
```

If you have [uv](https://docs.astral.sh/uv/) installed, `uv pip install .` is a drop-in, much faster replacement for the `pip install .` line above.

You can then `import uclchem` in any python script. You need to `pip install .` (or `uv pip install .`) whenever you change your network with Makerates. 

To see the contents of this python module, check our [Python API docs](https://uclchem.github.io/develop/api/index.html). To see some example notebooks, check the tutorial section of the docs or the notebooks in `Tutorials/`.

### Prerequisites
To build UCLCHEM, you'll need gfortran, make and python 3.11+. On MacOS, make sure to have xcode installed. See more detailed installation instructions [here](https://uclchem.github.io/develop/getting-started/installation.html).

## Citing UCLCHEM

If you use UCLCHEM for your research, please cite the following [paper](https://ui.adsabs.harvard.edu/abs/2026arXiv260620265V/abstract) for the framework
and any important papers mentioned therein that contributed to specific modules you used. 


## Changes and releases
You can check a broad description of changes in each of the [releases](https://github.com/uclchem/UCLCHEM/releases).


## Contributing
**Imposter syndrome disclaimer**: We want your help. No really, we do.

There might be a little voice inside that tells you you're not ready;
that you need to do one more tutorial, or learn another framework,
or write a few more blog posts before you can help with this project.

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

(This disclaimer was originally written by
[Adrienne Lowe](https://github.com/adriennefriend) for a
[PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo),
and was adapted from the README in the [yt project](https://github.com/yt-project/yt))

This is an open source science code for the community and are open to pull requests. We are also happy to work with you to produce a physics module if none of the models available in the python module `uclchem.model` suit the modeling work you wish to do. If you are contributing, please try to work with our current code style. Feel free to checkout the latest developments with `git fetch; git checkout develop` We have the following general guidelines:

### Editable development Setup
After cloning the repository, install the development dependencies and set up pre-commit hooks.

With [uv](https://docs.astral.sh/uv/) (recommended -- much faster, and `uv.lock` pins every dependency so everyone builds against the same versions):

```bash
uv sync --group build --extra dev # creates .venv with uclchem + dev tools, builds the Fortran extension
uv run pre-commit install # linting, type checking and formatting
```
Re-run `uv sync --group build --extra dev` whenever you change your network with Makerates, or after pulling changes that touch the Fortran sources. Use `uv run <command>` (e.g. `uv run pytest`) to run anything inside `.venv` without activating it first, or `source .venv/bin/activate` as usual.

Without uv, plain pip works the same way it always has:

```bash
pip install meson-python meson ninja # dependencies 
pip install -e ".[dev]" --no-build-isolation # uclchem
pre-commit install # linting, type checking and formatting
```

Either way, the pre-commit install step will automatically run linting and formatting checks before each commit and save you time with failed CI/CD.
This might introduce some weird behavior with fortran .mod files lingering if you quit compilation at any point, so beware!

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

# Developers
## Developer tools

After adding or renaming Fortran `PARAMETER` declarations, regenerate `src/uclchem/advanced/fortran_metadata.yaml` by running `uclchem-generate-metadata` (or `--dry-run` to preview changes). The CI workflow `check-fortran-metadata.yml` will fail on pull requests if the file is out of date.

## License
We use the [MIT License](https://github.com/uclchem/UCLCHEM/blob/main/LICENSE.txt), allowing the user liberties with how they use UCLCHEM. To foster our efforts and open science, we however kindly request the user to:
1. Cite the paper if you use our code.
2. Make public any flavors of UCLCHEM you derive and use for articles, in order of preference: as a pull request on this reposistory, as a fork, code uploaded somewhere.

