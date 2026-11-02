# UCLCHEM

[![arXiv](https://img.shields.io/badge/arXiv-2606.20265-gray?style=for-the-badge&logo=arxiv&logoSize=auto&labelColor=gray&color=%23B31B1B)](https://arxiv.org/abs/2606.20265)
[![GitHub Release](https://img.shields.io/github/v/release/uclchem/UCLCHEM?include_prereleases&style=for-the-badge)](https://github.com/uclchem/UCLCHEM/releases)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fuclchem%2FUCLCHEM%2Fmain%2Fpyproject.toml&style=for-the-badge&logo=python&logoSize=auto)
[![Read the Docs](https://img.shields.io/readthedocs/uclchem?style=for-the-badge&logo=readthedocs)](https://uclchem.github.io)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/uclchem/UCLCHEM/ci-lint.yml?style=for-the-badge&logo=precommit&label=pre-commit)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/uclchem/UCLCHEM/.github%2Fworkflows%2Fci-test-suite.yml?style=for-the-badge&logo=pytest&logoSize=auto)
[![MIT License](https://img.shields.io/github/license/uclchem/UCLCHEM?style=for-the-badge)](https://github.com/uclchem/UCLCHEM/LICENSE.txt)
[![ty](https://img.shields.io/badge/-ty-gray?style=for-the-badge&logo=ty&logoColor=%2346ebe1&labelColor=gray&color=%2330173d)](https://github.com/astral-sh/ty)
[![ruff](https://img.shields.io/badge/-ruff-gray?style=for-the-badge&logo=ruff&logoColor=%2346ebe1&labelColor=gray&color=%2330173d)](https://github.com/astral-sh/ruff)<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/github/all-contributors/uclchem/UCLCHEM?style=for-the-badge)](#contributors)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

<div align="center">
  <img src="https://raw.githubusercontent.com/uclchem/uclchem/main/images/uclchem.svg" width="200">
</div>

A gas-grain chemical code that propagates the abundances of chemical species through a network of user-defined reactions
according to the physical conditions of the gas. We provide several physical models to enable the modeling of different
astrophysical environments and a utility script `MakeRates` to help the user produce a chemical network from simple
lists of reactions and species.

> [!TIP]
>  The UCLCHEM 4.0 preprint is available [online](https://doi.org/10.48550/arXiv.2606.20265)

## Installation

### Prerequisites
To build UCLCHEM, you'll need gfortran, make and python 3.11+.
On MacOS, make sure to have xcode installed.
See more detailed installation instructions [here](https://uclchem.github.io/develop/getting-started/installation.html).

### Installation

UCLCHEM is intended to be used as a python module but must be installed from source rather than an online index such as PyPI.
This is because users are expected to modify the source code, at least by creating their own networks.
To obtain and install the code simply run:

```bash
git clone https://github.com/uclchem/UCLCHEM.git
cd UCLCHEM
pip install .
```

If you have [uv](https://docs.astral.sh/uv/) installed, `uv pip install .` is a drop-in, much faster replacement for the `pip install .` line above.

You can then `import uclchem` in any python script. You need to `pip install .` (or `uv pip install .`) whenever you change your network with Makerates. 

## Usage

Running a basic dark molecular cloud model:
```python
import uclchem
import matplotlib.pyplot as plt

if __name__ == "__main__":
    param_dict = {
        "initialDens": 1e5,  # Initial density of the cloud
        "initialTemp": 10,   # Temperature
        "finalTime": 1e6,    # Duration of the model in years
    }
    model = uclchem.model.Cloud(param_dict)
    model.check_error()

    # Plot the abundances of gas-phase and icy H2O and CO over lifetime of the cloud
    fig, axs = model.create_abundance_plot(species=["H2O", "CO", "$H2O", "$CO"])
    plt.show()
```

Full documentation is available from the website: [uclchem.github.io](https://uclchem.github.io)

To see the contents of this python module, check our [Python API docs](https://uclchem.github.io/develop/api/index.html).
To see some example notebooks, check the tutorial section of the docs or the notebooks in `Tutorials/`.

## Citing UCLCHEM

If you use UCLCHEM for your research, please cite the following [paper](https://ui.adsabs.harvard.edu/abs/2026arXiv260620265V/abstract)
for the framework and any important papers mentioned therein that contributed to specific modules you used. 


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

This is an open source science code for the community and are open to pull requests.
We are also happy to work with you to produce a physics module if none of the models available in the python module `uclchem.model`
suit the modeling work you wish to do. If you are contributing, please try to work with our current code style.
Feel free to checkout the latest developments with `git fetch; git checkout develop` We have the following general guidelines:

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

After adding or renaming Fortran `PARAMETER` declarations, regenerate `src/uclchem/advanced/fortran_metadata.yaml` by running `uclchem-generate-metadata` (or `--dry-run` to preview changes). The CI workflow `check-fortran-metadata.yml` will fail on pull requests if the file is out of date.

## License
We use the [MIT License](https://github.com/uclchem/UCLCHEM/blob/main/LICENSE.txt), allowing the user liberties with how they use UCLCHEM. To foster our efforts and open science, we however kindly request the user to:
1. Cite the paper if you use our code.
2. Make public any flavors of UCLCHEM you derive and use for articles, in order of preference: as a pull request on this repository, as a fork, code uploaded somewhere.


## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/en/reference/emoji-key/)):

<!-- spellchecker:off -->
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/jonholdship"><img src="https://avatars.githubusercontent.com/u/11349588?v=4&size=100" width="100px" height="100px" alt="Jon Holdship"/><br /><sub><b>Jon Holdship</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=jonholdship" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/GijsVermarien"><img src="https://avatars.githubusercontent.com/u/25642414?v=4&size=100" width="100px" height="100px" alt="Gijs Vermariën"/><br /><sub><b>Gijs Vermariën</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=GijsVermarien" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/TobiasDijkhuis"><img src="https://avatars.githubusercontent.com/u/99098052?v=4&size=100" width="100px" height="100px" alt="Tobias Dijkhuis"/><br /><sub><b>Tobias Dijkhuis</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=TobiasDijkhuis" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Marcus-Keil"><img src="https://avatars.githubusercontent.com/u/37833544?v=4&size=100" width="100px" height="100px" alt="Dr. Marcus Keil"/><br /><sub><b>Dr. Marcus Keil</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=Marcus-Keil" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/tomasjames"><img src="https://avatars.githubusercontent.com/u/5236734?v=4&size=100" width="100px" height="100px" alt="Tomas James"/><br /><sub><b>Tomas James</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=tomasjames" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Ross-ODonoghue"><img src="https://avatars.githubusercontent.com/u/59836299?v=4&size=100" width="100px" height="100px" alt="Ross O'Donoghue"/><br /><sub><b>Ross O'Donoghue</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=Ross-ODonoghue" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Bamash"><img src="https://avatars.githubusercontent.com/u/32938368?v=4&size=100" width="100px" height="100px" alt="Bamash"/><br /><sub><b>Bamash</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=Bamash" title="Code">💻</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/lengoctram"><img src="https://avatars.githubusercontent.com/u/91097424?v=4&size=100" width="100px" height="100px" alt="Le N. Tram"/><br /><sub><b>Le N. Tram</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=lengoctram" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://piyushsharda.wixsite.com/main"><img src="https://avatars.githubusercontent.com/u/34922596?v=4&size=100" width="100px" height="100px" alt="Piyush Sharda"/><br /><sub><b>Piyush Sharda</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=psharda" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/katarzynadutkowska"><img src="https://avatars.githubusercontent.com/u/73302015?v=4&size=100" width="100px" height="100px" alt="Katarzyna M. Dutkowska"/><br /><sub><b>Katarzyna M. Dutkowska</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=katarzynadutkowska" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Diederick-V"><img src="https://avatars.githubusercontent.com/u/79215711?v=4&size=100" width="100px" height="100px" alt="Diederick Vroom "/><br /><sub><b>Diederick Vroom </b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/commits?author=Diederick-V" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/JosephSalaris"><img src="https://avatars.githubusercontent.com/u/64165609?v=4&size=100" width="100px" height="100px" alt="JosephSalaris"/><br /><sub><b>JosephSalaris</b></sub></a><br /><a href="https://github.com/uclchem/UCLCHEM/issues?q=author%3AJosephSalaris" title="Bug reports">🐛</a></td>
    </tr>
  </tbody>
  <tfoot>
    <tr>
      <td align="center" size="13px" colspan="7">
        <img src="https://raw.githubusercontent.com/all-contributors/all-contributors-cli/1b8533af435da9854653492b1327a23a4dbd0a10/assets/logo-small.svg">
          <a href="https://all-contributors.js.org/docs/en/bot/usage">Add your contributions</a>
        </img>
      </td>
    </tr>
  </tfoot>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
<!-- spellchecker:on -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
