"""UCLCHEM MakeRates Module.

Chemical network builder for UCLCHEM.

MakeRates is UCLCHEM's system for building custom chemical networks from
reaction databases (UMIST, KIDA) and user-defined reactions. It generates
all files needed to compile UCLCHEM with your network.

**Key Components:**

- :class:`MakeratesConfig` - Configure network generation settings
- :class:`Network` - Chemical network container
- :class:`Reaction` - Individual reaction object
- :class:`Species` - Chemical species object
- :func:`build_network` - Build network from reaction databases
- :func:`run_makerates` - Execute the full MakeRates workflow

**Quick Example:**

    >>> from uclchem.makerates.config import MakeratesConfig
    >>> from uclchem.makerates import run_makerates
    >>>
    >>> from uclchem.utils import UCLCHEM_ROOT_DIR
    >>> makerates_dir = UCLCHEM_ROOT_DIR / "../../Makerates/data"
    >>>
    >>> # Configure network builder
    >>> config = MakeratesConfig(
    ...     species_file = makerates_dir / "default/default_species.csv",
    ...     database_reaction_file = makerates_dir / "databases/umist22.csv",
    ...     database_reaction_type = "UMIST",
    ...     custom_reaction_file = makerates_dir/"default/default_grain_network.csv",
    ...     custom_reaction_type = "UCL",
    ...     output_directory = "./network_files/",
    ... )
    >>>
    >>> # Generate network, write all necessary Fortran files
    >>> network = run_makerates(config)
    >>> print(f"Reactions: {len(network.get_reaction_list())}")
    Reactions: ...

**Workflow:**

1. **Configure**: Create :class:`MakeratesConfig` with input files
2. **Build**: Use :func:`build_network` or :func:`run_makerates`
3. **Validate**: Check for missing species, duplicate reactions
4. **Export**: Generate Fortran source files for compilation

**Network Databases:**

MakeRates can read reactions from:

- **UMIST Database**: Gas-phase astrochemical reactions
- **KIDA Database**: Kinetics Database for Astrochemistry
- **Custom CSV files**: User-defined reactions

**Reaction Format:**

Reactions are typically defined in CSV with columns:

- Reactants (R1, R2)
- Products (P1, P2, P3, P4)
- Reaction type code
- Rate coefficients (alpha, beta, gamma)

**Species Format:**

Species defined in CSV with:

- Name (chemical formula)
- Mass (amu)
- Binding energy (K)
- Enthalpy of formation (K)

**Advanced Features:**

- Three-phase chemistry (gas, surface, bulk ice)
- Surface reaction types (Langmuir-Hinshelwood, Eley-Rideal)
- Thermal and non-thermal desorption
- Freeze-out and grain surface reactions
- Network optimization and validation

**Example - Adding Custom Reactions:**

    >>> from uclchem.makerates.network import Network, load_network_from_csv
    >>> from uclchem.makerates.reaction import Reaction
    >>> from uclchem.utils import UCLCHEM_ROOT_DIR
    >>>
    >>> # Load existing network
    >>> network = load_network_from_csv(
    ...     UCLCHEM_ROOT_DIR / "species.csv",
    ...     UCLCHEM_ROOT_DIR / "reactions.csv",
    ... )
    >>>

    >>> # Add custom reaction
    >>> alpha = 1.0e-10
    >>> beta = 0.0
    >>> gamma = 0.0
    >>> templow = 0.0
    >>> temphigh = 10000.0
    >>> custom_reaction = Reaction(
    ...     ["H", "CO", "NAN", "HCO", "NAN", "NAN", "NAN", alpha, beta, gamma, templow, temphigh],
    ... )
    >>> network.add_reactions(custom_reaction)
    >>>
    >>> # Export modified network
    >>> io_functions.write_reactions("custom_reactions.csv", network.get_reaction_list())

**Configuration Files:**

MakeRates uses YAML configuration:

.. code-block:: yaml

    # user_settings.yaml
    reaction_files:
      - grain_reactions.csv
      - umist12_reactions.csv
    species_file: species.csv
    output_dir: ./output
    freeze_out_reactions: True
    three_phase: True

**See Also:**

- User guide for complete MakeRates documentation
- GitHub repository for example networks
- UMIST/KIDA database documentation

**Note:**

Each UCLCHEM installation compiles one network. To change networks,
rebuild UCLCHEM with different MakeRates output.
"""

from .config import MakeratesConfig as MakeratesConfig
from .makerates import run_makerates as run_makerates
from .network import BaseNetwork as BaseNetwork
from .network import MutableNetworkABC as MutableNetworkABC
from .network import Network as Network
from .network import NetworkABC as NetworkABC
from .network import build_network as build_network
from .network import create_network as create_network
from .network import load_network_from_csv as load_network_from_csv
from .reaction import Reaction as Reaction
from .species import Species as Species
