"""UCLCHEM MakeRates Module

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

.. code-block:: python

    from uclchem.makerates import MakeratesConfig, run_makerates

    # Configure network builder
    config = MakeratesConfig()
    config.reaction_files = [\"umist_reactions.csv\"]
    config.species_file = \"species.csv\"
    config.output_dir = \"./network_files\"

    # Generate network
    run_makerates(config)

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

.. code-block:: python

    from uclchem.makerates import Network, Reaction

    # Load existing network
    network = load_network_from_csv(\"reactions.csv\", \"species.csv\")

    # Add custom reaction
    custom_rxn = Reaction(
        reactants=[\"H\", \"CO\"],
        products=[\"HCO\"],
        alpha=1.0e-10,
        beta=0.0,
        gamma=0.0,
        reaction_type=1
    )
    network.add_reaction(custom_rxn)

    # Export modified network
    network.to_csv(\"custom_network.csv\")

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
from .network import (
    BaseNetwork as BaseNetwork,
)
from .network import (
    MutableNetworkABC as MutableNetworkABC,
)
from .network import (
    Network as Network,
)
from .network import (
    NetworkABC as NetworkABC,
)
from .network import (
    build_network as build_network,
)
from .network import (
    create_network as create_network,
)
from .network import (
    load_network_from_csv as load_network_from_csv,
)
from .reaction import Reaction as Reaction
from .species import Species as Species
