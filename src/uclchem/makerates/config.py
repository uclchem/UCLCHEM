"""
Pydantic-based configuration system for UCLCHEM Makerates.

This module provides validated configuration handling with clear defaults,
type checking, and automatic documentation generation.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


class MakeratesConfig(BaseModel):
    """
    Configuration for UCLCHEM Makerates chemical network generation.

    This class validates all configuration parameters and provides sensible
    defaults where appropriate. All file paths are resolved relative to the
    configuration file location.

    Example:
        >>> config = MakeratesConfig.from_yaml("user_settings.yaml")
        >>> print(config.species_file)
    """

    # ============================================================================
    # REQUIRED PARAMETERS
    # ============================================================================

    species_file: Path = Field(
        ...,
        description="Path to the species list file (CSV format with species names and properties)",
    )

    database_reaction_file: Union[Path, List[Path]] = Field(
        ...,
        description=(
            "Path(s) to database reaction files. Can be a single path or list of paths. "
            "Common databases: UMIST12, KIDA2014, etc."
        ),
    )

    database_reaction_type: Union[str, List[str]] = Field(
        ...,
        description=(
            "Type(s) of reaction database corresponding to database_reaction_file. "
            "Supported types: 'UMIST12', 'KIDA', 'UCL'. "
            "Must have same length as database_reaction_file if both are lists."
        ),
    )

    # ============================================================================
    # OPTIONAL PARAMETERS - Network Generation
    # ============================================================================

    custom_reaction_file: Optional[Union[Path, List[Path]]] = Field(
        default=None,
        description=(
            "Path(s) to custom reaction files to add to the network. "
            "Use this to supplement database reactions with custom chemistry."
        ),
    )

    custom_reaction_type: Optional[Union[str, List[str]]] = Field(
        default=None,
        description=(
            "Type(s) of custom reaction files. Must be provided if custom_reaction_file is set. "
            "Must have same length as custom_reaction_file if both are lists."
        ),
    )

    # ============================================================================
    # OPTIONAL PARAMETERS - Chemistry Features
    # ============================================================================

    add_crp_photo_to_grain: bool = Field(
        default=False,
        description=(
            "Add cosmic ray proton (CRP) and photon-induced reactions to grain surface species. "
            "Enables cosmic ray and UV photon chemistry on grain mantles."
        ),
    )

    gas_phase_extrapolation: bool = Field(
        default=False,
        description=(
            "Extrapolate gas-phase reactions to grain surface analogs. "
            "Creates surface versions of gas-phase reactions where applicable."
        ),
    )

    # ============================================================================
    # OPTIONAL PARAMETERS - Exothermicity & Heating
    # ============================================================================

    derive_reaction_exothermicity: Union[bool, str, List[str]] = Field(
        default=False,
        description=(
            "Calculate reaction exothermicity from species binding energies and formation enthalpies. "
            "Can be False (disabled), True (all reactions), or a string/list specifying reaction types e.g.: "
            "'GAS' (gas-phase only), 'TWOBODY' (Only two-body reactions), or ['CRP', 'TWOBODY'] (both cosmic rays and two-body)"
            "Automatically enables enable_rates_storage=True."
        ),
    )

    database_reaction_exothermicity: Optional[Union[Path, List[Path]]] = Field(
        default=None,
        description=(
            "Path(s) to files containing pre-calculated reaction exothermicity data. "
            "Can be a single file or list of files. CSV format with reaction indices and enthalpies. "
            "Automatically enables enable_rates_storage=True."
        ),
    )

    # ============================================================================
    # OPTIONAL PARAMETERS - Grain Surface Chemistry
    # ============================================================================

    grain_assisted_recombination_file: Optional[Path] = Field(
        default=None,
        description=(
            "Path to grain-assisted recombination (GAR) parameters file. "
            "Required if your network contains GAR reactions. "
            "Contains ion-specific recombination parameters."
        ),
    )

    # ============================================================================
    # OPTIONAL PARAMETERS - Output & Performance
    # ============================================================================

    enable_rates_storage: bool = Field(
        default=False,
        description=(
            "Store individual reaction rates during integration for post-processing analysis. "
            "Increases memory usage but enables detailed reaction pathway analysis. "
            "Automatically enabled when using exothermicity features."
        ),
    )

    output_directory: Optional[Path] = Field(
        default=None,
        description=(
            "Output directory for generated network files (network.f90, odes.f90, species.csv, etc.). "
            "If not specified, writes to src/fortran_src/ (default build location)."
        ),
    )
    # ==========================================================================
    # OPTIONAL PARAMETERS - Cooling / Coolants
    # ============================================================================

    coolants: Optional[List[Dict]] = Field(
        default=None,
        description=(
            "Optional inline list of coolant specifications. "
            "Each entry should be a dict with 'file' (filename) and 'name' (species label) keys. "
            "Optional keys: 'parent_species' (network species name for abundance lookup), "
            "'conversion_factor' (float, fraction of parent abundance to use). "
            "Example: [{'file': 'co.dat', 'name': 'CO'}, {'file': 'p-nh3.dat', 'name': 'p-NH3', "
            "'parent_species': 'NH3', 'conversion_factor': 0.5}]. "
            "If not specified, defaults to the 7 standard UCLCHEM coolants. "
            "Mutually exclusive with coolants_file."
        ),
    )

    coolants_file: Optional[Path] = Field(
        default=None,
        description=(
            "Optional path to a YAML file listing coolant specifications. "
            "Each entry should map 'file' -> filename and 'name' -> species label. "
            "If provided, this file will be loaded and used to define coolant constants. "
            "Mutually exclusive with coolants."
        ),
    )

    coolant_data_dir: Optional[str] = Field(
        default="",
        description=(
            "Optional directory path for collisional rate data files. "
            "If specified, this path will be written to f2py_constants.f90 as the default. "
            "Can be overridden at runtime via HeatingSettings.set_coolant_directory(). "
            "Default: empty string (path set by Python at runtime)."
        ),
    )

    # ============================================================================
    # DEPRECATED PARAMETERS
    # ============================================================================

    three_phase: bool = Field(
        default=True,
        description=(
            "DEPRECATED: Three-phase chemistry (gas/surface/bulk) is always enabled as of v3.5.0. "
            "Setting this to False will raise an error."
        ),
    )

    # ============================================================================
    # Internal fields (not set by user)
    # ============================================================================

    _config_dir: Optional[Path] = None  # Set during from_yaml, used for path resolution

    model_config = {
        "extra": "forbid",  # Catch typos in config files
        "validate_assignment": True,  # Validate on attribute changes
        "arbitrary_types_allowed": True,  # Allow Path objects
    }

    # ============================================================================
    # VALIDATORS
    # ============================================================================

    @field_validator(
        "database_reaction_file",
        "custom_reaction_file",
        "database_reaction_exothermicity",
        mode="before",
    )
    @classmethod
    def normalize_to_list(cls, v):
        """Convert single values to lists for consistent handling."""
        if v is None:
            return v
        if isinstance(v, (str, Path)):
            return [v]
        return v

    @field_validator("coolants_file", mode="before")
    @classmethod
    def normalize_coolants_file(cls, v):
        """Normalize a single coolant file path to a Path object."""
        if v is None:
            return v
        from pathlib import Path as _Path

        if isinstance(v, (str, _Path)):
            return _Path(v)
        raise ValueError("coolants_file must be a path to a YAML file listing coolants")

    @field_validator("coolants", mode="before")
    @classmethod
    def validate_coolants(cls, v):
        """Validate inline coolants format."""
        if v is None:
            return v
        if not isinstance(v, list):
            raise ValueError("coolants must be a list of dicts")

        from pathlib import Path as _Path

        validated = []
        for i, item in enumerate(v):
            if not isinstance(item, dict):
                raise ValueError(
                    f"coolants[{i}] must be a dict with 'file' and 'name' keys"
                )
            if "file" not in item or "name" not in item:
                raise ValueError(
                    f"coolants[{i}] must contain 'file' and 'name' keys. Got: {list(item.keys())}"
                )
            # Validate that file is a bare filename (no path)
            file_val = str(item["file"])
            if _Path(file_val).name != file_val or _Path(file_val).parent != _Path("."):
                raise ValueError(
                    f"coolants[{i}]['file'] must be a bare filename (no directories). Got: {file_val}"
                )
            entry = {"file": file_val, "name": str(item["name"])}
            if "parent_species" in item:
                entry["parent_species"] = str(item["parent_species"])
            if "conversion_factor" in item:
                entry["conversion_factor"] = float(item["conversion_factor"])
            validated.append(entry)
        return validated

    @field_validator("database_reaction_type", "custom_reaction_type", mode="before")
    @classmethod
    def normalize_type_to_list(cls, v):
        """Convert single type strings to lists for consistent handling."""
        if v is None:
            return v
        if isinstance(v, str):
            return [v]
        return v

    @model_validator(mode="after")
    def validate_reaction_files_and_types(self):
        """Ensure reaction files and types are consistent."""
        # Check database files and types match
        db_files = self.database_reaction_file
        db_types = self.database_reaction_type
        if isinstance(db_files, list) and isinstance(db_types, list):
            if len(db_files) != len(db_types):
                raise ValueError(
                    f"database_reaction_file has {len(db_files)} entries but "
                    f"database_reaction_type has {len(db_types)} entries. They must match."
                )

        # Check custom files and types match
        if self.custom_reaction_file is not None:
            if self.custom_reaction_type is None:
                raise ValueError(
                    "custom_reaction_type must be provided when custom_reaction_file is specified"
                )
            if isinstance(self.custom_reaction_file, list) and isinstance(
                self.custom_reaction_type, list
            ):
                if len(self.custom_reaction_file) != len(self.custom_reaction_type):
                    raise ValueError(
                        f"custom_reaction_file has {len(self.custom_reaction_file)} entries but "
                        f"custom_reaction_type has {len(self.custom_reaction_type)} entries. They must match."
                    )

        return self

    @model_validator(mode="after")
    def check_three_phase_deprecation(self):
        """Raise error if three_phase is explicitly set to False."""
        if self.three_phase is False:
            raise ValueError(
                "three_phase=False is deprecated as of UCLCHEM v3.5.0. "
                "Three-phase chemistry is now always enabled. "
                "Please remove 'three_phase: false' from your configuration."
            )
        return self

    @model_validator(mode="after")
    def auto_enable_rates_storage(self):
        """Automatically enable rates storage if needed for exothermicity."""
        if self.database_reaction_exothermicity or self.derive_reaction_exothermicity:
            if not self.enable_rates_storage:
                logging.warning(
                    "Exothermicity features require rate storage. "
                    "Automatically enabling enable_rates_storage=True. "
                    "Add 'enable_rates_storage: true' to your config to suppress this warning."
                )
                self.enable_rates_storage = True
        return self

    @model_validator(mode="after")
    def validate_coolants_mutual_exclusion(self):
        """Ensure coolants and coolants_file are mutually exclusive."""
        if self.coolants is not None and self.coolants_file is not None:
            raise ValueError(
                "Cannot specify both 'coolants' and 'coolants_file'. "
                "Use 'coolants' for inline specification or 'coolants_file' to reference an external file."
            )
        return self

    # ============================================================================
    # CLASS METHODS
    # ============================================================================

    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> "MakeratesConfig":
        """
        Load and validate configuration from a YAML file.

        All relative paths in the config file are resolved relative to the
        directory containing the YAML file.

        Args:
            yaml_path: Path to the YAML configuration file

        Returns:
            Validated MakeratesConfig instance

        Raises:
            ValidationError: If configuration is invalid
            FileNotFoundError: If config file doesn't exist
        """
        yaml_path = Path(yaml_path).resolve()

        if not yaml_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {yaml_path}")

        logging.info(f"Reading configuration from: {yaml_path}")
        logging.info(f"Configuration directory: {yaml_path.parent}")

        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f)

        # Create instance and store config directory for path resolution
        config = cls(**data)
        config._config_dir = yaml_path.parent

        # Coolants are no longer supplied inline via the configuration. To use
        # custom coolants, supply a `coolants_file` pointing to a YAML file with
        # mappings of 'file' and 'name' entries. If none is provided defaults will be used.

        return config

    @classmethod
    def generate_template(
        cls, output_path: Union[str, Path] = "user_settings_template.yaml"
    ):
        """
        Generate a template configuration file with all parameters documented.

        Args:
            output_path: Where to write the template file
        """
        output_path = Path(output_path)

        template = """# ============================================================================
# UCLCHEM Makerates Configuration Template
# ============================================================================
# This file configures chemical network generation for UCLCHEM.
# Lines starting with # are comments. Uncomment and modify as needed.
#
# For full documentation, see: https://uclchem.github.io
# ============================================================================

# ============================================================================
# REQUIRED PARAMETERS
# ============================================================================

# Path to species list file (relative to this config file)
species_file: "data/species.csv"

# Path(s) to reaction database files
database_reaction_file: "data/umist12_reactions.csv"
# Or use a list for multiple databases:
# database_reaction_file:
#   - "data/umist12_reactions.csv"
#   - "data/additional_reactions.csv"

# Type(s) of reaction databases (must match database_reaction_file length)
database_reaction_type: "UMIST12"
# Or as a list:
# database_reaction_type:
#   - "UMIST12"
#   - "UCL"

# ============================================================================
# OPTIONAL PARAMETERS - Network Generation
# ============================================================================

# Custom reactions to supplement the database (default: none)
# custom_reaction_file: "data/my_custom_reactions.csv"
# custom_reaction_type: "UCL"

# ============================================================================
# OPTIONAL PARAMETERS - Chemistry Features
# ============================================================================

# Add cosmic ray and photon-induced reactions to grain species (default: false)
# add_crp_photo_to_grain: false

# Extrapolate gas-phase reactions to grain surfaces (default: false)
# gas_phase_extrapolation: false

# ============================================================================
# OPTIONAL PARAMETERS - Cooling / Coolants
# ============================================================================

# Specify coolants directly inline (recommended):
# coolants:
#   - file: "co.dat"
#     name: "CO"
#   - file: "o-h2.dat"
#     name: "o-H2"
#   - file: "p-h2.dat"
#     name: "p-H2"

# OR use a separate YAML file (mutually exclusive with inline):
# coolants_file: "data/my_coolants.yaml"

# Optional: Specify the directory containing collisional rate data files
# This path will be written to f2py_constants.f90 as the default directory.
# If not specified (or empty string), the path is set by Python at runtime.
# coolant_data_dir: "/path/to/lamda/rates/"


# ============================================================================
# OPTIONAL PARAMETERS - Exothermicity & Heating
# ============================================================================

# Calculate reaction exothermicity from species data (default: false)
# Automatically enables enable_rates_storage=true
# derive_reaction_exothermicity: false

# Path(s) to pre-calculated reaction exothermicity files (default: none)
# Automatically enables enable_rates_storage=true
# database_reaction_exothermicity: "data/exothermicity.csv"
# Or as a list:
# database_reaction_exothermicity:
#   - "data/exothermicity1.csv"
#   - "data/exothermicity2.csv"

# ============================================================================
# OPTIONAL PARAMETERS - Grain Surface Chemistry
# ============================================================================

# Grain-assisted recombination parameters (required if using GAR reactions)
# grain_assisted_recombination_file: "data/gar_parameters.csv"

# ============================================================================
# OPTIONAL PARAMETERS - Output & Performance
# ============================================================================

# Store individual reaction rates for analysis (default: false)
# Increases memory usage. Automatically enabled by exothermicity features.
# enable_rates_storage: false

# Output directory for generated files (default: src/fortran_src/)
# output_directory: "output/"

# Optional: provide a custom list of coolant files (each with 'file' and 'name')
# coolants:
#   - file: "data/coolants/12co.dat"
#     name: "CO"
#   - file: "data/coolants/16o.dat"
#     name: "O"

# ============================================================================
# NOTES
# ============================================================================
# - All file paths are relative to this configuration file's directory
# - Absolute paths are also supported
# - After editing this file, run: python MakeRates.py
# - Then reinstall: pip install .
"""

        with open(output_path, "w") as f:
            f.write(template)

        print(f"✓ Template configuration written to: {output_path}")
        print(
            f"  Edit this file and use it with: python MakeRates.py --config {output_path}"
        )

    @classmethod
    def print_help(cls):
        """Print detailed help about all configuration parameters."""
        print("=" * 80)
        print("UCLCHEM MAKERATES CONFIGURATION PARAMETERS")
        print("=" * 80)
        print()

        # Get field info from the model
        for field_name, field_info in cls.model_fields.items():
            if field_name.startswith("_"):
                continue  # Skip internal fields

            # Determine if required
            is_required = field_info.is_required()
            req_text = (
                "REQUIRED" if is_required else f"Optional (default: {field_info.default})"
            )

            # Get type info
            type_text = str(field_info.annotation).replace("typing.", "")

            print(f"{field_name}")
            print(f"  Status: {req_text}")
            print(f"  Type: {type_text}")
            print(f"  Description: {field_info.description}")
            print()

        print("=" * 80)
        print("For a template configuration file, run:")
        print(
            "  python -c 'from uclchem.makerates.config import MakeratesConfig; MakeratesConfig.generate_template()'"
        )
        print("=" * 80)

    # ============================================================================
    # INSTANCE METHODS
    # ============================================================================

    def resolve_path(self, path: Union[str, Path]) -> Path:
        """
        Resolve a path relative to the configuration file directory.

        Args:
            path: Path to resolve (can be absolute or relative)

        Returns:
            Resolved absolute Path
        """
        path = Path(path)
        if path.is_absolute():
            return path
        elif self._config_dir:
            return (self._config_dir / path).resolve()
        else:
            return path.resolve()

    def get_all_reaction_files(self) -> List[Path]:
        """
        Get all reaction files (database + custom) as resolved paths.

        Returns:
            List of absolute paths to all reaction files
        """
        files = []

        # Add database files
        if isinstance(self.database_reaction_file, list):
            files.extend([self.resolve_path(f) for f in self.database_reaction_file])
        else:
            files.append(self.resolve_path(self.database_reaction_file))

        # Add custom files if present
        if self.custom_reaction_file:
            if isinstance(self.custom_reaction_file, list):
                files.extend([self.resolve_path(f) for f in self.custom_reaction_file])
            else:
                files.append(self.resolve_path(self.custom_reaction_file))

        return files

    def get_all_reaction_types(self) -> List[str]:
        """
        Get all reaction types (database + custom) in correct order.

        Returns:
            List of reaction type strings
        """
        types = []

        # Add database types
        if isinstance(self.database_reaction_type, list):
            types.extend(self.database_reaction_type)
        else:
            types.append(self.database_reaction_type)

        # Add custom types if present
        if self.custom_reaction_type:
            if isinstance(self.custom_reaction_type, list):
                types.extend(self.custom_reaction_type)
            else:
                types.append(self.custom_reaction_type)

        return types

    def log_configuration(self):
        """Log the current configuration for debugging."""
        logging.info("Configuration loaded successfully:")
        for field_name, field_info in self.__class__.model_fields.items():
            if field_name.startswith("_"):
                continue
            value = getattr(self, field_name)
            if value != field_info.default:  # Only log non-default values
                logging.info(f"  {field_name}: {value}")
