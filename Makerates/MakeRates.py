# ruff: ignore[undocumented-public-module]
####################################################################################################
# 				MakeRates
# 		Current version by Jonathan Holdship & Antonios Makrymallis. Original by Tom Bell.
# 		MakeRates reads in lists of species and reactions and produces the network files needed
# 		by UCLCHEM to run. It also performs basic cleaning and sanity checks on the network.
#
####################################################################################################
# All code that is run by this script resides in src/uclchem/makerates/makerates.py
try:
    from uclchem.makerates import run_makerates
    from uclchem.makerates.config import MakeratesConfig
except ModuleNotFoundError as err:
    msg = (
        "The uclchem module could not be found, please make sure it is "
        "installed.\nPlease refer to uclchem.github.io for installation "
        "instructions."
    )
    raise ModuleNotFoundError(msg) from err
import logging
import pathlib
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

_LOG_POSSIBLE_LEVELS = ("DEBUG", "INFO", "WARNING", "CRITICAL", "ERROR")


def get_args():  # ruff: ignore[missing-return-type-undocumented-public-function]
    """Get the parsed arguments.

    Allows for interacting with MakeRates.py via the command line.

    Returns
    -------
    Namespace
        Arguments passed via the CLI or their defaults

    Examples
    --------
    python3 MakeRates.py custom_settings.yaml --verbosity DEBUG
    python3 MakeRates.py --generate-template
    python3 MakeRates.py --help-config

    """
    parser = ArgumentParser(
        description="UCLCHEM Makerates: Generate chemical network files",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    # Main argument - config file path
    parser.add_argument(
        "settings_path",
        nargs="?",
        default="user_settings.yaml",
        type=pathlib.Path,
        help="Path to YAML configuration file",
        metavar="settings-path",
    )

    # Verbosity options
    parser.add_argument(
        "-v",
        "--verbosity-stdout",
        default="WARNING",
        type=str,
        help="Console output verbosity",
        choices=_LOG_POSSIBLE_LEVELS,
        metavar="LEVEL",
    )
    parser.add_argument(
        "-f",
        "--verbosity-file",
        default="INFO",
        help="Verbosity of output to 'makerates.log'",
        choices=_LOG_POSSIBLE_LEVELS,
        metavar="LEVEL",
    )

    # Helper options
    parser.add_argument(
        "--generate-template",
        action="store_true",
        help="Generate a template configuration file and exit",
    )
    parser.add_argument(
        "--help-config",
        action="store_true",
        help="Print detailed help about configuration parameters and exit",
    )

    # Options to help debugging the network
    parser.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        help="Only do a dry run, meaning network validation, without writing the files",
    )
    parser.add_argument(
        "-p",
        "--output-directory",
        type=pathlib.Path,
        required=False,
        help="Directory to write the output files. If not passed, write to the directory specified in one of the configuration files, or else the package source directory.",
    )

    return parser.parse_args()


def get_logger(verbosity_stdout: str, verbosity_file: str) -> None:
    """Define a logger that logs both to file and stdout.

    Parameters
    ----------
    verbosity_stdout : str
        stdout verbosity
    verbosity_file : str
        Verbosity to write to 'makerates.log'.

    """
    logger = logging.getLogger("uclchem")
    logger.setLevel("DEBUG")

    file_handler = logging.FileHandler(
        "makerates.log",
        mode="w",
    )
    file_handler.setLevel(verbosity_file)
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s %(levelname)s: %(message)s",
            datefmt="%m-%d %H:%M",
        )
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(verbosity_stdout)
    console_handler.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    logger.info(
        "Configured logging: file=%s, stdout=%s",
        logging.getLevelName(file_handler.level),
        logging.getLevelName(console_handler.level),
    )


if __name__ == "__main__":
    args = get_args()

    # Handle helper commands that exit immediately
    if args.help_config:
        MakeratesConfig.print_help()
        sys.exit(0)

    if args.generate_template:
        output_file = "user_settings_template.yaml"
        MakeratesConfig.generate_template(output_file)
        sys.exit(0)

    # Set up logging
    get_logger(args.verbosity_stdout, args.verbosity_file)

    # Run makerates with the specified config file
    run_makerates(
        args.settings_path,
        output_directory=args.output_directory,
        write_files=not args.dry_run,
    )
