#! /usr/bin/python
#
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
except ModuleNotFoundError as err:
    raise ModuleNotFoundError(
        "The uclchem module could not be found, please make sure it is installed\nPlease refer to uclchem.github.io for installation instructions."
    ) from err
import logging
from argparse import ArgumentParser
import pathlib


def get_args():
    """Allows for interacting with MakeRates.py via the command line i.e.
    ```
    python3 MakeRates.py /home/user/some_directory/custom_settings.yaml --verbosity DEBUG
    ```
    Returns:
        Namespace: Arguments passed via the CLI or their defaults
    """
    parser = ArgumentParser()
    parser.add_argument(
        "settings_path", nargs="?", default="user_settings.yaml", type=pathlib.Path
    )
    parser.add_argument("-v", "--verbosity_stdout", default="WARNING", type=str)
    parser.add_argument("-d", "--debug", action="store_true")
    return parser.parse_args()


def get_logger(verbosity_stdout: str, debug: bool):
    """Define a logger that logs both to file and stdout"""
    # TODO: fix that both verbosity for file and stdout are the same type, but it works for now.
    if debug:
        verbosity_file = logging.DEBUG
        verbosity_stdout = "DEBUG"
    else:
        verbosity_file = logging.INFO
        verbosity_stdout = verbosity_stdout
    # Make sure the verbosity to the file is always smaller than stdout to avoid confusion
    if verbosity_stdout.upper() == "DEBUG":
        verbosity_file = logging.DEBUG
    logging.basicConfig(
        level=verbosity_file,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%m-%d %H:%M",
        filename="makerates.log",
        filemode="w",
    )
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(verbosity_stdout)
    # set a format which is simpler for console use
    formatter = logging.Formatter(" %(levelname)-8s %(message)s")
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger("").addHandler(console)

    # Now, we can log to the root logger, or any other logger. First the root...
    logging.info(
        f"Configured the logging. Files verbosity is {logging.getLevelName(verbosity_file)} and stdout verbosity is {verbosity_stdout}"
    )


if __name__ == "__main__":
    args = get_args()
    get_logger(args.verbosity_stdout, args.debug)
    run_makerates(args.settings_path)
