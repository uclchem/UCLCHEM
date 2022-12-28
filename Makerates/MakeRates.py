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
    parser = ArgumentParser()
    parser.add_argument(
        "settings_path", nargs="?", default="user_settings.yaml", type=pathlib.Path
    )
    parser.add_argument("-v", "--verbosity", default="WARNING", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    logging.basicConfig(
        format="%(levelname)s: %(message)s",
        level=args.verbosity,
    )
    run_makerates(args.settings_path)
