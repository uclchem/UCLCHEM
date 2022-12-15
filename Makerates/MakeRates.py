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
    raise ModuleNotFoundError("The uclchem module could not be found, please make sure it is installed\nPlease refer to uclchem.github.io for installation instructions.") from err
    
if __name__=="__main__":
    run_makerates("user_settings.yaml", verbosity="INFO")
    