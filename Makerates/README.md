# MakeRates

MakeRates is a python script designed to turn simple lists of reactions and species into all the files UCLCHEM needs to run.

## Basic use
Whenever you want to change your network, you simply do the following from the main UCLCHEM directory.

```
cd Makerates
python MakeRates.py
cd ..
pip install .
```

Note the `pip install .`, you must always recompile UCLCHEM after changing the network. 

## User Options
Running the above will just build the default network. In order to tell MakeRates which input files to use, you should edit the file `Makerates/user_settings.yaml`. This has the following contents:

```yaml
#Your list of all species
species_file : inputFiles/default_species.csv

#core reactions from gas phase database
database_reaction_file : inputFiles/umist12-ucledit.csv
database_reaction_type : UMIST

#set of additional reactions: eg grain network
custom_reaction_file : inputFiles/default_grain_network.csv
custom_reaction_type : UCL

#whether to automatically expand to three phase network
three_phase : True

# Directory in which to store output. If not included we automatically move
# the files to the uclchem src folder. This is the default behaviour.
output_directory : "outputFiles"
```

Our [network documentation](https://uclchem.github.io/docs/network) explains what each of these parameters does as well the required format of each of the files. 

Note the output_directory is not specified by default. In this case, the output files are copied directly to `src/` and `src/fortran_src` so that the new network is incorporated on the next `pip install .`.