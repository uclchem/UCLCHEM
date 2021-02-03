<a name="uclchem"></a>
# uclchem

<a name="uclchem.run_model"></a>
#### run\_model

```python
run_model(param_dict)
```

Run UCLCHEM using variables taking from a dictionary of parameter values. Any parameter
not included in the dictionary will be taken from defaultparameters.f90.

**Arguments**:

- `param_dict`: A dictionary of parameters where keys are any of the variables in
defaultparameters.f90 and values are value for current run.

<a name="uclchem.run_model_for_abundances"></a>
#### run\_model\_for\_abundances

```python
run_model_for_abundances(param_dict)
```

Run UCLCHEM, returning the abundances of up to 50 species at the end of the run.

**Arguments**:

- `param_dict`: A dictionary of parameters where keys are any of the variables in
defaultparameters.f90 and values are value for current run.

<a name="uclchem.get_species_rates"></a>
#### get\_species\_rates

```python
get_species_rates(param_dict, input_abundances, reac_indxs)
```

Get the rate of up to 500 reactions from UCLCHEM for a given set of parameters and abundances.
Intended for use within the analysis script.

**Arguments**:

- `param_dict`: A dictionary of parameters where keys are any of the variables in
defaultparameters.f90 and values are value for current run.
- `input_abundances`: Abundance of every species in network
- `reac_indxs`: Index of reactions of interest in the network's reaction list.

<a name="uclchem.read_output_file"></a>
#### read\_output\_file

```python
read_output_file(output_file)
```

Read the output of a UCLCHEM run created with the outputFile parameter into a pandas DataFrame

**Arguments**:

- `output_file`: - path to file containing a full UCLCHEM output

<a name="uclchem.create_abundance_plot"></a>
#### create\_abundance\_plot

```python
create_abundance_plot(df, species, plot_file=None)
```

Produce a plot of the abundances of chosen species through time, returning the pyplot
figure and axis objects

**Arguments**:

- `df`: A dataframe created by :func:`read_output_file`
- `species`: A list of species names to be plotted
- `plot_file`: optional argument with path to file where the plot should be saved

<a name="uclchem.plot_species"></a>
#### plot\_species

```python
plot_species(ax, df, species)
```

Plot the abundance of several species through time onto an existing pyplot axis

**Arguments**:

- `ax`: pyplot axis on which to plot
- `df`: A dataframe created by :func:`read_output_file`
- `species`: A list of species names to be plotted

<a name="uclchem.param_dict_from_output"></a>
#### param\_dict\_from\_output

```python
param_dict_from_output(output_line)
```

Generate a parameter dictionary with enough variables to correctly estimate the rates of 
reactions.

<a name="uclchem.get_rates_of_change"></a>
#### get\_rates\_of\_change

```python
get_rates_of_change(rates, reactions, speciesList, species, row)
```

Calculate the terms in the rate of equation of a particular species using rates calculated using
get_species_rates() and a row from the full output of UCLCHEM.

