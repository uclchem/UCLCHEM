#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

speciesNamesLists=[["H","H2"],["CO","H2O","CH3OH","NH3","HCO+","CO2"],["#CO","#CO2","#CH3OH","#H2O"]]

for i,speciesNames in enumerate(speciesNamesLists):
	#call read_uclchem. 
	time,dens,temp,abundances=read_uclchem("output/full.dat",speciesNames)

	#plot species and save to test.png, alternatively send dens instead of time.
	axis,fig=plot_species(speciesNames,time,abundances)


	#since fig and axis were returned, optionally alter with plots.
	axis.set_xscale('log')
	if i != 0:
		axis.set_ylim(1e-15,1e-3)
	axis.set_xlim(1e1,6e6)
	#axis.set_title("This is a test plot")
	fig.savefig("output/test_abundances{0}.png".format(i+1))