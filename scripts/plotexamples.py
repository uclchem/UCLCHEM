#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

for folder in ["example-output/","test-output/"]:
	for model in ["phase1","phase2","static"]:
		#pick species, any number is fine
		speciesNames=["CO","CS","H2O","#CO","#CH3OH","NH3"]

		#call read_uclchem. 
		time,dens,temp,abundances=read_uclchem("examples/"+folder+model+"-full.dat",speciesNames)

		#plot species and save to test.png, alternatively send dens instead of time.
		axis,fig=plot_species(speciesNames,time,abundances,plotFile="examples/"+folder+model+".png")

		#plot species returns the axis so we can further edit
		axis.set(xscale='log',ylim=(1e-15,1e-3),xlim=(1e0,6e6))
		axis.set_title(model)
		fig.savefig("examples/"+folder+model+".png")

