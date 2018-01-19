#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["#CH4","#CH3","#H2O","#OH"]

#call read_uclchem. 
time,dens,temp,abundances=read_uclchem("output/full.dat",speciesNames)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances)

i=0
while time[i]<1.0e6:
	i+=1
axis.axhline(abundances[0][i]*0.01,ls="--")
axis.axvline(time[i])

while time[i]<3.0e6:
	i+=1
axis.axhline(abundances[0][i]*0.03,ls="--")
axis.axvline(time[i])

#since fig and axis were returned, optionally alter with plots.
axis.set_xlim(0,3e6)
axis.set_ylim(1e-12,1e-3)
axis.set_title("This is a test plot")
fig.savefig("output/abstractiontest.png")