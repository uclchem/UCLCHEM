#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["CO","CS","H2O","#CO","#CH3OH","NH3"]

#call read_uclchem. 
time,dens,temp,abundances=read_uclchem("output/full.dat",speciesNames)

#write out to columnated output,
write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances,plotFile="output/cloud.png")

#plot species returns the axis so we can further edit
axis.set_xscale('log')
axis.set_xlim(0.1,1e7)
axis.set_title('This is a Test Plot')
fig.savefig("output/cloud2.png")
