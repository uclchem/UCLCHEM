#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
#speciesNames=["CO","CS","H2O","#CO","#CH3OH","NH3"]
speciesNames=["N2H+","NH3","CO","HCO+"]

#call read_uclchem. 
time,dens,temp,abundances=read_uclchem("output/fullcloud.dat",speciesNames)

#write out to columnated output,
write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances,plotFile="output/test.png")

#plot species returns the axis so we can further edit
axis.set(xlim=(0.1,1e7),ylim=(1e-13,1e-3))
axis.set_title('This is a Test Plot')
fig.savefig("output/test.png")
