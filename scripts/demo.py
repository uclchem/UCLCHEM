#demonstration of plotfunctions. caleld from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["H2S","CO","CS","N2H+","HCO+"]

#call read_uclchem. 
time,dens,abundances=read_uclchem("output/full.dat",speciesNames)

#write out to columnated output,
write_cols("testcols.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances,plotFile="test.png")

#since fig and axis were returned, optionally alter with plots.
axis.set_xlim(0,3e5)
fig.savefig("test2.png")