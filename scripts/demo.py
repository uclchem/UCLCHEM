#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["CO","#CO","CS","NH3"]

#call read_uclchem. 
time,dens,temp,abundances=read_uclchem("output/full.dat",speciesNames)

#write out to columnated output,
write_cols("testcols.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances,plotFile="output/test.png")

#since fig and axis were returned, optionally alter with plots.
axis.set_xscale('log')
axis.set_title("This is a test plot")
fig.savefig("output/test2.png")