#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["CO","E-","HCO+"]
#speciesNames=["#H2O","#CO","CO","SIO","#SIH4","HCO+"]
input_file="output/TESTbJonnetwork.dat"
plot_file="output/testsv_2.png"
#call read_uclchem. 
time,dens,temp,abundances=read_uclchem(input_file,speciesNames)

#write out to columnated output,
#write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances,plotFile=plot_file)
ax2=axis.twinx()
ax2.plot(time,temp)
#plot species returns the axis so we can further edit
axis.set(xscale='log',xlim=(1,1e7),ylim=(5e-16,3e-3))
axis.set_title('This is a Test Plot')
fig.savefig(plot_file)
