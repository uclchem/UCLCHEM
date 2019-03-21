#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
#zspeciesNames=["CO","CS","H2O","#CO","#CH3OH","NH3"]
speciesNames=["#H2O","#CO","CO","SIO","#SIH4","HCO+"]
plot_file="output/test10.png"
#call read_uclchem. 
time,dens,temp,abundances=read_uclchem("output/fullsputter10.dat",speciesNames)

#write out to columnated output,
#write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
axis,fig=plot_species(speciesNames,time,abundances,plotFile=plot_file)
ax2=axis.twinx()
ax2.plot(time,temp)
#plot species returns the axis so we can further edit
axis.set(xlim=(0,50),ylim=(1e-20,1e-3))
axis.set_title('This is a Test Plot')
fig.savefig(plot_file)
