#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["H2","CO","H2O","CH3OH","#H2O","#CO","#CH3OH"]
input_file="output/full.dat"
plot_file="output/full.png"


#call read_uclchem. 
time,dens,temp,abundances=read_uclchem(input_file,speciesNames)

#write out to columnated output,
#write_cols("output/democolumns.dat",time,dens,abundances)

#plot species and save to test.png, alternatively send dens instead of time.
fig,ax=plt.subplots(figsize(16,9)) 

#plot_species will make an axis if you don't provide one
#it will save if plotFile is given otherwise just returns figure and axis
fig,ax=plot_species(speciesNames,time,abundances,ax=ax,plotFile=plot_file)


ax2=ax.twinx()
ax2.plot(time,temp)

#plot species returns the axis so we can further edit
ax.set(xscale='log',xlim=(1,1e7),ylim=(5e-16,3e-3))
ax.set_title('This is a Test Plot')

#overwrite our previous plot
fig.savefig(plot_file)
