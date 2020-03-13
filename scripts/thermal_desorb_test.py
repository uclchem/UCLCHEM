#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["CO","#CO","H2O","#H2O","H2"]
plot_file="output/therm_desorb_tests.png"
titles=["No thermal desorption","thermal desorption","No thermal @ 20K","thermal @ 20K"]
fig,axes=plt.subplots(2,2,figsize=(16,9))
axes=axes.flatten()
for i,input_file in enumerate(["output/therm_desorb_on.dat","output/therm_desorb_off.dat","output/therm_desorb_high_off.dat"
	,"output/therm_desorb_high.dat"]):
	#call read_uclchem. 
	time,dens,temp,abundances=read_uclchem(input_file,speciesNames)
	#write out to columnated output,
	#write_cols("output/democolumns.dat",time,dens,abundances)
	colours=make_colours(len(speciesNames))

	#plot species and save to test.png, alternatively send dens instead of time.
	for specIndx,specName in enumerate(speciesNames):
	    axes[i].plot(time,abundances[specIndx],color=next(colours),label=specName)

	axes[i].legend(loc=4,fontsize='small')

	axes[i].set_xlabel('Time / years')
	axes[i].set_ylabel("X$_{Species}$")

	axes[i].text(0.05,0.9,titles[i],transform=axes[i].transAxes)
	#plot species returns the axis so we can further edit
	axes[i].set(xscale='log',yscale='log',xlim=(1,1e7),ylim=(5e-16,1e-0))
fig.savefig(plot_file,bbox_inches="tight",dpi=300)
