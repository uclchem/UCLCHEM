#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

speciesNamesLists=[["H","H2"],["CO","H2O","CH3OH","NH3","HCO+","CO2"],["#CO","#CO2","#CH3OH","#H2O"]]
fig,axes=plt.subplots(2,2)
axes=axes.flatten()

for i,speciesNames in enumerate(speciesNamesLists):
	ax=axes[i]
	#call read_uclchem. 
	time,dens,temp,abundances=read_uclchem("output/full.dat",speciesNames)

	#write out to columnated output,
	write_cols("testcols.dat",time,dens,abundances)
	colours=make_colours(len(speciesNames))

	for specIndx,specName in enumerate(speciesNames):
	    ax.plot(time,abundances[specIndx],color=colours.next(),label=specName)

	ax.legend(loc=4,fontsize='small')

	ax.set_xlabel('Time / years')
	ax.set_ylabel("X$_{Species}$")

	ax.set_yscale('log')

	#since fig and axis were returned, optionally alter with plots.
	ax.set_xscale('log')
	if i !=0:
		ax.set_ylim(1e-15,1e-3)
	ax.set_xlim(1e1,6e6)
#axis.set_title("This is a test plot")
fig.savefig("output/test_abundances.png")