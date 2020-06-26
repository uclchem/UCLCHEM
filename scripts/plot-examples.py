#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

fig,axes=plt.subplots(2,3,figsize=(16,9))
axes=axes.flatten()
i=0
for folder in ["example-output/","test-output/"]:
	for model in ["phase1","phase2","static"]:
		#pick species, any number is fine
		speciesNames=["CO","CS","H2O","#CO","#CH3OH","NH3"]

		#call read_uclchem. 
		time,dens,temp,abundances=read_uclchem("examples/"+folder+model+"-full.dat",speciesNames)


		#plot species and save to test.png, alternatively send dens instead of time.
		axis=plot_species(speciesNames,time,abundances,ax=axes[i])

		#plot species returns the axis so we can further edit
		axis.set(xscale='log',ylim=(1e-15,1e-3),xlim=(1e0,6e6))
		axis.set_title(model.capitalize())
		i=i+1
axes[0].text(.02,0.98,"Example Row",horizontalalignment="left",verticalalignment="top",transform=axes[0].transAxes)
axes[3].text(.02,0.98,"Your Row",horizontalalignment="left",verticalalignment="top",transform=axes[3].transAxes)
fig.savefig("examples/example-comparisons.png",dpi=300)

