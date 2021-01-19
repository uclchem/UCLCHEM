####################################################################
# Example Plotting Script
# This is a simple script that demonstrates how to use the functions
# in plotfunctions.py to make plots from UCLCHEM outputs in 3 ways
###################################################################

from plotfunctions import *

#pick species, any number is fine
speciesNames=["#H2O","#CO","#CH3OH","@CO","@CH3OH","@H2O"]
input_file="output/full2.dat"


#call read_uclchem. 
model_data=read_uclchem(input_file)

#If we simply send a list of species and the data to plot_species
#we'll get an axis object back and can also save the file as is
ax=plot_species(speciesNames,model_data,plot_file="examples/example_plot.png")


#the returned object lets us make some edits
#lets plot the temperature on a second y axis 
ax2=ax.twinx()
ax2.plot(model_data["Time"],model_data["gasTemp"],color="black")
ax2.set(ylabel="Temperature / K")
#and make some slight adjustments to the plot before saving again
ax.set(xscale='log',xlim=(1,1e7),ylim=(9e-31,5e-4))
ax.set_title('This is a Test Plot')

#overwrite our previous plot
fig=plt.gcf()
fig.savefig("examples/improved_example_plot.png")


#finally, we can send axes to plot_species to completely control the plots
#lets plot four pairs of gas/grain species on four separate subplots
fig,axes=plt.subplots(2,2,figsize=(16,9))
axes=axes.flatten() #turn [2,2] array into [4]

for i,species in enumerate([["CO","#CO"],["H2O","#H2O"],["CH3OH","#CH3OH"],["CO2","#CO2"]]):
	axes[i]=plot_species(species,model_data,ax=axes[i])
fig.savefig("examples/multiplot_example.png")