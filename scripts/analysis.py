'''
analysis.py calculates the most important reactions causing the change of a species at each time in the model
this allows the user to understand the chemistry causing the abundances they get from the model.
Currently, the script only works for gas phase species as the grain surface reactions are not coded.
'''
from __future__ import print_function
from plotfunctions import *

################################################
#User Inputs Go Here
################################################

speciesName="CO"
resultFile="examples/example-output/phase1-full.dat"
reactionFile="src/reactions.csv"
speciesFile="src/species.csv"
parameterFile="src/defaultparameters.f90"


################################################################################################
#NO CHANGES REQUIRED BELOW THIS LINE
################################################################################################




species,masses=np.loadtxt(speciesFile,usecols=[0,1],dtype=str,skiprows=1,unpack=True,delimiter=',',comments="%")
species=list(species)
grains=[i for i,item in enumerate(species) if "#" in item]
network=getNetwork(reactionFile,speciesName)
cloud=getParameters(parameterFile)
model_output=read_uclchem(resultFile)

oldMostForms=[]
oldMostDestructs=[]
plotTimes=[]
oldTotalChange=0.0
destructions=[]
formations=[]
for i,row in model_output.iterrows():
	abundances=row[species].values
	cloud["density"]=row["Density"]	
	cloud["temp"]=row["gasTemp"]
	cloud['mantle']=sum(abundances[grains])
	cloud['av']=row["av"]
	changes,reacIndxs=getChanges(speciesName,species,masses,abundances,network,cloud)#speciesName,species,abundances,network,temp,dens

	A=zip(changes,reacIndxs)
	A=sorted(A)
	changes,reacIndxs=zip(*A)
	changes=np.asarray(changes)

	totalDestruct=sum(changes[np.where(changes<0)])
	destructions.append(-totalDestruct)
	totalProd=sum(changes[np.where(changes>0)])
	formations.append(totalProd)

	totalChange=sum(changes)
	mostForms=[]
	form=0.0
	i=-1
	while form < 0.99*totalProd:
		mostForms.append(reacIndxs[i])
		form+=changes[i]
		i-=1

	mostDestructs=[]	
	j=0
	destruct=0.0
	while abs(destruct) < 0.99*abs(totalDestruct):
		mostDestructs.append(reacIndxs[j])
		destruct+=changes[j]
		j+=1

	if set(oldMostDestructs)!=set(mostDestructs) or set(oldMostForms) !=set(mostForms):
		oldMostDestructs=mostDestructs[:]
		oldMostForms=mostForms[:]
		print("\n***************************\nNew Important Reactions At: {0:.2e} years\n".format(row["Time"]))
		print("Formation = {0:.2e} from:".format(totalProd))
		for k in range(-1,i,-1):
			outString="{x[0]} + {x[1]} -> {x[3]} + {x[4]}".format(x=network[reacIndxs[k]])
			outString+=": {0:.2f}%".format(float(changes[k]/totalProd)*100)
			print(outString)

		print("\nDestruction = {0:.2e} from:".format(totalDestruct))
		for k in range(0,j):
			outString="{x[0]} + {x[1]} -> {x[3]} + {x[4]}".format(x=network[reacIndxs[k]])
			outString+=": {0:.2f}%".format(float(changes[k]/totalDestruct)*100)
			print(outString)
		plotTimes.append(row["Time"])
		oldTotalChange=totalChange

fig,ax=plt.subplots()
ax.plot(model_output["Time"],model_output[speciesName],color="black")
ax.plot(model_output["Time"],destructions,color="red")
ax.plot(model_output["Time"],formations,color="green")
#for time in plotTimes:
	#ax.axvline(time)
ax.set_yscale('log')
ax.set_ylim(1e-35,)
plt.show()
