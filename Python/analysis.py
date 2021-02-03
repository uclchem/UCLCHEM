import uclchem
import numpy as np
import matplotlib.pyplot as plt
################################################
#User Inputs Go Here
################################################

speciesName="CO"
result_file="output/full.dat"
reaction_file="src/reactions.csv"
species_file="src/species.csv"


################################################

result_df=uclchem.read_output_file(result_file)
species,masses=np.loadtxt(species_file,usecols=[0,1],dtype=str,skiprows=1,unpack=True,delimiter=',',comments="%")
reactions=np.loadtxt(reaction_file,dtype=str,skiprows=1,delimiter=',',usecols=[0,1,2,3,4,5,6],comments="%")

fortran_reac_indxs=[i+1 for i,reaction in enumerate(reactions) if speciesName in reaction]
reac_indxs=[i for i,reaction in enumerate(reactions) if speciesName in reaction]


oldMostForms=[]
oldMostDestructs=[]
oldTotalChange=0.0
plotTimes=[]
destructions=[]
formations=[]
for i,row in result_df.iterrows():


    param_dict=uclchem.param_dict_from_output(row)
    rates=uclchem.get_species_rates(param_dict,row[species],fortran_reac_indxs)
    changes=uclchem.get_rates_of_change(rates,reactions[reac_indxs],species,speciesName,row)

    A=zip(changes,reac_indxs)
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
            outString="{x[0]} + {x[1]} -> {x[3]} + {x[4]}".format(x=reactions[reacIndxs[k]])
            outString+=": {0:.2f}%".format(float(changes[k]/totalProd)*100)
            print(outString)

        print("\nDestruction = {0:.2e} from:".format(totalDestruct))
        for k in range(0,j):
            outString="{x[0]} + {x[1]} -> {x[3]} + {x[4]}".format(x=reactions[reacIndxs[k]])
            outString+=": {0:.2f}%".format(float(changes[k]/totalDestruct)*100)
            print(outString)
        plotTimes.append(row["Time"])
        oldTotalChange=totalChange


fig,ax=plt.subplots()
ax.plot(result_df["Time"],result_df[speciesName],color="black")
ax.plot(result_df["Time"],destructions,color="red")
ax.plot(result_df["Time"],formations,color="green")
#for time in plotTimes:
	#ax.axvline(time)
ax.set_yscale('log')
ax.set_ylim(1e-35,)
plt.show()
