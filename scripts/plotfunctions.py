#A set of functions for working with UCLCHEM outputs
# adding "from plotfunctions import * to any python script in scripts/ will allow their use"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv

def read_uclchem(filename):
    f=open(filename)
    f.readline()
    bits=f.readline().split()
    radfield=bits[1]
    zeta=bits[3]
    data=pd.read_csv(f)
    data["zeta"]=zeta
    data["radfield"]=radfield
    data.columns=data.columns.str.strip()
    return data

#send a  list of species names and their abundances in  a list of lists
#with a list of times for each abundance point
#same as the species input and time/abundance output from read_uclchem
#optionally send an output filename to save the plot
#return ax,figure for further manipulation
def plot_species(species,df,ax=None,plot_file=None):
    if ax is None:
        fig,ax=plt.subplots()
    colours=make_colours(len(species))

    for specIndx,specName in enumerate(species):
        ax.plot(df["Time"],df[specName],color=next(colours),label=specName)

    ax.legend(loc=4,fontsize='small')

    ax.set_xlabel('Time / years')
    ax.set_ylabel("X$_{Species}$")

    ax.set_yscale('log')

    if plot_file is not None:
        fig.savefig(plot_file)
    return ax
    


def make_colours(n):
    return iter(cm.rainbow(np.linspace(0, 1, n)))

def formatSpecies(speciesName):
    speciesName=speciesName.upper()
    speciesName=speciesName.replace("+","$^+$")
    speciesName=speciesName.replace("2","$_2$")
    speciesName=speciesName.replace("3","$_3$")
    speciesName=speciesName.replace("4","$_4$")
    return speciesName


#########################################################################################################
#Analysis Functions
#########################################################################################################

def getRate(reactype,mass,alpha,beta,gamma,cloud):
    if reactype=='CRP':
        rate = alpha*cloud['zeta']

    elif reactype == 'PHOTON':
        #I ignored self-shielding
        #print alpha
        rate = alpha*np.exp(-gamma*cloud['av'])*cloud['radfield']/1.7

    elif reactype=='CRPHOT':
        rate=alpha*gamma*1.0/(1.0-cloud['omega'])*cloud['zeta']*(cloud['temp']/300.0)**beta
   
    elif reactype=='FREEZE':             
        if (cloud['evap'] != 0 or cloud['fr'] == 0.0):
            rate=0.0
        else:    
            if (beta==0.0 ): 
                rate=4.57e4*alpha*np.sqrt(cloud['temp']/mass)*cloud['grainArea']*cloud['fr']
            else:
                cion=1.0+16.71e-4/(cloud['radg']*cloud['temp'])
                rate=4.57e4*alpha*np.sqrt(cloud['temp']/mass)*cloud['grainArea']*cloud['fr']*cion
                    
    elif reactype=='DESOH2':
            if (cloud['desorb'] and cloud['h2desorb']):
                if (gamma < cloud['ebmaxh2'] and cloud['mantle']>1e-30):
                    rate = cloud['epsilon']*1.0e-17*np.sqrt(cloud['temp'])*cloud['h']*1.0/cloud['mantle']
                else:
                    rate=0.0
            else:
                rate = 0.0
            
    elif reactype=='DESCR':
            if (cloud['desorb'] and cloud['crdesorb']):
                if (cloud['mantle'] > 1e-30 and gamma < cloud['ebmaxcr']): 
                    rate = 4.0*np.pi*cloud['zeta']*1.64e-4*(cloud['grainArea'])*(1.0/cloud['mantle'])*cloud['phi']
                else:
                    rate=0.0
            else:
                rate = 0.0
            
    elif reactype=='DEUVCR':
        if (cloud['desorb'] and cloud['uvcr']):
            if( gamma < cloud['ebmaxuvcr'] and cloud['mantle'] > 1.0e-30):
                rate = cloud['grainArea']*cloud['uvy']*4.875e3*cloud['zeta']*(1.0/cloud['mantle'])
                rate = rate *(1.0+(cloud['radfield']/cloud['uvcreff'])*(1.0/cloud['zeta'])*np.exp(-1.8*cloud['av']))
            else:
                rate = 0.0
        else:
            rate=0.0
    else:
        rate = alpha*((cloud['temp']/300.)**beta)*np.exp(-gamma/cloud['temp'])

    return rate

def getChanges(speciesName,species,masses,abundances,network,cloud):
    reacIndxs=[]
    changes=[]
    for i in range(0,len(network)):
        reacIndxs.append(i)
        reacIndx=list(network[i]).index(speciesName)
        rate=getRate(network[i][1],float(masses[species.index(speciesName)]),network[i][-3],network[i][-2],network[i][-1],cloud)
        reactantCount=0
        reactants=network[i][0:3]
        for reactant in reactants:
            if reactant in species:
                indx=species.index(reactant)
                rate*=abundances[indx]
                reactantCount+=1
        if reactantCount>1:
            rate*=cloud['density']
        else:
            if reactants[1]=='FREEZE':
                rate*=cloud['density']
            elif reactants[1]=='DESOH2':
                rate*=cloud['density']
            elif reactants[1]=='electr':
                rate*=cloud['density']
        if reacIndx<3:
            changes.append(-rate)
        else:
            changes.append(rate)
    return changes,reacIndxs

def getNetwork(file,speciesName):
    reactions=np.loadtxt(file,dtype=str,skiprows=1,delimiter=',',usecols=[0,1,2,3,4,5,6],comments="%")
    alpha,beta,gamma=np.loadtxt(file,usecols=[7,8,9],unpack=True,skiprows=1,delimiter=',',comments="%")
    network=[]
    for i,reaction in enumerate(reactions):
        if speciesName in reaction:
            keep=list(reaction)
            keep.extend([alpha[i],beta[i],gamma[i]])
            network.append(keep)

    return network

def getParameters(file):
    cloud={'temp':10,
    'density':1e2,
    'zeta':1.0,
    'av':1.0,
    'radfield':1.0,
    'fr':1.0,
    'evap':0,
    'mantle':0.0,

    'desorb':True,
    'crdesorb':True,
    'uvcr':True,
    'h2desorb':True,

    'omega':0.5,
    'phi':1.0e5,
    'grainArea':2.4e-22,
    'radg':1.0e-5,
    'ebmaxcr':1.21e3,
    'ebmaxh2':1.21e3,
    'ebmaxcrf':1.21e3,
    'ebmaxuvcr':1.0e4,
    'uvcreff':1.0e-3,
    'epsilon':0.01,
    'uvy':0.1,
    'h':0.0}

    with open(file,"r") as inFile:
        for line in inFile.readlines():
            for bits in line.split(';'):
                #print bits
                vals=bits.split('=')
                if len(vals)>1:
                    #print vals[0].lower()
                    if vals[0] in cloud.keys():
                        if 'desorb' in vals[0] or vals[0]=='uvcr':
                            cloud[vals[0]]=(int(vals[1])==1)
                        else:
                            cloud[vals[0]]=float(vals[1].replace('d','e'))
    return cloud 