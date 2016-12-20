#A set of functions for working with UCLCHEM outputs
# adding "from plotfunctions import * to any python script in scripts/ will allow their use"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv

#to read from uclchem's full output send the filename of the output and a list of species names to read_uclchem()
#the list for example would be ["H","C","CO","#CO"]
#the return is a list of times and list of lists of abundances
#call time,abundances=read_uclchem("output-full",["H","C","CO","#CO"])
# abudnance[0] would be a list of "H" abundances
def read_uclchem(filename,species):
    a=open(filename).read()
    a=a.split('\n')
    #there are 68 lines per time step in the UCL_CHEM output file.
    lines=68
    timesteps=len(a)/lines

    abunds=[]
    for spec in species:
        abunds.append([])
    time=[]
    dens=[]
    #so now do the following until end of file
    with open(filename) as file:
        for line in file:
            bits=line.split()        
            #find time line
            if  'age' in bits:
                time.append(float(bits[-2].replace('D','E')))
            #read another line for dens
            if 'density' in bits:
                densi=float(bits[-2].replace('D','E'))
                if densi==0.0:
                    densi=1e-10
                dens.append(densi)
            #then read until we hit abundances
            if bits .count('=')>2:
                for specIndx,specName in enumerate(species):
                    if specName in bits:
                        abunds[specIndx].append(float(bits[2+bits.index(specName)].replace('D','E')))

    return time,dens,abunds

def write_cols(filename,times,dens,abundances):
    f=open(filename,"wb")
    for timeIndx,time in enumerate(times):
        outString="{0:.3e} {1:.3e}".format(time,dens[timeIndx])
        for i in range(0,len(abundances)):
            outString+=" {0:.3e}".format(abundances[i][timeIndx])
        outString+="\n"
        f.write(outString)
    f.close()

#send a  list of species names and their abundances in  a list of lists
#with a list of times for each abundance point
#same as the species input and time/abundance output from read_uclchem
#optionally send an output filename to save the plot
#return ax,figure for further manipulation
def plot_species(species,times,abundances,plotFile=None):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    colours=make_colours(len(species))

    for specIndx,specName in enumerate(species):
        ax.plot(times,abundances[specIndx],color=colours.next(),label=specName)

    ax.legend(loc=4,fontsize='small')

    ax.set_xlabel('Time / years')
    ax.set_ylabel("X$_{Species}$")

    ax.set_yscale('log')

    if plotFile is not None:
        fig.savefig(plotFile)
    return ax,fig
    


def make_colours(n):
    return iter(cm.rainbow(np.linspace(0, 1, n)))

