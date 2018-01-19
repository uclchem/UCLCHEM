#demonstration of plotfunctions. called from main UCLCHEM directory
#it reads full UCLCHEM output and saves a plot of the abudances of select species 

from plotfunctions import *

#pick species, any number is fine
speciesNames=["CH3OH","#CH3OH","H","#H","H2"]

fig,axes=plt.subplots(1,2,figsize=(16,6),dpi=300,sharey=True,sharex=True)

#call read_uclchem. 
time,dens,temp,abundances=read_uclchem("output/fullexplode.dat",speciesNames)
ax=axes[0]

time,dens,temp,Habunds=read_uclchem("output/fullexplode.dat",["H","#H","H2","#H2"])
Habunds=np.asarray(Habunds)
Habunds=Habunds[0,:]+Habunds[1,:]+2.0*(Habunds[2,:]+Habunds[3,:])
colours=make_colours(len(speciesNames))

for specIndx,specName in enumerate(speciesNames):
    ax.plot(time,abundances[specIndx],color=colours.next(),label=specName)
ax.plot(time,Habunds,color="black")
ax.legend(loc=4,fontsize='small')
ax.set_xlabel('Time / years')
ax.set_ylabel("X$_{Species}$")
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(10,6e6)
ax.set_ylim(1e-20,2)
ax.set_title('EXPLOSIONS')

#call read_uclchem. 
time,dens,temp,abundances=read_uclchem("output/full.dat",speciesNames)
time,dens,temp,Habunds=read_uclchem("output/full.dat",["H","#H","H2","#H2"])
Habunds=np.asarray(Habunds)
Habunds=Habunds[0,:]+Habunds[1,:]+2.0*(Habunds[2,:]+Habunds[3,:])
ax=axes[1]

colours=make_colours(len(speciesNames))
for specIndx,specName in enumerate(speciesNames):
    ax.plot(time,abundances[specIndx],color=colours.next(),label=specName)
ax.plot(time,Habunds,color="black")

ax.legend(loc=4,fontsize='small')
ax.set_xlabel('Time / years')
ax.set_ylabel("X$_{Species}$")
ax.set_yscale('log')
ax.set_title('Standard Cloud')

#plot species and save to test.png, alternatively send dens instead of time.
plt.savefig("output/audreycompare.png")
'''
time,density,temperature=np.loadtxt("fort.99",unpack=True)
H,H2,H2O,CH3OH=np.loadtxt("fort.89",unpack=True)
fig,[ax,ax2]=plt.subplots(1,2,figsize=(16,6),dpi=300)
ax.plot(time,H2O,label="H$_2$O")
ax.plot(time,CH3OH,label="CH$_3$OH")
ax.legend()
ax2.plot(time,density)
ax2.plot(time,temperature,color="red")
ax.set_yscale('log')
ax2.set_yscale('log')

ax.set_xscale('log')
ax2.set_xscale('log')

plt.savefig("output/explosioncycle.png")
'''
habund,mantle,totalatoms=np.loadtxt("fort.82",unpack=True)
#A=zip(mantle,habund,totalatoms)
#A.sort()
#mantle,habund,totalatoms=zip(*A)
mantle=np.asarray(mantle)
habund=np.asarray(habund)
totalatoms=np.asarray(totalatoms)
timeSteps=np.arange(0,len(mantle))
fig,ax=plt.subplots()
ax.plot(timeSteps,mantle,color='red')
ax.plot(timeSteps,habund/totalatoms,color='blue')
ax.axhline(0.05,lw=0.5,color='black')
#ax.legend()
plt.savefig("output/explosionlimit.png")
