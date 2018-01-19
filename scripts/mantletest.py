import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def make_colours(n):
    return iter(cm.rainbow(np.linspace(0, 1, n)))

atomsPerCore=1.0e8
jcrFactor=1.0/5.41061e8
coresAbund=1e-12


def myThresh(hAbund,mantle):
	return hAbund/(atomsPerCore*coresAbund+mantle)

def jcrThresh(hAbund,mantle):
	return hAbund/(atomsPerCore*jcrFactor+mantle)

def make_colours(n):
    return iter(cm.rainbow(np.linspace(0, 1, n)))

fig,axes=plt.subplots(1,2,figsize=(12,4.5),dpi=300,sharey=True)

hAbunds=10.0**np.arange(-4,-1,0.1)

colours=make_colours(5)
for mantle in np.arange(0,0.5,0.1):
	colour=colours.next()
	mantles=hAbunds+mantle
	axes[0].plot(hAbunds,myThresh(hAbunds,mantles),color=colour,label=mantle)
	axes[1].plot(hAbunds,jcrThresh(hAbunds,mantle),color=colour,label=mantle)

axes[0].axhline(0.05,color='black')
axes[1].axhline(0.05,color='black')
axes[0].set_yscale('log')
axes[1].set_yscale('log')

#axes[0].set_xlim(0,0.02)
#axes[1].set_xlim(0,0.02)

axes[1].legend()
plt.show()