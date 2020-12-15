import numpy as np
import matplotlib.pyplot as plt

specs={
		"H":1,
		"H2":2,
		"CO":12,
		"CH3OH":31}


prefac=4.0/(3.0*np.sqrt(np.pi))
prefac=prefac/((2*1.38e-16)**1.5)

temps=np.arange(1,1000,5)
# prefac=prefac/(temps**1.5)

fig,ax=plt.subplots()
# for species,mass in specs.items():
# 	sticking=prefac*((mass*1.67e-24)**1.5)
# 	ax.plot(temps,sticking,label=species)
sticking=1.0/(1.0+((temps/464.0)**1.5))
ax.plot(temps,sticking)
#ax.set(yscale="log")
ax.legend()
plt.show()