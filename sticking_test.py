import numpy as np
import matplotlib.pyplot as plt

specs={
		"H":1,
		"H2":2,
		"CO":12,
		"CH3OH":31}


amu=1.67e-24
kboltz=1.38e-16
temps=np.arange(1,1000,5)

stick_prefac=4.0/(3.0*np.sqrt(np.pi))
stick_prefac=stick_prefac/((2*kboltz)**1.5)
stick_prefac=stick_prefac/(temps**1.5)


freeze_prefac=np.sqrt(8.0*kboltz/(np.pi*amu))

fig,[ax,ax2]=plt.subplots(1,2)
sticking=1.0/(1.0+((temps/464.0)**1.5))
for species,mass in specs.items():
	freeze=freeze_prefac*np.sqrt(temps/mass)*8e-22
	#sticking=stick_prefac*((mass*amu)**1.5)
	ax.plot(temps,freeze*sticking,label=species)
	ax2.plot(temps,freeze,label=species)

 
# ax.plot(temps,sticking)
ax.set(yscale="log")
ax2.set(yscale="log")
ax.legend()
plt.show()