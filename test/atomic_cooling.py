import numpy as np
import matplotlib.pyplot as plt

elecDens=1.0e-4*1.0e2
hDens=0.25*1.0e2
hxDens=0.05*1.0e2
hexDens=0.05*1.0e2
heDens=0.05*1.0e2

h_coll_ex=[]
he_coll_ex=[]

h_coll_ion=[]
he_coll_ion=[]
hex_coll_ion=[]

hx_coll_recomb=[]
hex_coll_recomb=[]
dielec=[]
free_free=[]
temps=np.logspace(1,6,10)
for temp in temps:
    gauntFactor=1.1+(0.34*np.exp(-((5.5-np.log10(temp))**2.0)/3.0))
    rootT=np.sqrt(temp)
    collTFactor=(1.0/(1.0+np.sqrt(temp*1.0e-5)))
    invT=1.0/temp

    h_coll_ex.append(7.5e-19*collTFactor*np.exp(-118348.0*invT)*elecDens*hDens)
    he_coll_ex.append(5.54e-17*(temp**-0.397)*collTFactor*np.exp(-473638.0*invT)*elecDens*hexDens)

    h_coll_ion.append(1.27e-21*rootT*np.exp(-157809.1*invT)*elecDens*hDens*collTFactor)
    he_coll_ion.append(9.38e-22*rootT*np.exp(-285335.4*invT)*elecDens*heDens*collTFactor)
    hex_coll_ion.append(4.95e-22*rootT*np.exp(-631515.0*invT)*elecDens*hexDens*collTFactor)


    hx_coll_recomb.append(8.7e-27*rootT*((1.0e-3*temp)**-0.2)*elecDens*hxDens/(1.0+((0.1*temp*1.0e-5)**0.7)))
    hex_coll_recomb.append(1.55e-26*(temp**0.3647)*elecDens*hexDens)
    dielec.append(1.24e-13*(temp**-1.5)*np.exp(-470000.0*invT)*(1.0+0.3*np.exp(-94000.0*invT))*elecDens*hexDens)

    free_free.append(1.42e-27*rootT*elecDens*(hexDens+hxDens)*gauntFactor)

fig,ax=plt.subplots(figsize=(16,9))

labels=["h_coll_ex","he_coll_ex", "h_coll_ion","he_coll_ion","hex_coll_ion","hx_coll_recomb","hex_coll_recomb","dielec","free_free"]

for i,coolings in enumerate([h_coll_ex,he_coll_ex,h_coll_ion,he_coll_ion,hex_coll_ion,hx_coll_recomb,hex_coll_recomb,dielec,free_free]):
    ax.plot(temps,coolings,label=labels[i])
ax.set(yscale='log',xscale='log',ylim=(1e-40,1e-15))
ax.legend()
plt.show()
fig.savefig("atomic_cooling.png")