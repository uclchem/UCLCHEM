import numpy as np
​
velocities = np.linspace(30,60,7)
densities = np.logspace(3,6,4)
​
input = open("input.txt","w")
​
# Define the distance needed
# d=1e21 #1x10^21 cm
​
for v in velocities:
    for n in densities:
​
        # Define variables
	maxTemp=2000
	initialDens=1e2
        finalDens=n
	# finalTime=(d/(v*1e5))/(60*60*24*365) # Determine the time using velocity in cm/s and convert to years
        finalTime = 1e6
        fr=1.0
        switch=1
        first=1
        phase=1
        desorb=1
        h2desorb=1
        crdesorb=1
        uvcr=1
        evap=0
        vs=v
​
        input.write("%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n" % (maxTemp,initialDens,finalDens,finalTime,fr,switch,first,phase,desorb,h2desorb,crdesorb,uvcr,evap,vs))
​
        initialDens=n
        finalDens=n*2
        fr=0.0
        switch=0
        first=0
        phase=2
        desorb=1
        h2desorb=1
        crdesorb=1
        uvcr=1
        evap=1
​
        input.write("%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n" % (maxTemp,initialDens,finalDens,finalTime,fr,switch,first,phase,desorb,h2desorb,crdesorb,uvcr,evap,vs))