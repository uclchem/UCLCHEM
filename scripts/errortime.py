import os
import numpy as np
import subprocess as sp
import time
import string


results=open("results2.dat","w")
rule = string.maketrans('D', 'E')
times,dens,goodco,goodh2o,goodmco,goodnh3 = np.loadtxt('cloudform',unpack=True,\
                converters = {
                1: lambda val: float(val.translate(rule)),\
      			4: lambda val: float(val.translate(rule)),\
                2: lambda val: float(val.translate(rule)),\
                3: lambda val: float(val.translate(rule)),\
				5: lambda val: float(val.translate(rule)),\
				0: lambda val: float(val.translate(rule))})
sp.call('make',shell=True)
length=len(times)
#better write a shell script that does this....
#could have it change tolerances and outputfile and compile different output programs?
#then python runs the programs and times
j=0 		
for reltol in ['1d-5','1d-7','1d-9','1d-12','1d-15']:
	for abstol in ['1d-20','1d-15','1d-10']:
		with open("testparams","w") as w:
			w.write(abstol+","+reltol)
		start = time.time()
		sp.call('./main > /dev/null',shell=True)
		stop=time.time()
		dure=stop-start
		times,dens,co,h2o,mco,nh3 = np.loadtxt('testoutput',unpack=True,\
                converters = {
                1: lambda val: float(val.translate(rule)),\
      			4: lambda val: float(val.translate(rule)),\
                2: lambda val: float(val.translate(rule)),\
                3: lambda val: float(val.translate(rule)),\
				5: lambda val: float(val.translate(rule)),\
				0: lambda val: float(val.translate(rule))})
		xco=0;xh2o=0;xmco=0;xnh3=0
		if len(times) < length:
			xco=100000000
			xh2o=100000000
			xmco=10000000
			xnh3=10000000
		for i in range(len(times)):
			xco+=(((co[i]-goodco[i])/co[i])**2)
			xh2o+=(((h2o[i]-goodh2o[i])/h2o[i])**2)
			xmco+=(((mco[i]-goodmco[i])/mco[i])**2)
			xnh3+=(((nh3[i]-goodnh3[i])/nh3[i])**2)
		xco=np.sqrt(xco)/len(times)
		xh2o=np.sqrt(xh2o)/len(times)
		xmco=np.sqrt(xmco)/len(times)
		xnh3=np.sqrt(xnh3)/len(times)
		x=xco+xh2o+xmco+xnh3
		x=x/4
		j+=1
		results.write(str(reltol)+" "+str(abstol)+" "+str(dure)+" "+str(x)+" "+str(xco)+" "+str(xh2o)+" "+str(xmco)+" "+str(xnh3)+"\n")
		print str(j)+"th test complete of 15. legnth "+str(len(times))
