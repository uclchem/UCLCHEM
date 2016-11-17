import numpy as np
import csv
names,mass,type,be,monoFrac,volcFrac=np.loadtxt("inputFiles/species_latest_audrey.csv",dtype=str,delimiter=",",unpack=True,comments="NOT")
name2=np.loadtxt("inputFiles/basicspecies.csv",dtype=str,delimiter=",",unpack=True,comments="NOT")

f=open("inputFiles/uclspeciesbasic.csv","wb")
writer=csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
for i,name in enumerate(names):
	print name
	if name in name2:
		writer.writerow([names[i],mass[i],type[i],be[i],monoFrac[i],volcFrac[i]])
