from plotfunctions import *


cloud=getParameters("src/parameters.f90")
time,cloud,species,abunds=readTimestep("output/start.dat",1.00e7,cloud)

for element in ["O","C","H"]:
	print "\n**********\n"
	totalabund=0.0
	oxyAbunds=[]
	oxySpecs=[]
	for i in range(0,len(species)):
		nAtom=species[i].count(element)
		if element+"2" in species[i]:
			nAtom+=1
		if nAtom>0:
			totalabund+=(nAtom*abunds[i])
			oxyAbunds.append(nAtom*abunds[i])
			oxySpecs.append(species[i])

	print "{:.2e}".format(totalabund)
	A=zip(oxyAbunds,oxySpecs)
	A.sort(reverse=True)
	oxyAbunds,oxySpecs=zip(*A)
	for i in range(0,10):
		print "{0} {1:.3e}".format(oxySpecs[i],oxyAbunds[i])
