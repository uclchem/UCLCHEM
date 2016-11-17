#enter two reaction files to see all reactions that are one file but not the other

import math
import os
import string
import struct
import sys
import time
import fileinput
import itertools
from Functions import *


reactions1='inputFiles/umist12.csv'
reactions2='inputFiles/rate13.rates'
speciesFile = 'inputFiles/species_latest_audrey.csv'

nSpecies, speciesList, massList,evaptypes,bindener,monoevap,volcevap = read_species_file(speciesFile)

make_capitals(reactions1)
make_capitals(reactions2)
make_capitals(speciesFile)


print '\nReading reactions'
nReactions2, reactants2, products2, alpha2, beta2, gamma2, templow2, temphigh2 = read_reaction_file(reactions2,speciesList,'UMIST')
nReactions1, reactants1, products1, alpha1, beta1, gamma1, templow1, temphigh1 = read_reaction_file(reactions1, speciesList,'UMIST')

print "Reactions from file 1 not in file 2"
for i in range(0,len(reactants1)):
	matches=0
	for j in range(0,len(reactants2)):
		if reactants1[i]==reactants2[j] and products1[i]==products2[j]:
			matches+=1
	if matches==0:
		print reactants1[i],products1[i]

print "Reactions from file 2 not in file 1"
for i in range(0,len(reactants2)):
	matches=0
	for j in range(0,len(reactants1)):
		if reactants2[i]==reactants1[j] and products2[i]==products1[j]:
			matches+=1
	if matches==0:
		print reactants2[i],products2[i]
