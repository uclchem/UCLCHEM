#enter two species files to see all species that are one file but not the other

import math
import os
import string
import struct
import sys
import time
import fileinput
import itertools
from Functions import *



speciesFile = 'inputFiles/fullspecies.csv'
speciesFile2 = 'inputFiles/reducedspecies.csv'

nSpecies, speciesList = read_species_file(speciesFile)
nSpecies2, speciesList2 = read_species_file(speciesFile2)

make_capitals(speciesFile)
make_capitals(speciesFile2)

print "Species in ",speciesFile2," but not ",speciesFile
for species in speciesList:
	match=False
	for species2 in speciesList2:
		if species.name=species2.name:
			match=True
	if not match:
		print species.name

print "\n And vice-versa"
for species in speciesList2:
	match=False
	for species2 in speciesList:
		if species.name=species2.name:
			match=True
	if not match:
		print species.name