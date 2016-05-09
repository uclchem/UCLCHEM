#enter two grain files to see all reactions that are one file but not the other

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

nSpecies, speciesList, massList,evaptypes,bindener = read_species_file(speciesFile)
nSpecies2, speciesList2, massList2,evaptypes2,bindener2 = read_species_file(speciesFile2)

make_capitals(speciesFile)
make_capitals(speciesFile2)

for species in speciesList:
	if species not in speciesList2:
		print species