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
from numpy import abs

reactions1='inputFiles/uclpdr-network.csv'
reactions2='inputFiles/umist12-ucledit.csv'
speciesFile = 'inputFiles/uclpdrspecies.csv'

def same_reaction(reaction1,reaction2):
	same=False
	if set(reaction2.reactants)==set(reaction1.reactants):
		if set(reaction2.products)==set(reaction2.products):
			same=True
	return same

#differences are only relevant insofar as the missing reactions contain your species
nSpecies, speciesList = read_species_file(speciesFile)

make_capitals(reactions1)
make_capitals(reactions2)
make_capitals(speciesFile)


print('\nReading reactions')
nReactions2, reactions2,drops = read_reaction_file(reactions2,speciesList,'UMIST')
nReactions1, reactions1,drops = read_reaction_file(reactions1, speciesList,'UCL')

print("Reactions from file 1 not in file 2")
for reaction1 in reactions1:
	match=False
	for reaction2 in reactions2:
		if same_reaction(reaction1,reaction2):
				match=True
	if not match:
		print(reaction1.reactants,"-->",reaction1.products)

print("Reactions from file 2 not in file 1")
for reaction1 in reactions2:
	match=False
	for reaction2 in reactions1:
		if same_reaction(reaction1,reaction2):
				match=True
	if not match:
		print(reaction1.reactants,"-->",reaction1.products)

print("Different alphas")
for reaction1 in reactions1:
	for reaction2 in reactions2:
		if same_reaction(reaction1,reaction2):
			if abs(reaction1.alpha-reaction2.alpha) > 0.01*reaction1.alpha:
				diff_alpha=True
				for reaction3 in reactions2:
					if same_reaction(reaction1,reaction3):
						if abs(reaction1.alpha-reaction3.alpha) < 0.01*reaction1.alpha:
							diff_alpha=False
				if diff_alpha:
					print(f"{reaction1.reactants} -> {reaction1.products}")
					print(f"network 1 alpha = {reaction1.alpha}, network 2 alpha = {reaction2.alpha}\n")