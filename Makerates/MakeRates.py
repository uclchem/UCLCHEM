#! /usr/bin/python
#! /usr/fink/bin/python
#
####################################################################################################
#				MakeRates
#		Current version by Jonathan Holdship & Antonios Makrymallis. Original by Tom Bell.
#		MakeRates reads in lists of species and reactions and produces the network files needed
#		by UCLCHEM to run. It also performs basic cleaning and sanity checks on the network.		
#		
####################################################################################################

import math
import os
import string
import struct
import sys
import time
import fileinput
import itertools
from Functions import *

reactionFile = 'inputFiles/umist12.csv'
reactionFile_grain = 'inputFiles/uclgrainbasic.csv'
speciesFile = 'inputFiles/uclspeciesbasic.csv'

if not os.path.exists('outputFiles'):
    os.makedirs('outputFiles')

make_capitals(reactionFile)
make_capitals(reactionFile_grain)
make_capitals(speciesFile)

print '\n################################################'
print 'Reading and checking input'
print '################################################\n'

#read species names,masses and evaporation details from input speciesFile
nSpecies, speciesList = read_species_file(speciesFile)
speciesList=remove_duplicate_species(speciesList)

# Read the reactants, products, Arrhenius equation parameters and measurement labels for each reaction
# IF the reaction involves species in our Species List
nReactions1, reactions1 = read_reaction_file(reactionFile, speciesList,'UMIST')
nReactions2, reactions2 = read_reaction_file(reactionFile_grain,speciesList,'UCL')
reactionList=reactions1+reactions2
reactionList=add_desorb_reactions(speciesList,reactionList)

#Keep only the species that are involved in the final reaction list
print '\nRemoving unused species...'
speciesList = filter_species(speciesList,reactionList)

#TODO replace this with a atom counter
# Calculate the molecular mass and elemental constituents of each species
print 'Calculating molecular masses and elemental constituents...'
speciesList = find_constituents(speciesList)

#sort the species file according to mass
print 'Sorting species by mass...'
speciesList.sort(key=lambda x: int(x.mass))

speciesList.append(Species(["E-",0,0,0,0,0,0]))
speciesList[-1].n_atoms=1
#check reactions to see if there are potential problems
print "Checking reactions..."
reaction_check(speciesList,reactionList)

print '\n################################################'
print 'Checks complete, writing output files'
print '################################################\n'

#Create the species file
print '\nWriting final species file...'
filename = 'outputFiles/species.csv'
write_species(filename,speciesList)
print '\tFinal Species File:',filename

# Create the reaction file
print 'Writing final reaction file...'
filename = 'outputFiles/reactions.csv'
write_reactions(filename, reactionList)
print '\tFinal Reaction File:',filename

#TODO this doesn't work now I use species objects
# Write the ODEs in the appropriate language format
print 'Writing system of ODEs in F95 format...'
filename = 'outputFiles/odes.f90'
write_odes_f90(filename, speciesList, reactionList)
print '\tFinal ODE file:',filename

print 'Writing Network File...'
filename= 'outputFiles/network.f90'
write_network_file(filename,speciesList,reactionList)
print '\tFinal Network file:',filename

ngrain=0
for species in speciesList:
	if species.name[0]=='#':
		ngrain+=1

print '\nnspec= '+str(len(speciesList))
print 'nreac= '+str(len(reactionList))
print 'ngrain='+str(ngrain)