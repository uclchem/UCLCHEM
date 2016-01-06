#! /usr/bin/python
#! /usr/fink/bin/python
#
import math
import os
import string
import struct
import sys
import time
import fileinput
import itertools
from Functions import *

reactionFile = 'inputFiles/umist06.csv'
reactionFile_grain = 'inputFiles/uclgrainrates2.csv'
speciesFile = 'inputFiles/newspecies.csv'


make_capitals(reactionFile)
make_capitals(reactionFile_grain)
make_capitals(speciesFile)

print '\n########################\n########################\n'
# Read the species, abundances and masses in the specified species file
print '\nReading species file...'
nSpecies, speciesList, massList,evaptypes,bindener = read_species_file(speciesFile)


print '\n########################\n########################\n'
################################# GAS ##################################################################
# Read the reactants, products, Arrhenius equation parameters and measurement labels for each reaction
# IF the reaction involves species in our Species List
print '\nReading gas reaction file...'
nReactions1, reactants1, products1, alpha1, beta1, gamma1, templow1, temphigh1 = read_reaction_file(reactionFile, speciesList,'UMIST')

################################################# GRAIN ######################################################
# Read the reactants, products, Arrhenius equation parameters and measurement labels for each grain reaction
#if the reaction involves species in our species list
print '\nReading grain reaction file...'
nReactions2, reactants2, products2, alpha2, beta2, gamma2, templow2, temphigh2 = read_reaction_file(reactionFile_grain,speciesList,'UCL')


# join the reactions
print '\nJoin reactions...'
reactants = reactants1; reactants.extend(reactants2)
products = products1; products.extend(products2)
alpha = alpha1; alpha.extend(alpha2)
beta = beta1; beta.extend(beta2)
gamma = gamma1; gamma.extend(gamma2)
templow = templow1; templow.extend(templow2)
temphigh = temphigh1; temphigh.extend(temphigh2)

#Keep only the species that are involved in the final reaction list
print '\nGetting rid of unused species'
speciesList, massList,evaptypes,bindener = find_species(reactants,products,speciesList, massList,evaptypes,bindener)

# Calculate the molecular mass and elemental constituents of each species
print '\nCalculating molecular masses and elemental constituents...'
massList, constituentList, elementList = find_constituents(speciesList)

#sort the species file according to mass
print '\nSorting Species ...'
speciesList, massList, evaptypes, bindener = sortSpecies(speciesList, massList, evaptypes, bindener)

#Create the species file
print '\nWriting final species file '
filename = 'outputFiles/species.csv'
write_species(filename, speciesList, massList)
print 'Final Species File:',filename



# Create the reaction file
print 'Writing final reaction file '
filename = 'outputFiles/reactions.csv'
write_reactions(filename, reactants, products, alpha, beta, gamma, templow,temphigh)
print 'Final Reaction File:',filename

# Write the ODEs in the appropriate language format
print 'Writing system of ODEs in F95 format...'
filename = 'outputFiles/odes.f90'
write_odes_f90(filename, speciesList, constituentList, reactants, products)

print 'Writing Evaporation lists...'
filename= 'outputFiles/evaplists.csv'
evap_lists(filename,speciesList,evaptypes,bindener)

ngrain=0
for spec in speciesList:
	if spec[0]=='#':
		ngrain+=1
print '\nnspec= '+str(len(speciesList)+1)
print 'nreac= '+str(len(reactants))
print 'ngrain='+str(ngrain)