#! /usr/bin/python
#
####################################################################################################
#				MakeRates
#		Current version by Jonathan Holdship & Antonios Makrymallis. Original by Tom Bell.
#		MakeRates reads in lists of species and reactions and produces the network files needed
#		by UCLCHEM to run. It also performs basic cleaning and sanity checks on the network.		
#		
####################################################################################################
from Functions import *
import os


###################################################################################################
# USER OPTIONS
###################################################################################################

database_reaction_file = "inputFiles/umist12-ucledit.csv"
database_reaction_type="UMIST"

custom_reaction_file = 'inputFiles/default_grain_network.csv'


speciesFile = 'inputFiles/default_species.csv'

three_phase=True



#################################################################################################
if not os.path.exists('outputFiles'):
    os.makedirs('outputFiles')

make_capitals(database_reaction_file)
make_capitals(custom_reaction_file)
make_capitals(speciesFile)

print('\n################################################')
print('Reading and checking input')
print('################################################\n')

#read species names,masses and evaporation details from input speciesFile
nSpecies, speciesList = read_species_file(speciesFile)
speciesList=remove_duplicate_species(speciesList)
if three_phase:
	speciesList=create_bulk_species(speciesList)
	
# Read the reactants, products, Arrhenius equation parameters and measurement labels for each reaction
# IF the reaction involves species in our Species List
# Store user reactions (grain file) that are filtered out in list to write out
nReactions1, reactions1, dropped_reactions = read_reaction_file(database_reaction_file, speciesList,database_reaction_type)
nReactions2, reactions2, dropped_reactions = read_reaction_file(custom_reaction_file,speciesList,'UCL')
reactionList=reactions1+reactions2




#Need additional grain reactions including non-thermal desorption and chemically induced desorption
reactionList=add_desorb_reactions(speciesList,reactionList)
reactionList=add_chemdes_reactions(speciesList,reactionList)
if three_phase:
	reactionList=add_bulk_reactions(speciesList,reactionList)

#Keep only the species that are involved in the final reaction list
print('\nRemoving unused species...')
speciesList = check_and_filter_species(speciesList,reactionList)

#Print dropped reactions from grain file or write if many
if len(dropped_reactions)<6:
	print("Reactions dropped from grain file:\n")
	for reaction in dropped_reactions:
		print(reaction)
else:
	print("\nReactions dropped from grain file written to outputFiles/dropped_reactions.csv\n")
	f= open('outputFiles/dropped_reactions.csv','w')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')		
	for reaction in dropped_reactions:
		writer.writerow(reaction)
	f.close()

#sort the species file according to mass
print('Sorting species by mass...')
speciesList.sort(key=lambda x: int(x.mass))

speciesList.append(Species(["E-",0,0,0,0,0,0]))
speciesList[-1].n_atoms=1

#check reactions to see if there are potential problems
print("Checking reactions...")
reaction_check(speciesList,reactionList)


reactionList=sorted(reactionList,key=lambda x: (x.reac_type,x.reactants[0]))

print('\n################################################')
print('Checks complete, writing output files')
print('################################################\n')

#Create the species file
print('\nWriting final species file...')
filename = 'outputFiles/species.csv'
write_species(filename,speciesList)
print('\tFinal Species File:',filename)

# Create the reaction file
print('Writing final reaction file...')
filename = 'outputFiles/reactions.csv'
write_reactions(filename, reactionList)
print('\tFinal Reaction File:',filename)

# Write the ODEs in the appropriate language format
print('Writing system of ODEs in F95 format...')
filename = 'outputFiles/odes.f90'
write_odes_f90(filename, speciesList, reactionList,three_phase)
# filename = 'outputFiles/jacobian.f90'
# write_jacobian(filename,speciesList)
print('\tFinal ODE file:',filename)

print('Writing Network File...')
filename= 'outputFiles/network.f90'
write_network_file(filename,speciesList,reactionList,three_phase)
print('\tFinal Network file:',filename)

ngrain=0
for species in speciesList:
	if species.name[0]=='#':
		ngrain+=1

print('\nnspec= '+str(len(speciesList)))
print('nreac= '+str(len(reactionList)))
print('ngrain='+str(ngrain))