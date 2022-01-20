#enter two reaction files to see all reactions that are one file but not the other

from Functions import *


reactions1='../src/reactions.csv'
reactions2='/home/jon/Documents/zupcx4/Students/Damien/UCLCHEM/reactions.csv'
speciesFile = 'inputFiles/default_species.csv'

#differences are only relevant insofar as the missing reactions contain your species
nSpecies, speciesList = read_species_file(speciesFile)

make_capitals(reactions1)
make_capitals(reactions2)
make_capitals(speciesFile)


print('\nReading reactions')
nReactions2, reactions2,drops = read_reaction_file(reactions2,speciesList,'UCL')
nReactions1, reactions1,drops = read_reaction_file(reactions1, speciesList,'UCL')

print("Reactions from file 1 not in file 2")
for reaction1 in reactions1:
	match=False
	for reaction2 in reactions2:
		if reaction2==reaction1:
			match=True
			break
	if not match:
		print(reaction1.reactants,"-->",reaction1.products)

print("Reactions from file 2 not in file 1")
for reaction1 in reactions2:
	match=False
	for reaction2 in reactions1:
		if reaction2==reaction1:
			match=True
			break
	if not match:
		print(reaction1.reactants,"-->",reaction1.products)

print("Reactions with different coefficients")
for reaction1 in reactions1:
	for reation2 in reactions2:
		if reaction1==reaction2:
			if reaction1.alpha!=reaction2.alpha:
				print(reaction1.reactants,"-->",reaction1.products)
				print(f"alpha 1 = {reaction1.alpha}, alpha 2 = {reaction2.alpha}")
			if reaction1.beta!=reaction2.beta:
				print(reaction1.reactants,"-->",reaction1.products)
				print(f"beta 1 = {reaction1.beta}, beta 2 = {reaction2.beta}")
			if reaction1.gamma!=reaction2.gamma:
				print(reaction1.reactants,"-->",reaction1.products)
				print(f"gamma 1 = {reaction1.gamma}, gamma 2 = {reaction2.gamma}")