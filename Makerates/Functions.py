#! /usr/bin/python
#! /usr/fink/bin/python

import math
import os
import string
import struct
import sys
import time
import csv
import numpy

#create a file containing length of each list of moleculetypes and then the two lists (gas and grain) of species in each type
def evap_lists(filename,species,evaptype,bindener):
	colist=[];mcolist=[];intlist=[];mintlist=[];grainlist=[];mgrainlist=[]
	co2list=[];mco2list=[];int2list=[];mint2list=[];coener=[];co2ener=[];intener=[]
	for i in range(len(species)):
		if species[i][0]=='#':
			j=species.index(species[i][1:])
			#plus ones as fortran and python label arrays differently
			mgrainlist.append(i+1)
			grainlist.append(j+1)

			#a bunch of if statements better represented by a switch. put each grain species into the right TYPE of evap list
			if (evaptype[i] == 'CO1'):
				colist.append(j+1)
				mcolist.append(i+1)
				coener.append(bindener[i])
			if (evaptype[i] == 'INT1'):
				intlist.append(j+1)
				mintlist.append(i+1)
				intener.append(bindener[i])
			if (evaptype[i] == 'CO2'):
				co2list.append(j+1)
				mco2list.append(i+1)
				co2ener.append(bindener[i])
			if (evaptype[i] == 'INT2'):
				int2list.append(j+1)
				mint2list.append(i+1)
	f = open(filename, 'wb')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	f.write(str(len(colist))+'\n')
	writer.writerow(colist)
	writer.writerow(mcolist)
	writer.writerow(coener)
	f.write(str(len(co2list))+'\n')
	writer.writerow(co2list)
	writer.writerow(mco2list)
	writer.writerow(co2ener)
	f.write(str(len(intlist))+'\n')
	writer.writerow(intlist)
	writer.writerow(mintlist)
	writer.writerow(intener)
	f.write(str(len(int2list))+'\n')
	writer.writerow(int2list)
	writer.writerow(mint2list)
	f.write(str(len(grainlist))+'\n')
	writer.writerow(grainlist)
	writer.writerow(mgrainlist)


#capitalize files
def make_capitals(fileName):
	a=open(fileName).read()
	output = open(fileName, mode='w')
	output.write(a.upper())
	output.close()

# Read the entries in the specified species file
def read_species_file(fileName):
	f = open(fileName, 'rb')
	reader = csv.reader(f, delimiter=',', quotechar='|')

	species = [] ; mass = []; evaptype=[];bindener=[];c=0
	for row in reader:
		species.append(row[0])
		mass.append(row[1])
		evaptype.append(row[2])
		bindener.append(row[3])
	nSpecies = len(species)
	return nSpecies,species, mass,evaptype,bindener

def NANCheck(a):
	aa  = a if a else 'NAN'
	return aa


# Read the entries in the specified reaction file and keep the reactions that involve the species in our species list
def read_reaction_file(fileName, species, ftype):
	reactants = [] ; products = [] ; alpha = [] ; beta = [] ; gamma = [] ; templow = [] ;temphigh = []; keepList = []
	# keeplist includes the elements that ALL the reactions should be formed from 
	keepList.extend(['','NAN','#','E-','e-','ELECTR','PHOTON','CRP','CRPHOT','FREEZE','CRH','PHOTD','THERM','XRAY','XRSEC','XRLYA','XRPHOT','DESOH2','DESCR1','DESCR2','DEUVCR'])
	keepList.extend(species)			                                  
	if ftype == 'UMIST': # if it is a umist database file
		f = open(fileName, 'rb')
		reader = csv.reader(f, delimiter=':', quotechar='|')
		for row in reader:
			if all(x in keepList for x in [row[2],row[3],row[4],row[5],row[6],row[7]]): #if all the reaction elements belong to the keeplist
				reactants.append([row[2],row[3],'NAN'])
				products.append([row[4],NANCheck(row[5]),NANCheck(row[6]),NANCheck(row[7])])
				alpha.append(row[9])
				beta.append(row[10])
				gamma.append(row[11])
				templow.append(row[12])
				temphigh.append(row[13])
	if ftype == 'UCL':	# if it is a ucl made (grain?) reaction file
		f = open(fileName, 'rb')
		reader = csv.reader(f, delimiter=',', quotechar='|')
		for row in reader:
			if all(x in keepList for x in [row[1],row[2],row[3],row[4],row[5],row[6],row[7]]):	#if all the reaction elements belong to the keeplist
				reactants.append([row[1],row[2],NANCheck(row[3])])
				products.append([row[4],NANCheck(row[5]),NANCheck(row[6]),NANCheck(row[7])])
				alpha.append(float(eval(row[8])))
				beta.append(row[9])
				gamma.append(row[10])
				templow.append('NAN')
				temphigh.append('NAN')					
	nReactions = len(reactants)
	return nReactions, reactants, products, alpha, beta, gamma, templow, temphigh

#Get rid of species that are not involved in reactions
#we want to get rid of species that are not present in any reaction
def find_species(reactants,products,species, mass,evaptype,bindener):
	speciesList = [] ; keepList = []; extraSpecies = []
	keepList.extend(['','#','NAN','ELECTR','PHOTON','CRP','CRPHOT','FREEZE','CRH','PHOTD','THERM','XRAY','XRSEC','XRLYA','XRPHOT','DESOH2','DESCR1','DESCR2','DEUVCR'])
	speciesList.extend(reactants)
	speciesList.extend(products)
	speciesList = sum(speciesList,[])
	unique = set(speciesList)
	speciesList=list(unique)
	for n in keepList:
		try:
			speciesList.remove(n)
		except ValueError:
			pass 
		except AttributeError:
			pass	 
	extraSpecies.extend(list(set(species)-set(speciesList)))
	#extraSpecies.extend(list(set(speciesList)-set(species)))
	
	#Remove extra species
	for n in extraSpecies:
		ind = species.index(n)
		del species[ind]
		del bindener[ind]
		del evaptype[ind]
		del mass[ind]



	
	print '###################################################'
	print 'Differences between initial and final speciesList:' 
	print extraSpecies
	print '###################################################'
	return species, mass, evaptype, bindener





# Find the elemental constituents and molecular mass of each species in the supplied list
def find_constituents(speciesList):
    elementList = ['PAH','HE','LI','NA','MG','SI','CL','CA','FE','H','D','C','N','O','F','P','S','#','+','-']
    elementMass = [420.0,4.0,7.0,23.0,24.0,28.0,35.0,40.0,56.0,1.0,2.0,12.0,14.0,16.0,19.0,31.0,32.0,0,0,0]
    nElements = len(elementList)
    speciesConstituents = []
    speciesMass = []

    for species in speciesList:
        constituents = []
        for element in elementList:
            constituents.append(0)
            for n in range(species.count(element)):
                index = species.index(element)+len(element)
                if species[index:index+2].isdigit():
                    constituents[-1] += int(species[index:index+2])
                    species = species[:index-len(element)]+species[index+2:]
                elif species[index:index+1].isdigit():
                    constituents[-1] += int(species[index:index+1])
                    species = species[:index-len(element)]+species[index+1:]
                else:
                    constituents[-1] += 1
                    species = species[:index-len(element)]+species[index:]

        # Calculate the total molecular mass as the sum of the elemental masses of each constituent
        speciesMass.append(int(sum([float(constituents[i])*float(elementMass[i]) for i in range(nElements)])))

        # Sort the elements in the constituent list by their atomic mass
        zippedList = zip(elementMass, elementList, constituents)
        zippedList.sort()
        sortedMasses, sortedElements, constituents = zip(*zippedList)
        speciesConstituents.append(constituents)

    # Sort the list of elements by their atomic mass
    zippedList = zip(elementMass, elementList)
    zippedList.sort()
    sortedMasses, sortedElements = zip(*zippedList)
    return speciesMass, speciesConstituents, sortedElements

def sortSpecies(species, mass,evaptype,ener):
	sortedList = sorted(zip(mass,species,evaptype,ener))
	A= numpy.array(sortedList)
	species = list(A[:,1])
	mass = list(A[:,0])
	evaptype =list(A[:,2])
	ener = list(A[:,3])
	for i in range(len(species)):
		mass[i] = int(eval(mass[i]))
	return species, mass,evaptype,ener



# Write the reaction file in the desired format
def write_reactions(fileName, reactants, products, alpha, beta, gamma, templow, temphigh):
	f = open(fileName, 'wb')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	nReactions = len(reactants)
	f.write(str(nReactions)+'\n')
	for n in range(nReactions):
		#if statement changes beta for ion freeze out to 1. This is how ucl_chem recognises ions when calculating freeze out rate
		if (reactants[n][1]=='FREEZE' and reactants[n][0][-1]=='+'):
			beta[n]=1
		writer.writerow([reactants[n][0],reactants[n][1],reactants[n][2],products[n][0],products[n][1],products[n][2],products[n][3],alpha[n],beta[n],gamma[n],templow[n],temphigh[n]])

    
# Write the species file in the desired format
def write_species(fileName, speciesList, massList):
	f= open(fileName,'wb')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')		
	nSpecies = len(speciesList)
	f.write(str(nSpecies+1)+'\n')
	for n in range(nSpecies):
		writer.writerow([speciesList[n],massList[n]])

##############################################################################################################################
##############################################################################################################################################################

def write_odes_f90(fileName, speciesList, constituentList, reactants, products):
	nSpecies = len(speciesList)
	nReactions = len(reactants)
	output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
	if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
		xrayReactions = True
	else:
		xrayReactions = False

    # Prepare and write the electron conservation equation
    #output.write(conserve_species('e-', constituentList, codeFormat='F90'))
	output.write(electron_eq(speciesList,codeFormat='F90'))
    # Prepare and write the loss and formation terms for each ODE
	output.write('\n')
	for n in range(nSpecies):
		species = speciesList[n]
		lossString = '' ; formString = ''
		for i in range(nReactions):
#########################################################################################################       	
			if reactants[i].count(species) > 0:
				if is_H2_formation(reactants[i], products[i]):
					lossString += '-2*RATE('+str(i+1)+')*D'
					continue
				lossString += '-'+multiple(reactants[i].count(species))+'RATE('+str(i+1)+')'
				for reactant in speciesList:
					if reactant == species:
						for j in range(reactants[i].count(reactant)-1):
							lossString += '*Y('+str(speciesList.index(reactant)+1)+')'
						continue
					for j in range(reactants[i].count(reactant)):
						lossString += '*Y('+str(speciesList.index(reactant)+1)+')'
				for j in range(reactants[i].count('E-')):
					lossString += '*Y('+str(nSpecies+1)+')'
				if sum([speciesList.count(reactant) for reactant in reactants[i]]) > 1 or reactants[i].count('E-') > 0 or reactants[i].count('FREEZE') > 0 or reactants[i].count('DESOH2') > 0:
						lossString += '*D'	
##########################################################################################################
			if products[i].count(species) > 0:
				if is_H2_formation(reactants[i], products[i]):
					formString += '+RATE('+str(i+1)+')*Y('+str(speciesList.index('H')+1)+')*D'
					continue
				for k in range(products[i].count(species)):
					formString += '+RATE('+str(i+1)+')'
					for reactant in speciesList:
						for j in range(reactants[i].count(reactant)):
							formString += '*Y('+str(speciesList.index(reactant)+1)+')'
					for j in range(reactants[i].count('E-')):
						formString += '*Y('+str(nSpecies+1)+')'
					if sum([speciesList.count(reactant) for reactant in reactants[i]]) > 1 or reactants[i].count('E-') > 0 or reactants[i].count('FREEZE') > 0 or reactants[i].count('DESOH2') > 0:
						formString += '*D'
		if lossString != '':
			lossString = '      LOSS = '+lossString+'\n'
			lossString = truncate_line(lossString)
			output.write(lossString)
		if formString != '':
			formString = '      PROD = '+formString+'\n'
			formString = truncate_line(formString)
			output.write(formString)
		ydotString = '      YDOT('+str(n+1)+') = '
		if formString != '':
			ydotString += 'PROD'
			if lossString != '': ydotString += '+'
		if lossString != '':
			ydotString += 'Y('+str(n+1)+')*LOSS'
		ydotString += '\n'
		ydotString = truncate_line(ydotString)
		output.write(ydotString)
	output.close()

# Write the ODEs file in F77 language format
def write_odes_f77(fileName, speciesList, constituentList, reactants, products):
    nSpecies = len(speciesList)
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
    if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
        xrayReactions = True
    else:
        xrayReactions = False


    # Prepare and write the electron conservation equation
    output.write(electron_eq(speciesList, codeFormat='F77'))

    # Prepare and write the loss and formation terms for each ODE
    output.write('\n')
    for n in range(nSpecies):
        species = speciesList[n]
        lossString = '' ; formString = ''
        for i in range(nReactions):
            if reactants[i].count(species) > 0:
				if is_H2_formation(reactants[i], products[i]):
					lossString += '-2*K('+str(i+1)+')*D'
					continue
				lossString += '-'+multiple(reactants[i].count(species))+'K('+str(i+1)+')'
				for reactant in speciesList:
					if reactant == species:
						for j in range(reactants[i].count(reactant)-1):
							lossString += '*Y('+str(speciesList.index(reactant)+1)+')'
						continue
					for j in range(reactants[i].count(reactant)):
						lossString += '*Y('+str(speciesList.index(reactant)+1)+')'
				for j in range(reactants[i].count('E-')):
					#lossString += '*Y('+str(nSpecies+1)+')'
					lossString += '*X(1)'
				if sum([speciesList.count(reactant) for reactant in reactants[i]]) > 1 or reactants[i].count('E-') > 0 or reactants[i].count('FREEZE') > 0 or reactants[i].count('DESOH2') > 0:
						lossString += '*D'	
            if products[i].count(species) > 0:
                if is_H2_formation(reactants[i], products[i]):
                    formString += '+K('+str(i+1)+')*Y('+str(speciesList.index('H')+1)+')*D'
                    continue
                for k in range(products[i].count(species)):
                    formString += '+K('+str(i+1)+')'
                    for reactant in speciesList:
                        for j in range(reactants[i].count(reactant)):
                            formString += '*Y('+str(speciesList.index(reactant)+1)+')'
                    for j in range(reactants[i].count('E-')):
                        formString += '*X(1)'
                    if sum([speciesList.count(reactant) for reactant in reactants[i]]) > 1 or reactants[i].count('E-') > 0 or reactants[i].count('FREEZE') > 0 or reactants[i].count('DESOH2') > 0:
                        formString += '*D'
        if lossString != '':
            lossString = '      LOSS = '+lossString+'\n'
            lossString = truncate_line(lossString,codeFormat='F77')
            output.write(lossString)
        if formString != '':
            formString = '      PROD = '+formString+'\n'
            formString = truncate_line(formString,codeFormat='F77')
            output.write(formString)
        ydotString = '      YDOT('+str(n+1)+') = '
        if formString != '':
            ydotString += 'PROD'
            if lossString != '': ydotString += '+'
        if lossString != '':
            ydotString += 'Y('+str(n+1)+')*LOSS'
        ydotString += '\n'
        ydotString = truncate_line(ydotString,codeFormat='F77')
        output.write(ydotString)
    output.close()
    
def electron_eq(speciesList,codeFormat='F90'):
    elec_eq=''
    nSpecies=len(speciesList)
    for n in range(len(speciesList)):
        if speciesList[n][-1]=='+':
            if len(elec_eq)>0: elec_eq +='+'
            elec_eq+='Y('+str(n+1)+')'
        if speciesList[n][-1]=='-':
            if len(elec_eq)>0: elec_eq +='-'
            elec_eq+='+Y('+str(n+1)+')'        
    if len(elec_eq) > 0:
        if codeFormat == 'C':   elec_eq = '  x_e = '+elec_eq+';\n'
        if codeFormat == 'F90': elec_eq = '      Y('+str(nSpecies+1)+') = '+elec_eq+'\n'
        if codeFormat == 'F77': elec_eq = '      X(1)  = '+elec_eq+'\n'
    else:
        if codeFormat == 'C':   elec_eq = '  x_e = 0;\n'
        if codeFormat == 'F90': elec_eq = '      Y('+str(nSpecies+1)+') = 0\n'
        if codeFormat == 'F77': elec_eq = '      X(1)  = 0\n'
    if codeFormat == 'F77': elec_eq = truncate_line(elec_eq,codeFormat='F77')
    if codeFormat == 'F90': elec_eq = truncate_line(elec_eq)
    return elec_eq

# Create the conservation term for the desired species (an element, electron or dust grain)
def conserve_species(species, speciesConstituents, codeFormat='C'):
    elementList = ['#','+','-','H','D','HE','LI','C','N','O','F','NA','MG','SI','P','S','CL','CA','FE','PAH']
    nSpecies = len(speciesConstituents)
    conservationEquation = ''
    # Handle the special case of electrons (i.e., charge conservation with both anions and cations)
    if species == 'e-':
        indexPos = elementList.index('+')
        indexNeg = elementList.index('-')
        for n in range(nSpecies):
            if speciesConstituents[n][indexPos] > 0:
                if len(conservationEquation) > 0: conservationEquation += '+'
                if codeFormat == 'C':   conservationEquation += multiple(speciesConstituents[n][indexPos])+'x['+str(n)+']'
                if codeFormat == 'F90': conservationEquation += multiple(speciesConstituents[n][indexPos])+'Y('+str(n+1)+')'
                if codeFormat == 'F77': conservationEquation += multiple(speciesConstituents[n][indexPos])+'Y('+str(n+1)+')'
            if speciesConstituents[n][indexNeg] > 0:
                conservationEquation += '-'
                if codeFormat == 'C':   conservationEquation += multiple(speciesConstituents[n][indexNeg])+'x['+str(n)+']'
                if codeFormat == 'F90': conservationEquation += multiple(speciesConstituents[n][indexNeg])+'Y('+str(n+1)+')'
                if codeFormat == 'F77': conservationEquation += multiple(speciesConstituents[n][indexNeg])+'Y('+str(n+1)+')'
    else:
        index = elementList.index(species)
        for n in range(nSpecies):
            if speciesConstituents[n][index] > 0:
                if len(conservationEquation) > 0: conservationEquation += '+'
                if codeFormat == 'C':   conservationEquation += multiple(speciesConstituents[n][index])+'x['+str(n)+']'
                if codeFormat == 'F90': conservationEquation += multiple(speciesConstituents[n][index])+'Y('+str(n+1)+')'
                if codeFormat == 'F77': conservationEquation += multiple(speciesConstituents[n][index])+'Y('+str(n+1)+')'
    if len(conservationEquation) > 0:
        if codeFormat == 'C':   conservationEquation = '  x_e = '+conservationEquation+';\n'
        if codeFormat == 'F90': conservationEquation = '      Y('+str(nSpecies+1)+') = '+conservationEquation+'\n'
        if codeFormat == 'F77': conservationEquation = '      X(1)  = '+conservationEquation+'\n'
    else:
        if codeFormat == 'C':   conservationEquation = '  x_e = 0;\n'
        if codeFormat == 'F90': conservationEquation = '      Y('+str(nSpecies+1)+') = 0\n'
        if codeFormat == 'F77': conservationEquation = '      X(1)  = 0\n'
    if codeFormat == 'F77': conservationEquation = truncate_line(conservationEquation,codeFormat='F77')
    if codeFormat == 'F90': conservationEquation = truncate_line(conservationEquation)
    return conservationEquation

# Create the appropriate multiplication string for a given number
def multiple(number):
    if number == 1: return ''
    else: return str(number)+'*'

#check reactions to alert user of potential issues
def reaction_check(speciesList,reactants,products):
	print "\nDuplicate Species...."
	duplicates=0
	for i in range(0,len(speciesList)):
		for j in range(0,len(speciesList)):
			if speciesList[i]==speciesList[j]:
				if (j!=i):
					print str(spec)+" appears twice in input species list\n"
					duplicates+=1
	if duplicates==0:
		print"\nNo duplicate species"

	#first check for multiple freeze outs so user knows to do alphas
	print "\nSpecies with multiple freeze outs, check alphas:"
	for spec in speciesList:
		freezes=0
		for i in range(0,len(reactants)):
			if (reactants[i][0]==spec and reactants[i][1]=='FREEZE'):
				freezes+=1
		if (freezes>1):
			print spec+" freezes out through "+str(freezes)+" routes"
	#now check for duplicate reactions
	print "\nPossible duplicate reactions:"
	duplicates=0
	for i in range(0,len(reactants)):
		for j in range(0,len(reactants)):
			if (j != i):
				if (reactants[i][0]==reactants[j][0] and reactants[i][1]==reactants[j][1]):
					if (products[i][0]==products[j][0] and products[i][1]==products[j][1]): 
						print str(i+1), str(j+1), " reactants are ", reactants[i]
						duplicates+=1
					elif (products[i][1]==products[j][0] and products[i][0]==products[j][1]):
						print str(i+1), str(j+1), " reactants are ", reactants[i],products[i]
						duplicates+=1
				if (reactants[i][1]==reactants[j][0] and reactants[i][0]==reactants[j][1]):
					if (products[i][0]==products[j][0] and products[i][1]==products[j][1]): 
						print str(i+1), str(j+1), " reactants are ", reactants[i]
						duplicates+=1
					elif (products[i][1]==products[j][0] and products[i][0]==products[j][1]):
						print str(i+1), str(j+1), " reactants are ", reactants[i]
						duplicates+=1
	if (duplicates==0):
		print "None"


# Truncate long lines for use in fixed-format Fortran code
def truncate_line(input, codeFormat='F90', continuationCode=None):
    lineLength = 72
    maxlines=30
    lines=0
    result = ''
    index=input.index('=')
    lhs=input[:index]
    while len(input) > lineLength:
        #introduced a counter up to max continuation lines allow, if reached a new equation is started
        lines+=1
        if lines ==maxlines:
            index = max([input.rfind('+',0,lineLength),input.rfind('-',0,lineLength)])
            result += input[:index]+'\n'
            input=lhs+'='+lhs+input[index:]
            lines=0
        else:
            index = max([input.rfind('+',0,lineLength),input.rfind('-',0,lineLength),input.rfind('*',0,lineLength),input.rfind('/',0,lineLength)])
            if codeFormat == 'F90':
                if continuationCode != None: result += input[:index]+' '+continuationCode.strip()+'\n'
                else: result += input[:index]+' &\n'
            else:
                result += input[:index]+'\n'
            if continuationCode != None:
                input = continuationCode+input[index:]
            else:
                input = '     &       '+input[index:]
            

    result += input
    return result    


def is_H2_formation(reactants, products):
    nReactants = len([species for species in reactants if species != ''])
    nProducts  = len([species for species in products  if species != ''])
    if nReactants == 2 and nProducts == 1:
        if reactants[0] == 'H' and reactants[1] == 'H' and products[0] == 'H2': return True
    if nReactants == 3 and nProducts == 2:
        if reactants[0] == 'H' and reactants[1] == 'H' and reactants[2] == '#' and products[0] == 'H2' and products[1] == '#': return True
    return False