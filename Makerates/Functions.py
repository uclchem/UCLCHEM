from __future__ import print_function
import csv
import numpy

#functions including
#1. simple classes to store all the information about each species and reaction.
#2. Functions to read in the species and reaction file and check for sanity
#3. Functions to write out files necessary for UCLCHEM


##########################################################################################
#1. simple classes to store all the information about each species and reaction.
#largely just to make the other functions more readable.
##########################################################################################
class Species:
	def __init__(self,inputRow):
		self.name=inputRow[0]
		self.mass=inputRow[1]
		self.bindener=float(inputRow[2])
		self.solidFraction=float(inputRow[3])
		self.monoFraction=float(inputRow[4])
		self.volcFraction=float(inputRow[5])
		self.enthalpy=float(inputRow[6])

	def is_grain_species(self):
		return self.name[0]=='#'

	def is_ion(self):
		return (self.name[-1]=="+" or self.name[-1]=="-")

class Reaction:
	def __init__(self,inputRow):
		self.reactants=[inputRow[0],inputRow[1],self.NANCheck(inputRow[2])]
		self.products=[inputRow[3],self.NANCheck(inputRow[4]),self.NANCheck(inputRow[5]),self.NANCheck(inputRow[6])]
		self.alpha=float(inputRow[7])
		self.beta=float(inputRow[8])
		self.gamma=float(inputRow[9])
		if inputRow[10]=="":
			inputRow[10]=0.0
			inputRow[11]=40000.0
		self.templow=float(inputRow[10])
		self.temphigh=float(inputRow[11])
		self.duplicate=False
		if (len(inputRow)>12):
			try:
				self.energy=float(inputRow[12])
			except:
				pass

	def NANCheck(self,a):
		aa  = a if a else 'NAN'
		return aa

	def same_reaction(self,other):
		if set(self.reactants)==set(other.reactants):
			if set(self.products)==set(other.products):
				return True
		return False


reaction_types=['PHOTON','CRP','CRPHOT','FREEZE','CRH','PHOTD','THERM','XRAY','XRSEC','XRLYA','XRPHOT','DESOH2','DESCR','DEUVCR',"CHEMDES","DIFF"]
elementList=['H','D','HE','C','N','O','F','P','S','CL','LI','NA','MG','SI','PAH','15N']
elementMass=[1,2,4,12,14,16,19,31,32,35,3,23,24,28,420,15]
symbols=['#','+','-','(',')']

##########################################################################################
#2. Functions to read in the species and reaction file and check for sanity
##########################################################################################

# Read the entries in the specified species file
def read_species_file(fileName):
	speciesList=[]
	f = open(fileName, 'r')
	reader = csv.reader(f, delimiter=',', quotechar='|')
	for row in reader:
		if row[0]!="NAME" and "!" not in row[0] :
			speciesList.append(Species(row))
	nSpecies = len(speciesList)
	return nSpecies,speciesList

# Read the entries in the specified reaction file and keep the reactions that involve the species in our species list
def read_reaction_file(fileName, speciesList, ftype):
	print(fileName)
	reactions=[]
	dropped_reactions=[]
	keepList=['','NAN','#','E-','e-','ELECTR']
	keepList.extend(reaction_types)
	for species in speciesList:
		keepList.append(species.name)			                                  
	if ftype == 'UMIST': # if it is a umist database file
		f = open(fileName, 'r')
		reader = csv.reader(f, delimiter=':', quotechar='|')
		for row in reader:
			if all(x in keepList for x in [row[2],row[3],row[4],row[5],row[6],row[7]]): #if all the reaction elements belong to the keeplist
				#umist file doesn't have third reactant so add space and has a note for how reactions there are so remove that
				reactions.append(Reaction(row[2:4]+['']+row[4:8]+row[9:14]))
	if ftype == 'UCL':	# if it is a ucl made (grain?) reaction file
		f = open(fileName, 'r')
		reader = csv.reader(f, delimiter=',', quotechar='|')
		for row in reader:
			if all(x in keepList for x in row[0:7]):	#if all the reaction elements belong to the keeplist
				# row[10]=0.0
				# row[11]=10000.0
				reactions.append(Reaction(row))	
			else:
				dropped_reactions.append(row)

	nReactions = len(reactions)
	return nReactions, reactions, dropped_reactions

def remove_duplicate_species(speciesList):
		#check for duplicate species
	duplicates=0
	duplicate_list=[]
	for i in range(0,len(speciesList)):
		for j in range(0,len(speciesList)):
			if speciesList[i].name==speciesList[j].name:
				if (j!=i) and speciesList[i].name not in duplicate_list:
					print("\t {0} appears twice in input species list".format(speciesList[i].name))
					duplicate_list.append(speciesList[i].name)

	for duplicate in duplicate_list:
		removed=False
		i=0
		while not removed:
			if speciesList[i].name==duplicate:
				del speciesList[i]
				print("\tOne entry of {0} removed from list".format(duplicate))
				removed=True
			else:
				i+=1
	return speciesList

#Look for possibly incorrect parts of species list
def filter_species(speciesList,reactionList):
	#check for species not involved in any reactions
	lostSpecies=[]
	for species in speciesList:
		keepFlag=False
		for reaction in reactionList:
			if species.name in reaction.reactants or species.name in reaction.products:
				keepFlag=True
		if not keepFlag:
			lostSpecies.append(species.name)
			speciesList.remove(species)

	print('\tSpecies in input list that do not appear in final list:')
	print('\t',lostSpecies)
	print('\n')
	return speciesList

#All species should freeze out at least as themselves and all grain species should desorb according to their binding energy
#This function adds those reactions automatically to slim down the grain file
def add_desorb_reactions(speciesList,reactionList,therm_flag=False):
	if therm_flag:
		desorb_reacs=['DESOH2',"DESCR","DEUVCR","THERM"]
	else:
		desorb_reacs=['DESOH2',"DESCR","DEUVCR"]

	for species in speciesList:
		if species.is_grain_species():
			for reacType in desorb_reacs:
				newReaction=Reaction([species.name,reacType,'NAN',species.name[1:],'NAN','NAN','NAN',1,0,species.bindener,0.0,10000.0])
				reactionList.append(newReaction)
	return reactionList

#check reactions to alert user of potential issues including repeat reactions
#and multiple freeze out routes
def reaction_check(speciesList,reactionList,freeze_check=True):


	#first check for multiple freeze outs so user knows to do alphas
	print("\tSpecies with multiple freeze outs, check alphas:")
	for spec in speciesList:
		freezes=0
		for reaction in reactionList:
			if (spec.name in reaction.reactants and 'FREEZE' in reaction.reactants):
				freezes+=1
		if (freezes>1):
			print("\t{0} freezes out through {1} routes".format(spec.name,freezes))
		if freezes<1 and not spec.is_grain_species() and freeze_check:
			print("\t{0} does not freeze out".format(spec.name,freezes))
	
	#now check for duplicate reactions
	print("\n\tPossible duplicate reactions for manual removal:")
	duplicates=0
	for i, reaction1 in enumerate(reactionList):
		# if not reaction1.duplicate: #if a reaction is already marked then it's partner got printed already
			for j, reaction2 in enumerate(reactionList):
				if i!=j:
					if reaction1.same_reaction(reaction2):
						print("\tReactions {0} and {1} are possible duplicates".format(i+1,j+1))
						print("\t",str(i+1), reaction1.reactants, "-->", reaction1.products)
						print("\t",str(j+1), reaction1.reactants, "-->", reaction2.products)
						duplicates+=1
		
						#adjust temperatures so temperature ranges are adjacent
						if reaction1.temphigh > reaction2.temphigh:
							if reaction1.templow<reaction2.temphigh:
								print(f"\tReactions {i+1} and {j+1} have non-adjacent temperature ranges")
						reaction1.duplicate=True
						reaction2.duplicate=True
	
	if (duplicates==0):
		print("\tNone")

#capitalize files
def make_capitals(fileName):
	a=open(fileName).read()
	output = open(fileName, mode='w')
	output.write(a.upper())
	output.close()

def find_constituents(speciesList):
	for species in speciesList:
		speciesName=species.name
		i=0
		atoms=[]
		bracket=False
		bracketContent=[]
		#loop over characters in species name to work out what it is made of
		while i<len(speciesName):
			#if character isn't a #,+ or - then check it otherwise move on
			if speciesName[i] not in symbols:
				if i+1<len(speciesName):
					#if next two characters are (eg) 'MG' then atom is Mg not M and G
					if speciesName[i:i+3] in elementList:
						j=i+3
					elif speciesName[i:i+2] in elementList:
						j=i+2
					#otherwise work out which element it is
					elif speciesName[i] in elementList:
						j=i+1

				#if there aren't two characters left just try next one
				elif speciesName[i] in elementList:
					j=i+1
				#if we've found a new element check for numbers otherwise print error
				if j>i:
					if bracket:
						bracketContent.append(speciesName[i:j])
					else:
						atoms.append(speciesName[i:j])#add element to list
					if j<len(speciesName):
						if is_number(speciesName[j]):
							if int(speciesName[j])>1:
								for k in range(1,int(speciesName[j])):
									if bracket:
										bracketContent.append(speciesName[i:j])
									else:
										atoms.append(speciesName[i:j])
								i=j+1
							else:
								i=j
						else:
							i=j
					else:
						i=j
				else:
					print(speciesName[i])
					print("\t{0} contains elements not in element list:".format(speciesName))
					print(elementList)
			else:
				#if symbol is start of a bracketed part of molecule, keep track
				if (speciesName[i]=="("):
					bracket=True
					bracketContent=[]
					i+=1
				#if it's the end then add bracket contents to list
				elif speciesName[i]==")":
					if is_number(speciesName[i+1]):
						for k in range(0,int(speciesName[i+1])):
							atoms.extend(bracketContent)
						i+=2
					else:
						atoms.extend(bracketContent)
						i+=1
				#otherwise move on
				else:
					i+=1

		species.n_atoms=len(atoms)
		mass=0
		for atom in atoms:
			mass+=elementMass[elementList.index(atom)]
		if mass!=float(species.mass):
			species.mass=str(mass)
	return speciesList

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
##########################################################################################
#3. Functions to write out files necessary for UCLCHEM
##########################################################################################

# Write the species file in the desired format
def write_species(fileName, speciesList):
	f= open(fileName,'w')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')		
	nSpecies = len(speciesList)
	for species in speciesList:
		writer.writerow([species.name,species.mass,species.n_atoms])

# Write the reaction file in the desired format
def write_reactions(fileName, reactionList):
	f = open(fileName, 'w')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	nReactions = len(reactionList)
	for reaction in reactionList:
		#if statement changes beta for ion freeze out to 1. This is how ucl_chem recognises ions when calculating freeze out rate
		if ('FREEZE' in reaction.reactants and reaction.reactants[0][-1]=='+'):
			reaction.beta=1
		writer.writerow(reaction.reactants+reaction.products+[reaction.alpha,reaction.beta,reaction.gamma,reaction.templow,reaction.temphigh])

def write_odes_f90(fileName, speciesList, reactionList):
	output = open(fileName, mode='w')
	#go through every species and build two strings, one with eq for all destruction routes and one for all formation
	ydotString=build_ode_string(speciesList,reactionList)
	output.write(ydotString)
	output.close()    

def build_ode_string(speciesList, reactionList):
	odeString=""
	nSpecies = len(speciesList)
	nReactions = len(reactionList)

	for n,species in enumerate(speciesList):
		lossString = '' ; formString = ''
		#go through entire reaction list
		for i,reaction in enumerate(reactionList):
			
			
			#if species appear in reactants, reaction is a destruction route      	
			if species.name in reaction.reactants:
				bodyCount=0 #two or more bodies in a reaction mean we multiply rate by density so need to keep track
				#easy for h2 formation
				if is_H2_formation(reaction.reactants, reaction.products):
					lossString += '-2*RATE('+str(i+1)+')*D'
					continue
				#multiply string by number of time species appears in reaction. multiple() defined below
				#so far reaction string is rate(reaction_index) indexs are all +1 for fortran array indexing
				lossString += '-'+multiple(reaction.reactants.count(species.name))+'RATE('+str(i+1)+')'
				
				#now add *Y(species_index) to string for every reactant
				for reactant in set(reaction.reactants):
					n_appearances=reaction.reactants.count(reactant)
					#every species appears at least once in its loss reaction
					#so we multiply entire loss string by Y(species_index) at end
					#thus need one less Y(species_index) per reaction
					if reactant==species.name:
						bodyCount+=n_appearances
						if n_appearances > 1:
							for j,possibleReactants in enumerate(speciesList):
								if reactant == possibleReactants.name:
									for appearance in range(1,n_appearances):
										lossString += '*Y('+str(j+1)+')'
									continue
					else:
						#look through species list and find reactant
						for j,possibleReactants in enumerate(speciesList):
							if reactant == possibleReactants.name:
								for appearance in range(n_appearances):
									lossString += '*Y('+str(j+1)+')'
									bodyCount+=1
								continue
				#now string is rate(reac_index)*Y(species_index1)*Y(species_index2) may need *D if total rate is 
				#proportional to density
				if reaction.reactants.count('FREEZE') > 0 or reaction.reactants.count('DESOH2') > 0:
					lossString += '*D'
				for body in range(1,bodyCount):
					lossString+="*D"				


			#same process as above but rate is positive for reactions where species is positive
			if species.name in reaction.products:
				bodyCount=0 #two or more bodies in a reaction mean we multiply rate by density so need to keep track

				if is_H2_formation(reaction.reactants,reaction.products):
					#honestly H should be index 1 but lets check
					H_index=speciesList.index(next((x for x in speciesList if x.name=='H')))
					formString += '+RATE('+str(i+1)+')*Y('+str(H_index+1)+')*D'
					continue

				#multiply string by number of time species appears in reaction. multiple() defined below
				#so far reaction string is rate(reaction_index) indexs are all +1 for fortran array indexing
				formString += '+'+multiple(reaction.products.count(species.name))+'RATE('+str(i+1)+')'
				
				#now add *Y(species_index) to string for every reactant						
				for reactant in set(reaction.reactants):
					n_appearances=reaction.reactants.count(reactant)
					for j,possibleReactants in enumerate(speciesList):
						if reactant == possibleReactants.name:
							for appearance in range(n_appearances):
								formString += '*Y('+str(j+1)+')'
								bodyCount+=1
							continue

				#now string is rate(reac_index)*Y(species_index1)*Y(species_index2) may need *D if total rate is 
				#proportional to density
				if reaction.reactants.count('FREEZE') > 0 or reaction.reactants.count('DESOH2') > 0:
					formString += '*D'
				for body in range(1,bodyCount):
					formString+="*D"

		if lossString != '':
			lossString = '      LOSS = '+lossString+'\n'
			lossString = truncate_line(lossString)
			odeString+=lossString
		if formString != '':
			formString = '      PROD = '+formString+'\n'
			formString = truncate_line(formString)
			odeString+=formString
		
		#start with empty string and add production and loss terms if they exists
		ydotString=''
		if formString != '':
			ydotString += 'PROD'
			if lossString != '':
				ydotString += '+'
		if lossString != '':
			ydotString += 'Y('+str(n+1)+')*LOSS'

		#if we have prod and/or loss add ydotstring to odes
		if ydotString!='':
			ydotString = '      YDOT('+str(n+1)+') = '+ydotString+"\n"
			ydotString = truncate_line(ydotString)
			odeString+=ydotString
		else:
			ydotString = '      YDOT('+str(n+1)+') = 0.0\n'
			ydotString = truncate_line(ydotString)
			odeString+=ydotString
	return odeString

#create a file containing length of each list of moleculetypes and then the two lists (gas and grain) of species in each type
#as  well as fraction that evaporated in each type of event
def write_evap_lists(openFile,speciesList):
	grainlist=[];mgrainlist=[];solidList=[];monoList=[];volcList=[]
	bindEnergyList=[];enthalpyList=[]

	for i,species in enumerate(speciesList):
		if species.name[0]=='#':
			#find gas phase version of grain species. For #CO it looks for first species in list with just CO and then finds the index of that
			try:
				j=speciesList.index(next((x for x in speciesList if x.name==species.name[1:]))) 
			except:
				print("\n**************************************\nWARNING\n**************************************")
				print("{0} has no gas phase equivalent in network. Every species should at least freeze out and desorb.".format(species.name))
				print("ensure {0} is in the species list, and at least one reaction involving it exists and try again".format(species.name[1:]))
				print("Alternatively, provide the name of the gas phase species you would like {0} to evaporate as".format(species.name))
				input=raw_input("type x to quit Makerates or any species name to continue\n")
				if input.lower()=="x":
					exit()
				else:
					j=speciesList.index(next((x for x in speciesList if x.name==input.upper())))					

			#plus ones as fortran and python label arrays differently
			mgrainlist.append(i+1)
			grainlist.append(j+1)
			solidList.append(species.solidFraction)
			monoList.append(species.monoFraction)
			volcList.append(species.volcFraction)
			bindEnergyList.append(species.bindener)
			enthalpyList.append(species.enthalpy)

	openFile.write(array_to_string("gasGrainList",grainlist,type="int"))
	openFile.write(array_to_string("grainList",mgrainlist,type="int"))
	openFile.write(array_to_string("solidFractions",solidList,type="float"))
	openFile.write(array_to_string("monoFractions",monoList,type="float"))
	openFile.write(array_to_string("volcanicFractions",volcList,type="float"))
	openFile.write(array_to_string("bindingEnergy",bindEnergyList,type="float",parameter=False))
	openFile.write(array_to_string("formationEnthalpy",enthalpyList,type="float"))

	return len(grainlist)

# Create the appropriate multiplication string for a given number
def multiple(number):
    if number == 1: return ''
    else: return str(number)+'*'

def truncate_line(input, codeFormat='F90', continuationCode=None):
	lineLength = 72
	maxlines=300
	lines=0
	result = ''
	i=0
	j=0
	while i+j<len(input):
		j+=1
		if j>lineLength:
			#important not to break entries so split lines at ,s
			try:
				k=input[i+j-16:i+j].index(",")
			except:
				try:
					k=input[i+j-16:i+j].index("*")
				except:
					k=input[i+j-16:i+j].index(")")
			j=j-16+k
			result+=input[i:i+j]+"&\n    &"
			i=i+j
			j=0
	result+=input[i:i+j]
	return result    


def is_H2_formation(reactants, products):
    nReactants = len([species for species in reactants if species != ''])
    nProducts  = len([species for species in products  if species != ''])
    if nReactants == 2 and nProducts == 1:
        if reactants[0] == 'H' and reactants[1] == 'H' and products[0] == 'H2': return True
    if nReactants == 3 and nProducts == 2:
        if reactants[0] == 'H' and reactants[1] == 'H' and reactants[2] == '#' and products[0] == 'H2' and products[1] == '#': return True
    return False

def write_network_file(fileName,speciesList,reactionList,exotherm_reacs):
	openFile=open(fileName,"w")
	openFile.write("MODULE network\n    IMPLICIT NONE\n")
	openFile.write("    INTEGER, PARAMETER :: nSpec={0}, nReac={1}\n".format(len(speciesList),len(reactionList)))

	#write arrays of all species stuff
	names=[]
	atoms=[]
	masses=[]
	for species in speciesList:
		names.append(species.name)
		masses.append(float(species.mass))
		atoms.append(species.n_atoms)

	#store indices of important species for use in code
	speciesIndices=""
	for element in ["E-","C+","H+","H2","SI+","S+","CL+","CO","HE+","HCO+","H3O+","H3+","#H","#H2","#N","#O",'#OH']+elementList:
		try:
			species_index=names.index(element)+1
		except:
			print(element," not in network")
			species_index=9999
		name=element.lower().replace("+","x").replace("e-","elec").replace("#",'g')
		speciesIndices+="n{0}={1},".format(name,species_index)
	if len(speciesIndices)>50:
		speciesIndices=speciesIndices[:60]+"&\n&"+speciesIndices[60:]
	speciesIndices=speciesIndices[:-1]+"\n"
	openFile.write("    INTEGER, PARAMETER ::"+speciesIndices)


	#do the same for important reactions
	reactionIndices=""
	for i,reaction in enumerate(reactionList):
		if ("H2+" in reaction.reactants) and ("e-" in reaction.reactants):
			reactionIndices+="nR_H2x_e={0},".format(i+1)
		if ("H2+" in reaction.reactants) and ("H" in reaction.reactants):
			reactionIndices+="nR_H2x_H={0},".format(i+1)
		if ("H3+" in reaction.reactants) and ("e-" in reaction.reactants):
			reactionIndices+="nR_H3x_e={0},".format(i+1)
		if ("H3O+" in reaction.reactants) and ("e-" in reaction.reactants):
			reactionIndices+="nR_H3Ox_e={0},".format(i+1)
		if ("HCO+" in reaction.reactants) and ("e-" in reaction.reactants):
			reactionIndices+="nR_HCOx_e={0},".format(i+1)
		if ("HE+" in reaction.reactants) and ("e-" in reaction.reactants):
			reactionIndices+="nR_HEx_e={0},".format(i+1)
		if ("HE+" in reaction.reactants) and ("e-" in reaction.reactants):
			reactionIndices+="nR_HEx_e={0},".format(i+1)
		if ("C" in reaction.reactants) and ("PHOTON" in reaction.reactants):
			reactionIndices+="nR_C_hv={0},".format(i+1)

	if len(reactionIndices)>50:
		reactionIndices=reactionIndices[:60]+"&\n&"+reactionIndices[60:]
	reactionIndices=reactionIndices[:-1]+"\n"
	openFile.write("    INTEGER, PARAMETER ::"+reactionIndices)

	openFile.write(array_to_string("    specname",names,type="string"))
	openFile.write(array_to_string("    mass",masses,type="float"))
	openFile.write(array_to_string("    atomCounts",atoms,type="int"))


	#then write evaporation stuff
	nGrain=write_evap_lists(openFile,speciesList)

	#finally all reactions
	reactant1=[]
	reactant2=[]
	reactant3=[]
	prod1=[]
	prod2=[]
	prod3=[]
	prod4=[]
	alpha=[]
	beta=[]
	gama=[]
	reacTypes=[]
	duplicates=[]
	tmins=[]
	tmaxs=[]
	for i,reaction in enumerate(reactionList):
		reactant1.append(find_reactant(names,reaction.reactants[0]))
		reactant2.append(find_reactant(names,reaction.reactants[1]))
		reactant3.append(find_reactant(names,reaction.reactants[2]))
		prod1.append(find_reactant(names,reaction.products[0]))
		prod2.append(find_reactant(names,reaction.products[1]))
		prod3.append(find_reactant(names,reaction.products[2]))
		prod4.append(find_reactant(names,reaction.products[3]))
		alpha.append(reaction.alpha)
		beta.append(reaction.beta)
		gama.append(reaction.gamma)
		reacTypes.append(get_reaction_type(reaction.reactants[1],reaction.reactants[2]))
		if reaction.duplicate or reaction.gamma<-200.0:
			duplicates.append(i+1)
			tmaxs.append(reaction.temphigh)
			tmins.append(reaction.templow)

	openFile.write(array_to_string("\tre1",reactant1,type="int"))
	openFile.write(array_to_string("\tre2",reactant2,type="int"))
	openFile.write(array_to_string("\tre3",reactant3,type="int"))
	openFile.write(array_to_string("\tp1",prod1,type="int"))
	openFile.write(array_to_string("\tp2",prod2,type="int"))
	openFile.write(array_to_string("\tp3",prod3,type="int"))
	openFile.write(array_to_string("\tp4",prod4,type="int"))
	openFile.write(array_to_string("\talpha",alpha,type="float",parameter=False))
	openFile.write(array_to_string("\tbeta",beta,type="float",parameter=False))
	openFile.write(array_to_string("\tgama",gama,type="float",parameter=False))
	openFile.write(array_to_string("\treacType",reacTypes,type="string",parameter=True))

	openFile.write(array_to_string("\tduplicates",duplicates,type="int",parameter=True))
	openFile.write(array_to_string("\tminTemps",tmins,type="float",parameter=True))
	openFile.write(array_to_string("\tmaxTemps",tmaxs,type="float",parameter=True))


	exo_reactants1=[]
	exo_reactants2=[]
	exo_reac_indxs=[]
	exothermicities=[]
	for reaction in exotherm_reacs:
		for i,reaction2 in enumerate(reactionList):
			if reaction.same_reaction(reaction2):
				exo_reac_indxs.append(i+1)
				break
		exo_reactants1.append(names.index(reaction.reactants[0])+1)
		exo_reactants2.append(names.index(reaction.reactants[1])+1)
		exothermicities.append(reaction.energy)


	openFile.write(array_to_string("exoReactants1",exo_reactants1,type="int"))
	openFile.write(array_to_string("exoReactants2",exo_reactants2,type="int"))
	openFile.write(array_to_string("exoReacIdxs",exo_reac_indxs,type="int"))
	openFile.write(array_to_string("exothermicities",exothermicities,type="float"))


	openFile.write("END MODULE network")
	openFile.close()

def find_reactant(species_list,reactant):
	try:
		return species_list.index(reactant)+1
	except:
		return 9999

def get_reaction_type(reac2,reac3):
	if (reac3=="CHEMDES") or (reac3=="DIFF"):
		return reac3
	else:
		return reac2


def array_to_string(name,array,type="int",parameter=True):
	if parameter:
		outString=", parameter :: "+name+" ({0})=(/".format(len(array))
	else:
		outString=" :: "+name+" ({0})=(/".format(len(array))
	if type=="int":
		outString="integer"+outString
		for value in array:
			outString+="{0},".format(value)
	elif type=="float":
		outString="double precision"+outString
		for value in array:
			outString+="{0:.4e},".format(value)
	elif type=="string":
		strLength=len(max(array, key=len))
		outString="character(Len={0:.0f})".format(strLength)+outString
		for value in array:
			outString+="\""+value.ljust(strLength)+"\","
	else:
		print("Not a valid type for array to string")
	outString=outString[:-1]+"/)\n"
	outString=truncate_line(outString)
	return outString
