import csv
import numpy as np
from copy import deepcopy as copy
#functions including
#1. simple classes to store all the information about each species and reaction.
#2. Functions to read in the species and reaction file and check for sanity
#3. Functions to write out files necessary for UCLCHEM


##########################################################################################
#1. simple classes to store all the information about each species and reaction.
#largely just to make the other functions more readable.
##########################################################################################
reaction_types=['PHOTON','CRP','CRPHOT','FREEZE','THERM','DESOH2','DESCR','DEUVCR',
			"H2FORM","ER","ERDES","LH","LHDES","BULKSWAP","SURFSWAP"]
#these reaction types removed as UCLCHEM does not handle them. 'CRH','PHOTD','XRAY','XRSEC','XRLYA','XRPHOT'
elementList=['H','D','HE','C','N','O','F','P','S','CL','LI','NA','MG','SI','PAH','15N','13C','18O']
elementMass=[1,2,4,12,14,16,19,31,32,35,3,23,24,28,420,15,13,18]
symbols=['#','@','+','-','(',')']

class Species:
	def __init__(self,inputRow):
		self.name=inputRow[0]
		self.mass=inputRow[1]
		self.binding_energy=float(inputRow[2])
		self.solidFraction=float(inputRow[3])
		self.monoFraction=float(inputRow[4])
		self.volcFraction=float(inputRow[5])
		self.enthalpy=float(inputRow[6])
		self.n_atoms=0

	def is_grain_species(self):
		if self.name in ["BULK","SURFACE"]:
			return True
		else:
			return (self.name[0] in ['#','@'])

	def is_surface_species(self):
		return self.name[0]=="#"

	def is_bulk_species(self):
		return self.name[0]=="@"

	def is_ion(self):
		return (self.name[-1]=="+" or self.name[-1]=="-")

	def find_constituents(self):
		speciesName=self.name[:]
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

		self.n_atoms=len(atoms)
		mass=0
		for atom in atoms:
			mass+=elementMass[elementList.index(atom)]
		if mass!=float(self.mass):
			print(f"Input mass of {self.name} does not match calculated mass of constituents")
			print("using calculated mass")
			self.mass=str(mass)


class Reaction:
	def __init__(self,inputRow):
		self.reactants=[inputRow[0],inputRow[1],self.NANCheck(inputRow[2])]
		self.products=[inputRow[3],self.NANCheck(inputRow[4]),self.NANCheck(inputRow[5]),self.NANCheck(inputRow[6])]
		self.alpha=float(inputRow[7])
		self.beta=float(inputRow[8])
		self.gamma=float(inputRow[9])
		self.templow=float(inputRow[10])
		self.temphigh=float(inputRow[11])
		self.reac_type=self.get_reaction_type()
		self.duplicate=False

		#body_count is the number of factors of density to include in ODE
		#we drop a factor of density from both the LHS and RHS of ODES
		#So reactions with 1 body have no factors of density which we manage by counting from -1
		self.body_count=-1
		for reactant in self.reactants:
			if (reactant not in reaction_types) and reactant!="NAN":
				self.body_count+=1
			if reactant in ["DESOH2","FREEZE"]:
				self.body_count+=1
			if reactant in ["LH","LHDES"]:
				self.body_count-=1

	def NANCheck(self,a):
		aa  = a if a else 'NAN'
		return aa

	def get_reaction_type(self):
		if (self.reactants[2] in reaction_types):
			return self.reactants[2]
		else:
			if self.reactants[1] in reaction_types:
				return self.reactants[1]
			else:
				return "TWOBODY"

	def convert_to_bulk(self):
		for i in range(len(self.reactants)):
			self.reactants[i]=self.reactants[i].replace("#","@")
		for i in range(len(self.products)):
			self.products[i]=self.products[i].replace("#","@")

	def same_reaction(self,other):
		if set(self.reactants)==set(other.reactants):
			if set(self.products)==set(other.products):
				return True
		return False

	def changes_surface_count(self):
		"""
		This checks whether a grain reaction changes number of particles on the surface
		2 reactants to 2 products won't but two reactants combining to one will.
		"""
		if len([x for x in self.reactants if "#" in x]) != len([x for x in self.products if "#" in x]):
			return True
		if len([x for x in self.reactants if "@" in x]) != len([x for x in self.products if "@" in x]):
			return True
		return False


	def changes_total_mantle(self):
		#If it's not just a movement between ice phases
		if ("BULK" not in self.reactants[1]) and ("SWAP" not in self.reactants[1]):
			# if the number of ice species changes
			if self.changes_surface_count():
				return True
			else:
				return False
		else:
			return False

	def print(self):
		print(" + ".join(self.reactants), "->","+".join(self.products))


##########################################################################################
#2. Functions to read in the species and reaction file and check for sanity
##########################################################################################

# Read the entries in the specified species file
def read_species_file(fileName):
	species_list=[]
	f = open(fileName, 'r')
	reader = csv.reader(f, delimiter=',', quotechar='|')
	for row in reader:
		if row[0]!="NAME" and "!" not in row[0] :
			species_list.append(Species(row))
	nSpecies = len(species_list)
	return nSpecies,species_list

# Read the entries in the specified reaction file and keep the reactions that involve the species in our species list
def read_reaction_file(fileName, species_list, ftype):
	reactions=[]
	dropped_reactions=[]
	keepList=['','NAN','#','E-','e-','ELECTR']
	keepList.extend(reaction_types)

	for species in species_list:
		keepList.append(species.name)			                                  
	

	if ftype == 'UMIST': 
		f = open(fileName, 'r')
		reader = csv.reader(f, delimiter=':', quotechar='|')
		for row in reader:
			if all(x in keepList for x in [row[2],row[3],row[4],row[5],row[6],row[7]]): #if all the reaction elements belong to the keeplist
				#umist file doesn't have third reactant so add space and has a note for how reactions there are so remove that
				reactions.append(Reaction(row[2:4]+['']+row[4:8]+row[9:]))
	
	if ftype == 'UCL':
		f = open(fileName, 'r')
		reader = csv.reader(f, delimiter=',', quotechar='|')
		for row in reader:
			if all(x in keepList for x in row[0:7]):	#if all the reaction elements belong to the keeplist
				if row[10]=="":
					row[10]=0.0
					row[11]=10000.0
				reactions.append(Reaction(row))	
			else:
				dropped_reactions.append(row)

	if ftype == "KIDA":
		for row in kida_parser(fileName):
			if all(x in keepList for x in row[0:7]):
				reactions.append(Reaction(row))

	nReactions = len(reactions)
	return nReactions, reactions, dropped_reactions

def kida_parser(kida_file):
	"""
	KIDA used a fixed format file so we read each line in the chunks they specify
	and use python built in classes to convert to the necessary types.
	"""
	str_parse=lambda x: str(x).strip().upper()

	kida_contents=[
	    [3,{str_parse:11}],
	    [1,{"skip":1}],
	    [5,{str_parse:11}],
	    [1,{"skip":1}],
	    [3,{float:10,"skip":1}],
	    [1,{"skip":27}],
	    [2,{int:6,"skip":1}],
	    [1,{int:2}],
	    [1,{"skip":11}]
	]
	rows=[]
	with open(kida_file,"r") as f:
		f.readline()
		for line in f:
			row=[]
			for item in kida_contents:
				for i in range(item[0]):
					for func,count in item[1].items():
						if func!="skip":
							a=line[:count]
							print(func,a)
							row.append(func(a))
						else:
							print("skip",count)
						line=line[count:]
			#ignore the ionpol and 3 body reacs in KIDA
			if row[-1]<4:
				rows.append(row[:7]+row[8:-1])
	return rows

def remove_duplicate_species(species_list):
	#check for duplicate species
	duplicates=0
	duplicate_list=[]
	for i in range(0,len(species_list)):
		for j in range(0,len(species_list)):
			if species_list[i].name==species_list[j].name:
				if (j!=i) and species_list[i].name not in duplicate_list:
					print("\t {0} appears twice in input species list".format(species_list[i].name))
					duplicate_list.append(species_list[i].name)

	for duplicate in duplicate_list:
		removed=False
		i=0
		while not removed:
			if species_list[i].name==duplicate:
				del species_list[i]
				print("\tOne entry of {0} removed from list".format(duplicate))
				removed=True
			else:
				i+=1
	return species_list

#Look for possibly incorrect parts of species list
def check_and_filter_species(species_list,reaction_list):
	#check for species not involved in any reactions
	lostSpecies=[]
	for species in species_list:
		keepFlag=False
		for reaction in reaction_list:
			if species.name in reaction.reactants or species.name in reaction.products:
				keepFlag=True
		if not keepFlag:
			lostSpecies.append(species.name)
			species_list.remove(species)

	print('\tSpecies in input list that do not appear in final list:')
	print('\t',lostSpecies)
	print('\n')
	for species in species_list:
		species.find_constituents()

	#add in pseudo-species to track mantle
	mantle_specs=[]
	new_spec=[999]*7
	new_spec[0]="BULK"
	mantle_specs.append(Species(new_spec))
	new_spec[0]="SURFACE"
	mantle_specs.append(Species(new_spec))
	species_list=species_list+mantle_specs
	return species_list

def create_bulk_species(species_list):
	speciesNames=[species.name for species in species_list]
	new_species=[]
	try:
		h2o_binding_energy=speciesNames.index("#H2O")
		h2o_binding_energy=species_list[h2o_binding_energy].binding_energy
	except:
		print("You are trying to create a three phase model but #H2O is not in your network")
		print("This is likely an error so Makerates will not complete")
		print("Try adding #H2O or switching to three_phase=False in Makerates.py")
		quit()
	for species in species_list:
		if species.is_surface_species:
			if not species.name.replace("#","@") in speciesNames:
				new_spec=copy(species)
				new_spec.name=new_spec.name.replace("#","@")
				new_spec.binding_energy=h2o_binding_energy
				new_species.append(new_spec)
	return species_list+new_species


#All species should freeze out at least as themselves and all grain species should desorb according to their binding energy
#This function adds those reactions automatically to slim down the grain file
def add_desorb_reactions(species_list,reaction_list):
	desorb_reacs=['DESOH2',"DESCR","DEUVCR","THERM"]

	for species in species_list:
		if species.is_surface_species():
			for reacType in desorb_reacs:
				newReaction=Reaction([species.name,reacType,'NAN',species.name[1:],'NAN','NAN','NAN',1,0,species.binding_energy,0.0,10000.0])
				reaction_list.append(newReaction)
		if species.is_bulk_species():
			newReaction=Reaction([species.name,"THERM",'NAN',species.name[1:],'NAN','NAN','NAN',1,0,species.binding_energy,0.0,10000.0])
			reaction_list.append(newReaction)
	return reaction_list


def add_chemdes_reactions(species_list,reaction_list):
	new_reacs=[]
	for reaction in reaction_list:
		if reaction.reac_type in ["LH","ER"]:
			new_reac=copy(reaction)
			new_reac.reac_type=new_reac.reac_type+"DES"
			new_reac.reactants[2]=new_reac.reactants[2]+"DES"
			for i,product in enumerate(new_reac.products):
				if ("#" in product):
					new_reac.products[i]=new_reac.products[i][1:]
				else:
					if product!="NAN":
						print("All Langmuir-Hinshelwood and Eley-Rideal reactions should be input with products on grains only.")
						print("The fraction of products that enter the gas is dealt with by Makerates and UCLCHEM.")
						print("the following reaction caused this warning")
						reaction.print()
			new_reacs.append(new_reac)

	reaction_list=reaction_list+new_reacs
	return reaction_list

def add_bulk_reactions(species_list,reaction_list):
	lh_reactions=[x for x in reaction_list if "LH" in x.reactants]
	lh_reactions=lh_reactions+[x for x in reaction_list if "LHDES" in x.reactants]
	new_reactions=[]
	for reaction in lh_reactions:
		new_reac=copy(reaction)
		new_reac.convert_to_bulk()
		new_reactions.append(new_reac)

	bulk_species=[x for x in species_list if "@" in x.name]
	for species in bulk_species:
		#add individual swapping
		new_reac_list=[species.name,"BULKSWAP","NAN",species.name.replace("@","#")]
		new_reac_list=new_reac_list+["NAN","NAN","NAN",1,0,0,0,10000]
		new_reac=Reaction(new_reac_list)
		new_reactions.append(new_reac)

		#and the reverse
		new_reac_list[0]=species.name.replace("@","#")
		new_reac_list[1]="SURFSWAP"
		new_reac_list[3]=species.name
		new_reac=Reaction(new_reac_list)
		new_reactions.append(new_reac)

	return reaction_list+new_reactions

#check reactions to alert user of potential issues including repeat reactions
#and multiple freeze out routes
def reaction_check(species_list,reaction_list,freeze_check=True):


	#first check for multiple freeze outs so user knows to do alphas
	print("\tSpecies with multiple freeze outs, check alphas:")
	for spec in species_list:
		freezes=0
		for reaction in reaction_list:
			if (spec.name in reaction.reactants and 'FREEZE' in reaction.reactants):
				freezes+=1
		if (freezes>1):
			print("\t{0} freezes out through {1} routes".format(spec.name,freezes))
		if freezes<1 and not spec.is_grain_species() and freeze_check:
			print("\t{0} does not freeze out".format(spec.name,freezes))

	#now check for duplicate reactions
	duplicate_list=[]
	print("\n\tPossible duplicate reactions for manual removal:")
	duplicates=False
	for i, reaction1 in enumerate(reaction_list):
		#if i not in duplicate_list:
			for j, reaction2 in enumerate(reaction_list):
				if i!=j:
					if reaction1.same_reaction(reaction2):
						print("\tReactions {0} and {1} are possible duplicates".format(i+1,j+1))
						reaction1.print()
						reaction2.print()
						duplicates=True
						#adjust temperatures so temperature ranges are adjacent
						if reaction1.temphigh > reaction2.temphigh:
							if reaction1.templow<reaction2.temphigh:
								print(f"\tReactions {i+1} and {j+1} have non-adjacent temperature ranges")
						reaction1.duplicate=True
						reaction2.duplicate=True
	
	if (not duplicates):
		print("\tNone")

#capitalize files
def make_capitals(fileName):
	a=open(fileName).read()
	output = open(fileName, mode='w')
	output.write(a.upper())
	output.close()


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
def write_species(fileName, species_list):
	f= open(fileName,'w')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')		
	nSpecies = len(species_list)
	for species in species_list:
		writer.writerow([species.name,species.mass,species.n_atoms])

# Write the reaction file in the desired format
def write_reactions(fileName, reaction_list):
	f = open(fileName, 'w')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	nReactions = len(reaction_list)
	for reaction in reaction_list:
		#if statement changes beta for ion freeze out to 1. This is how ucl_chem recognises ions when calculating freeze out rate
		if ('FREEZE' in reaction.reactants and reaction.reactants[0][-1]=='+'):
			reaction.beta=1
		writer.writerow(reaction.reactants+reaction.products+[reaction.alpha,reaction.beta,reaction.gamma,reaction.templow,reaction.temphigh])

def write_odes_f90(file_name, species_list, reaction_list,three_phase):
	output = open(file_name, mode='w')
	#go through every species and build two strings, one with eq for all destruction routes and one for all formation
	ydotString=build_ode_string(species_list,reaction_list,three_phase)
	output.write(ydotString)
	output.close()    

def write_jacobian(file_name,species_list):
	output=open(file_name,"w")
	species_names=""
	for i,species in enumerate(species_list):
		species_names+=species.name
		losses=species.losses.split("+")
		gains=species.gains.split("+")
		for j in range(1,len(species_list)+1):
			if species.name=="SURFACE":
				di_dj=f"J({i+1},{j})=SUM(J(surfaceList,{j}))\n"
				output.write(di_dj)
			elif species.name=="BULK":
				if species_names.count("@")>0:
					di_dj=f"J({i+1},{j})=SUM(J(bulkList,{j}))\n"
					output.write(di_dj)
			else:
				#every time an ode bit has our species in it, we remove it (dy/dx=a for y=ax)
				di_dj=[f"-{x}".replace(f"*Y({j})","",1) for x in losses if f"*Y({j})" in x]
				di_dj+=[f"+{x}".replace(f"*Y({j})","",1) for x in gains if f"*Y({j})" in x]
				#of course there might be y=a*x*x so we only replace first instance and if there's still an instance
				#we put a factor of two in since dy/dx=2ax for y=a*x*x
				di_dj=[x+"*2" if f"*Y({j})" in x else x for x in di_dj]

				#safeMantle is a stand in for the surface so do it manually here
				# since it's divided by safemantle, derivative is negative so sign flips and we get another factor of 1/safeMantle
				if species_list[j-1].name=="SURFACE":
					di_dj=[f"+{x}/safeMantle" for x in losses if f"/safeMantle" in x]
					di_dj+=[f"-{x}/safeMantle" for x in gains if f"/safeMantle" in x]
				if len(di_dj)>0:
					di_dj=f"J({i+1},{j})="+"".join(di_dj)+"\n"
					output.write(di_dj)

		#tackle density separately.
		j=j+1
		if species.name=="SURFACE":
			di_dj=f"J({i+1},{j})=SUM(J(surfaceList,{j}))\n"
			output.write(di_dj)
		elif species.name=="BULK":
			if species_names.count("@")>0:
				di_dj=f"J({i+1},{j})=SUM(J(bulkList,{j}))\n"
				output.write(di_dj)
		else:		
			di_dj=[f"-{x}".replace(f"*D","",1) for x in losses if f"*D" in x]
			di_dj+=[f"+{x}".replace(f"*D","",1) for x in gains if f"*D" in x]
			di_dj=[x+"*2" if f"*D" in x else x for x in di_dj]
			if len(di_dj)>0:
				di_dj=f"J({i+1},{j})="+("".join(di_dj))+"\n"
				output.write(di_dj)
	i=i+2
	di_dj=f"J({i},{i})=ddensdensdot(D)\n"
	output.write(di_dj)

	output.close()




def build_ode_string(species_list, reaction_list,three_phase):
	species_names=[]
	for i, species in enumerate(species_list):
		species_names.append(species.name)
		species.losses=""
		species.gains=""

	bulk_index=species_names.index("BULK")
	surface_index=species_names.index("SURFACE")
	total_swap=""

	for i,reaction in enumerate(reaction_list):
		ODE_BIT=f"+RATE({i+1})"
		
		#every body after the first requires a factor of density
		for body in range(reaction.body_count):
			ODE_BIT=ODE_BIT+f"*D"
		
		#then bring in factors of abundances
		for species in reaction.reactants:
			if species in species_names:
				ODE_BIT+=f"*Y({species_names.index(species)+1})"
			elif species=="BULKSWAP":
				ODE_BIT+="*bulkLayersReciprocal"
			elif species=="SURFSWAP":
				ODE_BIT+="*totalSwap/safeMantle"
			elif species in ["DEUVCR","DESCR","DESOH2","ER","ERDES"]:
				ODE_BIT=ODE_BIT+f"/safeMantle"
				if species=="DESOH2":
					ODE_BIT=ODE_BIT+f"*Y({species_names.index('H')+1})"
			elif ((species in ["THERM"]) and not (three_phase)):
				ODE_BIT+=f"*D/safeMantle"
			if "H2FORM" in reaction.reactants:
				#only 1 factor of H abundance in Cazaux & Tielens 2004 H2 formation so stop looping after first iteration
				break

		if "LH" in reaction.reactants[2]:
			if "@" in reaction.reactants[0]:
				ODE_BIT+="*bulkLayersReciprocal"


		#now add to ydot strings for each species
		for species in reaction.reactants:
			if species in species_names:
				#Eley-Rideal reactions take a share of total freeze out rate which is already accounted for
				#so we add as a loss term to the frozen version of the species rather than the gas version
				if ("ER" in reaction.reactants) and (not species_list[species_names.index(species)].is_surface_species()):
					species_list[species_names.index("#"+species)].losses+=ODE_BIT
				else:
					species_list[species_names.index(species)].losses+=ODE_BIT
				if reaction.reactants[1]=="BULKSWAP":
					total_swap+=ODE_BIT
		for species in reaction.products:
			if species in species_names:
				species_list[species_names.index(species)].gains+=ODE_BIT
	
	ode_string=""
	if three_phase:
		ode_string+=truncate_line(f"totalSwap={total_swap[1:]}\n\n")
	#First get total rate of change of bulk and surface by adding ydots
	for n,species in enumerate(species_list):
		if species.name[0]=="@":
			species_list[bulk_index].gains+=f"+YDOT({n+1})"
		elif species.name[0]=="#":
			species_list[surface_index].gains+=f"+YDOT({n+1})"

	for n,species in enumerate(species_list):
		ydot_string=species_ode_string(n,species)
		ode_string+=ydot_string

	#now add bulk transfer to rate of change of surface species after they've already been calculated
	if three_phase:
		ode_string+="!Update surface species for bulk growth, replace surfaceCoverage with alpha_des\n"
		ode_string+="!Since ydot(surface_index) is negative, bulk is lost and surface forms\n"

		ode_string+=f"IF (YDOT({surface_index+1}) .lt. 0) THEN\n    surfaceCoverage = MIN(1.0,safeBulk/safeMantle)\n"
		for n,species in enumerate(species_list):
			if species.name[0]=="#":
				bulk_version=species_names.index(species.name.replace("#","@"))
				ode_string+=f"    YDOT({n+1})=YDOT({n+1})-YDOT({surface_index+1})*surfaceCoverage*Y({bulk_version+1})/safeBulk\n"
			if species.name[0]=="@":
				ode_string+=f"    YDOT({n+1})=YDOT({n+1})+YDOT({surface_index+1})*surfaceCoverage*Y({n+1})/safeBulk\n"
		ode_string+="ELSE\n"
		for n,species in enumerate(species_list):
			if species.name[0]=="@":
				surface_version=species_names.index(species.name.replace("@","#"))
				ode_string+=f"    YDOT({n+1})=YDOT({n+1})+YDOT({surface_index+1})*surfaceCoverage*Y({surface_version+1})\n"
			if species.name[0]=="#":
				ode_string+=f"    YDOT({n+1})=YDOT({n+1})-YDOT({surface_index+1})*surfaceCoverage*Y({n+1})\n"
		ode_string+="ENDIF\n"

		#once bulk transfer has been added, odes for bulk and surface must be updated to account for it
		ode_string+="!Update total rate of change of bulk and surface for bulk growth\n"
		ode_string+=species_ode_string(bulk_index,species_list[bulk_index])
		ode_string+=species_ode_string(surface_index,species_list[surface_index])

	return ode_string

def species_ode_string(n,species):
	ydot_string=""
	if species.losses != '':
		loss_string = '    LOSS = '+species.losses[1:]+'\n'
		ydot_string+=loss_string
	if species.gains != '':
		prod_string = '    PROD = '+species.gains[1:]+'\n'
		ydot_string+=prod_string
	
	if ydot_string!="":
		ydot_string+=f"    YDOT({n+1}) = "
		#start with empty string and add production and loss terms if they exists
		if species.gains != '':
			ydot_string += 'PROD'
		if species.losses != '':
			ydot_string += '-LOSS'
		ydot_string+="\n"
	else:
		ydot_string=f"    YDOT({n+1}) = {0.0}\n"

	ydot_string = truncate_line(ydot_string)
	return ydot_string



#create a file containing length of each list of moleculetypes and then the two lists (gas and grain) of species in each type
#as  well as fraction that evaporated in each type of event
def write_evap_lists(openFile,species_list):
	gasIceList=[];surfacelist=[];solidList=[];monoList=[];volcList=[]
	binding_energygyList=[];enthalpyList=[];bulkList=[];iceList=[]

	for i,species in enumerate(species_list):
		if species.name[0]=='#':
			#find gas phase version of grain species. For #CO it looks for first species in list with just CO and then finds the index of that
			try:
				j=species_list.index(next((x for x in species_list if x.name==species.name[1:]))) 
			except:
				print("\n**************************************\nWARNING\n**************************************")
				print("{0} has no gas phase equivalent in network. Every species should at least freeze out and desorb.".format(species.name))
				print("ensure {0} is in the species list, and at least one reaction involving it exists and try again".format(species.name[1:]))
				print("Alternatively, provide the name of the gas phase species you would like {0} to evaporate as".format(species.name))
				input=input("type x to quit Makerates or any species name to continue\n")
				if input.lower()=="x":
					exit()
				else:
					j=species_list.index(next((x for x in species_list if x.name==input.upper())))					

			#plus ones as fortran and python label arrays differently
			surfacelist.append(i+1)
			gasIceList.append(j+1)
			solidList.append(species.solidFraction)
			monoList.append(species.monoFraction)
			volcList.append(species.volcFraction)
			iceList.append(i+1)
			binding_energygyList.append(species.binding_energy)
			enthalpyList.append(species.enthalpy)
		elif species.name[0]=="@":
			j=species_list.index(next((x for x in species_list if x.name==species.name[1:]))) 
			gasIceList.append(j+1)
			bulkList.append(i+1)
			iceList.append(i+1)
			binding_energygyList.append(species.binding_energy)
			enthalpyList.append(species.enthalpy)

	openFile.write(array_to_string("surfaceList",surfacelist,type="int"))
	if len(bulkList)>0:
		openFile.write(array_to_string("bulkList",bulkList,type="int"))
	openFile.write(array_to_string("iceList",iceList,type="int"))
	openFile.write(array_to_string("gasIceList",gasIceList,type="int"))
	openFile.write(array_to_string("solidFractions",solidList,type="float"))
	openFile.write(array_to_string("monoFractions",monoList,type="float"))
	openFile.write(array_to_string("volcanicFractions",volcList,type="float"))
	openFile.write(array_to_string("bindingEnergy",binding_energygyList,type="float",parameter=False))
	openFile.write(array_to_string("formationEnthalpy",enthalpyList,type="float"))

# Create the appropriate multiplication string for a given number
def multiple(number):
	if number == 1: return ''
	else: return str(number)+'*'

def truncate_line(input_string, codeFormat='F90', continuationCode=None,lineLength = 72):
	result = ''
	i=0
	j=0
	splits=["*",")","+",","]
	while len(input_string[i:])>lineLength:
		j=i+lineLength
		while input_string[j] not in splits:
			j=j-1
		result+=input_string[i:j]+"&\n    &"
		i=j
	result+=input_string[i:]
	return result    


def write_network_file(fileName,species_list,reaction_list,three_phase):
	openFile=open(fileName,"w")
	openFile.write("MODULE network\nUSE constants\nIMPLICIT NONE\n")
	openFile.write("    INTEGER, PARAMETER :: nSpec={0}, nReac={1}\n".format(len(species_list),len(reaction_list)))

	#write arrays of all species stuff
	names=[]
	atoms=[]
	masses=[]
	for species in species_list:
		names.append(species.name)
		masses.append(float(species.mass))
		atoms.append(species.n_atoms)

	speciesIndices=""
	for element in ["E-","C+","H+","H2","SI+","S+","CL+","CO","HE+","#H","#H2","#N","#O",'#OH',
					"SURFACE","BULK"]+elementList:
		try:
			species_index=names.index(element)+1
	
		except:
			print(element," not in network, adding dummy index")
			species_index=len(species_list)+1
		name=element.lower().replace("+","x").replace("e-","elec").replace("#","g")
		speciesIndices+="n{0}={1},".format(name,species_index)
	if len(speciesIndices)>72:
		speciesIndices=truncate_line(speciesIndices)
	speciesIndices=speciesIndices[:-1]+"\n"
	openFile.write("    INTEGER, PARAMETER ::"+speciesIndices)
	if three_phase:
		openFile.write("    LOGICAL, PARAMETER :: THREE_PHASE = .TRUE.\n")
	else:
		openFile.write("    LOGICAL, PARAMETER :: THREE_PHASE = .FALSE.\n")
	openFile.write(array_to_string("    specname",names,type="string"))
	openFile.write(array_to_string("    mass",masses,type="float"))
	openFile.write(array_to_string("    atomCounts",atoms,type="int"))


	#then write evaporation stuff
	write_evap_lists(openFile,species_list)

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
	#store important reactions
	reactionIndices=""

	for i,reaction in enumerate(reaction_list):
		if ("CO" in reaction.reactants) and ("PHOTON" in reaction.reactants):
			if "O" in reaction.products and "C" in reaction.products:
				reactionIndices+="nR_CO_hv={0},".format(i+1)
		if ("C" in reaction.reactants) and ("PHOTON" in reaction.reactants):
			reactionIndices+="nR_C_hv={0},".format(i+1)
		if ("H2FORM" in reaction.reactants):
			reactionIndices+=f"nR_H2Form_CT={i+1},"
		if (("H" in reaction.reactants) and ("#H" in reaction.reactants)):
			if "H2" in reaction.products:
				reactionIndices+=f"nR_H2Form_ERDes={i+1},"
			elif "#H2" in reaction.products:
				reactionIndices+=f"nR_H2Form_ER={i+1},"
		if ((reaction.reactants.count("#H")==2) and ("LH" in reaction.reactants)):
			reactionIndices+=f"nR_H2Form_LH={i+1},"
		if ((reaction.reactants.count("#H")==2) and ("LHDES" in reaction.reactants)):
			reactionIndices+=f"nR_H2Form_LHDes={i+1},"
		if (("H" in reaction.reactants) and ("FREEZE" in reaction.reactants)):
			reactionIndices+=f"nR_HFreeze={i+1},"
		if (("E-" in reaction.reactants) and ("FREEZE" in reaction.reactants)):
			reactionIndices+=f"nR_EFreeze={i+1},"
	reactionIndices=reactionIndices[:-1]


	if len(reactionIndices)>60:
		reactionIndices=reactionIndices[:60]+"&\n&"+reactionIndices[60:]
	reactionIndices=reactionIndices[:-1]+"\n"
	openFile.write("    INTEGER, PARAMETER ::"+reactionIndices)

	for i,reaction in enumerate(reaction_list):
		if "CO" in reaction.reactants and "PHOTON" in reaction.reactants:
			if "O" in reaction.products and "C" in reaction.products:
				openFile.write(f"INTEGER, PARAMETER :: nrco={i+1}\n")
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
		if reaction.duplicate:
			duplicates.append(i+1)
			tmaxs.append(reaction.temphigh)
			tmins.append(reaction.templow)
		reacTypes.append(reaction.reac_type)
	if len(duplicates)==0:
		duplicates=[9999]
		tmaxs=[0]
		tmins=[0]

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
	openFile.write(array_to_string("\tduplicates",duplicates,type="int",parameter=True))
	openFile.write(array_to_string("\tminTemps",tmins,type="float",parameter=True))
	openFile.write(array_to_string("\tmaxTemps",tmaxs,type="float",parameter=True))

	reacTypes=np.asarray(reacTypes)

	partners=get_desorption_freeze_partners(reaction_list)
	openFile.write(array_to_string("\tfreezePartners",partners,type="int",parameter=True))

	
	for reaction_type in reaction_types+["TWOBODY"]:
		list_name=reaction_type.lower()+"Reacs"
		indices=np.where(reacTypes==reaction_type)[0]
		if len(indices>1):
			indices=[indices[0]+1,indices[-1]+1]
		else:
			#We still want a dummy array if the reaction type isn't in network
			indices=[99999,99999]
		openFile.write(array_to_string("\t"+list_name,indices,type="int",parameter=True))
	openFile.write("END MODULE network")
	openFile.close()

def find_reactant(species_list,reactant):
	try:
		return species_list.index(reactant)+1
	except:
		return 9999

def get_desorption_freeze_partners(reaction_list):
	freeze_species=[x.products[0] for x in reaction_list if x.reactants[1]=="DESCR"]
	partners=[]
	for spec in freeze_species:
		for i,reaction in enumerate(reaction_list):
			if reaction.reac_type=="FREEZE":
				if reaction.reactants[0]==spec:
					partners.append(i+1)
					break
	return partners


def array_to_string(name,array,type="int",parameter=True):
	if parameter:
		outString=", PARAMETER :: "+name+" ({0})=(/".format(len(array))
	else:
		outString=" :: "+name+" ({0})=(/".format(len(array))
	if type=="int":
		outString="INTEGER"+outString
		for value in array:
			outString+="{0},".format(value)
	elif type=="float":
		outString="REAL(dp)"+outString
		for value in array:
			outString+="{0:.4e},".format(value)
	elif type=="string":
		strLength=len(max(array, key=len))
		outString="CHARACTER(Len={0:.0f})".format(strLength)+outString
		for value in array:
			outString+="\""+value.ljust(strLength)+"\","
	else:
		print("Not a valid type for array to string")
	outString=outString[:-1]+"/)\n"
	outString=truncate_line(outString)
	return outString
