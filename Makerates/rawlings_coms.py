import csv

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def count_atoms(speciesName):
	atoms=[]
	i=0
	elementList=['H','D','HE','C','N','O','F','P','S','CL','LI','NA','MG','SI','PAH','15N',"E-","FE"]
	symbols=['#','+','-','(',')']
	bracket=False
	bracketContent=[]
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
					if is_number(speciesName[j:j+2]):
						if int(speciesName[j:j+2])>1:
							for k in range(1,int(speciesName[j])):
								if bracket:
									bracketContent.append(speciesName[i:j])
								else:
									atoms.append(speciesName[i:j])
							i=j+2
						else:
							i=j
					elif is_number(speciesName[j]):
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
				print(speciesName[i:])
				print("\t{0} contains elements not in element list:".format(speciesName))
				return 0
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
	return len(atoms)

f=open("inputFiles/umist12.csv","r")
out=open("jonathan_com_reacs.csv","w")
reader = csv.reader(f, delimiter=':', quotechar='|')
writer= csv.writer(out,delimiter=":",quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
for row in reader:
	com_destruct=False
	for species in row[2:4]:
		if count_atoms(species)>5:
			com_destruct=True
	if com_destruct:
		writer.writerow(row)