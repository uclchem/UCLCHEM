This is a disclaimer and help document for MakeRates and the networks used in UCLCHEM.

Disclaimer:
	UCLCHEM users are responsible for their own networks. A simplistic set of input files for MakeRates are provided as an example and skeleton from which to build. The gas-phase reactions in the UMIST database file are well referenced and experimentally verified. The grain file is user-defined and depends entirely on the assumptions made by the user. All users should take steps to verify their network and benchmark against established results.

############################
##HOW DO I MAKE A NETWORK?##
############################

MakeRates.py has, at the top, three file inputs. One each for a reaction file, a grain reaction file and a species file. These should be created according to the instructions below and then pointed to in MakeRates.py. After this, running "python MakeRates.py" from this directory will produce all necessary files for UCLCHEM in the output folder. Copy those into the UCLCHEM/src directory and make UCLCHEM to set up the code with the new network. It also produces reactions.csv and species.csv for inspection.

###INPUT FILES###

Reaction File:
	This is a UMIST database file by default. Any set of reactions formatted in the same way will work. This file contains the bulk of the network, a series of experimentally derived rate coefficients for gas phase reactions and so it is convenient to keep separate from user-defined reactions and rates.

Grain Reaction File:
	Formatted in exactly the same way as the gas reaction file, it is really just an extension of the same file. Here all the grain surface reactions can  be included. As a minimum, each species should have a freeze out reaction of the kind:
	
		CO,FREEZE,,#CO,,,,1,0,0,,,


	It can also be useful to define freeze out reactions where a species freezes as something other than itself. For example, we may expect CO hydrogenates quickly on the mantles and so choose to freeze it along multiple reactions to represent this:

		1762,CO,FREEZE,,#CO,,,,0.5,0,0,,,
		1857,CO,FREEZE,,#CH3OH,,,,0.1,0,0,,,
		1875,CO,FREEZE,,#HCO,,,,0.2,0,0,,,
		1875,CO,FREEZE,,#H2CO,,,,0.2,0,0,,,

	Note the alpha values are set to add up to 1. This ensures the total rate of CO freeze out remains the same and only the proportion in each form is set. For hydrogenation, we expect this to have no effect on the overall H abundance and so don't include it.

	Makerates automatically produces all the necessary code for the non-thermal desorption of any grain species in the network, these reactions do not need to be added to the grain file.


Species File:
	
	The species list is a list of every species used in the network. Any reactions with species not in this list are ignored. This list also gives the thermal desorption properties of each molecule. A sample entry looks like:

		#CO,28,1300,0.35,0.7,0.667,-27.2

	Where we give NAME,MASS,BINDING ENERGY,SOLID FRACTION,MONO FRACTION,VOLCANO FRACTION,ENTHALPY
. The first two are self-explantory. The latter five are desorption properties and only need to be set for grain species. 0,0,0,0,0 is acceptable for all gas species. The fractions are based on the work of Collings et al. 2004 and Viti et al. 2004. 


###OUTPUT FILES###
	Reactions:
		every reaction with its rate coefficients.

	Species:
		every species in the network. When referencing a species in the code (for output array for example) use the line number from this file -1 (the first line is the number of species).

	odes.f90:
		A fortran file with full set of ODEs for the network. Every species has it's own differential equation made up of reaction rate*abundances for every reaction the species is involved in. ODES.F90 IS COMPILED, NOT READ. You must make after changing odes.f90 and if you alter the filename, you must redirect chem.f90.f to the correct file.
	Network.f90:
		All of the reactions, species and their properties are written to arrays in network.f90. These are required by UCLCHEM to run.


#############################################################

There are a small number of useful scripts included in the directory. Comparenetworks, comparespecies and speciesnetwork print helpful information about input and output files to help manage networks. For example, comparenetworks prints the differences between two sets of reaction files and speciesnetwork prints every reaction that a species is involved in as two lists: formation and destruction.


 