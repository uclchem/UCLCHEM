SHELL = /bin/sh
FC =gfortran
f2pyFC=gnu95
#f2pyFC=intelem #this must match FC (ifort=intelem, gfortran=gnu95)
#FC=ifort
FFLAGS = -O2 -fPIC -ffree-line-length-0
PHYSICS=cloud.f90
EXT_SUFFIX = $(shell python3-config --extension-suffix)

##simple makefile for uclchem
##user must point the FC variable to their preferred fortran compiler
##builds ode solver, physics module and chemistry module before linking together for main

##User can also compile a python module which contains the subroutines in wrap.f90 as functions
##to do this run "python -m numpy.f2py -c --help-fcompiler" to find your fortran compiler
##and edit the python makef

##physics module selected by changing physics variable to chosen fortran file.
main: chemistry.o physics.o main.f90 dvode.o defaultparameters.f90 readparameters.f90
	${FC} ${FFLAGS} -o ../uclchem swap_function.o constants.o physics.o dvode.o chemistry.o network.o heating.o photoreactions.o coolant_module.o main.f90

python: chemistry.o physics.o wrap.f90 dvode.o defaultparameters.f90
	python3 -m numpy.f2py -c --fcompiler=${f2pyFC} swap_function.o constants.o physics.o dvode.o chemistry.o network.o heating.o photoreactions.o coolant_module.o wrap.f90 -m uclchem 
	mv *${EXT_SUFFIX} ../uclchem.so

chemistry.o: odes.f90 rates.f90  chemistry.f90 network.o physics.o dvode.o heating.o photoreactions.o
	${FC} ${FFLAGS} -c chemistry.f90

heating.o: swap_function.o coolant_module.o heating.f90
	${FC} ${FFLAGS} -c heating.f90

physics.o: ${PHYSICS} Makefile constants.o
	${FC} ${FFLAGS} -c ${PHYSICS} -o physics.o

%.o: %.f90 
	${FC} ${FFLAGS} -c $<

heating_test:
	python3 -m numpy.f2py -c --fcompiler=${f2pyFC}  -m heating_test test_heatings.f90
	mv heating_test*${EXT_SUFFIX} ../heating_tests/heating_tests.so

clean: 
	rm *.o *.mod ../uclchem ../uclchem.so