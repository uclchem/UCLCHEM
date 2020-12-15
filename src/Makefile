SHELL = /bin/sh
FC =gfortran
f2pyFC=gnu95 #this must match FC (ifort=intelem, gfortran=gnu95)
#FC=ifort
FFLAGS = -O2 -fPIC
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
	${FC} ${FFLAGS} -o ../uclchem physics.o dvode.o photoreactions.o surfacereactions.o chemistry.o network.o main.f90

python: chemistry.o physics.o wrap.f90 dvode.o defaultparameters.f90
	python3 -m numpy.f2py -c --fcompiler=${f2pyFC} physics.o dvode.o photoreactions.o surfacereactions.o chemistry.o network.o -m uclchem wrap.f90
	mv *${EXT_SUFFIX} ../uclchem.so

chemistry.o: odes.f90 constants.o photoreactions.o surfacereactions.o rates.f90 chemistry.f90 network.o physics.o dvode.o 
	${FC} ${FFLAGS} -c chemistry.f90

physics.o: ${PHYSICS} constants.o Makefile
	${FC} ${FFLAGS} -c ${PHYSICS} -o physics.o

%.o: %.f90 
	${FC} ${FFLAGS} -c $<

clean: 
	rm *.o *.mod ../uclchem ../uclchem.so