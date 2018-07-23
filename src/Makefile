SHELL = /bin/sh
FC =gfortran
#FC=ifort
FFLAGS = -O2 -fPIC
PHYSICS=cloud.f90
EXT_SUFFIX = $(shell python3-config --extension-suffix)

##simple makefile for uclchem
##user must point the compile variable to their preferred fortran compiler
##builds ode solver, physics module and chemistry module before linking together for main

##physics module selected by changing physics variable to chosen fortran file.
main: chem.o physics.o main.f90 dvode.o defaultparameters.f90 readparameters.f90
	${FC} ${FFLAGS} -o ../uclchem physics.o dvode.o chem.o main.f90

python: chem.o physics.o main.f90 dvode.o parameters.f90
	f2py -c --fcompiler=gnu95 physics.o dvode.o chem.o -m uclchem wrap.f90
	mv *${EXT_SUFFIX} ../uclchem.so

chem.o: odes.f90 chem.f90 physics.o dvode.o rates.f90
	${FC} ${FFLAGS} -c chem.f90

physics.o: ${PHYSICS} Makefile
	${FC} ${FFLAGS} -c ${PHYSICS} -o physics.o

dvode.o: dvode.f90
	${FC} ${FFLAGS} -c dvode.f90

clean: 
	rm *.o *.mod ../uclchem