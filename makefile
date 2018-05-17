SHELL = /bin/sh
FC =gfortran -fbacktrace -ffree-line-length-0
#FC=ifort
FFLAGS = -O2 -fPIC 
PHYSICS=cloud.f90

##simple makefile for uclchem
##user must point the compile variable to their preferred fortran compiler
##builds ode solver, physics module and chemistry module before linking together for main

##physics module selected by changing physics variable to chosen fortran file.
main: chem.o physics.o main.f90 dvode.o parameters.f90
	${FC} ${FFLAGS} -o uclchem physics.o dvode.o chem.o main.f90

chem.o: odes.f90 chem.f90 physics.o dvode.o rates.f90
	${FC} ${FFLAGS} -c chem.f90
 
physics.o: ${PHYSICS}
	${FC} ${FFLAGS} -c ${PHYSICS} -o physics.o

dvode.o: dvode.f
	${FC} ${FFLAGS} -c dvode.f

clean: 
	rm *.o *.mod uclchem