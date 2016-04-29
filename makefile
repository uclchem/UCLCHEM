#compile=/opt/intel/Compiler/11.1/046/bin/intel64/ifort -fp-stack-check -check all -g -traceback -O2
compile= /usr/bin/gfortran -ffree-line-length-0 -fbacktrace
#compile=/opt/intel/Compiler/11.1/046/bin/intel64/ifort
physics=hydro.f90


main: chem.o physics.o main.f90 dvode.o parameters.f90
	${compile} -o main physics.o dvode.o chem.o main.f90

chem.o: chem.f90 physics.o dvode.o rates.f90
	${compile} -c chem.f90
 
physics.o: ${physics}
	${compile} -c ${physics} -o physics.o

dvode.o: dvode.f
	${compile} -c dvode.f

clean: 
	rm *.o *.mod main
