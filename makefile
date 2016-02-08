#compile=/opt/intel/Compiler/11.1/046/bin/intel64/ifort -fp-stack-check -check all -g -traceback -O2
compile= /usr/bin/gfortran -ffree-line-length-0 -fbacktrace
#compile=/opt/intel/Compiler/11.1/046/bin/intel64/ifort
physics=cloud.f90


main: chem.o physics.mod main.f90 parameters.f90
	${compile} -o main physics.o DLSODE/opkdmain.o DLSODE/opkda1.o DLSODE/opkda2.o chem.o main.f90

chem.o: chem.f90 physics.o DLSODE/opkdmain.o DLSODE/opkda1.o DLSODE/opkda2.o rates.f90
	${compile} -c chem.f90
 
physics.o: ${physics}
	${compile} -c ${physics} -o physics.o

DLSODE/opkdmain.o: DLSODE/opkdmain.f
	${compile} -c DLSODE/opkdmain.f -o DLSODE/opkdmain.o
DLSODE/opkda1.o: DLSODE/opkda1.f
	${compile} -c DLSODE/opkda1.f -o DLSODE/opkda1.o 
DLSODE/opkda2.o: DLSODE/opkda2.f
	${compile} -c  DLSODE/opkda2.f -o DLSODE/opkda2.o
clean: 
	rm *.o *.mod DLSODE/*.o main
