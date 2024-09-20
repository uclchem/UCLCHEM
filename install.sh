# Compile the dvode library, but the binary here for now (easier to find this way)
gfortran -shared -O2 -o libdvode.so -fPIC fortran/src/dvode.f90
# Move to the fortran directory, build the wrapper there:
cd src/fortran_src
FFLAGS="-O3 -fPIC -ffree-line-length-0" python -m numpy.f2py --fcompiler=gfortran -I$PWD -L$PWD -ldvode \
-c defaultparameters.f90  constants.f90 f2py_constants.f90  network.f90 physics-core.f90 cloud.f90 hotcore.f90 \ 
postprocess.f90 surfacereactions.f90 sputtering.f90 cshock.f90 jshock.f90 collapse.f90 photoreactions.f90 \
 network.f90 chemistry.f90 rates.f90 io.f90 odes.f90 wrap.f90 --build-dir build -m uclchemwrap
# Copy the binary into the python source
cp uclchemwrap* ../uclchem
# Go back to the root dir and  pip install locally
cd ../../
pip install -e .
