###################### Script to run grid of models ####################
# Import statements
import numpy as np
import subprocess
import os

velocities = np.linspace(10,30,3)
densities = np.logspace(3,5,3)

# Create a single folder to house results
print('Making the results directory')
os.system('mkdir results')

# Create folders to store values and run UCLCHEM for the parameters
for v in velocities:
    v=int(v)
    for n in densities:
        n=int(n)

        print('Making: /results/v'+str(v)+'n'+str(n))
        os.system('mkdir /home/tjames/grid/results/v'+str(v)+'n'+str(n))

        # Compile code
        os.chdir('/home/tjames/grid/UCLCHEM/')
        print ('Compiling...')
        os.system('make')
        print('Compiling complete!')

        os.chdir('/home/tjames/grid/UCLCHEM/')
        print('Running UCLCHEM')
        os.system('./uclchem')

        os.system('cp /home/tjames/grid/UCLCHEM/output/full.dat /home/tjames/grid/results/'+'v'+str(v)+'n'+str(n)+'/full1.dat')
        os.system('cp /home/tjames/grid/UCLCHEM/output/column.dat /home/tjames/grid/results/'+'v'+str(v)+'n'+str(n)+'/column1.dat')
        os.system('cp /home/tjames/grid/UCLCHEM/output/startabund.dat /home/tjames/grid/results/'+'v'+str(v)+'n'+str(n)+'/startabund1.dat')

        os.chdir('/home/tjames/grid/UCLCHEM/')
        print('Running UCLCHEM')
        os.system('./uclchem')

        os.system('cp /home/tjames/grid/UCLCHEM/output/full.dat /home/tjames/grid/results/'+'v'+str(v)+'n'+str(n)+'/full2.dat')
        os.system('cp /home/tjames/grid/UCLCHEM/output/column.dat /home/tjames/grid/results/'+'v'+str(v)+'n'+str(n)+'/column2.dat')
        os.system('cp /home/tjames/grid/UCLCHEM/output/startabund.dat /home/tjames/grid/results/'+'v'+str(v)+'n'+str(n)+'/startabund2.dat')

        os.chdir('/home/tjames/grid/UCLCHEM/')
        os.system('make clean')
        print('Model run for v'+str(v)+'n'+str(n)+' complete!')
