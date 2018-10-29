from __future__ import print_function
import uclchem
import numpy as np
from multiprocessing import Pool

#uclchem(final_dens,max_temp,shock_vel,phase_flag,outFile,startFile)
def runuclchem(a):
    modelNo,initial_dens,final_dens,fin_time,temp,outFile,startFile=a
    uclchem.cloud(initial_dens,final_dens,fin_time,temp,outFile,startFile)

print(__name__)

#Set up grid by writing all parameter combinations to a file
#and creating a list of lists of parameters
#pass all lists to pool.map() to be done in parallel
#remember to do phase 1s first if doing multiple phases
if __name__ == '__main__':
    outputFolder="output/"
    f=open(outputFolder+"grid.dat","w")

    init_dens=1e2 
    densities = np.logspace(3,6,4)
    modelNo=1

    models=[]
    for fin_dens in densities:
        n=int(np.log(fin_dens))
        outFile=outputFolder+"data/phase1-n{0}.dat".format(n)
        startFile=outputFolder+"start/{0}.dat".format((n))     
        models.append([modelNo,init_dens,fin_dens,6e6,10.0,outFile,startFile])
        f.write("{0} {1} {2} {3} {4} {5} {6}\n".format(modelNo,init_dens,fin_dens,6e6,10.0,outFile,startFile))
        modelNo+=1

    f.close()
    pool=Pool()
    pool.map(runuclchem,models)
    pool.close()
    pool.join()
