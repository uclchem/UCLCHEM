import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd
from os import path,mkdir

av_fac=6.289E-22
col_dens=10.0/av_fac

for density in [1e4]:
	for cr in [1e5]:


		rout=col_dens/density
		params={
			"columnFile":f"test/sv_nh{density}_cr{cr}_low_met.dat",
			"initialTemp":10.0,
			"initialDens":density,
			"finalTime":1.0e6,
			"radfield":1700.0,
			"zeta":cr,
			"outSpecies":3,
			"h2col":0.5*col_dens,
			"ccol":1.0e+18,
			"coldens":col_dens,
			"rout":rout,
			"fc":1.42e-4,
			"fhe":7.5e-2,
			# "fh":0.4,
			 "fo":3.2e-4,
			"fn":6.5e-5,
			 "fs":1.43e-6,#1.318e-8,
			"fcl":1.1e-7,#3.162e-10,
			"fsi":8.2e-7,#1.0e-10,
			"fmg":5.10e-06,
			"ion":2,
			"fr":0,
			"metallicity":0.1,
			"heatingFlag":True,
			"heatWriteFlag":False,
			"gasDustMassRatio":10.0,
			"avFactor":6.289E-22}

		start=time.time()
		outSpecies="H2,C,H3O+,SO,#CO"
		uclchem.wrap.general(params,outSpecies,f"test/","1")
		end=time.time()

		print(f"Completed in {end-start} seconds")