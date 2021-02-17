##################################################################################
# Load up abundances from the last grid and pick a random one. Then run the model
# from that starting point.
##################################################################################
import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd
from os import path,mkdir
import random
import csv

limits=pd.read_csv("Dimensionality/component_ranges.csv").values
components=pd.read_csv("Dimensionality/component_vectors.csv").values
means=np.loadtxt("Dimensionality/species_means.csv")



def sample_set(id):
	av_factor=6.289E-22
	params=	{
		"finalTime":1.0e5,
		"heatingFlag":True,
		"heatWriteFlag":False,
		"avFactor":av_factor}	
	infile=open(f"samples/in{id}.csv","w")
	outfile=open(f"samples/out{id}.csv","w")
	input_writer = csv.writer(infile,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	output_writer = csv.writer(outfile,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')

	for i in range(1000):
		abundances=generate_abundances()
		params["initialTemp"]=random.uniform(10.0,1000.0)
		params["initialDens"]=10.0**random.uniform(2.0,6.0)
		params["radfield"]=10.0**random.uniform(0.0,5.0)
		params["zeta"]=10.0**random.uniform(0.0,5.0)
		params=generate_column_dens(params)
		#we're generating abundances but need metallicity for heating.
		#approximate as sum of abundances of major C species / solar C abundance
		params["metallicity"]=np.sum(abundances[[8,9,10,11,12,13,20,21,22,38,62,63,64,59,73,74,95,97]])/2.6e-4


		output_list=[params["initialDens"],params["initialTemp"],params["radfield"],params["zeta"],params["metallicity"]]
		output_list=output_list+[params["coldens"],params["h2col"],params["ccol"],*abundances]
		input_writer.writerow(output_list)
		abundances,dust_temp=run_model(i,params,abundances)

		output_list=[abundances[-2],dust_temp,*abundances[:-2]]
		output_writer.writerow(output_list)
	return 1

def run_model(i,params,abundances):
	outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
	abundances,dust_temp=uclchem.sampler(params,outSpecies,f"Dimensionality/raw/",f"{i:.0f}",abundances)
	return abundances,dust_temp

def generate_column_dens(param_dict):
	col_dens=random.uniform(0.0,24.0)
	col_dens=10.0**col_dens

	#after a certain av, c col is high.
	if col_dens>1e22:
		c_col=(10.0**random.uniform(16.7,19.7))
		ion=0
	else:
		c_col=(10.0**random.uniform(-15.0,-4.0))
		if c_col>1e-6:
			ion=1
		else:
			ion=2
		c_col=c_col*col_dens


	#same for H2
	if col_dens>1e23:
		h2_col=0.5
	else:
		h2_col=(10.0**random.uniform(-15.0,-0.3))

	param_dict["coldens"]=col_dens
	param_dict["ccol"]=c_col
	param_dict["h2col"]=h2_col*col_dens
	param_dict["rout"]=param_dict["coldens"]/param_dict["initialDens"]
	param_dict["ion"]=ion
	return param_dict


def generate_abundances():
	#first generate random point in PCA space
	abuns=np.zeros(15)
	for i in range(15):
		abuns[i]=np.random.uniform(*limits[i])
	#return to real space
	abuns=np.dot(abuns,components)	
	abuns=abuns+means
	#finally, very small abundances can come up negative so make them small again
	abuns=np.where(abuns<1.0e-20,1.0e-20,abuns)
	return abuns

if __name__ == '__main__':
	n_workers=24
	ids=range(n_workers)
	start=time.time()
	pool=Pool(n_workers)
	print("mapping...")
	rel=pool.map(sample_set,ids)
	print(rel)
	pool.close()
	pool.join()
	end=time.time()
	print(start,end,end-start)
