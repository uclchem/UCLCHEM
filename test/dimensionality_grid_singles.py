import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd
from os import path,mkdir
from pyDOE import lhs
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob

def run_model(a):
	i,params=a
	outSpecies=""
	uclchem.wrap.general(params,outSpecies,f"test/",f"{i:.0f}")
	return 1

def generate_elements(param_dict):
	metallicity=param_dict["metallicity"]
	param_dict["fc"]=metallicity*2.6e-4
	param_dict["fo"]=metallicity*4.6e-4
	param_dict["fn"]=metallicity*6.1e-5
	param_dict["fs"]=metallicity*1.3e-5
	param_dict["fsi"]=metallicity*1.0e-7
	param_dict["fcl"]=metallicity*3.162e-7
	param_dict["fmg"]=metallicity*3.98e-5
	return param_dict



sample_df=pd.read_csv("Dimensionality/grid_first_run.csv")

av_factor=6.289E-22
file_number=101
idx=sample_df["outputFile"]==f"Dimensionality/raw/{file_number:.0f}.csv"
run_df=sample_df.loc[idx]

for i, row in run_df.iterrows():

	params=	{
		"finalTime":1.0e6,
		"heatingFlag":True,
		"heatWriteFlag":False,
		"avFactor":av_factor}
	print(row.index)		
	for key in row.index:
		print(key,row[key])
		params[key]=row[key]
	params["ion"]=1
	params["gasDustMassRatio"]=100.0
	params["grainRadius"]=1.0e-5
	params=generate_elements(params)

	params["outputFile"]=f"test/{file_number:.0f}.csv"
	start=time.time()
	run_model((file_number,params))
	end=time.time()
	print(start,end,end-start)
