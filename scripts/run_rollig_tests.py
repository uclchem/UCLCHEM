import pandas as pd
import numpy as np
import sys
sys.path.insert(0,'Benchmarking')
import uclchem
import time
from multiprocessing import get_context

data_file="rollig_abundances_mod.dat"


classic_params={
	"collapse":0,
	"finalTime":1.0e6,
	"heatingFlag":True,
	"heatWriteFlag":False,
	"zeta":5.0/1.3,
	"metallicity":1.0,
	"avFactor":6.289E-22,
	"fr":0,
	"ion":2,
	"fh":0.4,
	"fhe":0.1,
	"fc":1.0e-4,
	"fo":3.0e-4,
	"fn":0.0,
	"fs":0.0,
	"fcl":0.0,
	"fsi":0.0,
	"fmg":5.00e-06,
	"fp":0.0,
	"ff":0.0,
	"pahAbund":6.0e-7,
	"gasDustMassRatio":1.0e2,
	"grainRadius":1.0e-5,
	"outSpecies":12}



def run_uclchem(params):
	outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
	abunds,gas_temp=uclchem.wrap.run_model_for_abundances(params,outSpecies)
	return gas_temp

if __name__=="__main__":

	columns=["x", "AV", "density", "temp","H2","H", "O+","O", "O2",
			"OH","C+", "C", "CH", "CO", "H2O", "OHPlus", "H2O+", "H3O+", "CO+","O2+","HCO+", "UV"]
	model_df=pd.read_csv(f"Benchmarking/haworth/{data_file}",header=None)
	model_df=model_df.rename(dict(zip(range(len(columns)),columns)),axis=1)
	print(model_df.head())

	#x in cm
	model_df["x"]=model_df["x"].values*1.0e10
	model_df["x"]=np.abs(model_df["x"].max()-model_df["x"].values)
	model_df=model_df.sort_values("x").reset_index(drop=True)
	print(model_df.head(10))

	#density in h nuclei
	model_df["H_nuclei_density"]=model_df["density"].values/1.67372346e-24

	model_df.loc[1:,"delta_x"]=(model_df.iloc[1:]["x"].values-model_df.iloc[0:-1]["x"].values)
	#can calulate change in column density between each particle as well
	model_df["delta_h2_col"]=model_df["delta_x"].values*model_df["H2"].values*model_df["H_nuclei_density"].values
	model_df["delta_c_col"]=model_df["delta_x"].values*model_df["C"].values*model_df["H_nuclei_density"].values
	model_df["delta_col_dens"]=model_df["delta_x"].values*model_df["H_nuclei_density"].values



	#and then integrate
	model_df["total_h2_col"]=model_df["delta_h2_col"].cumsum()
	model_df["total_c_col"]=model_df["delta_c_col"].cumsum()
	model_df["total_col_dens"]=model_df["delta_col_dens"].cumsum()
	model_df=model_df.fillna(0.0)
	
	models=np.linspace(0,len(model_df),24).astype(int)
	models=[0,10]
	print(model_df.head())
	param_sets=[]

	if "classic" in data_file:
		basic_params=classic_params

	for i in models:
		row=model_df.iloc[i]
		params=basic_params.copy()
		params["coldens"]=row["total_col_dens"]
		params["h2col"]=row["total_h2_col"]
		params["ccol"]=row["total_c_col"]
		params["rout"]=row["x"]
		params["initialTemp"]=row["temp"]
		params["radfield"]=1.7e5
		param_sets.append(params)

	result_df=pd.DataFrame({"Av":model_df.loc[models,"AV"],"Haworth Temp":model_df.loc[models,"temp"]})

	start=time.time()
	with get_context("spawn").Pool(2) as pool:
		result=pool.map(run_uclchem,param_sets)
		pool.close()
		pool.join()
	result=np.asarray(result)
	result_df["UCLCHEM Temp"]=result
	print(result)
	result_df.to_csv("Benchmarking/haworth/rollig_classic.csv",index=False)