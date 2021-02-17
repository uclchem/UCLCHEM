import uclchem
import numpy as np
import time
from multiprocessing import Pool
import pandas as pd
from os import path,mkdir

def run_model(a):
	i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field,init_temp,dens,heating_flag,av_factor,zeta,metallicity=a	
	params={
		"columnFile":f"test/{model_type}_col.dat",
		"outputFile":f"test/{model_type}_full.dat",
		"initialTemp":init_temp,
		"initialDens":dens,
		"finalTime":2.0e6,
		"radfield":rad_field,
		"outSpecies":12,
		"h2col":h2_col,
		"ccol":c_col,
		"coldens":col_dens,
		"rout":cloud_size,
		"fc":1.0e-4,
		"fhe":0.1,
		"fh":0.4,
		"fo":3.0e-4,
		"fn":0.0,
		"fs":0.0,
		"fcl":0.0,
		"fsi":0.0,
		"fmg":5.00e-06,
		"ion":2,
		"zeta":zeta,
		"fr":0.0,
		"metallicity":metallicity,
		"heatingFlag":heating_flag,
		"avFactor":av_factor}	

	outSpecies="H2O,H2,H,H+,HE+,C,C+,O,O+,CO,CO+,E-"
	uclchem.general(params,outSpecies,f"test/","1")
	return 1


av_converts={"10_1e3":6.289E-22,
			"10_1e5.5":6.289E-22,
			"1e5_1e3":6.289E-22,
			"1e5_1e5.5":6.289E-22,
			"fixed_cooling":6.289E-22,
			"low_metallicity":6.289E-22,
			"10_linear_increase":6.289E-22,
			"10_linear_decrease":6.289E-22,
			"low_rad":1.6e-21,
			"low_rad_fixed":1.6e-21,
			"square":6.289E-22,
			"sine":6.289E-22,
			"high_cr":1.6e-21}

zetas={"10_1e3":3.84,
	"10_1e5.5":3.84,
	"1e5_1e3":3.84,
	"1e5_1e5.5":3.84,
	"low_metallicity":3.84,
	"fixed_cooling":3.84,
	"10_linear_increase":3.84,
	"10_linear_decrease":3.84,
	"low_rad":1.23,
	"low_rad_fixed":1.23,
	"square":3.84,
	"sine":3.84,
	"high_cr":123.0
}


start=time.time()


model_type="1e5_1e5.5"


fixed_temp=("low_rad_fixed" in model_type)
fixed_cooling=("fixed_cooling" in model_type)
fixed_heating=("fixed_heating" in model_type)
low_metallicity= (model_type=="low_metallicity")

model_df=pd.read_csv(f"Benchmarking/grid_inputs/{model_type}.csv")
i=0
row=model_df.iloc[i]
cloud_size=row["size"]/3.086e18
h2_col=row["total_h2_col"]
c_col=row["total_c_col"]
col_dens=row["total_col_dens"]
dens=row["n_H"]
rad_field=row["FUV"]*1.7

if fixed_temp:
	initialTemp=row["T_g"]
	heatingFlag=False

else:
	initialTemp=row["T_g"]
	heatingFlag=True

initialTemp=10.0
if fixed_cooling:
	cooling=row["totalCooling"]
else:
	cooling=0.0

if fixed_heating:
	heating=row["Total"]
else:
	heating=0.0


if low_metallicity:
	metallicity=0.2
	print(metallicity)
else:
	metallicity=1.0

# av conversion factor CHANGED TO MATCH UCL_PDR, should find most agreed value after benchmarking
av=col_dens*av_converts[model_type]
zeta=zetas[model_type]
success=run_model([i,cloud_size,h2_col,c_col,col_dens,model_type,rad_field,initialTemp,dens,heatingFlag,av_converts[model_type],zeta,metallicity])
end=time.time()

print(f"success flag {success} in {end-start} seconds")