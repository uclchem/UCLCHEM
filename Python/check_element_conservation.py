import uclchem
import numpy as np
input_file="output/test.csv"
df=uclchem.read_output_file(input_file)

print("Variance in total abundance as a fraction of initial elemental abundance.")
for element in ["H","C","N"]:
    total_element=uclchem.check_abunds(element,df)
    print(f"{element} {total_element.var()/total_element.values[0]:.2e}")

    sums=np.where(df.columns.str.contains(element),1,0)
    for i in range(2,10):
        sums+=np.where(df.columns.str.contains(element+f"{i:.0f}"),i-1,0)
    for variable in ['Time', 'Density', 'gasTemp', 'av', 'point',"SURFACE","BULK"]:
        sums=np.where(df.columns==variable,0,sums)
    print(df.mul(sums,axis=1).max(axis=0).sort_values(ascending=False).head(1))

