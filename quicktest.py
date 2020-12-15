import pandas as pd
df=pd.read_csv("examples/test-output/phase1-full.dat",skiprows=2)
df.columns=df.columns.str.strip()
print(df[[x for x in df.columns if "#" in x]].sum(axis=1))