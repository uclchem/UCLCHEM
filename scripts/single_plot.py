import pandas as pd
import matplotlib.pyplot as plt

fig,[ax,ax2]=plt.subplots(2,1)
df=pd.read_csv("test/1e5_1e5.5_full.dat")
df=df.rename(lambda x: x.strip(),axis=1)
species=["CO","C+","C","N","E-","H","H+","H2","H2O","#CO","#H2O"]

ax.set(xscale="log",yscale="log")

for spec in species:
	ax.plot(df["Time"],df[spec],label=spec)
ax.legend()
ax.set(xscale="log",yscale="log",ylim=(1e-12,9.99e-1))
ax2.plot(df["Time"],df["gasTemp"])
ax2.set(xscale="log")
plt.show()
