import numpy as np
import matplotlib.pyplot as plt
from seaborn import color_palette
masses=np.asarray([1,5, 10, 15, 25,60])
tempa=np.log10([1.927e-1,4.8560e-2,7.8470e-3,9.6966e-4,1.706e-4,4.74e-7])
tempb=[0.5339,0.6255,0.8395,1.085,1.289,1.98]
solidtemp=[20.0,19.6,19.45,19.3,19.5,20.35]
volctemp=[84.0,86.3,88.2,89.5,90.4,92.2]
codestemp=[95.0,97.5,99.4,100.8,101.6,103.4]

fig,ax=plt.subplots()
colors=color_palette(n_colors=5)
for i,y in enumerate([tempa,tempb,solidtemp,volctemp,codestemp]):
    poly_fit=np.polynomial.polynomial.Polynomial.fit(masses,y,3)
    poly_line=poly_fit(masses)
    print(y,np.sum(np.abs(poly_line-y)))
    ax.plot(masses,y,label=y,marker="*",ls="",color=colors[i])
    ax.plot(masses,poly_line,color=colors[i])
plt.show()