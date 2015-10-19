import numpy as np
import matplotlib.pyplot as plt
import csv

#function to read the output of UCL_CHEM, give it the output file
#and 5 species names as strings.
def read_uclchem(filename,s1,s2,s3,s4,s5):
    a=open(filename).read()
    a=a.split('\n')
    #there are 68 lines per time step in the UCL_CHEM output file.
    lines=68
    timesteps=len(a)/lines

    time=[]
    dens=[]
    spec1=[]
    spec2=[]
    spec3=[]
    spec4=[]
    spec5=[]
    #so now do the following until end of file
    for i in np.arange(0,timesteps):
        j=i*lines
        b=a[j].split()
        b=float(b[-2].replace("D", "E"))
        time.append(b)
        b=a[j+1].split()
        b=float(b[-2].replace("D", "E"))
        dens.append(b)
        for k in np.arange(j+11,j+1+lines):
            b=a[k].split()
            for element in b:
                if element==s1:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec1.append(b)
                elif element==s2:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec2.append(b)
                elif element==s3:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec3.append(b)
                elif element==s4:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec4.append(b)      
                elif element==s5:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec5.append(b)
    time=np.log10(time)
    dens=np.log10(dens)
    spec1=np.log10(spec1)
    spec2=np.log10(spec1)
    spec3=np.log10(spec1)
    spec4=np.log10(spec1)
    spec5=np.log10(spec1)
    return time,dens,spec1,spec2,spec3,spec4,spec5

def read_olduclchem(filename,s1,s2,s3,s4,s5):
    a=open(filename).read()
    a=a.split('\n')
    #there are 68 lines per time step in the UCL_CHEM output file.
    lines=72
    timesteps=(len(a)-1)/lines
    time=[]
    dens=[]
    spec1=[]
    spec2=[]
    spec3=[]
    spec4=[]
    spec5=[]
    #so now do the following until end of file
    for i in np.arange(0,timesteps):
        j=(i*lines)+1
        b=a[j].split()
        b=float(b[-2].replace("D", "E"))
        dens.append(b)
        b=a[j+11].split()
        b=float(b[-4].replace("D", "E"))
        time.append(b)
        for k in np.arange(j+14,j+lines+1):
            b=a[k].split()
            for element in b:
                if element==s1:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec1.append(b)
                elif element==s2:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec2.append(b)
                elif element==s3:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec3.append(b)
                elif element==s4:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec4.append(b)      
                elif element==s5:
                    pos=b.index(element)
                    b=float(b[pos+2].replace("D", "E"))
                    spec5.append(b)
    time=np.log10(time)
    dens=np.log10(dens)
    spec1=np.log10(spec1)
    spec2=np.log10(spec1)
    spec3=np.log10(spec1)
    spec4=np.log10(spec1)
    spec5=np.log10(spec1)
    return time,dens,spec1,spec2,spec3,spec4,spec5


def write_plot(filename,time,dens,s1,s2,s3,s4,s5):
    outfile=open(filename,'wb')
    writer = csv.writer(outfile,delimiter=' ')
    n=len(time)
    for i in range(n):
        writer.writerow([time[i],dens[i],s1[i],s2[i],s3[i],s4[i],s5[i]])











def gen_plot(x,y1,y2,y3,y4,y5,outfile):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(x,y1,ls='-',color='black')
    # ax.text(0.9,0.9,'CO',transform=ax.transAxes)
    ax.plot(x,y2,ls='-',color='purple')
    # ax.text(0.9,0.85,'HCO+',color='purple',transform=ax.transAxes)
    ax.plot(x,y3,ls='-',color='red')
    # ax.text(0.9,0.8,'CS',color='red',transform=ax.transAxes)
    ax.plot(x,y4,ls='-',color='green')
    # ax.text(0.9,0.75,'NH$_3$',color='green',transform=ax.transAxes)
    ax.plot(x,y5,ls='-',color='brown')
    # ax.text(0.9,0.7,'N$_2$H+',color='brown',transform=ax.transAxes)
    ax.set_xlabel("log(time / yr)")
    ax.set_ylabel("log(X$_{species}$)")
    ax.minorticks_on()
    plt.savefig(outfile,bbox_inches='tight',pad_inches=0.1,dpi=300)
    plt.close()
