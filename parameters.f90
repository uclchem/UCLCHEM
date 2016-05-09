!Initial physics variables
inittemp=10.0;initdens=1.00d2;dfin=1.00d6;t0=0.0;tfin=1.00d7
fr=1.0;radfield=1.0;zeta=1.00;avic=2.0
rout=0.05;rin=0;bc=1.0;maxtemp=300
points=1

!Behavioural switches
!switch (0/1) -> (tfin/dfin)
!evap (0/1/2) -> (none/temp dependent/ instantaneous)
!In phase 2, temp profile depends on mass of star
!Tempindx selects mass: 1=5Msol,2=10M,3=15M,4=25M,5=60M
!other switches are on/off (1/0)
switch=0;collapse=1;first=1;desorb=1;startr=.true.
h2desorb=1;crdesorb=1;crdesorb2=1;uvcr=1;evap=0;ion=2
phase=1;tempindx=1

!initial fractional abundances (from Asplund et al. 2009 ARAA table 1 -SOLAR)
fh=0.0;fhe = 0.085;fc  = 2.692d-04;fo  = 4.898d-04;fn  = 6.761d-05
fs  = 1.318d-05;fmg = 3.981d-05;fsi = 3.236d-05;fcl = 3.162d-07

!output species
outindx=(/66,35,61,30/);writestep=10

!DVODE SETTINGS        
ISTATE=1;MF=22;ITOL=1;ITASK=1;IOPT=1;MESFLG=1
abstol=1e-25;reltol=1e-7;MXSTEP=10000

!open files for reading=writing
open(1,file='odeoutput',status='unknown')
open(2,file='reactions.csv',status='old')
open(3,file='species.csv',status='old')
open(4,file='testoutput',status='unknown')
open(7,file='startabund',status='unknown')
open(8,file='evaplists.csv',status='old')
open(79,file='debuglog',status='unknown')
open(88,file='analysis',status='unknown')

!Chemistry variables
!dopw = doppler width (in s-1) of a typical transition
!(assuming turbulent broadening with beta=3e5cms-1)
!radw = radiative line width of typ. transition (in s-1)
!fosc = oscillator strength of a typical transition
ebmaxh2=1.21d3;epsilon=0.01;ebmaxcrf=1.21d3;uvcreff=1.0d-3
ebmaxcr=1.21d3;phi=1.0d5;ebmaxuvcr=1.0d4; uvy=0.1
omega=0.5;grain=1.1d-17;radg=1.d-5
dopw=3.0e10;radw=8.0e07;xl=1000.0;fosc  = 1.0d-2

!CO self-shielding
corates =reshape((/0.000d+00, -1.408d-02, -1.099d-01, -4.400d-01,&
     &  -1.154d+00, -1.888d+00, -2.760d+00,&
     &  -8.539d-02, -1.015d-01, -2.104d-01, -5.608d-01,&
     &  -1.272d+00, -1.973d+00, -2.818d+00,&
     &  -1.451d-01, -1.612d-01, -2.708d-01, -6.273d-01,&
     &  -1.355d+00, -2.057d+00, -2.902d+00,&
     &  -4.559d-01, -4.666d-01, -5.432d-01, -8.665d-01,&
     &  -1.602d+00, -2.303d+00, -3.146d+00,&
     &  -1.303d+00, -1.312d+00, -1.367d+00, -1.676d+00,&
     &  -2.305d+00, -3.034d+00, -3.758d+00,&
     &  -3.883d+00, -3.888d+00, -3.936d+00, -4.197d+00,&
     &  -4.739d+00, -5.165d+00, -5.441d+00 /),shape(corates))

ncogr =(/12.0d+00, 13.0d+00, 14.0d+00, 15.0d+00,&
      &16.0d+00, 17.0d+00, 18.0d+00 /)
nh2gr=(/18.0d+00, 19.0d+00, 20.0d+00, 21.0d+00,&
       &22.0d+00, 23.0d+00 /)

!Arrays for phase 2 temp profiles. Parameters for equation chosen by index
!arrays go [5 Msun, 10, 15, 25,60]
tempa=(/4.8560d-2,7.8470d-3,9.6966d-4,1.706d-4,4.74d-7/)
tempb=(/0.6255,0.8395,1.085,1.289,1.98/)
solidtemp=(/19.6,19.45,19.3,19.5,20.35/)
volctemp=(/86.3,88.2,89.5,90.4,92.2/)
codestemp=(/97.5,99.4,100.8,101.6,103.4/)

dimco=7; dimh2=6
startr=.true.
         
