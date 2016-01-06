!Initial physics variables
temp=30.0;d0=1.00d6;dfin=1.00d4;t0=0.0;tfin=1.00d6
fr=1.0;radfield=1.0;zeta=1.000;avic=2.0
rout=0.3;rin=0;oldtemp=temp;bc=1.0;maxt=300
tempa=0.1927;tempb=0.5339;points=3

!Behavioural switches
!switch (0/1) -> (tfin/dfin)
!evap (0/1/2) -> (none/temp dependent/ instantaneous)
!other switches are on/off (1/0)
switch=0;collapse=0;first=1;desorb=1;startr=.true.
h2desorb=1;crdesorb=1;crdesorb2=1;uvcr=1;evap=0;ion=2
phase=1

!initial fractional abundances (from Asplund et al. 2009 ARAA table 1 -SOLAR)
fh=0.0;fhe = 0.085;fc  = 2.692d-04;fo  = 4.898d-04;fn  = 6.761d-05
fs  = 1.318d-05;fmg = 3.981d-05;fsi = 3.236d-05;fcl = 3.162d-07

!output species
outindx=(/1,2,3,4,5,6/);writestep=10

!open files for reading=writing
open(1,file='output',status='unknown')
open(2,file='reactions.csv',status='old')
open(3,file='species.csv',status='old')
open(4,file='outspecies',status='unknown')
open(7,file='end_step.d',status='unknown')
open(8,file='evaplists.csv',status='old')
open(78,file='debuglog',status='unknown')
!open(88,file='analysis',status='unknown')
open(79,file='jondebug',status='unknown')

!Chemistry variables
!dopw = doppler width (in s-1) of a typical transition
!(assuming turbulent broadening with beta=3e5cms-1)
!radw = radiative line width of typ. transition (in s-1)
!fosc = oscillator strength of a typical transition
ebmaxh2=1.21d3;epsilon=0.01;ebmaxcrf=1.21d3;uvcreff=1.0d-3
ebmaxcr=1.21d3;phi=1.0d5;ebmaxuvcr=1.0d4; uvy=0.1
omega=0.5;grain=1.1d-17;radg=1.d-5
dopw=3.0e10;radw=8.0e07;xl=1000.0;fosc  = 1.0d-2


!DLSODE SETTINGS        
RWORK=0.0;IWORK=0.0
ITOL=1;ITASK=1;ISTATE=1;IOPT=1;MESFLG=1
LUNIT=6;LRW=100000;LIW=500;MXSTEP=10000

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
dimco=7; dimh2=6
startr=.true.
         
