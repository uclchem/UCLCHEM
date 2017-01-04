!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Physical Conditions and Initial Abundances
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Initial physics variables and final  values. for temp, density and time
initialTemp=10.0;;maxTemp=300;initialDens=1.00d3;finalDens=1.00d7;t0=0.0;finalTime=1.00d6
!radfield in habing, cosmic ray ionisation rates as multiple of standard
radfield=1.0;zeta=1.0
fr=1.0;
!size of cloud set by inner and outer radii (rin and rout). used to calculate extinction.
!baseAv is extinction at cloud edge
!points is number of parcels to run model for. spaced  evenly between rin and rout
rout=0.05;rin=0;baseAv=2.0;points=1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Behavioural switches
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!switch (0/1) -> finish model run at (finalTime/finalDens)
switch=0

!collapse (0/1/2/3/4) 1/0 are on/off for standard free-fall. 2/3/4 are different collapse modes noted in cloud.f90
!collape=0/1 ONLY if not using cloud.f90
!In all cases collapse=1 lets chem.f90 know it should call  densdot in the physics module to get time derivative of density and include it in ODES
!Any other values tells it to use density value as set by physics module
collapse=1
!for collapse=1 can introduce factor bc to slow freefall
bc=1.0

!first chooses whether first run (So write final abudances) or second phase run (So read abudances from previous phase)
!so file 7 will be written to at end (first=0) or read from to initialise abundances (first=1)
!phase chooses behaviour. ie. heating in phase2 for cloud model.
!you may choose to run phase1 physics twice with the second run building from the first so first and phase are separated
first=1;phase=1;

!non-thermal Desorption. Turn it all on/off. Turn off h2, cosmic ray induced and uv induced off separately too
desorb=1;
h2desorb=1;crdesorb=1;uvcr=1;
!evap sets thermal desorption  (0/1/2) -> (none/temp dependent/ instantaneous)
evap=0;


!In phase 2, temp profile depends on mass of star
!Tempindx selects mass: 1=5Msol,2=10M,3=15M,4=25M,5=60M
tempindx=5

!ion sets ionisatoin fraction of carbon. See chem.f90:initialise
ion=2

!cshock specific variable, uncomment or comment as  needed
!vs=40.0

!initial fractional abundances of elements(from Asplund et al. 2009 ARAA table 1 -SOLAR)
fh=0.0;fhe = 0.085;fc  = 2.692d-04;fo  = 4.898d-04;fn  = 6.761d-05
fs  = 1.318d-05;fmg = 3.981d-05;fsi = 3.236d-05;fcl = 3.162d-07;
fp=2.57d-09 ; ff = 3.6d-08 !fp depleted 1/100 of solar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Input and output Files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!output species. Numbers are position of species in species.csv starting from 1.
!Due to line 1 being the number of species, the line number of a  species is 1 higher than it's index for this list. ie. H is the first species and is on line 2
!file 4 will print columnated time,dens,temp,abudances for all species listed in outindx every writestep timesteps.
!length of this array set by nout in chem.f90
outindx=(/73,260,262,220,219,274/);writestep=1

!open files for reading=writing
!output files
open(10,file='output-full.dat',status='unknown') !full output
open(11,file='output-column.dat',status='unknown')!columnated output based  on outindx
open(12,file='analysis.dat',status='unknown') !analysis file showing main reacction and formation routes of outindx species.

!input files
open(21,file='species.csv',status='old')         !species file
open(22,file='reactions.csv',status='old')       !reaction file
open(23,file='evaplists.csv',status='old')       !lists of species to evaporate in different thermal desorption events
open(7,file='startabund.dat',status='unknown')      !initialise abundance file. saved to at end or loaded from at start depending on first=(0/1)

open(79,file='debuglog',status='unknown')       !debug file.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!More complicated parameters that affect core code below. Do not alter without reading articles associated  !
!with each process. eg chemistry variables below are reference in rates.f90                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Chemistry variables
!Description and use found in rate.f90
!Desorption treatment is described in Roberts et al. 2007, these are rates and efficiencies of processes that cause desorption
ebmaxh2=1.21d3;epsilon=0.01;ebmaxcrf=1.21d3;uvcreff=1.0d-3
ebmaxcr=1.21d3;phi=1.0d5;ebmaxuvcr=1.0d4; uvy=0.1
omega=0.5;grainArea=2.4d-22;radg=1.d-5

!dopw = doppler width (in s-1) of a typical transition
!(assuming turbulent broadening with beta=3e5cms-1)
!radw = radiative line width of typ. transition (in s-1)
!fosc = oscillator strength of a typical transition
dopw=3.0e10;radw=8.0e07;xl=1000.0;fosc  = 1.0d-2

!DVODE SETTINGS
!You may wish to change abstol and reltol. Larger values run faster but can lose accuracy. In extreme cases, the model will crash.     
ISTATE=1;MF=22;ITOL=1;ITASK=1;IOPT=1;MESFLG=1
abstol=1e-20;reltol=1e-15;MXSTEP=10000

!Arrays for phase 2 temp profiles. Parameters for equation chosen by index
!arrays go [5 Msun, 10, 15, 25,60]
tempa=(/4.8560d-2,7.8470d-3,9.6966d-4,1.706d-4,4.74d-7/)
tempb=(/0.6255,0.8395,1.085,1.289,1.98/)
solidtemp=(/19.6,19.45,19.3,19.5,20.35/)
volctemp=(/86.3,88.2,89.5,90.4,92.2/)
codestemp=(/97.5,99.4,100.8,101.6,103.4/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CO and H2 self-shielding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
startr=.true.
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
