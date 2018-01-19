!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Physical Conditions and Initial Abundances
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Initial physics variables and final  values. for temp, density and time
initialTemp=10.0;maxTemp=300;initialDens=1.00d4;finalDens=1.00d4;currentTime=0.0;finalTime=6.00d6
!radfield in habing, cosmic ray ionisation rates as multiple of standard
radfield=1.0;zeta=1.0
fr=1.0;
!Size of cloud set by inner and outer radii (rin and rout). used to calculate extinction.
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
collapse=0
!for collapse=1 can introduce factor bc to slow freefall
bc=1.0

!First chooses whether first run (So write final abudances) or second phase run (So read abudances from previous phase)
!First=1 starts from elemental abundances and writes a file (file 7) at the end, first=0 reads a file (file 7) to get initial abundances
!phase chooses behaviour. ie. heating in phase2 for cloud models
!you may choose to run phase1 physics twice with the second run building from the first so first and phase are separated
first=1;phase=1;

!non-thermal Desorption. Turn it all on/off. Turn off h2, cosmic ray induced and uv induced off separately too
desorb=1;
h2desorb=1;crdesorb=1;uvcr=1;
!evap sets thermal desorption  (0/1/2) -> (none/temp dependent/ instantaneous)
evap=0;

!ion sets ionisatoin fraction of carbon. See chem.f90:initialise
ion=2


!cloud module specific variable for phase 2, temp profile depends on mass of star
!Tempindx selects mass: 1=1Msol,2=5,3=10M,4=15M,5=25M,6=60M
tempindx=5


!cshock module specific variable, uncomment or comment as  needed
!vs=40.0

!initial fractional abundances of elements(from Asplund et al. 2009 ARAA table 1 -SOLAR)
fh=0.0;fhe = 0.1;fc  = 2.6d-04;fo  = 4.6d-04;fn  = 6.1d-05
fs  = 1.0d-07;fmg = 3.981d-05;fsi = 1.0d-07;fcl = 3.162d-07;
fp=2.57d-09 ; ff = 3.6d-08 !fp depleted 1/100 of solar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Input and output Files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!A full output of abundances is written by default. Additionally, name species here for 
!a columnated output of time,density,temperature and abundances of those species
!writeStep sets how often this is written out. Columns written every n steps for writeStep=n.
!Fortran will reject this array if species with shorter names are not padded with spaaces at the end.
outSpecies=(/'CO ','H2S','OCS'/);writeStep=1

!open files for reading=writing
!output files
open(10,file='output/fullexplode.dat',status='unknown') !full output
open(11,file='output/columnexplode.dat',status='unknown')!columnated output based  on outindx

!input files
open(21,file='species.csv',status='old')         !species file
open(22,file='reactions.csv',status='old')       !reaction file
open(23,file='evaplists.csv',status='old')       !lists of species to evaporate in different thermal desorption events
open(7,file='output/startabund.dat',status='unknown')      !initialise abundance file. saved to at end or loaded from at start depending on first=(0/1)

open(79,file='output/debuglog',status='unknown')       !debug file.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!More complicated parameters that affect core code below. Do not alter without reading articles associated  !
!with each process. eg chemistry variables below are reference in rates.f90                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Chemistry variables
!Description and use found in rate.f90
!Desorption treatment is described in Roberts et al. 2007, these are rates and efficiencies of processes that cause desorption
!Careful, grainArea is actually grain surface area per hydrogen atom. grainRadius is actual radius of grains in cm.
ebmaxh2=1.21d3;epsilon=0.01;ebmaxcrf=1.21d3;uvcreff=1.0d-3
ebmaxcr=1.21d3;phi=1.0d5;ebmaxuvcr=1.0d4; uvy=0.1
omega=0.5;
!dopw = doppler width (in s-1) of a typical transition
!(assuming turbulent broadening with beta=3e5cms-1)
!radw = radiative line width of typ. transition (in s-1)
!fosc = oscillator strength of a typical transition
dopw=3.0e10;radw=8.0e07;xl=1000.0;fosc  = 1.0d-2