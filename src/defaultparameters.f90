!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Physical Conditions and Initial Abundances
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Initial physics variables and final  values. for temp, density and time
initialTemp=10.0;maxTemp=300
initialDens=1.00d2;finalDens=1.00d6
currentTime=0.0;finalTime=1.00d7

!radfield in habing, cosmic ray ionisation rates as multiple of standard
radfield=1.0;zeta=1.0

!Scale freeze out efficiency by an aribtrary value. 
fr=1.0;

!Size of cloud set by inner and outer radii (rin and rout). used to calculate extinction.
!baseAv is extinction at cloud edge
!points is number of parcels to run model for. spaced  evenly between rin and rout
rout=0.051847;rin=0;baseAv=0.0;points=1


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

!If you want to read abundances from a file (abundFile) set to 1
!If you set to zero, final abundances are written to abundFile
readAbunds=0
!phase chooses behaviour. Basic clouds in phase 1. Hot core /cshock in phase 2 (depending on physics module)
phase=1;

!non-thermal Desorption. Turn it all on/off. Turn off h2, cosmic ray induced and uv induced off separately too
desorb=1;
h2desorb=1;crdesorb=1;uvcr=1;

!evap sets thermal desorption  (0/1/2) -> (none/temp dependent/ instantaneous)
evap=0

!ion sets ionization fraction of carbon. See chem.f90:initialise
ion=2


!cloud module specific variable for phase 2, temp profile depends on mass of star
!Tempindx selects mass: 1=1Msol,2=5,3=10M,4=15M,5=25M,6=60M
tempindx=3


!cshock module specific variable, uncomment or comment as  needed
vs=20.0

!initial fractional abundances of elements(from Asplund et al. 2009 ARAA table 1 -SOLAR)
fh=0.0;fhe = 0.1;fc  = 2.6d-04;fo  = 4.6d-04;fn  = 6.1d-05
fs  = 1.318d-05;fmg = 3.981d-05;fsi = 1.0d-07;fcl = 3.162d-07;
fp=2.57d-09 ; ff = 3.6d-08 !fp depleted 1/100 of solar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Input and output Files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!A full output of abundances is written by default. Additionally, name species here for 
!a columnated output of time,density,temperature and abundances of those species
!Fortran will reject this array if species with shorter names are not padded with spaces at the end.
!array commented so it does not override array in input file. 
!If no array is passed in input, no column file is written
!outSpecies=(/'CO ','H2S','OCS','CS '/)

!writeStep sets how often columns written out. Columns written every n steps for writeStep=n.
writeStep=1

!If readAbund=1, starting abundances are read from abundFile
!If readAbund=0, final abundances are written to abundFile
abundFile="output/startcloud.dat"
!Full output written to outputFlie
outputFile="output/fullcloud.dat"
!columnated output of time,dens,temp and outSpecies written to column file
columnFile='output/columncloud.dat'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!More complicated parameters that affect core code below. Do not alter without reading articles associated  !
!with each process. eg chemistry variables below are reference in rates.f90                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Chemistry variables
!Description and use found in rate.f90
!Desorption treatment is described in Roberts et al. 2007, these are rates and efficiencies of processes that cause desorption
ebmaxh2=1.21d3;epsilon=0.01;ebmaxcrf=1.21d3;uvcreff=1.0d-3
ebmaxcr=1.21d3;phi=1.0d5;ebmaxuvcr=1.0d4; uvy=0.1
omega=0.5;
