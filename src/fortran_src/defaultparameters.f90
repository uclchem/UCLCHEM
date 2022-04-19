!!This file is used to autogenerate the docs. So please ignore the mess!
!!Double !! lines do not show up in docs, single ones do. 
!!If you add a parameter, please take the time to add a useful descriptor comment on the same line
!!and then re-run utils/generate_param_docs.py to update the docs.
!!note the resuting md file needs manually adding to the website.


!---
!id: parameters
!title: Model Parameters
!---
!UCLCHEM will default to these values unless they are overridden by user. Users can override these by adding the variable name as written here in the param_dict argument of any UCLCHEM model function. param_dict is not case sensitive.
!
!## Physical Variables
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
initialTemp=10.0 !Initial gas temperature for all gas parcels in model
initialDens=1.00d2 !Initial gas density for all gas parcels in model
finalDens=1.00d5 !Final gas density achieved via freefall
currentTime=0.0 !Time at start of model
finalTime=5.0d6 !Time at end of model, see switch below
radfield=1.0 !Interstellar radiatin field in Habing
zeta=1.0 !Cosmic ray ionisation rate as multiple of 1.3e-17 s-1
rout=0.05 !Inner radius of cloud being modelled
rin=0.0 !Minimum radial distance from cloud centre to consider
baseAv=2.0 !Extinction at cloud edge, Av of a parcel at rout
points=1 !Number of gas parcels equally spaced between rin to rout to consider
!
!## Behavioural Controls
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
fr=1.0 !Modify freeze out rate of gas parcels by factor fr
endAtFinalDensity=.False. !Choose to end model at final density, otherwise end at final time
freefall=.False. !Controls whether models density increaes following freefall equation.
bc=1.0 !Modify freefall rate by factor bc
desorb=.True. !Toggles all non-thermal desoprtion processes on or off
h2desorb=.True. !Individually turn on and off non-thermal desorption due to H2 formation
crdesorb=.True. !Individually turn on and off non-thermal desorption due to cosmic rays
uvdesorb=.True. !Individually turn on and off non-thermal desorption due to uv photons
thermdesorb=.True. !continuous thermal desorption - best have on for three phase and off for two phase
instantSublimation=.False. !Instantaneous sublimation of gas parcels at t=0
cosmicRayAttentuation=.False. !Use column density to attenuate cosmic ray ionisation rate following [Padovani et al. 2018](https://arxiv.org/abs/1803.09348)
ionModel='L' !L/H model for cosmic ray attenuation [Padovani et al. 2018](https://arxiv.org/abs/1803.09348)
improvedH2CRPDissociation=.False. !Use h2 CRP dissociation rate from [Padovani et al. 2018b](https://arxiv.org/abs/1809.04168)
!
!## Input and Output
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
outputFile="output/full.dat" !File to write full output of UCLCHEM
columnFile="output/column.dat" !File to write specific species abundances, see outSpecies
writeStep=1 !Writing to column file only happens every writeStep timesteps
abundSaveFile="output/start.dat" !File to store final abundances at the end of the model so future models can use them as the initial abundances. If not provided, no file will be produced.
!|abundLoadFile |None| File from which to load initial abundances for the model, created through `abundSaveFile`. If not provided, the model starts from elemental gas.
!|outSpecies|None| A space separated list of species to output to columnFile. Supplied as a separate list argument to most python functions
!
!### Initial Abundances
!*Unless otherwise specified, we take all abundances from Jenkins et al. 2009, using the heavily depleted case from Table 4.*
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
metallicity=1.0 !Scale the abundances of all elements heavier than He by this factor
ion=2 !Sets how much elemental C is initially atomic (0= all atomic/1=50:50/2=fully ionized).
fh=0.5 !Total abundance of H is always 1 but set how much to put in atomic H here, rest is H2
fhe = 0.1 !Total elemental abundance of He
fc=1.77d-04 !Total elemental abundance of C
fo  = 3.34d-04 !Total elemental abundance of O
fn  = 6.18d-05 !Total elemental abundance of N
fs  = 3.51d-6 !Total elemental abundance of S
fmg = 2.256d-07 !Total elemental abundance of Mg
fsi = 1.78d-06 !Total elemental abundance of Si
fcl = 3.39d-08 !Total elemental abundance of Cl
fp=7.78d-08 !Total elemental abundance of P
ffe=2.01d-7!Total elemental abundance of Fe
ff = 3.6d-08 !fp depleted 1/100 of solar from Asplund 2009
fd=0.0 ! The following elements are not typically used. We do not recommend any particular value
fli=0.0 !Total elemental abundance of Li
fna=0.0 !Total elemental abundance of Na
fpah=0.0 !Total initial abundance of PAHs
f15n=0.0 !Total initial abundance of 15N
f13c=0.0 !Total initial abundance of 13C
f18O=0.0 !Total initial abundance of 18O
!!
!! We used to use Asplund et al. 2009,kept here for reference
!! !initial fractional abundances of elements(from Asplund et al. 2009 ARAA table 1 -SOLAR)
!! !note fh is fraction of H initially in H atoms. Total H is always 1.
!! !fh=0.5;fhe = 0.1;fc  = 2.6d-04;fo  = 4.6d-04;fn  = 6.1d-05
!! fs  = 1.318d-05;fmg = 3.981d-05;fsi = 1.0d-07;fcl = 3.162d-07;
!! fp=2.57d-09 ; ff = 3.6d-08 !fp depleted 1/100 of solar
!
!## Integration Controls
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
reltol=1d-8 !Relative tolerance for integration, see [integration docs](/docs/integration) for advice
abstol_factor=1.0d-14 !Absolute tolerance for integration is calculated by multiplying species abundance by this factor
abstol_min=1.0d-25 !Minimum value absolute tolerances can take
MXSTEP=10000 !Maximum steps allowed in integration before warning is thrown
!
!## Here be Dragons
!*These are not recommended to be changed unless you know what you are doing*
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
ebmaxh2=1.21d3 !
epsilon=0.01 !
ebmaxcrf=1.21d3 !
uvcreff=1.0d-3 !
ebmaxcr=1.21d3 !
phi=1.0d5 !
ebmaxuvcr=1.0d4 !
uv_yield=0.1 !
omega=0.5 !