!!This file is used to autogenerate the docs. So please ignore the mess!
!!Double !! lines do not show up in docs, single ones do. 
!!If you add a parameter, please take the time to add a useful descriptor comment on the same line
!!and then re-run utils/generate_param_docs.py to update the docs.
!!note the resuting md file needs manually adding to the website.
MODULE DEFAULTPARAMETERS
USE constants
!---  
!id: parameters
!title: Model Parameters
!---
!UCLCHEM will default to these values unless they are overridden by user. Users can override these by adding the variable name as written here in the param_dict argument of any UCLCHEM model function. param_dict is not case sensitive.
!
!## Physical Variables
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
IMPLICIT NONE
REAL(dp) :: initialTemp=10.0 !Initial gas temperature in Kelvin for all gas parcels in model.
REAL(dp) :: initialDens=1.00d2 !Initial gas density in H nuclei per cm$^{-3}$ for all gas parcels in model.
REAL(dp) :: finalDens=1.00d5 !Final gas density achieved via freefall. 
REAL(dp) :: currentTime=0.0 !Time at start of model in years.
REAL(dp) :: finalTime=5.0d6 !Time to stop model in years, if not using `endAtFinalDensity` below.
REAL(dp) :: radfield=1.0 !Interstellar radiation field in Habing
REAL(dp) :: zeta=1.0 !Cosmic ray ionisation rate as multiple of $1.3 10^{-17} s^{-1}$
REAL(dp) :: rout=0.05 !Outer radius of cloud being modelled in pc.
REAL(dp) :: rin=0.0 !Minimum radial distance from cloud centre to consider.
REAL(dp) :: baseAv=2.0 !Extinction at cloud edge, Av of a parcel at rout.
INTEGER(dp) :: points=1 !Number of gas parcels equally spaced between rin to rout to consider
REAL(dp) :: bm0=1.0 !magnetic parameter [microgauss]: B0 = bm0*sqrt(initialDens)
!
!## Behavioural Controls
!*The following parameters generally turn on or off features of the model. If a parameter is set to `True`, then it is turned on. If it is set to `False`, then it is turned off.*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
REAL(dp) :: freezeFactor=1.0 !Modify freeze out rate of gas parcels by this factor.
LOGICAL :: endAtFinalDensity=.False. !Choose to end model at final density, otherwise end at final time.
LOGICAL :: freefall=.False. !Controls whether models density increaes following freefall equation.
REAL(dp) :: freefallFactor=1.0 !Modify freefall rate by factor, usually to slow it.
LOGICAL :: desorb=.True. !Toggles all non-thermal desoprtion processes on or off.
LOGICAL :: h2desorb=.True. !Individually toggle non-thermal desorption due to H2 formation.
LOGICAL :: crdesorb=.True. !Individually toggle non-thermal desorption due to cosmic rays.
LOGICAL :: uvdesorb=.True. !Individually toggle non-thermal desorption due to uv photons.
LOGICAL :: thermdesorb=.True. !Toggle continuous thermal desorption.
LOGICAL :: instantSublimation=.False. !Toggle instantaneous sublimation of the ices at t=0
LOGICAL :: cosmicRayAttenuation=.False. !Use column density to attenuate cosmic ray ionisation rate following [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
CHARACTER :: ionModel='L' !L/H model for cosmic ray attenuation [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
LOGICAL :: improvedH2CRPDissociation=.False. !Use H2 CRP dissociation rate from [Padovani et al. 2018b](https://arxiv.org/abs/1809.04168).
LOGICAL :: heatingFlag=.false. !If True, heating is applied to the gas parcels.
LOGICAL :: enforceChargeConservation = .false. ! Enforce the chrage by keeping track of charged ions.
!
!## Input and Output
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
CHARACTER(256) :: outputFile="" !File to write full output of UCLCHEM. This includes physical parameter values and all abundances at every time step.
CHARACTER(256) :: columnFile="" !File to write specific species abundances, see outSpecies.
CHARACTER(256) :: rateConstantFile="" !File to write rate 'constants' at each timestep. This includes physical parameter values.
CHARACTER(256) :: ratesFile="" !File to write reaction rates (flux) at each timestep. This includes physical parameter values.
CHARACTER(256) :: heatingFile="" !File to write heating and cooling rates at each timestep.
INTEGER :: writeStep=1 !Writing to columnFile only happens every writeStep timesteps.
CHARACTER(256) :: abundSaveFile="" ! The file to save the abundances to at the end of the model.
CHARACTER(256) :: abundLoadFile="" ! The file to load the abundances from at the start of the model.
CHARACTER(256) :: coolantDataDir="" ! Directory where the collisional rate files are stored.
!|abundSaveFile |None| File to store final abundances at the end of the model so future models can use them as the initial abundances. If not provided, no file will be produced.
!|abundLoadFile |None| File from which to load initial abundances for the model, created through `abundSaveFile`. If not provided, the model starts from elemental gas.
!|outSpecies|None| A space separated list of species to output to columnFile. Supplied as a separate list argument to most python functions, see python API docs.
!
!## Initial Abundances
!*Unless otherwise specified, we take all abundances from Jenkins et al. 2009, using the heavily depleted case from Table 4.*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
REAL(dp) :: metallicity=1.0 !Scale the abundances of all elements heavier than He by this factor.
INTEGER(dp) :: ion=2 !Sets how much elemental C is initially atomic (0= all atomic/1=50:50/2=fully ionized).
REAL(dp) :: fh=0.5 !Total elemental abundance of H is always 1 by definition because abundances are relative to number of H nuclei. Use fh to set how much to initially put in atomic H, the rest goes to H2.
REAL(dp) :: fhe = 0.1 !Total elemental abundance of He.
REAL(dp) :: fc=1.77d-04 !Total elemental abundance of C.
REAL(dp) :: fo  = 3.34d-04 !Total elemental abundance of O.
REAL(dp) :: fn  = 6.18d-05 !Total elemental abundance of N.
REAL(dp) :: fs  = 3.51d-6 !Total elemental abundance of S.
REAL(dp) :: fmg = 2.256d-06 !Total elemental abundance of Mg.
REAL(dp) :: fsi = 1.78d-06 !Total elemental abundance of Si.
REAL(dp) :: fcl = 3.39d-08 !Total elemental abundance of Cl.
REAL(dp) :: fp =7.78d-08 !Total elemental abundance of P.
REAL(dp) :: ffe =2.01d-7!Total elemental abundance of Fe.
REAL(dp) :: ff = 3.6d-08 !fp depleted 1/100 of solar from Asplund 2009.
REAL(dp) :: fd=0.0 ! The following elements are not typically used. We do not recommend any particular value.
REAL(dp) :: fli=0.0 !Total elemental abundance of Li.
REAL(dp) :: fna=0.0 !Total elemental abundance of Na.
REAL(dp) :: fpah=0.0 !Total initial abundance of PAHs.
REAL(dp) :: f15n=0.0 !Total initial abundance of 15N.
REAL(dp) :: f13c=0.0 !Total initial abundance of 13C.
REAL(dp) :: f18O=0.0 !Total initial abundance of 18O.
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
REAL(dp) :: reltol=1d-8 !Relative tolerance for integration, see [integration docs](/docs/trouble-integration) for advice.
REAL(dp) :: abstol_factor=1.0d-14 !Absolute tolerance for integration is calculated by multiplying species abundance by this factor.
REAL(dp) :: abstol_min=1.0d-25 !Minimum value absolute tolerances can take.
INTEGER :: MXSTEP=10000 !Maximum steps allowed in integration before warning is thrown. ! HAS TO BE INT4 instead of INT8
!
!## Here be Dragons
!*These are not recommended to be changed unless you know what you are doing*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
REAL(dp) :: ebmaxh2=1.21d3 ! Maximum binding energy of species desorbed by H2 formation.
REAL(dp) :: ebmaxcr=1.21d3 ! Maximum binding energy of species desorbed by cosmic ray ionisation.
REAL(dp) :: ebmaxuvcr=1.0d4 ! Maximum binding energy of species desorbed by UV photons.
REAL(dp) :: epsilon=0.01 !Number of molecules desorbed per H2 formation.
REAL(dp) :: uv_yield=0.03 !Number of molecules desorbed per UV photon. The yield is extrapolated from Oberg et al. 2009
REAL(dp) :: phi=1.0d5 !Number of molecules desorbed per cosmic ray ionisation.
REAL(dp) :: uvcreff=1.0d-3 !Ratio of CR induced UV photons to ISRF UV photons.
REAL(dp) :: omega=0.5 !Dust grain albedo.
!|alpha|{1:0.0,2:0.0}| Set alpha coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how alpha is used for each reaction type.|
!|beta|{1:0.0,2:0.0}| Set beta coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how beta is used for each reaction type.|
!|gama|{1:0.0,2:0.0}| Set gama coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how gama is used for each reaction type.|
CONTAINS
! Add a dummy subroutine to help f2py compile: https://github.com/numpy/numpy/issues/27167
SUBROUTINE DUMMY_TWO(dummy_two_output)
        integer, intent(out) :: dummy_two_output
        dummy_two_output = 2
    END SUBROUTINE DUMMY_TWO
END MODULE DEFAULTPARAMETERS