!!This file is used to autogenerate the docs. So please ignore the mess!
!!Double !! lines do not show up in docs, single ones do.
!!If you add a parameter, please take the time to add a useful descriptor comment on the same line
!!and then re-run utils/generate_param_docs.py to update the docs.
!!note the resulting md file needs manually adding to the website.
module DEFAULTPARAMETERS
    use constants, only: dp
    !---
    !id: parameters
    !title: Model Parameters
    !---
    !UCLCHEM will default to these values unless they are overridden by user. Users can override these by adding the variable name as written here in the param_dict argument of any UCLCHEM model function. param_dict is not case sensitive.
    !
    !## Physical Variables
    !|Parameter|Default Value |Description|
    !| ----- | ------| ------ |
    implicit none

    public

    real(dp) :: initialTemp=10.0  !Initial gas temperature in Kelvin for all gas parcels in model.
    real(dp) :: initialDens=1.00e2_dp  !Initial gas density in H nuclei per cm$^{-3}$ for all gas parcels in model.
    real(dp) :: finalDens=1.00e5_dp  !Final gas density achieved via freefall.
    real(dp) :: currentTime=0.0  !Time at start of model in years (matches finalTime units).
    real(dp) :: finalTime=5.0e6_dp  !Time to stop model in years, if not using `endAtFinalDensity` below.
    real(dp) :: radfield=1.0  !Interstellar radiation field in Habing
    real(dp) :: zeta=1.0  !Cosmic ray ionization rate as multiple of $1.3 10^{-17} s^{-1}$
    real(dp) :: rout=0.05_dp  !Outer radius of cloud being modeled in pc.
    real(dp) :: rin=0.0  !Minimum radial distance from cloud center to consider.
    real(dp) :: baseAv=2.0  !Extinction at cloud edge, Av of a parcel at rout.
    integer :: points=1  !Number of gas parcels equally spaced between rin to rout to consider
    real(dp) :: bm0=1.0  !magnetic parameter [microgauss]: B0 = bm0*sqrt(initialDens)
    !Physical profiles for 1D model with pre-described gas density
    real(dp) :: density_scale_radius=0.05_dp  !unit of pc, distance below which the gas volume density is constant, and above which the gas density drops as n ~ r^{-a}
    real(dp) :: density_power_index=2.0  !Power-law index for density profile: n(r) = n0/(1 + (r/density_scale_radius)^density_power_index)
    !Luminosity source for hotcore in 1D model
    real(dp) :: lum_star=1.00e6_dp  !unit of Lsun, bolometric luminosity of the central source
    real(dp) :: temp_star=4.50e4_dp  !unit of K, temperature of the central source
    !
    !## Behavioral Controls
    !*The following parameters generally turn on or off features of the model. If a parameter is set to `True`, then it is turned on. If it is set to `False`, then it is turned off.*
    !
    !|Parameter|Default Value |Description|
    !| ----- | ------| ------ |
    real(dp) :: freezeFactor=1.0  !Modify freeze out rate of gas parcels by this factor.
    logical :: endAtFinalDensity=.false.  !Choose to end model at final density, otherwise end at final time.
    logical :: freefall=.false.  !Controls whether models density increases following freefall equation.
    real(dp) :: freefallFactor=1.0  !Modify freefall rate by factor, usually to slow it.
    logical :: desorb=.true.  !Toggles all non-thermal desoprtion processes on or off.
    logical :: h2desorb=.false.  !Individually toggle non-thermal desorption due to H2 formation.
    logical :: crdesorb=.true.  !Individually toggle non-thermal desorption due to cosmic rays.
    logical :: uvdesorb=.true.  !Individually toggle non-thermal desorption due to uv photons.
    logical :: chemdesorb=.true.  !Individually toggle non-thermal desorption due to chemical reactions.
    logical :: thermdesorb=.true.  !Toggle continuous thermal desorption.

    logical :: instantSublimation=.false.  !Toggle instantaneous sublimation of the ices at t=0
    logical :: cosmicRayAttenuation=.false.  !Use column density to attenuate cosmic ray ionization rate following [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
    character :: ionModel="L"  !L/H model for cosmic ray attenuation [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
    logical :: improvedH2CRPDissociation=.false.  !Use H2 CRP dissociation rate from [Padovani et al. 2018b](https://arxiv.org/abs/1809.04168).
    real(dp) :: diffToBindRatio=0.5  !Ratio of diffusion barrier to binding energy of all species
    logical :: h2EncounterDesorption=.true.  !Encounter desorption mechanism of Hincelin et al 2015 (H2 on H2)
    logical :: hEncounterDesorption=.false.  !Encounter desorption mechanism of Hincelin et al 2015 (H on H2)
    real(dp) :: EDEndothermicityFactor=0.0  !Account for endothermicity of moving off of H2O onto H2 by fraction of diff in binding energies
    logical :: h2StickingCoeffByh2Coverage=.false.  !Decrease sticking coeff of H2 by H2 coverage of surface
    logical :: hStickingCoeffByh2Coverage=.false.  !Decrease sticking coeff of H by H2 coverage of surface
    real(dp) :: HdiffusionBarrier=-1.0  !Diffusion barrier for atomic H on grain surface (K).
    !!This is later corrected to diffToBind*Ebind(#H) if no other value is input
    logical :: useCustomDiffusionBarriers=.true.  !Use custom diffusion barriers, instead of assuming they're a fraction of the binding energy
    logical :: separateDiffAndDesorbPrefactor=.true.  !Calculate different prefactors for diffusion and desorption
    logical :: useTSTprefactors=.false.  !Calculate diffusion and desorption prefactors using TST. Otherwise, use Hasegawa-Herbst equation.
    logical :: useCustomPrefactors=.false.  !Use custom diffusion and desorption prefactors, instead of TST or Hasegawa-Herbst values.
    logical :: useMinissaleIceChemdesEfficiency=.false.  !Use Minissale 2016 efficiency for chemical desorption on ices. If False, use Fredon 2021
    logical :: heatingFlag=.false.  !If True, heating is applied to the gas parcels.
    logical :: enforceChargeConservation = .false.  ! Enforce the chrage by keeping track of charged ions.
    logical :: enable_radiative_transfer=.false.  !Enable 1D radiative transfer calculations for spatial models (points>1).
    integer :: parcelStoppingMode=0  !Controls when parcels stop evolving in 1D freefall models: 0=never stop (default), 1=stop when outermost parcel reaches max density, 2=stop each parcel individually at max density.
    logical :: useGrainActivatedRecombination=.true.  !Whether to use grain-activated recombination (GAR) reactions
    !
    !## Input and Output
    !|Parameter|Default Value |Description|
    !| ----- | ------| ------ |
    character(256) :: outputFile=""  !File to write full output of UCLCHEM. This includes physical parameter values and all abundances at every time step.
    character(256) :: columnFile=""  !File to write specific species abundances, see outSpecies.
    character(256) :: rateConstantFile=""  !File to write rate 'constants' at each timestep. This includes physical parameter values.
    character(256) :: ratesFile=""  !File to write reaction rates (flux) at each timestep. This includes physical parameter values.
    character(256) :: heatingFile=""  !File to write heating and cooling rates at each timestep.
    integer :: writeStep=1  !Writing to columnFile only happens every writeStep timesteps.
    logical :: writeTimestepInfo=.false.  !If True, print timestep progress (current time, final time, next timestep goal) each time the target time is set.
    character(256) :: abundSaveFile=""  ! The file to save the abundances to at the end of the model.
    character(256) :: abundLoadFile=""  ! The file to load the abundances from at the start of the model.

    !## Coolant / Validation tolerances
    !|Parameter|Default Value|Description|
    !| ----- | ------| ------ |
    real(dp) :: freq_rel_tol = 1.0e-1_dp  ! Relative tolerance (fraction) for comparing file vs calculated frequencies. Default 10%; overridden by Python layer with makerates-computed value when available.
    real(dp) :: pop_rel_tol  = 1.0e-1_dp  ! Relative tolerance (fraction) for checking LTE population consistency. Can be adjusted at runtime via Generalsettings (tutorial 6).

    !## DVODE Solver Mode
    !|Parameter|Default Value|Description|
    !| ----- | ------| ------ |
    integer  :: solverMode         = 2      !DVODE ISTATE strategy: 0=always restart (ISTATE=1), 1=always continue (ISTATE=2), 2=adaptive (default)
    real(dp) :: logChangeThreshold = 1.0_dp  !log10 per-step abundance change that triggers forced ISTATE=1 restart in adaptive mode (solver_mode=2)

    !|abundSaveFile |None| File to store final abundances at the end of the model so future models can use them as the initial abundances. If not provided, no file will be produced.
    !|abundLoadFile |None| File from which to load initial abundances for the model, created through `abundSaveFile`. If not provided, the model starts from elemental gas.
    !|outSpecies|None| A space separated list of species to output to columnFile. Supplied as a separate list argument to most python functions, see python API docs.
    !
    !## Initial Abundances
    !*Unless otherwise specified, we take all abundances from Jenkins et al. 2009, using the heavily depleted case from Table 4.*
    !
    !|Parameter|Default Value |Description|
    !| ----- | ------| ------ |
    real(dp) :: metallicity=1.0  !Scale the abundances of all elements heavier than He by this factor.
    integer :: ion=2  !Sets how much elemental C is initially atomic (0= all atomic/1=50:50/2=fully ionized).
    real(dp) :: fh=0.5  !Total elemental abundance of H is always 1 by definition because abundances are relative to number of H nuclei. Use fh to set how much to initially put in atomic H, the rest goes to H2.
    real(dp) :: fhe = 0.1_dp  !Total elemental abundance of He.
    real(dp) :: fc=1.77e-04_dp  !Total elemental abundance of C.
    real(dp) :: fo  = 3.34e-04_dp  !Total elemental abundance of O.
    real(dp) :: fn  = 6.18e-05_dp  !Total elemental abundance of N.
    real(dp) :: fs  = 3.51e-6_dp  !Total elemental abundance of S.
    real(dp) :: fmg = 2.256e-06_dp  !Total elemental abundance of Mg.
    real(dp) :: fsi = 1.78e-06_dp  !Total elemental abundance of Si.
    real(dp) :: fcl = 3.39e-08_dp  !Total elemental abundance of Cl.
    real(dp) :: fp =7.78e-08_dp  !Total elemental abundance of P.
    real(dp) :: ffe =2.01e-7_dp  !Total elemental abundance of Fe.
    real(dp) :: ff = 3.6e-08_dp  !fp depleted 1/100 of solar from Asplund 2009.
    real(dp) :: fd=0.0_dp  ! The following elements are not typically used. We do not recommend any particular value.
    real(dp) :: fli=0.0_dp  !Total elemental abundance of Li.
    real(dp) :: fna=0.0_dp  !Total elemental abundance of Na.
    real(dp) :: fpah=0.0_dp  !Total initial abundance of PAHs.
    real(dp) :: f15n=0.0_dp  !Total initial abundance of 15N.
    real(dp) :: f13c=0.0_dp  !Total initial abundance of 13C.
    real(dp) :: f18O=0.0_dp  !Total initial abundance of 18O.
    !!
    !! We used to use Asplund et al. 2009,kept here for reference
    !! !initial fractional abundances of elements(from Asplund et al. 2009 ARAA table 1 -SOLAR)
    !! !note fh is fraction of H initially in H atoms. Total H is always 1.
    !! !fh=0.5;fhe = 0.1;fc  = 2.6e-04_dp;fo  = 4.6e-04_dp;fn  = 6.1e-05_dp
    !! fs  = 1.318e-05_dp;fmg = 3.981e-05_dp;fsi = 1.0e-07_dp;fcl = 3.162e-07_dp;
    !! fp=2.57e-09_dp ; ff = 3.6e-08_dp !fp depleted 1/100 of solar
    !
    !## Integration Controls
    !|Parameter|Default Value |Description|
    !| ----- | ------| ------ |
    real(dp) :: reltol=1e-8_dp  !Relative tolerance for integration, see [integration docs](/docs/trouble-integration) for advice.
    real(dp) :: abstol_factor=1.0e-14_dp  !Absolute tolerance for integration is calculated by multiplying species abundance by this factor.
    real(dp) :: abstol_min=1.0e-25_dp  !Minimum value absolute tolerances can take.
    real(dp) :: abstol_ice_factor=1.0e-10_dp  !Absolute tolerance factor for ice (grain surface + bulk) species; looser than gas to reduce stiffness from ice intermediates.
    real(dp) :: abstol_ice_min=1.0e-20_dp  !Minimum absolute tolerance for ice species.
    real(dp) :: negative_abundance_tol=1.0e-10_dp  !Abundances in (-negative_abundance_tol, 0) are solver noise clamped to 1e-30 after integration; more negative triggers NEGATIVE_ABUNDANCE_ERROR.
    real(dp) :: runtime_conservation_tolerance=0.01_dp  !Fractional tolerance for runtime element conservation check (1% by default). Set negative to disable.
    real(dp) :: reltol_phys=1.0e-4_dp  !Relative tolerance for physical variables (temperature, density) in integration.
    real(dp) :: abstol_phys_factor=1.0e-4_dp  !Absolute tolerance factor for physical variables (temperature, density).
    real(dp) :: abstol_T_min=1.0e-2_dp  !Minimum absolute tolerance for gas temperature (K).
    real(dp) :: abstol_nH_min=1.0_dp  !Minimum absolute tolerance for gas density (cm^-3).
    integer :: MXSTEP=10000  !Maximum steps allowed in integration before warning is thrown. ! HAS TO BE INT4 instead of INT8
    !
    !## Here be Dragons
    !*These are not recommended to be changed unless you know what you are doing*
    !
    !|Parameter|Default Value |Description|
    !| ----- | ------| ------ |
    real(dp) :: ebmaxh2=1.21e3_dp  ! Maximum binding energy of species desorbed by H2 formation.
    real(dp) :: ebmaxcr=1.21e3_dp  ! Maximum binding energy of species desorbed by cosmic ray ionization.
    real(dp) :: ebmaxuvcr=1.0e4_dp  ! Maximum binding energy of species desorbed by UV photons.
    real(dp) :: epsilon=0.01_dp  !Number of molecules desorbed per H2 formation.
    real(dp) :: uv_yield=0.03_dp  !Number of molecules desorbed per UV photon. The yield is extrapolated from Oberg et al. 2009
    real(dp) :: phi=1.0e5_dp  !Number of molecules desorbed per cosmic ray ionization.
    real(dp) :: uvcreff=1.0e-3_dp  !Ratio of CR induced UV photons to ISRF UV photons.
    real(dp) :: omega=0.5  !Dust grain albedo.
    real(dp) :: lower_limit_gastemp=10.0  !Lower limit for gas temperature in K when heating is enabled.
    real(dp) :: upper_limit_gastemp=1.0e4_dp  !Upper limit for gas temperature in K when heating is enabled.
    real(dp) :: heating_temp_abstol=1.0_dp  !Absolute temperature change threshold (K) for recalculating heating/cooling in ODE function.
    real(dp) :: heating_temp_reltol=1.0e-2_dp  !Relative temperature change threshold for recalculating heating/cooling in ODE function.
    real(dp) :: lower_limit_dusttemp=10.0  !Lower limit for dust temperature in K when heating is enabled.
    real(dp) :: upper_limit_dusttemp=1.0e3_dp  !Upper limit for dust temperature in K when heating is enabled.
    real(dp) :: maxGrainTemp=150.0  !Dust temperature (K) above which grain surface chemistry is disabled and H2 formation is parameterized.
    integer :: parameterizeH2Form=2  !H2 formation mode: 0=always off, 1=always on (parameterized), 2=explicit LH/ER below maxGrainTemp, parameterized above (default).
    real(dp) :: min_desorption_rate_constant = 1.0e-60_dp  ! Floor on desorption rate constants k (s^-1): k in (0, min_desorption_rate) are zeroed to avoid underflow. 0 disables.
    real(dp) :: max_desorption_rate_constant_factor = 10.0_dp  ! Dynamic cap on thermal desorption: effective cap = clamp(factor/(targetTime-currentTime), min_constant_cap, max_constant_cap) [s^-1]. 0 disables.
    real(dp) :: min_desorption_rate_constant_cap = 1.0_dp  ! Lower bound on the dynamic cap, in yr^-1 (timescale 1 yr): k slower than this are never capped.
    real(dp) :: max_desorption_rate_constant_cap = 3.16e7_dp  ! Upper bound on the dynamic cap, in yr^-1 (= 1 s^-1): k faster than 1 s are always capped regardless of timestep.
    !|alpha|{1:0.0,2:0.0}| Set alpha coefficients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how alpha is used for each reaction type.|
    !|beta|{1:0.0,2:0.0}| Set beta coefficients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how beta is used for each reaction type.|
    !|gama|{1:0.0,2:0.0}| Set gama coefficients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how gama is used for each reaction type.|
contains
    ! Add a dummy subroutine to help f2py compile: https://github.com/numpy/numpy/issues/27167
    subroutine DUMMY_TWO(dummy_two_output)
        integer, intent(out) :: dummy_two_output
        dummy_two_output = 2
    end subroutine DUMMY_TWO
end module DEFAULTPARAMETERS
