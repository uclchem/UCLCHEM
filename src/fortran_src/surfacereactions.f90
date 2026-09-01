module SurfaceReactions
  use constants, only: dp, K_BOLTZ, K_BOLTZ_SI, PI, AMU, HP_SI, N_AVOGADRO, &
      REDUCED_PLANCK, KCAL_TO_JOULE
  use DEFAULTPARAMETERS
  use f2py_constants, only: nSpec, nReac, NO_REACTANT_OR_PRODUCT, &
      MISSING_VALUE_FLOAT, MISSING_VALUE_INTEGER
  use network, only: iceList, bulkList, surfaceList, gasIceList, bindingEnergy, diffusionBarrier, &
      re1, re2, p1, p2, p3, p4, alpha, beta, gama, reducedMasses, &
      mass, specName, atomCounts, surfSwapReacs, bulkSwapReacs, &
      moleculeIsLinear, formationEnthalpy, LHDEScorrespondingLHreacs, &
      customVdes, customVdiff, THREE_PHASE, inertiaProducts, &
      ngh2, nh2, ngn, ngo, ngoh, nh, ngh
  !f2py INTEGER, parameter :: dp
  use numerics, only: is_equal, linear_interp

  implicit none

  public

  real(dp) :: surfaceCoverage,totalSwap,bulkLayersReciprocal
  real(dp) :: safeMantle,safeBulk

  !Silicate grain properties for H2 Formation
  real(dp), parameter :: SILICATE_MU=0.005_dp  ! Fraction of newly formed H2 that stays on the grain surface
  real(dp), parameter :: SILICATE_E_S=110.0_dp  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
  real(dp), parameter :: SILICATE_E_H2=320.0_dp  ! Desorption energy of H2 molecules (K)
  real(dp), parameter :: SILICATE_E_HP=450.0_dp  ! Desorption energy of physisorbed H atoms (K)
  real(dp), parameter :: SILICATE_E_HC=3.0e4_dp   ! Desorption energy of chemisorbed H atoms (K)
  real(dp), parameter :: SILICATE_NU_H2=3.0e12_dp  ! Vibrational frequency of H2 molecules in surface sites (s^-1)
  real(dp), parameter :: SILICATE_NU_HC=1.3e13_dp  ! Vibrational frequency of H atoms in their surface sites (s^-1)
  real(dp), parameter :: SILICATE_CROSS_SECTION=8.473e-22_dp  !*CROSS_SECTION_SCALE ! Silicate grain cross section per H nucleus (cm^-2/nucleus)

  !Graphite grain properties for H2 Formation
  real(dp), parameter :: GRAPHITE_MU=0.005_dp   ! Fraction of newly formed H2 that stays on the grain surface
  real(dp), parameter :: GRAPHITE_E_S=260.0_dp  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
  real(dp), parameter :: GRAPHITE_E_H2=520.0_dp  ! Desorption energy of H2 molecules (K)
  real(dp), parameter :: GRAPHITE_E_HP=800.0_dp  ! Desorption energy of physisorbed H atoms (K)
  real(dp), parameter :: GRAPHITE_E_HC=3.0e4_dp   ! Desorption energy of chemisorbed H atoms (K)
  real(dp), parameter :: GRAPHITE_NU_H2=3.0e12_dp  ! Vibrational frequency of H2 molecules in surface sites (s^-1)
  real(dp), parameter :: GRAPHITE_NU_HC=1.3e13_dp  ! Vibrational frequency of H atoms in their surface sites (s^-1)
  real(dp), parameter :: GRAPHITE_CROSS_SECTION=7.908e-22_dp  !*CROSS_SECTION_SCALE ! Graphite grain cross section per H nucleus (cm^-2/nucleus)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Grain surface parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), parameter :: GAS_DUST_MASS_RATIO=100.0_dp,GRAIN_RADIUS=1.0e-5_dp, GRAIN_DENSITY = 3.0_dp  ! Mass density of a dust grain
  real(dp), parameter :: THERMAL_VEL= sqrt(8.0_dp*K_BOLTZ/(PI*AMU))  !Thermal velocity without the factor of SQRT(T/m) where m is moelcular mass in amu

  !reciprocal of fractional abundance of dust grains (we only divide by number density so better to store reciprocal)
  real(dp), parameter :: GAS_DUST_DENSITY_RATIO = (4.0_dp*PI*(GRAIN_RADIUS**3)*GRAIN_DENSITY*GAS_DUST_MASS_RATIO)/(3.0_dp * AMU)
  !Grain area per h nuclei, values taken from Cazaux & Tielens 2004 via UCL-PDR to match H2 formation rate
  real(dp), parameter :: GRAIN_CROSSSECTION_PER_H=0.5_dp*(7.908e-22_dp+8.473e-22_dp)
  real(dp), parameter :: GRAIN_SURFACEAREA_PER_H=4.0_dp*GRAIN_CROSSSECTION_PER_H  !2.0*4.0*PI*GRAIN_RADIUS*GRAIN_RADIUS/GAS_DUST_DENSITY_RATIO

  !Below are values for grain surface reactions
  logical, parameter :: DIFFUSE_REACT_COMPETITION=.true., GRAINS_HAVE_ICE=.true.
  real(dp), parameter :: NUM_MONOLAYERS_IS_SURFACE=2.0_dp  ! Number of monolayers to count as surface
  logical, parameter :: useGarrod2011Transfer=.true.  ! Use Garrod 2011 transfer upon net desorption
  logical, parameter :: useCustomReducedMass=.true.  ! Use custom predicted reduced mass for tunneling
  real(dp), parameter :: CHEMICAL_BARRIER_THICKNESS = 1.40e-8_dp  ! Parameter used to compute the probability for a surface reaction with
  !! activation energy to occur through quantum tunneling (Hasegawa et al. Eq 6 (1992).)
  real(dp), parameter :: SURFACE_SITE_DENSITY = 1.5e15_dp  ! site density on one grain [cm-2]
  real(dp), parameter :: VDIFF_PREFACTOR=2.0_dp*K_BOLTZ*SURFACE_SITE_DENSITY/PI/PI/AMU
  real(dp), parameter :: NUM_SITES_PER_GRAIN = SURFACE_SITE_DENSITY * (4.0_dp*PI*GRAIN_RADIUS**2)
  real(dp), parameter :: EFFECTIVE_SURFACE_MASS = 120.0_dp  ! Effective mass of grain surface in Dalton

  ! TST prefactor constants
  real(dp), parameter :: HH_VDES_PREFACTOR=2.0_dp*K_BOLTZ_SI*SURFACE_SITE_DENSITY*1.0e4_dp/(PI*PI*AMU)
  real(dp), parameter :: TST_VDES_PREFACTOR = 2.0_dp*PI*K_BOLTZ_SI**2*AMU/ &
    (SURFACE_SITE_DENSITY*1.0e4_dp*(HP_SI**3)*1.0e3_dp)

  ! Encounter desorption constants
  real(dp), parameter :: H2_ON_H2_BINDING_ENERGY=23.0_dp  ! K
  real(dp), parameter :: H_ON_H2_BINDING_ENERGY=45.0_dp  ! K

  real(dp), parameter :: MIN_SURFACE_ABUND=1.0e-20_dp

  ! Desorption fraction arrays for LHDES/ERDES reactions (pre-calculated at initialization)
  real(dp), dimension(nReac) :: desorptionFractionsBare, desorptionFractionsFullCoverage

  real(dp), allocatable :: vdiff(:), vdes(:)
contains
  !=======================================================================
  !
  !  Calculate the rate of molecular hydrogen (H2) formation on grains
  !  using the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
  !  Cazaux & Tielens (2004, ApJ, 604, 222).
  !
  !-----------------------------------------------------------------------
  pure function h2FormEfficiency(gasTemp,dustTemp) result(rate)
    real(dp) :: rate
    real(dp), intent(in) :: gasTemp,dustTemp

    real(dp) :: THERMAL_VELOCITY,STICKING_COEFFICIENT
    real(dp) :: FLUX,FACTOR1,FACTOR2,EPSILON
    real(dp) :: SILICATE_FORMATION_EFFICIENCY,GRAPHITE_FORMATION_EFFICIENCY
    !  Mean thermal velocity of hydrogen atoms (cm s^-1)
    THERMAL_VELOCITY=1.45e5_dp*sqrt(gasTemp/1.0e2_dp)

    !  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
    !  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
    STICKING_COEFFICIENT=1.0_dp/(1.0_dp+0.04_dp*sqrt(gasTemp+dustTemp) &
                    & + 0.2_dp*(gasTemp/1.0e2_dp)+0.08_dp*(gasTemp/1.0e2_dp)**2)

    FLUX=1.0e-10_dp  ! Flux of H atoms in monolayers per second (mLy s^-1)

    FACTOR1=SILICATE_MU*FLUX/(2*SILICATE_NU_H2*exp(-SILICATE_E_H2/dustTemp))

   FACTOR2=1.0_dp*(1.0_dp+sqrt((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2 &
        & /4.0_dp*exp(-SILICATE_E_S/dustTemp)

   EPSILON=1.0_dp/(1.0_dp+SILICATE_NU_HC/(2*FLUX)*exp(-1.5_dp*SILICATE_E_HC/dustTemp) &
              & *(1.0_dp+sqrt((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2)

   SILICATE_FORMATION_EFFICIENCY=1.0_dp/(1.0_dp+FACTOR1+FACTOR2)*EPSILON


   FACTOR1=GRAPHITE_MU*FLUX/(2*GRAPHITE_NU_H2*exp(-GRAPHITE_E_H2/dustTemp))

   FACTOR2=1.0_dp*(1.0_dp+sqrt((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2 &
        & /4.0_dp*exp(-GRAPHITE_E_S/dustTemp)

   EPSILON=1.0_dp/(1.0_dp+GRAPHITE_NU_HC/(2*FLUX)*exp(-1.5_dp*GRAPHITE_E_HC/dustTemp) &
              & *(1.0_dp+sqrt((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2)

   GRAPHITE_FORMATION_EFFICIENCY=1.0_dp/(1.0_dp+FACTOR1+FACTOR2)*EPSILON

!  Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222)
   rate=0.5_dp*THERMAL_VELOCITY*(SILICATE_CROSS_SECTION*SILICATE_FORMATION_EFFICIENCY &
    & + GRAPHITE_CROSS_SECTION*GRAPHITE_FORMATION_EFFICIENCY)*STICKING_COEFFICIENT

  end function h2FormEfficiency

  pure subroutine bulkSurfaceExchangeReactions(rate,dustTemperature)
    real(dp), intent(inout), dimension(nReac) :: rate
    real(dp), intent(in) :: dustTemperature
    if (THREE_PHASE) then
      call bulkToSurfaceSwappingRateConstants(rate,bulkswapReacs(1),bulkswapReacs(2),dustTemperature)
      rate(surfSwapReacs(1):surfSwapReacs(2))=surfaceToBulkSwappingRateConstants(dustTemperature)
    end if
  end subroutine bulkSurfaceExchangeReactions

  !surface abundance multiplied by this value gives fraction of surface covered by material
  pure function bulkGainFromMantleBuildUp() result(rate)
    real(dp) :: rate
    rate=0.5_dp*GAS_DUST_DENSITY_RATIO/NUM_SITES_PER_GRAIN
  end function bulkGainFromMantleBuildUp

  pure function surfaceToBulkSwappingRateConstants(dustTemperature) result(rate)
    real(dp), intent(in) :: dustTemperature
    real(dp) :: rate

    if ((dustTemperature > maxGrainTemp) .or. (safeMantle < MIN_SURFACE_ABUND)) then
              rate = 0.0_dp
    else
        rate = 1.0_dp
    end if
  end function surfaceToBulkSwappingRateConstants


  pure subroutine bulkToSurfaceSwappingRateConstants(rate,idx1,idx2,dustTemperature)
    real(dp), intent(inout), dimension(nReac) :: rate
    real(dp), intent(in) :: dustTemperature
    integer, intent(in) :: idx1, idx2

    integer :: i, j

    if ((dustTemperature > maxGrainTemp) .or. (safeMantle < MIN_SURFACE_ABUND)) then
        rate(idx1:idx2) = 0.0_dp
    else
        do i=idx1,idx2
            do j=lbound(iceList,1),ubound(iceList,1)
                if (iceList(j) == re1(i)) then
                  rate(i)=vdiff(j)*exp(-bindingEnergy(j)/dustTemperature)
                end if
            end do
        end do
    end if
  end subroutine bulkToSurfaceSwappingRateConstants

  !----------------------------------------------------------------------------------------------------
!Reactions on the surface treated by evaluating diffusion rates across grains and accounting
!For competition with chemical desorption. Products remain bound ('DIFF') or are desorbed ('CHEMDES')
!Units of s-1.
!David Quenard 2017 Arxiv:1711.05184
!----------------------------------------------------------------------------------------------------
  pure function getDiffusionBarrier(iceListIndex) result(diffusionBarrierValue)
    !! Calculate the diffusion barrier of a species in ice.
    ! For all species, except hydrogen, this is assumed to be a fraction of the binding energy.
    ! For hydrogen, it is this by default, but can be set differently by an option in the input
    integer, intent(in) :: iceListIndex

    real(dp) :: diffusionBarrierValue

    if ((useCustomDiffusionBarriers) .and. &
        (.not. is_equal(diffusionBarrier(iceListIndex), MISSING_VALUE_FLOAT))) then
        diffusionBarrierValue = diffusionBarrier(iceListIndex)
    else if (iceList(iceListIndex) == ngh) then
        diffusionBarrierValue = HdiffusionBarrier
    else
        diffusionBarrierValue = diffToBindRatio*bindingEnergy(iceListIndex)
    end if
  end function getDiffusionBarrier

  pure function getReactionProbability(reacIndx, index1, index2, dustTemperature) result(reacProb)
    integer, intent(in) :: reacIndx, index1, index2
    real(dp), intent(in) :: dustTemperature

    real(dp) :: reducedMass
    real(dp) :: reacProb, tunnelProb
    !Calculate classical activation energy barrier exponent
    reacProb = gama(reacIndx)/dustTemperature
    !Calculate quantum activation energy barrier exponent
    reducedMass = reducedMasses(reacIndx)
    if ((.NOT. useCustomReducedMass) .OR. (is_equal(reducedMass, MISSING_VALUE_FLOAT))) then
        ! reducedMasses(reacIndx) should never be MISSING_VALUE_FLOAT,
        ! because it is set by MakeRateConstants, but just in case we calculate it here.
        reducedMass = mass(iceList(index1)) * mass(iceList(index2)) / (mass(iceList(index1)) + mass(iceList(index2)))
    end if
    tunnelProb = 2.0_dp *CHEMICAL_BARRIER_THICKNESS/REDUCED_PLANCK * sqrt(2.0_dp*AMU*reducedMass*K_BOLTZ*gama(reacIndx))

    !Choose fastest between classical and tunneling
    if (reacProb>tunnelProb) reacProb=tunnelProb
  end function getReactionProbability

  pure function getDiffusionReactionRateConstant(reacIndx,dustTemperature) result(diffusionReactionRateConstant)
    integer, intent(in) :: reacIndx
    real(dp), intent(in) :: dustTemperature

    real(dp) :: diffusionReactionRateConstant

    real(dp) :: diffuseProb,desorbProb,reactionProb
    integer :: index1,index2,i


    !want position of species in the grain array but gas phase species aren't in there
    !so store species index
    index1=re1(reacIndx)
    index2=re2(reacIndx)

    !then try to overwrite with position in grain array
    do i=lbound(iceList,1),ubound(iceList,1)
        if (iceList(i) == index1) index1 = i
        if (iceList(i) == index2) index2 = i
    end do

    !Hasegawa 1992 diffusion rate. RateConstant that two species diffuse and meet on grain surface
    diffuseProb = vdiff(index1)*exp(-getDiffusionBarrier(index1)/dustTemperature)
    diffuseProb = diffuseProb+ (vdiff(index2)*exp(-getDiffusionBarrier(index2)/dustTemperature))

    !probability a reactant will just desorb
    desorbProb = vdes(index1)*exp(-bindingEnergy(index1)/dustTemperature)
    desorbProb = desorbProb + vdes(index2)*exp(-bindingEnergy(index2)/dustTemperature)

    !Overall reaction probability is chance of reaction occurring on meeting * diffusion rate
    reactionProb = max(vdiff(index1),vdiff(index2)) * exp(-getReactionProbability(reacIndx, index1, index2, dustTemperature))


    ! Keff from Garrod & Pauly 2011 and Ruaud+2016
    ! Actual reaction probability is Preac/(Preac+Pevap+Pdiffuse), accounting for the other possible processes
    if(DIFFUSE_REACT_COMPETITION) then
       reactionProb = reactionProb/(reactionProb + desorbProb + diffuseProb)
    else
        reactionProb = 1.0_dp
    end if

    !see Eq A1 of Quenard et al. 2018
    !NUM_SITES_PER_GRAIN should be multiplied by n_dust as in A1
    !n_dust=density/GAS_DUST_DENSITY_RATIO so we use the 1/density to cancel the density in odes.f90 and drop it here
    diffusionReactionRateConstant=alpha(reacIndx) *reactionProb* diffuseProb*GAS_DUST_DENSITY_RATIO/NUM_SITES_PER_GRAIN
  end function getDiffusionReactionRateConstant

! ---------------------------------------------------------------------
!  Chemical Reactive Desorption (CRD)
! David Quenard 2017 Arxiv:1711.05184
! From Minissalle+ 2016 and Vasyunin+ 2016
!
! Modern implementation uses separate functions for bare grains and
! ice-covered grains with ice-coverage-dependent interpolation:
! - getDesorptionFractionBare: Minissale+ 2016 formulation
! - getDesorptionFractionFullCoverage: Fredon+ 2021, Furuya+ 2022
! - desorptionFractionIncludingIce: interpolates between the two
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Get bare grain desorption fraction (Minissale+ 2016)
! ---------------------------------------------------------------------
  real(dp) function getDesorptionFractionBare(reacIndx, LHDESindex) result(desorptionFractionBare)
    integer, intent(in) :: reacIndx, LHDESindex

    integer :: reactIndex1,reactIndex2,degreesOfFreedom,i,j
    integer :: desorbingIndex, desorbingOnGrainIndex, desorbingIceListIndex
    integer :: productIndex(4)

    real(dp) :: deltaEnthalpy,epsilonCd,productEnthalpy
    real(dp) :: bindingEnergyDesorbingSpec, chi
    logical :: twoProductReaction

    reactIndex1=0
    reactIndex2=0
    productIndex=0


     if (.NOT.(ANY(iceList == re1(reacIndx)) .OR. (ANY(iceList == re2(reacIndx))))) then
        ! Gasphase reactions do not need to be calculated, should be 0
        desorptionFractionBare = 0.0_dp
        return
     end if
     if (ANY(bulkList == re1(reacIndx))) then
         ! No chemical desorption from bulk ice allowed
         desorptionFractionBare = 0.0_dp
         return
     end if

    !Get indices of grain surface version of products products
    productIndex=0
    !Arrays like binding energy and formation enthalpy are indexed by position in iceList
    !rather than species list. Need to find position in grain list where every reactant and product appears
    !bearing in mind that Eley-Rideal reactions can have reactants in gas phase and CHEMDES has products in gas
    do i=lbound(iceList,1),ubound(iceList,1)
        !check grain lists for reactants
        if (iceList(i) == re1(reacIndx)) reactIndex1 = i
        if (gasIceList(i) == re1(reacIndx)) reactIndex1 = i
        !check equivalent gas list in case of ER reaction.
        if (iceList(i) == re2(reacIndx)) reactIndex2 = i
        if (gasIceList(i) == re2(reacIndx)) reactIndex2 = i

        if (iceList(i) == p1(reacIndx)) productIndex(1) = i
        if (iceList(i) == p2(reacIndx)) productIndex(2) = i
        if (iceList(i) == p3(reacIndx)) productIndex(3) = i
        if (iceList(i) == p4(reacIndx)) productIndex(4) = i

        if (gasIceList(i) == p1(reacIndx)) productIndex(1) = i
        if (gasIceList(i) == p2(reacIndx)) productIndex(2) = i
        if (gasIceList(i) == p3(reacIndx)) productIndex(3) = i
        if (gasIceList(i) == p4(reacIndx)) productIndex(4) = i
    end do

    if ((reactIndex1 == 0) .or. (reactIndex2 == 0)) then
        write(*, *) "COULD NOT DETERMINE INDEX OF REACTANT FOR REACTION:"
        write(*,*) specName(re1(reacIndx)), specName(re2(reacIndx)), "->", &
            specName(p1(reacIndx)), specName(p2(reacIndx)), specName(p3(reacIndx))
        stop 1
    end if

    if (p2(reacIndx) == NO_REACTANT_OR_PRODUCT) then
        ! Only one product, and so that one product is desorbing
        desorbingIndex = 1
        desorbingOnGrainIndex = p1(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .false.
    else if (p1(LHDEScorrespondingLHreacs(LHDESindex)) /= p1(reacIndx)) then  ! p1 is desorbing
        desorbingIndex = 1
        desorbingOnGrainIndex = p1(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .true.
    else if (p2(LHDEScorrespondingLHreacs(LHDESindex)) /= p2(reacIndx)) then  ! p2 is desorbing
        desorbingIndex = 2
        desorbingOnGrainIndex = p2(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .true.
    else if (p3(LHDEScorrespondingLHreacs(LHDESindex)) /= p3(reacIndx)) then  ! p3 is desorbing
        desorbingIndex = 3
        desorbingOnGrainIndex = p3(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .true.
    else
        write(*,*) "COULD NOT DETERMINE DESORBING PRODUCT INDEX OF REACTION:"
        write(*,*) specName(re1(reacIndx)), specName(re2(reacIndx)), "->", &
            specName(p1(reacIndx)), specName(p2(reacIndx)), specName(p3(reacIndx))
        write(*,*) "LHDES INDEX:", LHDESindex
        write(*,*) "REAC INDEX:", reacIndx
        write(*,*) "CORRESPONDING LH INDEX:", LHDEScorrespondingLHreacs(LHDESindex)
        write(*,*) "CORRESPONDING LH REACTION:"
        write(*,*) specName(re1(LHDEScorrespondingLHreacs(LHDESindex))), &
            specName(re2(LHDEScorrespondingLHreacs(LHDESindex))), "->", &
            specName(p1(LHDEScorrespondingLHreacs(LHDESindex))), &
            specName(p2(LHDEScorrespondingLHreacs(LHDESindex))), &
            specName(p3(LHDEScorrespondingLHreacs(LHDESindex)))
        stop 1
    end if

    ! Now we know which product desorbs, we just have to calculate bare desorption prob using Minissale et al 2016.

    desorbingIceListIndex = 0
    productEnthalpy = 0.0_dp
    do i = 1,4
        if (productIndex(i) /= 0) then
            if (i == desorbingIndex) then
                do j = LBOUND(iceList, 1), UBOUND(iceList, 1)
                    if (iceList(j) == desorbingOnGrainIndex) then
                        desorbingIceListIndex = j
                        productEnthalpy = productEnthalpy + formationEnthalpy(j)
                    end if
                end do
            else
                productEnthalpy = productEnthalpy + formationEnthalpy(productIndex(i))
            end if
        end if
    end do

    deltaEnthalpy = productEnthalpy - (formationEnthalpy(reactIndex1) + formationEnthalpy(reactIndex2))
    ! If deltaEnthalpy > 0: endothermic
    ! If deltaEnthalpy < 0: exothermic, energy released to environment

    if (deltaEnthalpy > 0.0_dp) then
        ! Endothermic reactions do not induce chemical desorption
        desorptionFractionBare = 0.0_dp
        return
    end if

    ! Now we use deltaEnthalpy as a measure of exothermicity, i.e. the amount of energy released
    deltaEnthalpy = -deltaEnthalpy

    !Convert from kcal to J, from J to K and from moles-1 to reactions-1
    deltaEnthalpy = deltaEnthalpy*KCAL_TO_JOULE/(K_BOLTZ_SI*N_AVOGADRO)

    bindingEnergyDesorbingSpec = bindingEnergy(desorbingIceListIndex)
    if (deltaEnthalpy < bindingEnergyDesorbingSpec) then
        desorptionFractionBare = 0.0_dp
        return
    end if

    epsilonCd = mass(desorbingOnGrainIndex)
    !epsilonCd is the fraction of kinetic energy kept my the product when it collides with grain surface
    epsilonCd = ((epsilonCd - EFFECTIVE_SURFACE_MASS) / (epsilonCd + EFFECTIVE_SURFACE_MASS))**2

    if (.NOT. twoProductReaction) then
        chi = 1.0_dp
    else
        ! Distribute energy in case of two product reaction
        ! chi_i = m_j/(m_i+m_j)
        if (desorbingIndex == 1) then
            chi = mass(p2(reacIndx)) / (mass(p1(reacIndx))+mass(p2(reacIndx)))
        else if (desorbingIndex == 2) then
            chi = mass(p1(reacIndx)) / (mass(p1(reacIndx))+mass(p2(reacIndx)))
        else
            write(*,*) "MINISSALE 2016 METHOD FOR CHEMICAL DESORPTION IS NOT VALID FOR DESORBINDEX > 2"
            stop 1
        end if
    end if

    epsilonCd = epsilonCd * chi

    if (epsilonCd * deltaEnthalpy < bindingEnergyDesorbingSpec) then
        desorptionFractionBare = 0.0_dp
        return
    end if

    degreesOfFreedom = 3 * atomCounts(desorbingOnGrainIndex)
    desorptionFractionBare = exp((-bindingEnergyDesorbingSpec*REAL(degreesOfFreedom)) / (epsilonCd * deltaEnthalpy))
  end function getDesorptionFractionBare

! ---------------------------------------------------------------------
! Get full ice coverage desorption fraction (Fredon+ 2021, Furuya+ 2022)
! ---------------------------------------------------------------------
  function getDesorptionFractionFullCoverage(reacIndx, LHDESindex) result (desorptionFractionFullCoverage)
    integer, intent(in) :: reacIndx, LHDESindex

    integer :: reactIndex1,reactIndex2,i,j
    integer :: productIndex(4)

    real(dp) :: deltaEnthalpy,productEnthalpy
    real(dp) :: bindingEnergyDesorbingSpec, chi
    integer :: desorbingIndex, desorbingOnGrainIndex, desorbingIceListIndex
    logical :: twoProductReaction

    real(dp) :: desorptionFractionFullCoverage

    if (.NOT.(ANY(iceList == re1(reacIndx)) .OR. (ANY(iceList == re2(reacIndx))))) then
       ! Gasphase reactions do not need to b calculated, should be 0
       desorptionFractionFullCoverage = 0.0_dp
       return
    end if
    if (ANY(bulkList == re1(reacIndx))) then
        ! No chemical desorption from bulk ice allowed
        desorptionFractionFullCoverage = 0.0_dp
        return
    end if

    if (useMinissaleIceChemdesEfficiency) then
        desorptionFractionFullCoverage = desorptionFractionsBare(reacIndx)/10.0_dp    !< See Minisalle et al. 2016 for icy grain surface.
        ! Special case of OH+H, O+H, N+N on ices, see same paper
        if (re1(reacIndx)==ngn.and.re2(reacIndx)==ngn) desorptionFractionFullCoverage = 0.5_dp
        if ((re1(reacIndx)==ngo.and.re2(reacIndx)==nh) &
            &.or. (re1(reacIndx)== nh.and.re2(reacIndx)==ngo)) desorptionFractionFullCoverage = 0.3_dp
        if ((re1(reacIndx)==ngoh.and.re2(reacIndx)==nh) &
            &.or. (re1(reacIndx)==nh.and.re2(reacIndx)==ngoh)) desorptionFractionFullCoverage = 0.25_dp
        return
    end if

    !Get indices of grain surface version of products products
    productIndex = 0
    reactIndex1 = 0
    reactIndex2 = 0
    desorbingIceListIndex = 0
    !Arrays like binding energy and formation enthalpy are indexed by position in iceList
    !rather than species list. Need to find position in grain list where every reactant and product appears
    !bearing in mind that Eley-Rideal reactions can have reactants in gas phase and CHEMDES has products in gas
    do i=lbound(iceList,1),ubound(iceList,1)
        !check grain lists for reactants
        if (iceList(i) == re1(reacIndx)) reactIndex1 = i
        if (gasIceList(i) == re1(reacIndx)) reactIndex1 = i
        !check equivalent gas list in case of ER reaction.
        if (iceList(i) == re2(reacIndx)) reactIndex2 = i
        if (gasIceList(i) == re2(reacIndx)) reactIndex2 = i

        if (iceList(i) == p1(reacIndx)) productIndex(1) = i
        if (iceList(i) == p2(reacIndx)) productIndex(2) = i
        if (iceList(i) == p3(reacIndx)) productIndex(3) = i
        if (iceList(i) == p4(reacIndx)) productIndex(4) = i

        if (gasIceList(i) == p1(reacIndx)) productIndex(1) = i
        if (gasIceList(i) == p2(reacIndx)) productIndex(2) = i
        if (gasIceList(i) == p3(reacIndx)) productIndex(3) = i
        if (gasIceList(i) == p4(reacIndx)) productIndex(4) = i
    end do

    if (ALL(productIndex == 0)) then
        write(*,*) "ERROR Could not determine productIndex for LHDES reaction"
        write(*,*) "ERROR", specname(re1(reacIndx)), "+", specname(re2(reacIndx)), "->", &
            & specname(p1(reacIndx)), "+", specname(p2(reacIndx))
        stop
    end if

    if (p2(reacIndx) == NO_REACTANT_OR_PRODUCT) then
        ! Only one product, and so that one product is desorbing
        desorbingIndex = 1
        desorbingOnGrainIndex = p1(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .false.
    else if (p1(LHDEScorrespondingLHreacs(LHDESindex)) /= p1(reacIndx)) then  ! p1 is desorbing
        desorbingIndex = 1
        desorbingOnGrainIndex = p1(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .true.
    else if (p2(LHDEScorrespondingLHreacs(LHDESindex)) /= p2(reacIndx)) then  ! p2 is desorbing
        desorbingIndex = 2
        desorbingOnGrainIndex = p2(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .true.
    else if (p3(LHDEScorrespondingLHreacs(LHDESindex)) /= p3(reacIndx)) then  ! p3 is desorbing
        desorbingIndex = 3
        desorbingOnGrainIndex = p3(LHDEScorrespondingLHreacs(LHDESindex))
        twoProductReaction = .true.
    else
        write(*,*) "COULD NOT DETERMINE DESORBING PRODUCT INDEX OF REACTION:"
        write(*,*) specName(re1(reacIndx)), specName(re2(reacIndx)), "->", &
            specName(p1(reacIndx)), specName(p2(reacIndx)), specName(p3(reacIndx))
        stop
    end if


    ! Now we know which product desorbs, we just have to calculate bare desorption prob using Minissale et al 2016.
    productEnthalpy = 0.0_dp
    do i = 1,4
        if (productIndex(i) /= 0) then
            if (i == desorbingIndex) then
                do j = LBOUND(iceList, 1), UBOUND(iceList, 1)
                    if (iceList(j) == desorbingOnGrainIndex) then
                        desorbingIceListIndex = j
                        productEnthalpy = productEnthalpy + formationEnthalpy(j)
                    end if
                end do
            else
                productEnthalpy = productEnthalpy + formationEnthalpy(productIndex(i))
            end if
        end if
    end do

    if (reactIndex1 == 0 .OR. reactIndex2 == 0) then
        write(*,*) "ERROR getDesorptionFractionFullCoverage: reactIndex not set, returning 0 for reacIndx", &
            & reacIndx, " reactIndex1=", reactIndex1, " reactIndex2=", reactIndex2
        stop
    end if

    deltaEnthalpy = productEnthalpy - (formationEnthalpy(reactIndex1) + formationEnthalpy(reactIndex2))
    ! If deltaEnthalpy > 0: endothermic
    ! If deltaEnthalpy < 0: exothermic, energy released to environment

    if (deltaEnthalpy > 0.0_dp) then
        ! Endothermic reactions do not induce chemical desorption
        desorptionFractionFullCoverage = 0.0_dp
        return
    end if

    ! Now we use deltaEnthalpy as a measure of exothermicity, i.e. the amount of energy released
    deltaEnthalpy = -deltaEnthalpy

    !Convert from kcal to J, from J to K and from moles-1 to reactions-1
    deltaEnthalpy = deltaEnthalpy*KCAL_TO_JOULE/(K_BOLTZ_SI*N_AVOGADRO)

    if (desorbingIceListIndex == 0) then
        write(*,*) "ERROR getDesorptionFractionFullCoverage: desorbingIceListIndex is 0; reacIndx=", reacIndx
        desorptionFractionFullCoverage = 0.0_dp
        write(*,*) "ERROR", specname(re1(reacIndx)), "+", specname(re2(reacIndx)), "->", specname(p1(reacIndx)), "+", &
            & specname(p2(reacIndx))
        write(*,*) "ERROR Corresponding LH reaction:", specname(re1(LHDEScorrespondingLHreacs(LHDESindex))), &
            & specname(re2(LHDEScorrespondingLHreacs(LHDESindex))), specname(p1(LHDEScorrespondingLHreacs(LHDESindex)))
        stop
    end if
    if (desorbingIceListIndex < LBOUND(bindingEnergy,1) .OR. desorbingIceListIndex > UBOUND(bindingEnergy,1)) then
        write(*,*) "ERROR getDesorptionFractionFullCoverage: desorbingIceListIndex out of bounds; reacIndx=", &
            & reacIndx, " idx=", desorbingIceListIndex
        desorptionFractionFullCoverage = 0.0_dp
        write(*,*) "ERROR", re1(reacIndx), "+", re2(reacIndx), "->", p1(reacIndx), "->", p2(reacIndx)
        stop
    end if

    bindingEnergyDesorbingSpec = bindingEnergy(desorbingIceListIndex)
    if (deltaEnthalpy < bindingEnergyDesorbingSpec) then
        desorptionFractionFullCoverage = 0.0_dp
        return
    end if

    if (.NOT. twoProductReaction) then
        chi = 0.07_dp  ! chi_1 approx 0.07, Furuya et al, 2022
    else
        if (desorbingIndex == 1) then
            chi = 0.2_dp * mass(p2(reacIndx)) / (mass(p1(reacIndx))+mass(p2(reacIndx)))  ! Assume that chi_2 = 0.2
        else if (desorbingIndex == 2) then
            chi = 0.2_dp * mass(p1(reacIndx)) / (mass(p1(reacIndx))+mass(p2(reacIndx)))
        else
            write(*,*) "FREDON 2021 METHOD FOR CHEMICAL DESORPTION IS NOT VALID FOR DESORBINDEX > 2"
            stop
        end if
    end if
    if (chi * deltaEnthalpy - bindingEnergyDesorbingSpec < 0.0_dp) then
        desorptionFractionFullCoverage = 0.0_dp
    else
        desorptionFractionFullCoverage = 0.5_dp*(1.0_dp-exp(-(chi * deltaEnthalpy - &
            bindingEnergyDesorbingSpec)/(3.0_dp*bindingEnergyDesorbingSpec)))
    end if
  end function getDesorptionFractionFullCoverage

! ---------------------------------------------------------------------
! Get ice-coverage-dependent desorption fraction
! Interpolates between bare grain and full ice coverage
! ---------------------------------------------------------------------
  function getDesorptionFractionIncludingIce(reacIndx, numMonolayers) result(desorptionFraction)
    integer, intent(in) :: reacIndx
    real(dp), intent(in) :: numMonolayers
    real(dp) :: desorptionFraction


    real(dp) :: desorptionFractionBare
    real(dp) :: desorptionFractionFullCoverage

    desorptionFractionBare = desorptionFractionsBare(reacIndx)
    if (.NOT. GRAINS_HAVE_ICE) then
        ! If we do not simulate with ice, return bare grain desorption efficiency
        desorptionFraction = desorptionFractionBare
        return
    end if

    desorptionFractionFullCoverage = desorptionFractionsFullCoverage(reacIndx)

    desorptionFraction = linear_interp(MIN(1.0_dp, numMonolayers), &
        desorptionFractionBare, desorptionFractionFullCoverage)
  end function getDesorptionFractionIncludingIce

! ---------------------------------------------------------------------
! Get number of monolayers from abundance
! ---------------------------------------------------------------------
  pure function getNumberMonolayers(abundance) result(numberMonolayers)
    real(dp), intent(in) :: abundance
    real(dp) :: numberMonolayers

    if (.NOT. GRAINS_HAVE_ICE) then
        numberMonolayers = 0.0_dp
        return
    end if

    numberMonolayers = abundance * GAS_DUST_DENSITY_RATIO / NUM_SITES_PER_GRAIN
  end function getNumberMonolayers

! ---------------------------------------------------------------------
! Update diffusion and desorption rates with TST or Hasegawa-Herbst
! ---------------------------------------------------------------------
  pure subroutine updateVdiffAndVdes(dustTemp, nIce, vdiff, vdes)
    real(dp), intent(in) :: dustTemp
    integer, intent(in) :: nIce
    real(dp), intent(out) :: vdiff(nIce), vdes(nIce)

    ! inertiaProducts are scaled up by 1e50 in Makerates to avoid numerical issues
    ! We need to scale them down again here
    real(dp), parameter :: scaleFactor = 1e-50_dp

    real(dp) :: estimatedInertiaProduct
    integer :: i, j

    if (.NOT. useTSTprefactors) then
        do i=1,nIce
            j = iceList(i)
            ! Original treatment by Hasegawa et al, 1992
            vdes(i) = sqrt(HH_VDES_PREFACTOR*bindingEnergy(i)*mass(j))
            vdiff(i) = vdes(i)
        end do
    else
        ! TST treatment - use dust temperature for desorption
        vdes(:) = TST_VDES_PREFACTOR * dustTemp * dustTemp

        do  i=1,nIce
            j=iceList(i)
            if (atomCounts(j) == 1) then
                ! Atomic species, no rotational partition function
                vdes(i) = vdes(i) * mass(j)
            else if (.not. is_equal(inertiaProducts(i), MISSING_VALUE_FLOAT)) then
                ! Custom supplied 1/sigma*sqrt(Ix*Iy*Iz)
                if (moleculeIsLinear(i)) then
                    ! Linear molecule (H2, OH, CO2, etc)
                    vdes(i) = vdes(i) * mass(j)*scaleFactor * &
                        sqrt(PI) / (HP_SI**2)*(8.0_dp*PI**2*K_BOLTZ_SI*dustTemp)*&
                        inertiaProducts(i)
                else
                    ! Nonlinear molecule
                    vdes(i) = vdes(i) * mass(j) *scaleFactor* &
                        sqrt(PI) / (HP_SI**3)*(8.0_dp*PI**2*K_BOLTZ_SI*dustTemp)**(3.0_dp/2.0_dp)*&
                        inertiaProducts(i)
                end if
            else
                ! No inertia data available - estimate for polyatomic molecules
                if (atomCounts(j) >= 3) then
                    ! Fitted function to estimate inertia product for nonlinear molecules
                    estimatedInertiaProduct = 2.35425621e-21_dp*exp(1.04448580e-01_dp*mass(j))
                    vdes(i) = vdes(i) * mass(j) *scaleFactor* &
                        sqrt(PI) / (HP_SI**3)*(8.0_dp*PI**2*K_BOLTZ_SI*dustTemp)** &
                        (3.0_dp/2.0_dp)*estimatedInertiaProduct
                else
                    ! Diatomic - use HH as fallback
                    vdes(i) = sqrt(HH_VDES_PREFACTOR*bindingEnergy(i)*mass(j))
                end if
            end if
        end do

        if (separateDiffAndDesorbPrefactor) then
            ! Under stationary adsorbate assumption, q^TS = q^RS, so vdiff = kB*T/h
            vdiff = K_BOLTZ_SI * dustTemp / HP_SI
        else
            vdiff = vdes
        end if

    end if


    if (useCustomPrefactors) then
        do i=LBOUND(iceList, 1), UBOUND(iceList, 1)
            if (.not. is_equal(customVdiff(i), MISSING_VALUE_FLOAT)) vdiff(i) = customVdiff(i)
            if (.not. is_equal(customVdes(i), MISSING_VALUE_FLOAT)) vdes(i) = customVdes(i)
        end do
    end if
  end subroutine updateVdiffAndVdes

! ---------------------------------------------------------------------
! Encounter Desorption for H and H2 on H2-covered surfaces
! Hincelin et al. 2015
! ---------------------------------------------------------------------
  pure function getEncounterDesorptionRateConstant(reacIndx,dustTemperature) result(encounterDesorptionRateConstant)
    real(dp), intent(in) :: dustTemperature
    integer, intent(in) :: reacIndx
    real(dp) :: encounterDesorptionRateConstant

    real(dp) :: meetProb,desorbProb,diffuseProb

    integer :: index1,index2,i

    ! Get position of reactant in grain array
    index1=re1(reacIndx)

    do i=lbound(iceList,1),ubound(iceList,1)
        if (iceList(i) == index1) index1 = i
        if (iceList(i) == ngh2) index2 = i
    end do

    ! Diffusion rate that species meets H2 on grain surface
    meetProb = vdiff(index1)*exp(-getDiffusionBarrier(index1)/dustTemperature)
    meetProb = meetProb + (vdiff(index2)*exp(-getDiffusionBarrier(index2)/dustTemperature))

    ! Adjust for energy required to move from H2O onto H2
    if (.not. is_equal(EDEndothermicityFactor, 0.0_dp)) then
        meetProb = meetProb &
        & * exp(-EDEndothermicityFactor*(bindingEnergy(index1)-H2_ON_H2_BINDING_ENERGY)/dustTemperature)
    end if

    ! RateConstant of diffusion of index1 specie off of H2
    diffuseProb = vdiff(index1)*exp(-diffToBindRatio*H2_ON_H2_BINDING_ENERGY/dustTemperature)

    ! RateConstant of desorption of index1 specie off of H2
    desorbProb = vdiff(index1)*exp(-H2_ON_H2_BINDING_ENERGY/dustTemperature)

    ! Probability of desorbing when on top of H2
    desorbProb = desorbProb / (desorbProb + diffuseProb)

    ! Actual rate of EncounterDesorption mechanism
    encounterDesorptionRateConstant = 0.5_dp*desorbProb* meetProb*GAS_DUST_DENSITY_RATIO/NUM_SITES_PER_GRAIN
  end function getEncounterDesorptionRateConstant

end module SurfaceReactions
