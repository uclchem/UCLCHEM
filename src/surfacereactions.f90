MODULE SurfaceReactions
  USE constants
  USE network
  REAL(dp) :: bulkGain,bulkLoss,surfaceGain,surfaceLoss,totalSwap
  REAL(dp) :: safeMantle,safeBulk

  !Silicate grain properties for H2 Formation
  REAL(dp),PARAMETER :: SILICATE_MU=0.005D0 ! Fraction of newly formed H2 that stays on the grain surface
  REAL(dp),PARAMETER :: SILICATE_E_S=110.0D0 ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
  REAL(dp),PARAMETER :: SILICATE_E_H2=320.0D0 ! Desorption energy of H2 molecules (K)
  REAL(dp),PARAMETER :: SILICATE_E_HP=450.0D0 ! Desorption energy of physisorbed H atoms (K)
  REAL(dp),PARAMETER :: SILICATE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
  REAL(dp),PARAMETER :: SILICATE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
  REAL(dp),PARAMETER :: SILICATE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)
  REAL(dp),PARAMETER :: SILICATE_CROSS_SECTION=8.473D-22!*CROSS_SECTION_SCALE ! Silicate grain cross section per H nucleus (cm^-2/nucleus)

  !Graphite grain properties for H2 Formation
  REAL(dp),PARAMETER :: GRAPHITE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
  REAL(dp),PARAMETER :: GRAPHITE_E_S=260.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
  REAL(dp),PARAMETER :: GRAPHITE_E_H2=520.0D0 ! Desorption energy of H2 molecules (K)
  REAL(dp),PARAMETER :: GRAPHITE_E_HP=800.0D0 ! Desorption energy of physisorbed H atoms (K)
  REAL(dp),PARAMETER :: GRAPHITE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
  REAL(dp),PARAMETER :: GRAPHITE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
  REAL(dp),PARAMETER :: GRAPHITE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)
  REAL(dp),PARAMETER :: GRAPHITE_CROSS_SECTION=7.908D-22!*CROSS_SECTION_SCALE ! Graphite grain cross section per H nucleus (cm^-2/nucleus)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Grain surface parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp), PARAMETER :: GAS_DUST_MASS_RATIO=100.0,GRAIN_RADIUS=1.d-5, GRAIN_DENSITY = 3.0 ! Mass density of a dust grain
  REAL(dp), PARAMETER :: THERMAL_VEL= SQRT(8.0d0*K_BOLTZ/(PI*AMU)) !Thermal velocity without the factor of SQRT(T/m) where m is moelcular mass in amu

  !reciprocal of fractional abundance of dust grains (we only divide by number density so better to store reciprocal)
  REAL(dp), PARAMETER :: GAS_DUST_DENSITY_RATIO = (4.0*PI*(GRAIN_RADIUS**3)*GRAIN_DENSITY*GAS_DUST_MASS_RATIO)/(3.0 * AMU)
  !Grain area per h nuclei, values taken from Cazaux & Tielens 2004 via UCL-PDR to match H2 formation rate
  REAL(dp), PARAMETER :: GRAIN_CROSSSECTION_PER_H=0.5*(7.908D-22+8.473D-22)
  REAL(dp), PARAMETER :: GRAIN_SURFACEAREA_PER_H=4.0*GRAIN_CROSSSECTION_PER_H!2.0*4.0*PI*GRAIN_RADIUS*GRAIN_RADIUS/GAS_DUST_DENSITY_RATIO

  !Below are values for grain surface reactions
  LOGICAL, PARAMETER :: DIFFUSE_REACT_COMPETITION=.True., GRAINS_HAVE_ICE=.True.
  REAL(dp), PARAMETER :: CHEMICAL_BARRIER_THICKNESS = 1.40d-8  !gre Parameter used to compute the probability for a surface reaction with 
  !! activation energy to occur through quantum tunneling (Hasegawa et al. Eq 6 (1992).)
  REAL(dp), PARAMETER :: SURFACE_SITE_DENSITY = 1.5d15 ! site density on one grain [cm-2]
  REAL(dp), PARAMETER :: VDIFF_PREFACTOR=2.0*K_BOLTZ*SURFACE_SITE_DENSITY/PI/PI/AMU
  REAL(dp), PARAMETER :: NUM_SITES_PER_GRAIN = GRAIN_RADIUS*GRAIN_RADIUS*SURFACE_SITE_DENSITY*4.0*PI

  REAL(dp), ALLOCATABLE ::vdiff(:)
CONTAINS
  !=======================================================================
  !
  !  Calculate the rate of molecular hydrogen (H2) formation on grains
  !  using the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
  !  Cazaux & Tielens (2004, ApJ, 604, 222).
  !
  !-----------------------------------------------------------------------
  FUNCTION h2FormEfficiency(gasTemp,dustTemp) RESULT(rate)
    REAL(dp) :: rate
    REAL(dp), INTENT(IN) :: gasTemp,dustTemp

    REAL(dp) :: THERMAL_VELOCITY,STICKING_COEFFICIENT,CROSS_SECTION_SCALE
    REAL(dp) :: FLUX,FACTOR1,FACTOR2,EPSILON
    REAL(dp) :: SILICATE_FORMATION_EFFICIENCY,GRAPHITE_FORMATION_EFFICIENCY
    !  Mean thermal velocity of hydrogen atoms (cm s^-1)
    THERMAL_VELOCITY=1.45D5*SQRT(gasTemp/1.0D2)

    !  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
    !  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
    STICKING_COEFFICIENT=1.0D0/(1.0D0+0.04D0*SQRT(gasTemp+dustTemp) &
                    & + 0.2D0*(gasTemp/1.0D2)+0.08D0*(gasTemp/1.0D2)**2)

    FLUX=1.0D-10 ! Flux of H atoms in monolayers per second (mLy s^-1)

    FACTOR1=SILICATE_MU*FLUX/(2*SILICATE_NU_H2*EXP(-SILICATE_E_H2/dustTemp))

   FACTOR2=1.0D0*(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2 &
        & /4.0D0*EXP(-SILICATE_E_S/dustTemp)

   EPSILON=1.0D0/(1.0D0+SILICATE_NU_HC/(2*FLUX)*EXP(-1.5*SILICATE_E_HC/dustTemp) &
              & *(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2)

   SILICATE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON


   FACTOR1=GRAPHITE_MU*FLUX/(2*GRAPHITE_NU_H2*EXP(-GRAPHITE_E_H2/dustTemp))

   FACTOR2=1.0D0*(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2 &
        & /4.0D0*EXP(-GRAPHITE_E_S/dustTemp)

   EPSILON=1.0D0/(1.0D0+GRAPHITE_NU_HC/(2*FLUX)*EXP(-1.5*GRAPHITE_E_HC/dustTemp) &
              & *(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2)

   GRAPHITE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!  Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222)
   rate=0.5D0*THERMAL_VELOCITY*(SILICATE_CROSS_SECTION*SILICATE_FORMATION_EFFICIENCY &
    & + GRAPHITE_CROSS_SECTION*GRAPHITE_FORMATION_EFFICIENCY)*STICKING_COEFFICIENT

     RETURN
  END FUNCTION h2FormEfficiency


  FUNCTION bulkGainFromMantleBuildUp() RESULT(rate)
    REAL(dp) :: rate
    IF (safeMantle .lt. 1e-20) THEN
        rate = 0.0
    ELSE
        rate=0.5*GAS_DUST_DENSITY_RATIO/NUM_SITES_PER_GRAIN
    END IF
  END FUNCTION bulkGainFromMantleBuildUp


  FUNCTION bulkLossFromMantleLoss() RESULT(rate)
    REAL(dp) :: rate
    IF (safeBulk .lt. 1e-20) THEN
        rate = 0.0
    ELSE
        rate=1.0
    END IF
  END FUNCTION bulkLossFromMantleLoss


  FUNCTION surfaceToBulkSwappingRates(gasTemperature) RESULT(rate)
    REAL(dp) ::rate,gasTemperature
    IF ((safeMantle .lt. 1e-20) .or. (gasTemperature .gt. 100)) THEN
        rate = 0.0
    ELSE
        rate = 1.0
    END IF
  END FUNCTION surfaceToBulkSwappingRates


  SUBROUTINE bulkToSurfaceSwappingRates(rate,idx1,idx2,gasTemperature)
    REAL(dp), INTENT(INOUT) :: rate(*)
    REAL(dp) :: gasTemperature
    INTEGER :: idx1,idx2
    IF ((safeMantle .lt. 1e-20) .or. (gasTemperature .gt. 100)) THEN
        rate(idx1:idx2) = 0.0
    ELSE
        DO i=idx1,idx2
            DO j=lbound(iceList,1),ubound(iceList,1)
                IF (iceList(j) .eq. re1(i)) rate(i)=vdiff(j)*DEXP(-bindingEnergy(j)/gasTemperature)
            END DO
        END DO
    END IF
  END SUBROUTINE bulkToSurfaceSwappingRates
END MODULE SurfaceReactions