SUBROUTINE calculateReactionRates
!Assuming the user has temperature changes or uses the desorption features of phase 1, these need working out on a timestep by time step basis
    cion=1.0+16.71d-4/(grainRadius*gasTemp(dstep))
    DO j=1,nreac
        !This case structure looks at the reaction type. species-species happens in default.
        !Other cases are special reactions, particularly desorption events (photons, CRs etc)
        SELECT CASE (reacType(j))
        !Cosmic ray reactions            
        CASE ('CRP')
            rate(j) = alpha(j)*zeta
        !UV photons, radfield has (factor of 1.7 conversion from habing to Draine)
        CASE ('PHOTON')
            rate(j) = alpha(j)*dexp(-gama(j)*av(dstep))*radfield/1.7
            !co photodissoction number is stored as nrco
            IF (re1(j).eq.nco) THEN
                IF(p1(j).eq.no .and. p2(j).eq.nc) nrco=j
                IF(p1(j).eq.nc .and. p2(j).eq.no) nrco=j
            ENDIF
        !cosmic ray induced photon
        CASE ('CRPHOT')
            rate(j)=alpha(j)*gama(j)*1.0/(1.0-omega)*zeta*(gasTemp(dstep)/300)**beta(j)
        !freeze out only happens if fr>0 and depending on evap choice 
        !freeze out only happens if fr>0 and depending on evap choice 
        CASE ('FREEZE')             
            IF (fr .eq. 0.0 .or. gasTemp(dstep) .gt. 30.0) then
                rate(j)=0.0
            ELSE
                IF (re1(j).eq.nelec) THEN
                    rate(j)=THERMAL_VEL*alpha(j)*GRAIN_AREA_PER_H*fr*cion
                ELSE
                    !taken from Rawlings et al. 1992
                    !0.25 factor converts surface area of a spherical grain to cross-sectional area
                    !Saves recalculating pi.r^2/n_H when we have 4pi.r^2/n_H for other parts of the code
                    rate(j)=alpha(j)*THERMAL_VEL*dsqrt(gasTemp(dstep)/mass(re1(j)))*0.25*GRAIN_AREA_PER_H*fr
                    IF (beta(j).ne.0.0 ) rate(j)=rate(j)*cion
                END IF
            ENDIF
        !The below desorption mechanisms are from Roberts et al. 2007 MNRAS with
        !the addition of direct UV photodesorption. DESOH2,DESCR1,DEUVCR
        CASE ('DESOH2')
            IF (desorb .eq. 1 .and. h2desorb .eq. 1&
            & .and. gama(j) .le. ebmaxh2 .and.&
            &  mantle(dstep) .ge. 1.0d-30) THEN
                !Epsilon is efficieny of this process, number of molecules removed per event
                !h2form is formation rate of h2, dependent on hydrogen abundance. 
                rate(j) = epsilon*h2form*abund(nh,dstep)*1.0/mantle(dstep)
            ELSE
                rate(j) = 0.0
            ENDIF
        CASE ('DESCR')
            IF (desorb .eq. 1 .and. crdesorb .eq. 1&
            &.and.mantle(dstep).ge. 1d-30&
            &.and. gama(j) .le. ebmaxcr) THEN
                !4*pi*zeta = total CR flux. 1.64d-4 is iron to proton ratio of CR
                !as iron nuclei are main cause of CR heating.
                !GRAIN_AREA is the total area per hydrogen atom. ie total grain area per cubic cm when multiplied by density.
                !phi is efficieny of this reaction, number of molecules removed per event.
                rate(j) = 4.0*pi*zeta*1.64d-4*(GRAIN_AREA_PER_H)*&
                      &(1.0/mantle(dstep))*phi
            ELSE
                rate(j) = 0.0
            ENDIF
        CASE ('DEUVCR')
            IF (desorb .eq. 1 .and. uvcr .eq. 1 &
             &.and. gama(j) .le. ebmaxuvcr .and. mantle(dstep) .ge. 1.0d-30) THEN
                !4.875d3 = photon flux, Checchi-Pestellini & Aiello (1992) via Roberts et al. (2007)
                !UVY is yield per photon.
                rate(j) = 0.25*GRAIN_AREA_PER_H*uv_yield*4.875d3*zeta*(1.0/mantle(dstep))
                !additional factor accounting for UV desorption from ISRF. UVCREFF is ratio of 
                !CR induced UV to ISRF UV.
                rate(j) = rate(j) * (1+(radfield/uvcreff)*(1.0/zeta)*dexp(-1.8*av(dstep)))
            ELSE
                rate(j) = 0.0
            ENDIF
        CASE('THERM')
            IF ((thermdesorb .eq.1) .and. (mantle(dstep) .gt. 1.0d-20)) THEN
                !then try to overwrite with position in grain array
                DO i=lbound(grainList,1),ubound(grainList,1)
                    !See Cuppen, Walsh et al. 2017 review (section 4.1)
                    IF (grainList(i) .eq. re1(j)) THEN
                        !Basic rate at which thermal desorption occurs
                        rate(j)=vdiff(i)*exp(-gama(j)/dustTemp(dstep))
                        !factor of 2.0 adjusts for fact only top two monolayers (Eq 8)
                        !becayse GRAIN_AREA_PER_H is per H nuclei, multiplying it by density gives area/cm-3
                        !that is roughly sigma_g.n_g from cuppen et al. 2017 but using surface instead of cross-sectional
                        !area seems more correct for this process.
                        rate(j)=rate(j)*(2.0/mantle(dstep))*SURFACE_SITE_DENSITY*GRAIN_AREA_PER_H*density(dstep)
                    END IF
                END DO
            ELSE
                rate(j)=0.0
            END IF

        !Reactions on surface can be treated considering diffusion of reactants
        !See work of David Quenard 2017 Arxiv:1711.05184
        !abstracted to functions below for ease of reading
        CASE ('DIFF')
            rate(j)=diffusionReactionRate()
        case ('CHEMDES')
            rate(j) = diffusionReactionRate()
        CASE DEFAULT
            rate(j) = alpha(j)*((gasTemp(dstep)/300.)**beta(j))*dexp(-gama(j)/gasTemp(dstep))
        END SELECT
        IF (ISNAN(rate(j))) write(*,*) "NAN RATE:", j,gasTemp(dstep)
    END DO

    !this catches the rates with a large negative gamma
    where(rate(duplicates).gt.HUGE(rate(duplicates))) rate(duplicates)=0.0d0
    !this multiplies rate by 0 or 1 depending on whether gastemp>mintemp of a reaction
    rate(duplicates)=rate(duplicates)*min(real(floor(gasTemp(dstep)/minTemps)),1.0)
    !and this multiplies by 0,1 if gastemp>max temp
    rate(duplicates)=rate(duplicates)*min(real(floor(maxTemps/gasTemp(dstep))),1.0)

    !Photoreactions for which we have a more detailed treatment
    h2dis=H2PhotoDissRate(h2Col,radField,av(dstep),turbVel) !H2 photodissociation
    rate(nrco)=COPhotoDissRate(h2Col,coCol,radField,av(dstep)) !CO photodissociation
    rate(nR_C_hv)=cIonizationRate(alpha(nR_C_hv),gama(nR_C_hv),gasTemp(dstep),ccol,h2col,av(dstep),radfield) !C photoionization

END SUBROUTINE calculateReactionRates


!----------------------------------------------------------------------------------------------------
!Reactions on the surface treated by evaluating diffusion rates across grains and accounting
!For competition with chemical desorption. Products remain bound ('DIFF') or are desorbed ('CHEMDES')
!Assuming Eb = 0.5 Ed. Units of s-1. 
!David Quenard 2017 Arxiv:1711.05184
!----------------------------------------------------------------------------------------------------
double precision FUNCTION diffusionReactionRate()
    double precision :: diffuseRate,activationBarrier,reducedMass,tunnelProb
    double precision :: diffuseProb,desorbProb
    integer :: index1,index2


    !want position of species in the grain array but gas phase species aren't in there
    !so store species index
    index1=re1(j)
    index2=re2(j)

    !then try to overwrite with position in grain array
    DO i=lbound(grainList,1),ubound(grainList,1)
        IF (grainList(i) .eq. index1) index1 = i
        IF (grainList(i) .eq. index2) index2 = i
    END DO

    !Hasegawa 1992 diffusion rate. Rate that two species diffuse and meet on grain surface
    diffuseRate = vdiff(index1)*dexp(-0.5*bindingEnergy(index1)/dustTemp(dstep))
    diffuseRate = (diffuseRate+ (vdiff(index2)*dexp(-0.5*bindingEnergy(index2)/dustTemp(dstep))))/NUM_SITES_PER_GRAIN

    !Calculate classical activation energy barrier exponent
    activationBarrier = gama(j)/gasTemp(dstep)

    !Calculate quantum activation energy barrier exponent
    reducedMass = mass(grainList(index1)) * mass(grainList(index2)) / (mass(grainList(index1)) + mass(grainList(index2)))
    tunnelProb = 2.0d0 *CHEMICAL_BARRIER_THICKNESS/REDUCED_PLANCK * dsqrt(2.0d0*AMU*reducedMass*K_BOLTZ*gama(j))

    !Choose fastest between classical and tunnelling
    IF (activationBarrier.GT.tunnelProb) activationBarrier=tunnelProb
    !set activationBarrier to probability of reaction Ruaud+2016
    activationBarrier=dexp(-activationBarrier)

    ! Keff from Garrod & Pauly 2011 and Ruaud+2016
    ! Actual reaction probability is Preac/(Preac+Pevap+Pdiffuse), accounting for the other possible processes
    IF(DIFFUSE_REACT_COMPETITION) THEN
       activationBarrier = max(vdiff(index1),vdiff(index2)) * activationBarrier       
       !probability a reactant will just desorb
       desorbProb = vdiff(index1)*dexp(-bindingEnergy(index1)/dustTemp(dstep))
       desorbProb = desorbProb + vdiff(index2)*dexp(-bindingEnergy(index2)/dustTemp(dstep)) 
       !Probability reactants wills diffuse onwards when they meet
       diffuseProb = (diffuseRate) * NUM_SITES_PER_GRAIN
       !Therefore, total probability the reactants that meet will react
       activationBarrier = activationBarrier/(activationBarrier + desorbProb + diffuseProb)
    END IF
    
    diffusionReactionRate=alpha(j) * diffuseRate * activationBarrier* GAS_DUST_DENSITY_RATIO / density(dstep)

    !Now adjust for fraction of this reaction's products that will desorb due to energy released
    IF (reacType(j).eq.'DIFF') THEN
        diffusionReactionRate = diffusionReactionRate * (1.0-desorptionFraction(j,index1,index2))
    ELSE IF(reacType(j).eq.'CHEMDES') THEN
        diffusionReactionRate = diffusionReactionRate * desorptionFraction(j,index1,index2)
    ENDIF
END FUNCTION diffusionReactionRate

! ---------------------------------------------------------------------
!  Chemical Reactive Desorption (CRD)
! David Quenard 2017 Arxiv:1711.05184
! From Minissalle+ 2016 and Vasyunin+ 2016
! ---------------------------------------------------------------------
double precision FUNCTION desorptionFraction(j,reactIndex1,reactIndex2)
    integer :: j,reactIndex1,reactIndex2,degreesOfFreedom
    integer :: productIndex(4)

    double precision :: deltaEnthalpy,maxBindingEnergy,epsilonCd,productEnthalpy
    double precision, parameter :: EFFECTIVE_SURFACE_MASS = 120.0

    
    !Get indices of grain surface version of products products 
    productIndex=0
    productIndex = 0.0
    maxBindingEnergy=0.0
    productEnthalpy=0.0

    IF (reacType(j).eq.'DIFF') THEN
        DO i=lbound(grainList,1),ubound(grainList,1)
            IF (grainList(i) .eq. p1(j)) productIndex(1)=grainList(i)
            !Go through grain list and try to find species in product list
            !If it has a binding energy larger than largest energy found so far, update maxBindingEnergy
            IF (grainList(i) .eq. p1(j)) THEN
                productIndex(1) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

            IF (grainList(i) .eq. p2(j)) THEN
                productIndex(2) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

            IF (grainList(i) .eq. p3(j)) THEN
                productIndex(3) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

            IF (grainList(i) .eq. p4(j)) THEN
                productIndex(4) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF

        END DO
    ELSE IF (reacType(j).eq.'CHEMDES') THEN
        DO i=lbound(gasGrainList,1),ubound(gasGrainList,1)
            IF (gasGrainList(i) .eq. p1(j)) THEN
                productIndex(1) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
            IF (gasGrainList(i) .eq. p2(j)) THEN
                productIndex(2) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
            IF (gasGrainList(i) .eq. p3(j)) THEN
                productIndex(3) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
            IF (gasGrainList(i) .eq. p4(j)) THEN
                productIndex(4) = grainList(i)
                productEnthalpy=productEnthalpy+formationEnthalpy(i)
                if (bindingEnergy(i) .ge. maxBindingEnergy) maxBindingEnergy=bindingEnergy(i)
            END IF
        END DO
    ENDIF
    
    !
    !epsilonCd is the fraction of kinetic energy kept my the product when it collides with grain surface
    epsilonCd = mass(productIndex(1)) + mass(productIndex(2)) + mass(productIndex(3)) + mass(productIndex(4))
    epsilonCd = ((epsilonCd - EFFECTIVE_SURFACE_MASS) / (epsilonCd + EFFECTIVE_SURFACE_MASS))**2
    
    !Now calculate the change in enthalpy of the reaction.
    deltaEnthalpy= formationEnthalpy(reactIndex1)+formationEnthalpy(reactIndex2)-productEnthalpy
    
    !Convert from kcal to J, from J to K and from moles-1 to reactions-1
    deltaEnthalpy = deltaEnthalpy*4.184d03/(1.38054D-23*6.02214129d23)
    ! Total energy change includes activation energy of the reaction !
    deltaEnthalpy = deltaEnthalpy + gama(j)


    IF (deltaEnthalpy.eq.0.00) deltaEnthalpy = 1e-30 

    !Degrees of freedom = 3 * number of atoms in the molecule
    degreesOfFreedom = atomCounts(productIndex(1))
    if (productIndex(2).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(2)))
    if (productIndex(3).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(3)))
    if (productIndex(4).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(4)))                    
    degreesOfFreedom = 3 * degreesOfFreedom
        
    desorptionFraction = dexp((-maxBindingEnergy*real(degreesOfFreedom)) / (epsilonCd * deltaEnthalpy))
    
   IF (deltaEnthalpy.lt.0.d0) THEN        !< If reaction is endothermic, no CRD
        desorptionFraction = 0.d0
    END IF
    
    IF (GRAINS_HAVE_ICE) THEN
        desorptionFraction = desorptionFraction/10    !< See Minisalle et al. 2016 for icy grain surface.
        ! Special case of OH+H, O+H, N+N on ices, see same paper
        if (re1(j).eq.ngn.and.re2(j).eq.ngn) desorptionFraction = 0.5
        if ((re1(j).eq.ngo.and.re2(j).eq.nh) .or. (re1(j).eq. nh.and.re2(j).eq.ngo)) desorptionFraction = 0.3
        if ((re1(j).eq.ngoh.and.re2(j).eq.nh) .or. (re1(j).eq.nh.and.re2(j).eq.ngoh)) desorptionFraction = 0.25
    ENDIF
END FUNCTION desorptionFraction



!=======================================================================
!
!  Calculate the rate of molecular hydrogen (H2) formation on grains
!  using the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222).
!
!-----------------------------------------------------------------------
FUNCTION h2FormRate(GAS_TEMPERATURE,GRAIN_TEMPERATURE) RESULT(RATE)
   IMPLICIT NONE

    REAL(dp) :: RATE
    REAL(dp), INTENT(IN) :: GAS_TEMPERATURE,GRAIN_TEMPERATURE

    REAL(dp) :: THERMAL_VELOCITY,STICKING_COEFFICIENT,TOTAL_CROSS_SECTION,CROSS_SECTION_SCALE
    REAL(dp) :: FLUX,FACTOR1,FACTOR2,EPSILON
    REAL(dp) :: SILICATE_FORMATION_EFFICIENCY,GRAPHITE_FORMATION_EFFICIENCY
    REAL(dp) :: SILICATE_CROSS_SECTION,SILICATE_MU,SILICATE_E_S,SILICATE_E_H2,SILICATE_E_HP,SILICATE_E_HC,SILICATE_NU_H2,SILICATE_NU_HC
    REAL(dp) :: GRAPHITE_CROSS_SECTION,GRAPHITE_MU,GRAPHITE_E_S,GRAPHITE_E_H2,GRAPHITE_E_HP,GRAPHITE_E_HC,GRAPHITE_NU_H2,GRAPHITE_NU_HC

    !  Mean thermal velocity of hydrogen atoms (cm s^-1)
    THERMAL_VELOCITY=1.45D5*SQRT(GAS_TEMPERATURE/1.0D2)

    !  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
    !  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
    STICKING_COEFFICIENT=1.0D0/(1.0D0+0.04D0*SQRT(GAS_TEMPERATURE+GRAIN_TEMPERATURE) &
                    & + 0.2D0*(GAS_TEMPERATURE/1.0D2)+0.08D0*(GAS_TEMPERATURE/1.0D2)**2)

    FLUX=1.0D-10 ! Flux of H atoms in monolayers per second (mLy s^-1)

    !Our cross-sectional area per H is different to Cazaux and Tielens so we scale theirs
    !by the fraction of our standard value. 
    ! CROSS_SECTION_SCALE=GRAIN_AREA_PER_H/1.660539E-021

    TOTAL_CROSS_SECTION=6.273D-22!*CROSS_SECTION_SCALE ! Total mixed grain cross section per H nucleus (cm^-2/nucleus)
    SILICATE_CROSS_SECTION=8.473D-22!*CROSS_SECTION_SCALE ! Silicate grain cross section per H nucleus (cm^-2/nucleus)
    GRAPHITE_CROSS_SECTION=7.908D-22!*CROSS_SECTION_SCALE ! Graphite grain cross section per H nucleus (cm^-2/nucleus)

    !  Silicate grain properties
    SILICATE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
    SILICATE_E_S=110.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
    SILICATE_E_H2=320.0D0 ! Desorption energy of H2 molecules (K)
    SILICATE_E_HP=450.0D0 ! Desorption energy of physisorbed H atoms (K)
    SILICATE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
    SILICATE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
    SILICATE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

    FACTOR1=SILICATE_MU*FLUX/(2*SILICATE_NU_H2*EXP(-SILICATE_E_H2/GRAIN_TEMPERATURE))

   FACTOR2=1.0D0*(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2 &
        & /4.0D0*EXP(-SILICATE_E_S/GRAIN_TEMPERATURE)

   EPSILON=1.0D0/(1.0D0+SILICATE_NU_HC/(2*FLUX)*EXP(-1.5*SILICATE_E_HC/GRAIN_TEMPERATURE) &
              & *(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2)

   SILICATE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!  Graphite grain properties
   GRAPHITE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
   GRAPHITE_E_S=260.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
   GRAPHITE_E_H2=520.0D0 ! Desorption energy of H2 molecules (K)
   GRAPHITE_E_HP=800.0D0 ! Desorption energy of physisorbed H atoms (K)
   GRAPHITE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
   GRAPHITE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
   GRAPHITE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

   FACTOR1=GRAPHITE_MU*FLUX/(2*GRAPHITE_NU_H2*EXP(-GRAPHITE_E_H2/GRAIN_TEMPERATURE))

   FACTOR2=1.0D0*(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2 &
        & /4.0D0*EXP(-GRAPHITE_E_S/GRAIN_TEMPERATURE)

   EPSILON=1.0D0/(1.0D0+GRAPHITE_NU_HC/(2*FLUX)*EXP(-1.5*GRAPHITE_E_HC/GRAIN_TEMPERATURE) &
              & *(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2)

   GRAPHITE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!!$!  Use the tradional rate, with a simple temperature dependence based on the
!!$!  thermal velocity of the H atoms in the gas and neglecting any temperature
!!$!  dependency of the formation and sticking efficiencies
!!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)

!!$!  Use the treatment of de Jong (1977, A&A, 55, 137, p140 right column).
!!$!  The second exponential dependence on the gas temperature reduces the
!!$!  efficiency at high temperatures and so prevents runaway H2 formation
!!$!  heating at high temperatures:
!!$!
!!$!  k_H2 = 3E-18 * T^0.5 * exp(-T/1000)   [cm3/s]
!!$!
!!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)*EXP(-(GAS_TEMPERATURE/1.0D3))

!!$!  Use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722, eqn 4)
!!$   RATE=0.5D0*THERMAL_VELOCITY*TOTAL_CROSS_SECTION*STICKING_COEFFICIENT

!  Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222)
   RATE=0.5D0*THERMAL_VELOCITY*(SILICATE_CROSS_SECTION*SILICATE_FORMATION_EFFICIENCY &
    & + GRAPHITE_CROSS_SECTION*GRAPHITE_FORMATION_EFFICIENCY)*STICKING_COEFFICIENT
!!$!  Use the expression given by Markus Rollig during the February 2012 Leiden workshop
!!$   RATE=0.5D0*THERMAL_VELOCITY &
!!$      & *(SILICATE_CROSS_SECTION/((1.0D0 + 6.998D24/EXP(1.5*SILICATE_E_HC/GRAIN_TEMPERATURE)) &
!!$      & *(1.0D0 + 1.0D0/(EXP(SILICATE_E_HP/GRAIN_TEMPERATURE) &
!!$      & *(0.427D0/EXP((SILICATE_E_HP-SILICATE_E_S)/GRAIN_TEMPERATURE) + 2.5336D-14*SQRT(GRAIN_TEMPERATURE))))) &
!!$      & + GRAPHITE_CROSS_SECTION/((1.0D0 + 4.610D24/EXP(1.5*GRAPHITE_E_HC/GRAIN_TEMPERATURE)) &
!!$      & *(1.0D0 + 1.0D0/(EXP(GRAPHITE_E_HP/GRAIN_TEMPERATURE) &
!!$      & *(0.539D0/EXP((GRAPHITE_E_HP-GRAPHITE_E_S)/GRAIN_TEMPERATURE) + 5.6334D-14*SQRT(GRAIN_TEMPERATURE)))))) &
!!$      & *STICKING_COEFFICIENT

   RETURN
END FUNCTION h2FormRate

!=======================================================================