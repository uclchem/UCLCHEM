SUBROUTINE calculateReactionRates
    INTEGER:: idx1,idx2
    REAL(dp) ::h2Efficiency
    !Calculate all reaction rates
    !Assuming the user has temperature changes or uses the desorption features of phase 1,
    !these need to be recalculated every time step.


    idx1=crpReacs(1)
    idx2=crpReacs(2)
    rate(idx1:idx2)=alpha(idx1:idx2)*zeta

    !UV photons, radfield has (factor of 1.7 conversion from habing to Draine)
    idx1=photonReacs(1)
    idx2=photonReacs(2)
    rate(idx1:idx2) = alpha(idx1:idx2)*dexp(-gama(idx1:idx2)*av(dstep))*radfield/1.7

    !Reactions involving cosmic ray induced photon
    idx1=crphotReacs(1)
    idx2=crphotReacs(2)
    rate(idx1:idx2)=alpha(idx1:idx2)*gama(idx1:idx2)*1.0/(1.0-omega)*zeta*(gasTemp(dstep)/300)**beta(idx1:idx2)


    !freeze out only happens if fr>0 and depending on evap choice 
    idx1=freezeReacs(1)
    idx2=freezeReacs(2)
    cion=1.0+16.71d-4/(GRAIN_RADIUS*gasTemp(dstep))
    IF (fr .eq. 0.0 .or. gasTemp(dstep) .gt. 30.0) then
        rate(idx1:idx2)=0.0
    ELSE
        DO j=idx1,idx2

            IF (re1(j).eq.nelec) THEN
                rate(j)=THERMAL_VEL*alpha(j)*GRAIN_CROSSSECTION_PER_H*fr*cion
            ELSE
                !taken from Rawlings et al. 1992
                rate(j)=alpha(j)*THERMAL_VEL*dsqrt(gasTemp(dstep)/mass(re1(j)))*GRAIN_CROSSSECTION_PER_H*fr
                IF (beta(j).ne.0.0 ) rate(j)=rate(j)*cion
            END IF
        END DO
    END IF



    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !The below desorption mechanisms are from Roberts et al. 2007 MNRAS with
    !the addition of direct UV photodesorption. DESOH2,DESCR1,DEUVCR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Desorption due to energy released by H2 Formations
    idx1=desoh2Reacs(1)
    idx2=desoh2Reacs(2)
    IF (desorb .eq. 1 .and. h2desorb .eq. 1 .and. gama(j) .le. ebmaxh2 .and.&
    &  mantle(dstep) .ge. 1.0d-30) THEN
        !Epsilon is efficieny of this process, number of molecules removed per event
        !h2form is formation rate of h2, dependent on hydrogen abundance. 
        rate(idx1:idx2) = epsilon*h2FormEfficiency(gasTemp(dstep),dustTemp(dstep))*abund(nh,dstep)/mantle(dstep)
    ELSE
        rate(idx1:idx2) = 0.0
    ENDIF
    !Desorption due to energy from cosmic rays
    idx1=descrReacs(1)
    idx2=descrReacs(2)
    IF (desorb .eq. 1 .and. crdesorb .eq. 1 .and.mantle(dstep).ge. 1d-30&
        &.and. gama(j) .le. ebmaxcr) THEN
        !4*pi*zeta = total CR flux. 1.64d-4 is iron to proton ratio of CR
        !as iron nuclei are main cause of CR heating.
        !GRAIN_SURFACEAREA_PER_H is the total surfaace area per hydrogen atom. ie total grain area per cubic cm when multiplied by density.
        !phi is efficieny of this reaction, number of molecules removed per event.
        rate(idx1:idx2) = 4.0*pi*zeta*1.64d-4*(GRAIN_SURFACEAREA_PER_H)*&
              &(1.0/mantle(dstep))*phi
    ELSE
        rate(idx1:idx2) = 0.0
    ENDIF

    !Desorption due to UV, partially from ISRF and partially from CR creating photons
    idx1=deuvcrReacs(1)
    idx2=deuvcrReacs(2)
    IF (desorb .eq. 1 .and. uvdesorb .eq. 1 .and. gama(j) .le. ebmaxuvcr &
        &.and. mantle(dstep) .ge. 1.0d-20) THEN
        !4.875d3 = photon flux, Checchi-Pestellini & Aiello (1992) via Roberts et al. (2007)
        !UVY is yield per photon.
        rate(idx1:idx2) = GRAIN_CROSSSECTION_PER_H*uv_yield*4.875d3*zeta*(1.0/mantle(dstep))
        !additional factor accounting for UV desorption from ISRF. UVCREFF is ratio of 
        !CR induced UV to ISRF UV.
        rate(idx1:idx2) = rate(idx1:idx2) * (1+(radfield/uvcreff)*(1.0/zeta)*dexp(-1.8*av(dstep)))
    ELSE
        rate(idx1:idx2) = 0.0
    ENDIF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Continuous Thermal Desorption. Reactions can be generated through a flag in Makerates
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    idx1=thermReacs(1)
    idx2=thermReacs(2)
    IF (thermdesorb .eq.1 .and. mantle(dstep) .ge. 1.0d-20) THEN
        DO j=idx1,idx2
            !then try to overwrite with position in grain array
            DO i=lbound(grainList,1),ubound(grainList,1)
                !See Cuppen, Walsh et al. 2017 review (section 4.1)
                IF (grainList(i) .eq. re1(j)) THEN
                    !Basic rate at which thermal desorption occurs
                    rate(j)=vdiff(i)*exp(-gama(j)/gasTemp(dstep))
                    !factor of 2.0 adjusts for fact only top two monolayers (Eq 8)
                    !becayse GRAIN_SURFACEAREA_PER_H is per H nuclei, multiplying it by density gives area/cm-3
                    !that is roughly sigma_g.n_g from cuppen et al. 2017 but using surface instead of cross-sectional
                    !area seems more correct for this process.
                    rate(j)=rate(j)*(2.0/mantle(dstep))*SURFACE_SITE_DENSITY*GRAIN_SURFACEAREA_PER_H*density(dstep)
                END IF
            END DO
        END DO
    ELSE
        rate(idx1:idx2)=0.0
    END IF


    !Reactions on surface can be treated considering diffusion of reactants
    !as in Langmuir-Hinshelwood mechanism
    !See work of David Quenard 2017 Arxiv:1711.05184
    !First calculate rate of the diffusion reaction
    idx1=lhReacs(1)
    idx2=lhReacs(2)
    DO j=idx1,idx2
        rate(j)=diffusionReactionRate(j)
    END DO
    !two routes for every diffusion reaction: products to gas or products remain on surface
    rate(lhdesReacs(1):lhdesReacs(2))=rate(idx1:idx2)

    !calculate fraction of reaction that goes down desorption route
    idx1=lhdesReacs(1)
    idx2=lhdesReacs(2)
    DO j=idx1,idx2
        rate(j)=desorptionFraction(j)*rate(j)
    END DO
    !remove that fraction from total rate of the diffusion route
    rate(lhReacs(1):lhReacs(2))=rate(lhReacs(1):lhReacs(2))-rate(idx1:idx2)


    !Account for Eley-Rideal reactions in a similar way.
    !First calculate overall rate and then split between desorption and sticking
    idx1=erReacs(1)
    idx2=erReacs(2)
    IF (mantle(dstep) .gt. 1e-20) THEN
        rate(idx1:idx2)=alpha(idx1:idx2)*THERMAL_VEL*dsqrt(gasTemp(dstep)/mass(re1(idx1:idx2)))*GRAIN_CROSSSECTION_PER_H*fr
        rate(idx1:idx2)=rate(idx1:idx2)*dexp(-gama(idx1:idx2)/gasTemp(dstep))/mantle(dstep)
    END IF
    rate(erdesReacs(1):erdesReacs(2))=rate(idx1:idx2)
    !calculate fraction of reaction that goes down desorption route
    idx1=erdesReacs(1)
    idx2=erdesReacs(2)
    DO j=idx1,idx2
        rate(j)=desorptionFraction(j)*rate(j)
    END DO
    !remove that fraction from total rate of the diffusion route
    rate(erReacs(1):erReacs(2))=rate(erReacs(1):erReacs(2))-rate(idx1:idx2)


    IF (PARAMETERIZE_H2FORM) THEN
        rate(nR_H2Form_CT)=h2FormEfficiency(gasTemp(dstep),dustTemp(dstep))
        ! rate(nR_H2Form_LH)=0.0
         rate(nR_H2Form_ER)=0.0
        ! rate(nR_H2Form_LHDes)=0.0
         rate(nR_H2Form_ERDes)=0.0
    ELSE
        rate(nR_H2Form_CT)= 0.0
    END IF


    !Basic gas phase reactions 
    !They only change if temperature has so we can save time with an if statement
    idx1=twobodyReacs(1)
    idx2=twobodyReacs(2)
    IF (lastTemp .ne. gasTemp(dstep)) THEN
        rate(idx1:idx2) = alpha(idx1:idx2)*((gasTemp(dstep)/300.)**beta(idx1:idx2))*dexp(-gama(idx1:idx2)/gasTemp(dstep)) 
    END IF


    lastTemp=gasTemp(dstep)
    IF (duplicates(1) .ne. 9999) THEN
        !this multiplies rate by 0 or 1 depending on whether gastemp>mintemp of a reaction
        rate(duplicates)=rate(duplicates)*min(real(floor(gasTemp(dstep)/minTemps)),1.0)
        !and this multiplies by 0,1 if gastemp>max temp
        rate(duplicates)=rate(duplicates)*min(real(floor(maxTemps/gasTemp(dstep))),1.0)
    END IF  

    !Reactions for which we have a more detailed photorecation treatment

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
double precision FUNCTION diffusionReactionRate(reacIndx)
    double precision :: diffuseRate,activationBarrier,reducedMass,tunnelProb
    double precision :: diffuseProb,desorbProb,reacProb,n_dust
    integer :: ns,index1,index2,reacIndx
    character(len=10) :: reactant1,reactant2


    !want position of species in the grain array but gas phase species aren't in there
    !so store species index
    index1=re1(reacIndx)
    index2=re2(reacIndx)

    !then try to overwrite with position in grain array
    DO i=lbound(grainList,1),ubound(grainList,1)
        IF (grainList(i) .eq. index1) index1 = i
        IF (grainList(i) .eq. index2) index2 = i
    END DO



    !Hasegawa 1992 diffusion rate. Rate that two species diffuse and meet on grain surface
    diffuseProb = vdiff(index1)*dexp(-0.5*bindingEnergy(index1)/gasTemp(dstep))
    diffuseProb = diffuseProb+ (vdiff(index2)*dexp(-0.5*bindingEnergy(index2)/gasTemp(dstep)))

    !probability a reactant will just desorb
    desorbProb = vdiff(index1)*dexp(-bindingEnergy(index1)/gasTemp(dstep))
    desorbProb = desorbProb + vdiff(index2)*dexp(-bindingEnergy(index2)/gasTemp(dstep)) 

    !Calculate classical activation energy barrier exponent
    reacProb = gama(reacIndx)/gasTemp(dstep)
    !Calculate quantum activation energy barrier exponent
    reducedMass = mass(grainList(index1)) * mass(grainList(index2)) / (mass(grainList(index1)) + mass(grainList(index2)))
    tunnelProb = 2.0d0 *CHEMICAL_BARRIER_THICKNESS/REDUCED_PLANCK * dsqrt(2.0d0*AMU*reducedMass*K_BOLTZ*gama(reacIndx))

    !Choose fastest between classical and tunnelling
    IF (reacProb.GT.tunnelProb) reacProb=tunnelProb
    !set activationBarrier to probability of reaction Ruaud+2016
    reacProb=dexp(-reacProb)

    !Overall reaction probability is chance of reaction occuring on meeting * diffusion rate
    reacProb = max(vdiff(index1),vdiff(index2)) * reacProb       


    ! Keff from Garrod & Pauly 2011 and Ruaud+2016
    ! Actual reaction probability is Preac/(Preac+Pevap+Pdiffuse), accounting for the other possible processes
    IF(DIFFUSE_REACT_COMPETITION) THEN
       reacProb = reacProb/(reacProb + desorbProb + diffuseProb)
    END IF
    
    !see Eq A1 of Quenard et al. 2018
    !GAS_DUST_DENSITY_RATIO = 1/ frac abundance of dust so ratio/density is 1/number density of dust
    n_dust=density(dstep)/GAS_DUST_DENSITY_RATIO
    diffusionReactionRate=alpha(reacIndx) *reacProb* diffuseProb/(NUM_SITES_PER_GRAIN*n_dust)
END FUNCTION diffusionReactionRate

! ---------------------------------------------------------------------
!  Chemical Reactive Desorption (CRD)
! David Quenard 2017 Arxiv:1711.05184
! From Minissalle+ 2016 and Vasyunin+ 2016
! ---------------------------------------------------------------------
double precision FUNCTION desorptionFraction(reacIndx)
    integer :: reacIndx,reactIndex1,reactIndex2,degreesOfFreedom
    integer :: productIndex(4)
    LOGICAL :: isDiff

    double precision :: deltaEnthalpy,maxBindingEnergy,epsilonCd,productEnthalpy
    double precision, parameter :: EFFECTIVE_SURFACE_MASS = 120.0

    
    !Get indices of grain surface version of products products 
    productIndex=0
    !Arrays like binding energy and formation enthalpy are indexed by position in grainList
    !rather than species list. Need to find position in grain list where every reactant and product appears
    !bearing in mind that Eley-Rideal reactions can have reactants in gas phase and CHEMDES has products in gas
    DO i=lbound(grainList,1),ubound(grainList,1)
        !check grain lists for reactants
        IF (grainList(i) .eq. re1(reacIndx)) reactIndex1 = i
        IF (gasGrainList(i) .eq. re1(reacIndx)) reactIndex1 = i
        !check equivalent gas list in case of ER reaction.
        IF (grainList(i) .eq. re2(reacIndx)) reactIndex2 = i
        IF (gasGrainList(i) .eq. re2(reacIndx)) reactIndex2 = i

        IF (grainList(i) .eq. p1(reacIndx)) productIndex(1) = i
        IF (grainList(i) .eq. p2(reacIndx)) productIndex(2) = i
        IF (grainList(i) .eq. p3(reacIndx)) productIndex(3) = i
        IF (grainList(i) .eq. p4(reacIndx)) productIndex(4) = i

        IF (gasGrainList(i) .eq. p1(reacIndx)) productIndex(1) = i
        IF (gasGrainList(i) .eq. p2(reacIndx)) productIndex(2) = i
        IF (gasGrainList(i) .eq. p3(reacIndx)) productIndex(3) = i
        IF (gasGrainList(i) .eq. p4(reacIndx)) productIndex(4) = i
    END DO

    maxBindingEnergy=0.0
    productEnthalpy=0.0
    DO i=1,4
        maxBindingEnergy=MAX(maxBindingEnergy,bindingEnergy(productIndex(i)))
        productEnthalpy=productEnthalpy+formationEnthalpy(productIndex(i))
    END DO

    !epsilonCd is the fraction of kinetic energy kept my the product when it collides with grain surface
    epsilonCd = mass(productIndex(1)) + mass(productIndex(2)) + mass(productIndex(3)) + mass(productIndex(4))
    epsilonCd = ((epsilonCd - EFFECTIVE_SURFACE_MASS) / (epsilonCd + EFFECTIVE_SURFACE_MASS))**2
    
    !Now calculate the change in enthalpy of the reaction.
    deltaEnthalpy= formationEnthalpy(reactIndex1)+formationEnthalpy(reactIndex2)-productEnthalpy
    
    !Convert from kcal to J, from J to K and from moles-1 to reactions-1
    deltaEnthalpy = deltaEnthalpy*4.184d03/(1.38054D-23*6.02214129d23)
    ! Total energy change includes activation energy of the reaction !
    deltaEnthalpy = deltaEnthalpy + gama(reacIndx)


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
        if (re1(reacIndx).eq.ngn.and.re2(reacIndx).eq.ngn) desorptionFraction = 0.5
        if ((re1(reacIndx).eq.ngo.and.re2(reacIndx).eq.nh) &
            &.or. (re1(reacIndx).eq. nh.and.re2(reacIndx).eq.ngo)) desorptionFraction = 0.3
        if ((re1(reacIndx).eq.ngoh.and.re2(reacIndx).eq.nh) &
            &.or. (re1(reacIndx).eq.nh.and.re2(reacIndx).eq.ngoh)) desorptionFraction = 0.25
    ENDIF
END FUNCTION desorptionFraction