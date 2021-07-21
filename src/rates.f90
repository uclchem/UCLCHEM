SUBROUTINE calculateReactionRates
    INTEGER:: idx1,idx2,k
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
    rate(idx1:idx2)=freezeOutRate(idx1,idx2)
    !freeze out rate uses thermal velocity but mass of E is 0 giving us infinite rates
    !just assume it's same as H
    rate(nR_EFreeze)=rate(nR_HFreeze)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !The below desorption mechanisms are from Roberts et al. 2007 MNRAS with
    !the addition of direct UV photodesorption. DESOH2,DESCR1,DEUVCR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Desorption due to energy released by H2 Formations
    idx1=desoh2Reacs(1)
    idx2=desoh2Reacs(2)
    IF (desorb .eq. 1 .and. h2desorb .eq. 1 .and. gama(j) .le. ebmaxh2) THEN
        !Epsilon is efficieny of this process, number of molecules removed per event
        !h2form is formation rate of h2, dependent on hydrogen abundance. 
        rate(idx1:idx2) = epsilon*h2FormEfficiency(gasTemp(dstep),dustTemp(dstep))
    ELSE
        rate(idx1:idx2) = 0.0
    ENDIF
    !Desorption due to energy from cosmic rays
    idx1=descrReacs(1)
    idx2=descrReacs(2)
    IF ((desorb .eq. 1) .and. (crdesorb .eq. 1) &
        &.and. (gama(j) .le. ebmaxcr)) THEN
        !4*pi*zeta = total CR flux. 1.64d-4 is iron to proton ratio of CR
        !as iron nuclei are main cause of CR heating.
        !GRAIN_SURFACEAREA_PER_H is the total surfaace area per hydrogen atom. ie total grain area per cubic cm when multiplied by density.
        !phi is efficieny of this reaction, number of molecules removed per event.
        rate(idx1:idx2) = 4.0*pi*zeta*1.64d-4*(GRAIN_SURFACEAREA_PER_H)*phi
    ELSE
        rate(idx1:idx2) = 0.0
    ENDIF

    !Desorption due to UV, partially from ISRF and partially from CR creating photons
    idx1=deuvcrReacs(1)
    idx2=deuvcrReacs(2)
    IF ((desorb .eq. 1) .and. (uvdesorb .eq. 1) .and. (gama(j) .le. ebmaxuvcr)) THEN
        !4.875d3 = photon flux, Checchi-Pestellini & Aiello (1992) via Roberts et al. (2007)
        !UVY is yield per photon.
        rate(idx1:idx2) = GRAIN_CROSSSECTION_PER_H*uv_yield*4.875d3*zeta
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
    IF (thermdesorb .eq.1) THEN
        DO j=idx1,idx2
            !then try to overwrite with position in grain array
            DO i=lbound(iceList,1),ubound(iceList,1)
                !See Cuppen, Walsh et al. 2017 review (section 4.1)
                IF (iceList(i) .eq. re1(j)) THEN
                    !Basic rate at which thermal desorption occurs
                    rate(j)=vdiff(i)*exp(-gama(j)/gasTemp(dstep))
                    !At some point, rate is so fast that there's no point freezing out any more
                    !Save some numerical trouble by setting freeze out to zero.
                    IF (rate(j) .ge. 1.0) THEN
                        do k=freezeReacs(1),freezeReacs(2)
                            IF (p1(j)==re1(k)) rate(k)=0.0
                        END DO
                    END IF
                    !factor of 2.0 adjusts for fact only top two monolayers (Eq 8)
                    !becayse GRAIN_SURFACEAREA_PER_H is per H nuclei, multiplying it by density gives area/cm-3
                    !that is roughly sigma_g.n_g from cuppen et al. 2017 but using surface instead of cross-sectional
                    !area seems more correct for this process.
                    IF (.NOT. THREE_PHASE) rate(j)=rate(j)*2.0*SURFACE_SITE_DENSITY*GRAIN_SURFACEAREA_PER_H
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

    if ((gasTemp(dstep) .lt. MAX_GRAIN_TEMP)) THEN
        DO j=idx1,idx2
            rate(j)=diffusionReactionRate(j,gasTemp(dstep))
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
    ELSE
        rate(idx1:idx2)=0.0
        rate(lhdesReacs(1):lhdesReacs(2))=0.0
    END IF


    !Account for Eley-Rideal reactions in a similar way.
    !First calculate overall rate and then split between desorption and sticking
    idx1=erReacs(1)
    idx2=erReacs(2)
    rate(idx1:idx2)=freezeOutRate(idx1,idx2)
    rate(idx1:idx2)=rate(idx1:idx2)*dexp(-gama(idx1:idx2)/gasTemp(dstep))
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
        !rate(nR_H2Form_LH)=0.0
        rate(nR_H2Form_ER)=0.0
        !rate(nR_H2Form_LHDes)=0.0
        rate(nR_H2Form_ERDes)=0.0
    ELSE
        rate(nR_H2Form_CT)= 0.0
    END IF

    CALL bulkSurfaceExchangeReactions(rate,gasTemp(dstep))

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
    rate(887)=0.0
END SUBROUTINE calculateReactionRates



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Freeze out determined by rate of collisions with grain
!No sticking coefficient is used because typical values are >0.95 below 150 K
! eg Le Bourlot et al. 2013, Molpeceres et al. 2020
!Above 150 K, thermal desorption will completely remove grain species
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION freezeOutRate(idx1,idx2) RESULT(freezeRates)
    REAL(dp) :: freezeRates(idx2-idx1+1)
    INTEGER :: idx1,idx2
    
    !additional factor for ions (beta=0 for neutrals)
    freezeRates=1.0+beta(idx1:idx2)*16.71d-4/(GRAIN_RADIUS*gasTemp(dstep))
    IF (fr .eq. 0.0 .or. gasTemp(dstep) .gt. 50) then
        freezeRates=0.0
    ELSE
        freezeRates=freezeRates*alpha(idx1:idx2)*THERMAL_VEL*dsqrt(gasTemp(dstep)/mass(re1(idx1:idx2)))*GRAIN_CROSSSECTION_PER_H
    END IF

END FUNCTION freezeOutRate

