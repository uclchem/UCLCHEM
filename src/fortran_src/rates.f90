module RATES
    use constants
    use DEFAULTPARAMETERS
    use f2py_constants
    !f2py INTEGER, parameter :: dp
    use network
    use photoreactions, only: H2PhotoDissRate, COPhotoDissRate, cIonizationRate, ICE_GAS_PHOTO_CROSSSECTION_RATIO
    use physicscore
    use SurfaceReactions
    implicit none

    !Variables controlling chemistry:
    real(dp) :: grainArea,cion,h2dis,lastTemp=0.0

    ! Controlling ice chemistry
    real(dp), parameter :: h2StickingZero=0.87d0,hStickingZero=1.0d0, h2StickingTemp=87.0d0,hStickingTemp=52.0d0
    !Flags to control desorption processes
    real(dp) :: turbVel=1.0  !unit? km/s or cm/s
    ! TODO: integrate into makerates and put it in network.f90

    ! Pre-split LH and ER base rates, saved in calculateReactionRates and used
    ! by the F callback to re-apply the chemdes split at the current ice thickness.
    real(dp) :: rate_lh_unsplit(nreac) = 0.0d0
    real(dp) :: rate_er_unsplit(nreac) = 0.0d0

contains
    subroutine calculateReactionRates(abund, safemantle,  h2col, cocol, ccol, rate)
        real(dp), intent(in) :: abund(:, :), safemantle, h2col, cocol, ccol
        real(dp), intent(inout) :: rate(:)
        integer :: idx1,idx2
        real(dp) :: vA,vB
        integer :: i,j
        integer :: k
        real(dp) :: numMonolayers
        real(dp) :: dynamic_cap, effective_cap, min_cap_s, max_cap_s
        ! REAL(dp) :: vdiff(:)

        !Calculate all reaction rates
        !Assuming the user has temperature changes or uses the desorption features of phase 1,
        !these need to be recalculated every time step.

        ! CRP
        idx1=crpReacs(1)
        idx2=crpReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            rate(idx1:idx2) = alpha(idx1:idx2)*zeta
        end if
        if (improvedH2CRPDissociation) then
            rate(nR_H2_CRP)=h2CRPRate
        end if

        idx1=photonReacs(1)
        idx2=photonReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            rate(idx1:idx2) = alpha(idx1:idx2) * ( &
                radfield              * dexp(-gama(idx1:idx2)*av(dstep))          + &
                radfield_internal(dstep) * dexp(-gama(idx1:idx2)*av_internal(dstep)) &
                ) / 1.7
            ! For all solid species, decrease rate by 0.3 (Kalvans 2018)
            ! For bulk species, also decrease rate by (1-Pabs)**(Bs+0.5*Bb) (Kalvans 2014)
            do j=idx1,idx2
                if (ANY(bulkList==re1(j))) then
                    rate(j) = rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO * (1.0-0.007)**(1.0+0.5/bulkLayersReciprocal)
                else if (ANY(surfaceList==re1(j))) then
                    rate(j) = rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO
                end if
            end do
        end if

        !Reactions involving cosmic ray induced photon
        idx1=crphotReacs(1)
        idx2=crphotReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            rate(idx1:idx2)=alpha(idx1:idx2)*gama(idx1:idx2)*1.0/(1.0-omega)*zeta*(gasTemp(dstep)/300)**beta(idx1:idx2)
            ! For all solid species, decrease rate by 0.3 (Kalvans 2018)
            ! For bulk species, also decrease rate by (1-Pabs)**(Bs+0.5*Bb) (Kalvans 2014)
            do j=idx1,idx2
                if (ANY(bulkList==re1(j))) then
                    rate(j) = rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO * (1-0.007)**(1+0.5/bulkLayersReciprocal)
                else if (ANY(surfaceList==re1(j))) then
                    rate(j) = rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO
                end if
            end do
        end if

        !freeze out only happens if freezeFactor>0 and depending on evap choice
        idx1=freezeReacs(1)
        idx2=freezeReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            rate(idx1:idx2)=freezeOutRate(idx1,idx2)
            !freeze out rate uses thermal velocity but mass of E is 0 giving us infinite rates
            !just assume it's same as H
            rate(nR_EFreeze)=rate(nR_HFreeze)

            rate(nR_H2Freeze)=stickingCoefficient(h2StickingZero,h2StickingTemp,gasTemp(dstep))*rate(nR_H2Freeze)
            if (h2StickingCoeffByh2Coverage) then
                ! If all surface is H2, (i.e. x_#H2 = safeMantle), assume no H2 sticks
                ! and so set the sticking coeff to 0. Linearly interpolate according to chance it will hit a H2 molecule on surface
                rate(nR_H2Freeze)=rate(nR_H2Freeze)*(1.0D0-abund(ngh2, dstep)/safeMantle)
            end if

            rate(nR_HFreeze)=stickingCoefficient(hStickingZero,hStickingTemp,gasTemp(dstep))*rate(nR_HFreeze)
            if (hStickingCoeffByh2Coverage) then
                ! If all surface is H2, (i.e. x_#H2 = safeMantle), assume no H sticks
                ! and so set the sticking coeff to 0. Linearly interpolate according to chance it will hit a H2 molecule on surface
                rate(nR_HFreeze)=rate(nR_HFreeze)*(1.0D0-abund(ngh2, dstep)/safeMantle)
            end if
        end if
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !The below desorption mechanisms are from Roberts et al. 2007 MNRAS with
        !the addition of direct UV photodesorption. DESOH2,DESCR1,DEUVCR
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Desorption due to energy released by H2 Formations
        idx1=desoh2Reacs(1)
        idx2=desoh2Reacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            if ((desorb) .and. (h2desorb) .and. (safeMantle > MIN_SURFACE_ABUND)) then
                !Epsilon is efficiency of this process, number of molecules removed per event
                !h2form is formation rate of h2, dependent on hydrogen abundance.
                rate(idx1:idx2) = epsilon*h2FormEfficiency(gasTemp(dstep),dustTemp(dstep))
                !alpha is a branching ratio (default 1.0; use <1.0 for isomer desorption channels)
                rate(idx1:idx2) = alpha(idx1:idx2)*rate(idx1:idx2)

                !Don't remove species with binding energy > max BE removed by this process
                where(gama(idx1:idx2) > ebmaxh2) rate(idx1:idx2)=0.0
            else
                rate(idx1:idx2) = 0.0
            end if
            !turn off freeze out if desorption due to H2 formation is much faster
            !both rates combine with density to get rate of change so drop that factor
            ! Commented out: let DVODE handle the stiff competition rather than zeroing at step boundaries,
            ! which makes the solution step-count-dependent.
            ! WHERE((rate(freezePartners)*abund(re1(freezePartners),dstep))<&
            ! &MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
        end if
        !Desorption due to energy from cosmic rays
        idx1=descrReacs(1)
        idx2=descrReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            if ((desorb) .and. (crdesorb) .and. (safeMantle > MIN_SURFACE_ABUND)) then
                !4*pi*zeta = total CR flux. 1.64d-4 is iron to proton ratio of CR
                !as iron nuclei are main cause of CR heating.
                !GRAIN_SURFACEAREA_PER_H is the total surface area per hydrogen atom. ie total grain area per cubic cm when multiplied by density.
                !phi is efficiency of this reaction, number of molecules removed per event.
                rate(idx1:idx2) = 4.0*pi*zeta*1.64d-4*(GRAIN_SURFACEAREA_PER_H)*phi
                !alpha is a branching ratio (default 1.0; use <1.0 for isomer desorption channels)
                rate(idx1:idx2) = alpha(idx1:idx2)*rate(idx1:idx2)

                !Don't remove species with binding energy > max BE removed by this process
                where(gama(idx1:idx2) > ebmaxcr) rate(idx1:idx2)=0.0

            else
                rate(idx1:idx2) = 0.0
            end if
            !turn off freeze out if desorption due to CR formation is much faster
            ! WHERE((rate(freezePartners)*abund(re1(freezePartners),dstep)*density(dstep))&
            ! <MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
        end if

        !Desorption due to UV, partially from ISRF and partially from CR creating photons
        idx1=deuvcrReacs(1)
        idx2=deuvcrReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            if ((desorb) .and. (uvdesorb) .and. (safeMantle > MIN_SURFACE_ABUND)&
                    &.and.(zeta > 0)) then
                !4.875d3 = photon flux, Checchi-Pestellini & Aiello (1992) via Roberts et al. (2007)
                !UVY is yield per photon.
                rate(idx1:idx2) = GRAIN_CROSSSECTION_PER_H*uv_yield*4.875d3*zeta
                !additional factor accounting for UV desorption from ISRF. UVCREFF is ratio of
                !CR induced UV to ISRF UV.
                rate(idx1:idx2) = rate(idx1:idx2) &
                    & * (1 + (radfield/uvcreff)*(1.0/zeta)*dexp(-1.8*av(dstep)) &
                    & + (radfield_internal(dstep)/uvcreff)*(1.0/zeta) &
                    & * dexp(-1.8*av_internal(dstep)))
                !alpha is a branching ratio (default 1.0; use <1.0 for isomer desorption channels)
                rate(idx1:idx2) = alpha(idx1:idx2)*rate(idx1:idx2)

                !Don't remove species with binding energy > max BE removed by this process
                where(gama(idx1:idx2) > ebmaxuvcr) rate(idx1:idx2)=0.0
            else
                rate(idx1:idx2) = 0.0
            end if
            !turn off freeze out if desorption due to UV is much faster
            ! WHERE((rate(freezePartners)*abund(re1(freezePartners),dstep)*density(dstep))&
            ! &<MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
        end if

        !CRS reactions represent the production of excited species from cosmic ray bombardment
        !rate equations from Shingledecker et. al. 2018
        idx1=crsReacs(1)
        idx2=crsReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            !8.6 is the Spitzer-Tomasko cosmic ray flux in cm^-2 s^-1
            !1.3 converts to: ionization rate/10^-17
            rate(idx1:idx2)=alpha(idx1:idx2)*(beta(idx1:idx2)*(gama(idx1:idx2)/100)*(8.6*zeta*1.3))
        end if

        !EXRELAX, relaxation reactions for each excited species
        idx1=exrelaxReacs(1)
        idx2=exrelaxReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            do j=idx1,idx2
                do i=lbound(iceList,1),ubound(iceList,1)
                    if (iceList(i) == re1(j)) then
                        vA=vdiff(i)
                    end if
                end do
                rate(j)=vA
            end do
        end if

        !EXSOLID reactions represent the reactions of excited species on the grain
        idx1=exsolidReacs(1)
        idx2=exsolidReacs(2)

        if (idx1 /= REAC_NOT_PRESENT) then
            !reaction rates calculated outside of UCLCHEM as per Shingledecker et al. 2018 and included in grain network
            !alpha are branching ratios and beta is reaction rate
            do j=idx1,idx2
                do i=lbound(iceList,1),ubound(iceList,1)
                    if (iceList(i) == re1(j)) then
                        vA = vdiff(i)
                    end if
                    if (iceList(i) == re2(j)) then
                        vB = vdiff(i)
                    end if
                end do
                rate(j) = (vB + vA)/(SURFACE_SITE_DENSITY*1.8d-8)
                rate(j) = alpha(j) * rate(j)
            end do
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Continuous Thermal Desorption. Reactions can be generated through a flag in Makerates
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        idx1=thermReacs(1)
        idx2=thermReacs(2)
        if (idx1 /= REAC_NOT_PRESENT) then
            if (thermdesorb) then
                do j=idx1,idx2
                    !then try to overwrite with position in grain array
                    do i=lbound(iceList,1),ubound(iceList,1)
                        !See Cuppen, Walsh et al. 2017 review (section 4.1)
                        if (iceList(i) == re1(j)) then
                            !Basic rate at which thermal desorption occurs
                            rate(j)=vdiff(i)*exp(-gama(j)/dustTemp(dstep))
                            !alpha is a branching ratio (default 1.0; use <1.0 for isomer channels)
                            rate(j)=alpha(j)*rate(j)
                            !factor of 2.0 adjusts for fact only top two monolayers (Eq 8)
                            !because GRAIN_SURFACEAREA_PER_H is per H nuclei, multiplying it by density gives area/cm-3
                            !that is roughly sigma_g.n_g from cuppen et al. 2017 but using surface instead of cross-sectional
                            !area seems more correct for this process.
                            if (.NOT. THREE_PHASE) rate(j)=rate(j)*2.0*SURFACE_SITE_DENSITY*GRAIN_SURFACEAREA_PER_H
                        end if
                    end do
                end do
                !At some point, rate is so fast that there's no point freezing out any more
                !Save the integrator some trouble and turn freeze out off
                ! Compare against surface thermal desorption only (first SIZE(freezePartners)
                ! reactions of thermReacs).
                ! (surface # first, then bulk @); freeze-out competes with surface only.
                ! WHERE(rate(freezePartners)*abund(re1(freezePartners),dstep)*density(dstep)&
                ! &<MIN_SURFACE_ABUND*rate(idx1:idx1+SIZE(freezePartners)-1)) rate(freezePartners)=0.0
                if (safeMantle < MIN_SURFACE_ABUND) rate(idx1:idx2)=0.0
            else
                rate(idx1:idx2)=0.0
            end if
        end if


    !Reactions on surface can be treated considering diffusion of reactants
    !as in Langmuir-Hinshelwood mechanism
    !See work of David Quenard 2017 Arxiv:1711.05184
    !First calculate rate of the diffusion reaction

    ! Calculate ice coverage in monolayers for chemical desorption
    numMonolayers = getNumberMonolayers(safeMantle + safeBulk)

    idx1=lhReacs(1)
    idx2=lhReacs(2)
    if (idx1 /= REAC_NOT_PRESENT) then
        if ((dustTemp(dstep) < maxGrainTemp) .and. (safeMantle > MIN_SURFACE_ABUND)) then
            do j=idx1,idx2
                rate(j)=diffusionReactionRate(j,dustTemp(dstep))
            end do

            ! Save unsplit LH rates for dynamic re-split inside the F callback
            rate_lh_unsplit(lhReacs(1):lhReacs(2)) = rate(lhReacs(1):lhReacs(2))

            if ((desorb) .and. (chemdesorb)) then
                !two routes for every diffusion reaction: products to gas or products remain on surface
                !calculate fraction of reaction that goes down desorption route
                idx1=lhdesReacs(1)
                idx2=lhdesReacs(2)
                k = 0
                do i=idx1, idx2
                    k = k + 1
                    rate(i)=desorptionFractionIncludingIce(i, numMonolayers)*rate(LHDEScorrespondingLHreacs(k))
                    if (ANY(bulkList==re1(i))) rate(i)=0.0  ! Bulk species are not able to chemically desorb
                end do

                !remove that fraction from total rate of the diffusion route
                k = 0
                do i = idx1, idx2
                    k = k + 1
                    rate(LHDEScorrespondingLHreacs(k)) = rate(LHDEScorrespondingLHreacs(k)) - rate(i)
                end do
            else
                rate(idx1:idx2)=0.0
                rate(lhdesReacs(1):lhdesReacs(2))=0.0
            end if
        end if
    end if

    !Account for Eley-Rideal reactions in a similar way.
    !First calculate overall rate and then split between desorption and sticking
    idx1=erReacs(1)
    idx2=erReacs(2)
    if (idx1 /= REAC_NOT_PRESENT) then
        rate(idx1:idx2)=freezeOutRate(idx1,idx2)
        rate(idx1:idx2)=rate(idx1:idx2)*dexp(-gama(idx1:idx2)/dustTemp(dstep))

        ! Save unsplit ER rates for dynamic re-split inside the F callback
        rate_er_unsplit(erReacs(1):erReacs(2)) = rate(erReacs(1):erReacs(2))

        if ((desorb) .and. (chemdesorb)) then
            !calculate fraction of reaction that goes down desorption route
            idx1 = erdesReacs(1)
            idx2 = erdesReacs(2)
            k = 0
                do i=idx1, idx2
                    k = k + 1
                    rate(i)=desorptionFractionIncludingIce(i, numMonolayers)*rate(ERDEScorrespondingERreacs(k))
                if (ANY(bulkList==re1(i))) rate(i)=0.0  ! Bulk species are not able to chemically desorb
            end do

            !remove that fraction from total rate of the diffusion route
            k = 0
            do i = idx1, idx2
                k = k + 1
                rate(ERDEScorrespondingERreacs(k)) = rate(ERDEScorrespondingERreacs(k)) - rate(i)
            end do
        else
            rate(idx1:idx2)=0.0
            rate(erdesReacs(1):erdesReacs(2))=0.0
        end if
    end if

    select case (parameterizeH2Form)
    case (0)  ! Always disable parameterized H2 formation
        rate(nR_H2Form_CT) = 0.0
    case (1)  ! Always enable parameterized H2 formation
        rate(nR_H2Form_CT) = h2FormEfficiency(gasTemp(dstep), dustTemp(dstep))
        rate(nR_H2Form_ER) = 0.0
        rate(nR_H2Form_ERDes) = 0.0
    case (2)  ! Temperature-dependent: explicit LH/ER below maxGrainTemp, parameterized above
        if (dustTemp(dstep) >= maxGrainTemp) then
            rate(nR_H2Form_CT) = h2FormEfficiency(gasTemp(dstep), dustTemp(dstep))
            rate(nR_H2Form_ER) = 0.0
            rate(nR_H2Form_ERDes) = 0.0
        else
            rate(nR_H2Form_CT) = 0.0
        end if
    end select

    call bulkSurfaceExchangeReactions(rate,dustTemp(dstep))

    !Basic gas phase reactions
    !They only change if temperature has so we can save time with an if statement
    idx1=twobodyReacs(1)
    idx2=twobodyReacs(2)
    if (lastTemp /= gasTemp(dstep)) then
        rate(idx1:idx2) = alpha(idx1:idx2)*((gasTemp(dstep)/300.0)**beta(idx1:idx2))*dexp(-gama(idx1:idx2)/gasTemp(dstep))
    end if

    idx1=ionopol1Reacs(1)
    idx2=ionopol1Reacs(2)
    if (idx1 /= REAC_NOT_PRESENT) then
        !This formula including the magic numbers come from KIDA help page.
        rate(idx1:idx2)=alpha(idx1:idx2)*beta(idx1:idx2)*(0.62d0+0.4767d0*gama(idx1:idx2)*dsqrt(300.0d0/gasTemp(dstep)))
    end if

    idx1=ionopol2Reacs(1)
    idx2=ionopol2Reacs(2)
    if (idx1 /= REAC_NOT_PRESENT) then
        !This formula including the magic numbers come from KIDA help page.
        rate(idx1:idx2)=alpha(idx1:idx2)*beta(idx1:idx2)*(1.0d0+0.0967d0*gama(idx1:idx2)&
        &*dsqrt(300.0d0/gasTemp(dstep))+gama(idx1:idx2)*gama(idx1:idx2)*300.0/(10.526*gasTemp(dstep)))
    end if
    lastTemp=gasTemp(dstep)

    idx1=garReacs(1)
    idx2=garReacs(2)
    ! Adapted from NEATH (Priestley et al 2023)
    ! grain-assisted recombination stuff from Weingartner & Draine (2001)
    ! https://ui.adsabs.harvard.edu/abs/2001ApJ...563..842W/abstract
    ! We use the 0.6 factor as provided in Gong et al 2017 (DOI:10.3847/1538-4357/aa7561)
    phi = (radfield*exp(-2.5*av(dstep)) + radfield_internal(dstep)*exp(-2.5*av_internal(dstep))) &
    & * sqrt(gasTemp(dstep)) / (abund(nspec+1,dstep)*abund(nelec,dstep))  ! phi = G T^0.5 / n_e

    ! Ensure phi is within the 1e2 to 1e6 range from the paper:
    phi = min(max(phi,1e2), 1e6)

    if (idx1 /= REAC_NOT_PRESENT) then
        rate(idx1:idx2)= 0.6 * alpha(idx1:idx2) * garParams(:,1) / (1.0 + garParams(:,2) *&
        &phi**garParams(:,3) * (1.0 + garParams(:,4) * gasTemp(dstep)**garParams(:,5) *&
        &phi**(-garParams(:,6)-garParams(:,7)*log(gasTemp(dstep)))))
    end if

    !turn off reactions outside their temperature range
    where((.not. ExtrapolateRates) .and. (gasTemp(dstep) < minTemps)) rate=0.0

    where((.not. ExtrapolateRates) .and. (gasTemp(dstep) > maxTemps)) rate=0.0

    !Overwrite reactions for which we have a more detailed photoreaction treatment
    rate(nR_H2_hv)=H2PhotoDissRate(h2Col,radField,av(dstep),turbVel) &
               + H2PhotoDissRate(h2Col,radfield_internal(dstep),av_internal(dstep),turbVel)  !H2 photodissociation
    rate(nR_CO_hv)=COPhotoDissRate(h2Col,coCol,radField,av(dstep)) &
               + COPhotoDissRate(h2Col,coCol,radfield_internal(dstep),av_internal(dstep))  !CO photodissociation
    rate(nR_C_hv)=cIonizationRate(alpha(nR_C_hv),gama(nR_C_hv),gasTemp(dstep),ccol,h2col,av(dstep),radfield) &
               + cIonizationRate(alpha(nR_C_hv),gama(nR_C_hv),gasTemp(dstep),ccol,h2col, &
                                &av_internal(dstep),radfield_internal(dstep))  !C photoionization

    ! Encounter Desorption mechanism (Hincelin et al. 2015)
    ! Species diffuse onto H2-covered surfaces and can desorb upon encountering H2
    if ((h2EncounterDesorption) .and. (safeMantle > MIN_SURFACE_ABUND)) then
        rate(nR_H2_ED)=EncounterDesorptionRate(nR_H2_ED, dustTemp(dstep))  !H2 Encounter Desorption
    else
        rate(nR_H2_ED)=0.0D0
    end if

    if ((hEncounterDesorption) .and. (safeMantle > MIN_SURFACE_ABUND)) then
        ! H atom encounter desorption on H2-covered surfaces
        rate(nR_H_ED)=EncounterDesorptionRate(nR_H_ED, dustTemp(dstep))  !H Encounter Desorption
    else
        rate(nR_H_ED)=0.0D0
    end if

    ! Min floor: zero desorption rate constants k below numerical threshold to eliminate
    ! near-zero Arrhenius terms (e.g. strongly-bound species at low T) that waste solver work.
    if (min_desorption_rate > 0.0d0) then
        if (thermReacs(1)  /= REAC_NOT_PRESENT) then
          where(rate(thermReacs(1):thermReacs(2)) > 0.0d0 .AND. &
                  rate(thermReacs(1):thermReacs(2)) < min_desorption_rate) &
                rate(thermReacs(1):thermReacs(2)) = 0.0d0
        end if
        if (desoh2Reacs(1) /= REAC_NOT_PRESENT) then
          where(rate(desoh2Reacs(1):desoh2Reacs(2)) > 0.0d0 .AND. &
                  rate(desoh2Reacs(1):desoh2Reacs(2)) < min_desorption_rate) &
                rate(desoh2Reacs(1):desoh2Reacs(2)) = 0.0d0
        end if
        if (descrReacs(1)  /= REAC_NOT_PRESENT) then
          where(rate(descrReacs(1):descrReacs(2)) > 0.0d0 .AND. &
                  rate(descrReacs(1):descrReacs(2)) < min_desorption_rate) &
                rate(descrReacs(1):descrReacs(2)) = 0.0d0
        end if
        if (deuvcrReacs(1) /= REAC_NOT_PRESENT) then
          where(rate(deuvcrReacs(1):deuvcrReacs(2)) > 0.0d0 .AND. &
                  rate(deuvcrReacs(1):deuvcrReacs(2)) < min_desorption_rate) &
                rate(deuvcrReacs(1):deuvcrReacs(2)) = 0.0d0
        end if
        if (lhdesReacs(1)  /= REAC_NOT_PRESENT) then
          where(rate(lhdesReacs(1):lhdesReacs(2)) > 0.0d0 .AND. &
                  rate(lhdesReacs(1):lhdesReacs(2)) < min_desorption_rate) &
                rate(lhdesReacs(1):lhdesReacs(2)) = 0.0d0
        end if
        if (erdesReacs(1)  /= REAC_NOT_PRESENT) then
          where(rate(erdesReacs(1):erdesReacs(2)) > 0.0d0 .AND. &
                  rate(erdesReacs(1):erdesReacs(2)) < min_desorption_rate) &
                rate(erdesReacs(1):erdesReacs(2)) = 0.0d0
        end if
        if (rate(nR_H2_ED) > 0.0d0 .AND. rate(nR_H2_ED) < min_desorption_rate) rate(nR_H2_ED) = 0.0d0
        if (rate(nR_H_ED)  > 0.0d0 .AND. rate(nR_H_ED)  < min_desorption_rate) rate(nR_H_ED)  = 0.0d0
    end if

    ! Dynamic max cap: clamp thermal desorption k to prevent DVODE stiffness.
    ! Three-regime: effective_cap = clamp(factor/Dt_outer, min_cap, max_cap)
    !   k < min_cap_s -> always kept;  k > max_cap_s -> always capped;  in between -> dynamic.
    ! Cap bounds in yr^-1 are converted to s^-1; targetTime and currentTime are in seconds.
    if (max_desorption_rate_factor > 0.0d0 .AND. thermReacs(1) /= REAC_NOT_PRESENT) then
        min_cap_s   = min_desorption_rate_cap / SECONDS_PER_YEAR
        max_cap_s   = max_desorption_rate_cap / SECONDS_PER_YEAR
        dynamic_cap = max_desorption_rate_factor / MAX(1.0d-300, targetTime - currentTime)
        effective_cap = MIN(MAX(dynamic_cap, min_cap_s), max_cap_s)
        where(rate(thermReacs(1):thermReacs(2)) > effective_cap) &
            rate(thermReacs(1):thermReacs(2)) = effective_cap
    end if

    end subroutine calculateReactionRates



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Freeze out determined by rate of collisions with grain
!No sticking coefficient is used because typical values are >0.95 below 150 K
! eg Le Bourlot et al. 2013, Molpeceres et al. 2020
!Above 150 K, thermal desorption will completely remove grain species
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function freezeOutRate(idx1,idx2) result(freezeRates)
    real(dp) :: freezeRates(idx2-idx1+1)
    integer :: idx1,idx2

    !additional factor for ions (beta=0 for neutrals)
    freezeRates=1.0+beta(idx1:idx2)*16.71d-4/(GRAIN_RADIUS*gasTemp(dstep))
    if ((freezeFactor == 0.0) .or. (dustTemp(dstep) > maxGrainTemp)) then
        freezeRates=0.0
    else
        freezeRates=freezeRates*freezeFactor*alpha(idx1:idx2)*THERMAL_VEL&
        &*dsqrt(gasTemp(dstep)/mass(re1(idx1:idx2)))*GRAIN_CROSSSECTION_PER_H
    end if

    end function freezeOutRate


    function stickingCoefficient(stickingZero,criticalTemp,gasTemp) result(stickingCoeff)
        !Sticking coefficient for freeze out taken from Chaabouni et al. 2012 A&A 538 Equation 1
        real(dp) :: stickingCoeff
        real(dp) :: stickingZero,criticalTemp,gasTemp,tempRatio
        real(dp), parameter :: beta=2.5d0

        tempRatio=gasTemp/criticalTemp
        stickingCoeff=stickingZero*(1.0d0+beta*tempRatio)/((1.0d0+tempRatio)**beta)
    end function stickingCoefficient

end module RATES
