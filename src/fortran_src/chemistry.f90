! Chemistry module of UCL_CHEM.                                                               !
! Contains all the core machinery of the code, not really intended to be altered in standard  !
! use. Use a (custom) physics module to alter temp/density behaviour etc.                     !
!                                                                                             !
! chemistry module contains rates.f90, a series of subroutines to calculate all reaction rates!
! when updateChemistry is called from main, these rates are calculated, the ODEs are solved   !
! from currentTime to targetTime to get abundances at targetTime and then all abundances are  !
! written to the fullOutput file.                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE chemistry
USE constants
USE DEFAULTPARAMETERS
USE rates, only: lastTemp
!f2py INTEGER, parameter :: dp
USE physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
& zeta, currentTime, targetTime, timeinyears, freefall, density, ion, densdot, gasTemp, dustTemp, av, colDens
USE DVODE_F90_M !dvode_f90_m
USE network
USE photoreactions
USE surfacereactions
use f2py_constants, only: nspec, nreac
USE postprocess_mod, only: lusecoldens,usepostprocess,tstep,lnh,lnh2,lnco,lnc
USE rates
USE odes
USE heating
IMPLICIT NONE
    !f2py integer, intent(aux) :: points
    !These integers store the array index of important species and reactions, x is for ions    
    !loop counters    
    INTEGER :: i,j,l,writeCounter=0,loopCounter,failedIntegrationCounter
    INTEGER, PARAMETER :: maxLoops=10,maxConsecutiveFailures=10

    !Array to store reaction rates
    REAL(dp) :: rate(nreac)

    !DLSODE variables    
    INTEGER :: ITASK,ISTATE,NEQ
    REAL(dp), ALLOCATABLE :: abstol(:)
    REAL(dp), ALLOCATABLE :: reltol_vec(:)
    ! TYPE(VODE_OPTS) :: OPTIONS
    !initial fractional elemental abudances and arrays to store abundances
    REAL(dp) :: h2col,cocol,ccol,h2colToCell,cocolToCell,ccolToCell
    REAL(dp), ALLOCATABLE :: abund(:,:)
    REAL(dp) :: numMonolayers,ratioSurfaceToBulk
    
    REAL(dp) :: MIN_ABUND = 1.0d-30 !Minimum abundance allowed

    INTEGER :: nion,ionlist(nspec)

    REAL(dp) :: tempDot, oldTemp=0.0d0
    REAL(dp) :: h2form

    REAL(dp)::lastGasTemp,lastDustTemp

    !DVODE solver statistics (populated by integrateODESystem, read by output)
    REAL(dp) :: dvode_rstats(22)
    INTEGER :: dvode_istats(31)
    INTEGER :: dvode_istate_out
    REAL(dp) :: dvode_cpu_start, dvode_cpu_end, dvode_cpu_time

    !Solver statistics counter - tracks all DVODE calls including retries
    INTEGER :: solver_stats_counter

    ! Error code set inside the F callback; cannot use successFlag directly due to fixed DVODE signature.
    INTEGER :: f_callback_error = 0

CONTAINS
    SUBROUTINE initializeChemistry(readAbunds, successFlag)
        LOGICAL, INTENT(IN) :: readAbunds
        INTEGER, INTENT(INOUT) :: successFlag
        !f2py integer, intent(aux) :: points

        ! Sets variables at the start of every run.
        ! Since python module persists, it's not enough to set initial
        ! values in module definitions above. Reset here.
        NEQ=nspec+2
        IF (ALLOCATED(abund)) DEALLOCATE(abund,vdiff,vdes)
        IF (ALLOCATED(reltol_vec)) DEALLOCATE(reltol_vec)
        ALLOCATE(abund(NEQ,points),vdiff(SIZE(iceList)),vdes(SIZE(iceList)))
        !Set abundances to initial elemental if not reading them in.
        IF (.NOT. readAbunds) THEN
            !ensure abund is initially zero
            ! abund= MIN_ABUND
            abund(1:nspec,:)=MIN_ABUND

            !Start by filling all metallicity scaling elements
            !neutral atoms  
            abund(no,:) = fo  
            abund(nn,:) = fn               
            abund(nmg,:) = fmg
            abund(np,:) = fp
            abund(nf,:) = ff
            !abund(nfe,:) = ffe
            abund(nna,:) = fna
            abund(nli,:) = fli
            abund(npah,:) = fpah
            !default to ions
            abund(nsx,:) = fs
            abund(nsix,:) = fsi                
            abund(nclx,:) = fcl 
            !Decide how much carbon is initiall ionized using parameters.f90
            SELECT CASE (ion)
                CASE(0)
                    abund(nc,:)=fc
                    abund(ncx,:)=1.d-10
                CASE(1)
                    abund(nc,:)=fc*0.5
                    abund(ncx,:)=fc*0.5
                CASE(2)
                    abund(nc,:)=1.d-10
                    abund(ncx,:)=fc
            END SELECT

            !isotopes
            abund(n18o,:) = f18o  
            abund(n15n,:) = f15n           
            abund(n13c,:) = f13c    

            abund(nelec,:)=abund(ncx,:)+abund(nsix,:)+abund(nsx,:)+abund(nclx,:)+abund(nmgx,:)

            abund=abund*metallicity

            !Total H nuclei is always 1 so put fh into H and whatever is left over in H2
            abund(nh,:) = fh
            abund(nh2,:) = 0.5*(1.0e0-fh) 
            abund(nd,:)=fd

            abund(nhe,:) = fhe  
        ENDIF
        abund(nspec+2,:)=density      !Gas density
        abund(nspec+1,:)=gasTemp    !Gas temperature
        !Initial calculations of diffusion and desorption frequencies
        !Uses updateVdiffAndVdes which supports both HH1992 and TST treatments
        CALL updateVdiffAndVdes(gasTemp(1), dustTemp(1), SIZE(iceList), vdiff, vdes)

        ! get list of positive-charged species to conserve charge later
        nion = 0
        do i=1,nspec
           if (index(specname(i),'+') .ne. 0) then
              nion = nion + 1
              ionlist(nion) = i
           end if
        end do
        
        !DVODE SETTINGS
        ISTATE=1
        ITASK=1

        !set integration counts
        loopCounter=0
        failedIntegrationCounter=0
        solver_stats_counter=0  !Reset solver statistics counter

        IF (.NOT. ALLOCATED(abstol)) THEN
            ALLOCATE(abstol(NEQ))
        END IF
        IF (.NOT. ALLOCATED(reltol_vec)) THEN
            ALLOCATE(reltol_vec(NEQ))
        END IF
        !OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR=reltol,USER_SUPPLIED_JACOBIAN=.FALSE.)
        
        IF (heatingFlag) THEN
            !Initializing heating.f90 --> get coolants
            CALL initializeHeating(gasTemp(dstep),density(dstep),abund(:,1),colDens(dstep),cloudSize,successFlag)
            IF (successFlag .lt. 0) RETURN
        END IF

        !Set rates to zero to ensure they don't hold previous values or random ones if we don't set them in calculateReactionRates
        rate=0.0

        !We typically don't recalculate rates that only depend on temperature if the temp hasn't changed
        !use arbitrarily high value to make sure they are calculated at least once.
        lastGasTemp=99.0d99
        lastDustTemp=99.0d99
        lastTemp=99.0d99  ! Reset rates module lastTemp to force recalculation

        ! If the hydrogen diffusion energy is still its default value (-1.0 in default_parameters.f90),
        ! i.e. no custom value was set in the input dictionary, set it to the correct value
        ! according to the ratio of the ratio of diffusion energy to binding energy.
        IF (HdiffusionBarrier .eq. -1.0) THEN
            DO i = LBOUND(iceList, 1), UBOUND(iceList, 1)
                IF (iceList(i) .eq. ngh) HdiffusionBarrier = diffToBindRatio*bindingEnergy(i)
            END DO
        END IF
        
        ! Pre-calculate desorption fractions for LHDES and ERDES reactions
        desorptionFractionsBare = 0.0D0
        desorptionFractionsFullCoverage = 0.0D0
        DO j = lhdesReacs(1), lhdesReacs(2)
            desorptionFractionsBare(j) = getDesorptionFractionBare(j, j-lhdesReacs(1)+1)
            desorptionFractionsFullCoverage(j) = getDesorptionFractionFullCoverage(j, j-lhdesReacs(1)+1)
        END DO
        DO j = erdesReacs(1), erdesReacs(2)
            desorptionFractionsBare(j) = getDesorptionFractionBare(j, j-erdesReacs(1)+1)
            desorptionFractionsFullCoverage(j) = getDesorptionFractionFullCoverage(j, j-erdesReacs(1)+1)
        END DO
        
    END SUBROUTINE initializeChemistry

    SUBROUTINE updateChemistry(successFlag, statsarray, statsarray_size, dtime)
    !Updates the abundances for the next time step, first updating chemical variables and reaction rates,
    !then by solving the ODE system to obtain new abundances.
    !Solving ODEs is complex so we have two checks to try to automatically overcome difficulties and end stalled models
    !Firstly, the integration subroutine is called up to maxLoops times whilst adjusting variables to help integration converge.
    !If it succeeds before maxLoops, we continue as normal, otherwise we'll call it a fail.
    !Secondly, we check for stalls caused by the the solver loop reducing the targetTime to overcome difficulties.
    !That reduction the possibility of the code "succeeding" by integrating tiny target times. We have a counter that resets each time
    !the code integrates to the planned targetTime rather than a reduced one. If the counter reaches maxConsecutiveFailures, we end the code.
        !f2py integer, intent(aux) :: points
        INTEGER, INTENT(OUT) :: successFlag
        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL, DIMENSION(:,:,:) :: statsarray
        INTEGER, INTENT(IN), OPTIONAL :: statsarray_size
        INTEGER, INTENT(IN), OPTIONAL :: dtime
        real(dp) :: originalTargetTime !targetTime can be altered by integrator but we'd like to know if it was changed
        real(dp) :: surfaceCoverage
        real(dp) :: h2form_CT_vol, h2form_LH_vol, h2form_ER_vol  ! per-mechanism volumetric H2 formation rates [cm^-3 s^-1]
        real(dp) :: h2form_heat  ! mechanism-weighted H2 formation heating [erg cm^-3 s^-1]
        real(dp) :: h2heatfac, h2_denom  ! H&M79 eq. 6.45 thermalization efficiency factor

        ! write(*,*) "update Chemistry ..."
        !Integration can fail in a way that we can manage. Allow maxLoops tries before giving up.
        loopCounter=0
        successFlag=0
        originalTargetTime=targetTime
        DO WHILE((currentTime .lt. targetTime) .and. (loopCounter .lt. maxLoops)) 
            !allow option for dens to have been changed elsewhere.
            IF (.not. freefall) abund(nspec+2,dstep)=density(dstep)

            !First sum the total column density over all points further towards edge of cloud
            IF (dstep.gt.1) THEN
                h2ColToCell=(sum(abund(nh2,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                coColToCell=(sum(abund(nco,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                cColToCell=(sum(abund(nc,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
            ELSE
                h2ColToCell=0.0
                coColToCell=0.0
                cColToCell=0.0
            ENDIF
            !then add half the column density of the current point to get average in this "cell"
            h2Col=h2ColToCell+0.5*abund(nh2,dstep)*density(dstep)*(cloudSize/real(points))
            coCol=coColToCell+0.5*abund(nco,dstep)*density(dstep)*(cloudSize/real(points))
            cCol=cColToCell+0.5*abund(nc,dstep)*density(dstep)*(cloudSize/real(points))

            ! Postprocessed tracers have column densities provided
            if (lusecoldens) then
               h2col = lnh2(tstep)
               cocol = lnco(tstep)
               ! ccol = lnc(dstep, tstep) ! TODO enable C column density support
               ccol = lnh(tstep) * abund(nc,dstep) ! No C column densities yet...
            end if

            !Reset surface and bulk values in case of integration error or sputtering
            abund(nBulk,dstep)=sum(abund(bulkList,dstep))
            abund(nSurface,dstep)=sum(abund(surfaceList,dstep))
            !recalculate coefficients for ice processes
            safeMantle=MAX(1d-30,abund(nSurface,dstep))
            safeBulk=MAX(1d-30,abund(nBulk,dstep))

            if (refractoryList(1) .gt. 0) safeBulk=safeBulk-SUM(abund(refractoryList,dstep))

            ratioSurfaceToBulk=MIN(1.0D0, safeMantle/safeBulk)
            bulkLayersReciprocal=MIN(1.0,NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
            surfaceCoverage=bulkGainFromMantleBuildUp()

            IF ((.NOT. dustTemp(dstep) .eq. lastDustTemp) .OR. &
                (.NOT. gasTemp(dstep) .eq. lastGasTemp)) THEN
                CALL updateVdiffAndVdes(gasTemp(dstep), dustTemp(dstep), SIZE(icelist), vdiff, vdes)
            END IF

            CALL calculateReactionRates(abund,safeMantle, h2col, cocol, ccol, rate)
            if (heatingFlag) then
                ! TODO: check that the local and global radiation fields are correct.
                dustTemp(dstep)=calculateDustTemp(radfield*EXP(-UV_FAC*av(dstep)),radfield,av(dstep),zeta)


                ! Per-mechanism volumetric H2 formation rates [cm^-3 s^-1] 
                ! Only accounting for H2 ending up in the gas phase.
                h2form_CT_vol = rate(nR_H2Form_CT) * abund(nspec+2,dstep)**2 * abund(nh,dstep)
                h2form_LH_vol = (rate(nR_H2Form_LHDes)) &
                              &  * abund(ngh,dstep)**2 * abund(nspec+2,dstep)
                h2form_ER_vol = (rate(nR_H2Form_ERDes)) &
                              &  * abund(nspec+2,dstep)**2 * abund(nh,dstep) * abund(ngh,dstep) / safeMantle
                ! H&M79 eq. 6.45: critical density for H2 thermalization
                ! (18100 coefficient for consistency with h2FUVPumpHeating in heating.f90)
                h2_denom = 1.6d0*abund(nh,dstep)*EXP(-((400.0d0/gasTemp(dstep))**2)) &
                         &+ 1.4d0*abund(nh2,dstep)*EXP(-(18100.0d0/(gasTemp(dstep)+1200.0d0)))
                IF (h2_denom > 0.0d0) THEN
                    h2heatfac = 1.0d0 / (1.0d0 + 1.0d6/(SQRT(gasTemp(dstep))*h2_denom*abund(nspec+2,dstep)))
                ELSE
                    h2heatfac = 0.0d0
                END IF
                ! H&M79 eq. 6.43: LH gives 0.1 eV kinetic + 4.2 eV vibrational (fraction h2heatfac goes to gas)
                ! ER: 0.6 eV (Bourlot et al. 2012), thermalization-corrected
                ! CT: 1.5 eV (Hollenbach & Tielens 1999), no thermalization correction
                h2form_heat = eV * (1.5d0*h2form_CT_vol &
                            &+ (0.1d0 + 4.2d0*h2heatfac)*h2form_LH_vol &
                            &+ 0.6d0*h2heatfac*h2form_ER_vol)
                tempDot= getTempDot(&
                                &    timeinyears, &                       ! time
                                &    abund(nspec+1,dstep), &              ! gas temperature
                                &    abund(nspec+2,dstep), &              ! gas density
                                &    colDens(dstep), &                    ! gas column density
                                &    radfield*EXP(-UV_FAC*av(dstep)), &   ! attenuated radiation field
                                &    abund(:,dstep), &                    ! full abundance vector
                                &    rate(nR_H2_hv), &                     ! H2 dissociation rate
                                &    h2form_heat, &                       ! mechanism-weighted H2 formation heating [erg cm^-3 s^-1]
                                &    zeta, &                              ! cosmic ray ionization rate
                                &    rate(nR_C_hv), &                     ! C-photo rate
                                &    1.0/GAS_DUST_DENSITY_RATIO, &        ! dust-to-gas ratio
                                &    grain_Radius, &                      ! grain radius
                                &    metallicity, &                       ! metallicity
                                ! &    heatWriteFlag, &                     ! write flag
                                &    dusttemp(dstep), &                   ! dust temperature
                                &    turbVel &                            ! turbulence velocity
                                )
            end if

            !Integrate chemistry, and return fail if unrecoverable error was reached
            IF (PRESENT(statsarray) .AND. PRESENT(statsarray_size) .AND. PRESENT(dtime)) THEN
                CALL integrateODESystem(successFlag, statsarray, statsarray_size, dtime)
            ELSE
                CALL integrateODESystem(successFlag)
            END IF
            IF (successFlag .lt. 0) THEN
                write(*,*) "Integration failed, exiting"
                RETURN
            END IF

            !1.d-30 stops numbers getting too small for fortran.
            ! WHERE(abund<MIN_ABUND) abund=MIN_ABUND
            WHERE(abund(1:nspec,:)<MIN_ABUND) abund(1:nspec,:)=MIN_ABUND
            gasTemp(dstep)=abund(nspec+1,dstep)
            density(dstep)=abund(nspec+2,dstep)
            ! IF (gasTemp(dstep) .lt. 10) gasTemp(dstep)=10.0
            IF (gasTemp(dstep) .lt. lower_limit_gastemp) gasTemp(dstep)=lower_limit_gastemp
            loopCounter=loopCounter+1

            ! For postprocessing, force solver to try and reach original target time
            if (usepostprocess) targettime = originaltargettime
        END DO

        ! Postprocessing needs to reach next timestep whatever the cost
        if (.not. usepostprocess) then
        IF (loopCounter .eq. maxLoops) successFlag=INT_TOO_MANY_FAILS_ERROR

        !Since targetTime can be altered, eventually leading to "successful" integration we want to
        !check if integrator ever just reaches the planned target time. If it doesn't for many attempts,
        !we will call the run a failure. This stops the target being constantly reduced to tiny increments
        !so that the code all but stalls as the time is increased by seconds each integraiton.
        IF (ABS(originalTargetTime- targetTime) .lt. 0.001*originalTargetTime) THEN
            failedIntegrationCounter=0
        ELSE
            failedIntegrationCounter=failedIntegrationCounter+1
        END IF
        IF (failedIntegrationCounter .gt. maxConsecutiveFailures)&
             &successFlag=INT_TOO_MANY_FAILS_ERROR
        end if

    END SUBROUTINE updateChemistry

    SUBROUTINE integrateODESystem(successFlag, statsarray, statsarray_size, dtime)
        INTEGER, INTENT(OUT) :: successFlag
        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL, DIMENSION(:,:,:) :: statsarray
        INTEGER, INTENT(IN), OPTIONAL :: statsarray_size
        INTEGER, INTENT(IN), OPTIONAL :: dtime
        TYPE(VODE_OPTS) :: OPTIONS
        ! Non-negativity constraints for all chemical species: DVODE predictor/Jacobian
        ! steps can push species (especially surface ones starting at 1e-30) to small
        ! negative values at high densities, causing /safeMantle explosions in F.
        INTEGER :: species_indices(nspec)
        REAL(dp) :: species_clower(nspec), species_cupper(nspec)
        successFlag=0
        f_callback_error=0

    !This subroutine calls DVODE (3rd party ODE solver) until it can reach targetTime with acceptable errors (reltol/abstol)
        !reset parameters for DVODE
        ITASK=1 !try to integrate to targetTime
        ISTATE=1 !pretend every step is the first
        !Alternative: if (ISTATE .lt. 0) ISTATE=1  !only reset on error, allows Jacobian reuse
        !Gas-phase species: absolute tolerances scaled by abundance
        abstol=abstol_factor*abund(:,dstep)
        WHERE(abstol<abstol_min) abstol=abstol_min
        !Ice species (surface + bulk): separate, looser absolute tolerances
        abstol(iceList) = abstol_ice_factor*abund(iceList,dstep)
        WHERE(abstol(iceList)<abstol_ice_min) abstol(iceList)=abstol_ice_min
        !Physical variables: separate tolerance heuristic (T and nH need looser tolerances)
        abstol(nspec+1) = MAX(abstol_phys_factor * ABS(abund(nspec+1,dstep)), abstol_T_min)
        abstol(nspec+2) = MAX(abstol_phys_factor * ABS(abund(nspec+2,dstep)), abstol_nH_min)
        !Per-component relative tolerances: tight for chemistry, relaxed for physics
        reltol_vec(1:nspec) = reltol
        reltol_vec(nspec+1) = reltol_phys
        reltol_vec(nspec+2) = reltol_phys
        !Call the integrator with ITOL=4 (vector reltol + vector abstol).
        species_indices = [(i, i=1, nspec)]
        species_clower  = 0.0_dp
        species_cupper  = HUGE(1.0_dp)
        OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR_VECTOR=reltol_vec, &
                           USER_SUPPLIED_JACOBIAN=.False., MXSTEP=MXSTEP, &
                           CONSTRAINED=species_indices, CLOWER=species_clower, CUPPER=species_cupper)
        CALL CPU_TIME(dvode_cpu_start)
        CALL DVODE_F90(F,NEQ,abund(:,dstep),currentTime,targetTime,ITASK,ISTATE,OPTIONS)
        CALL CPU_TIME(dvode_cpu_end)
        dvode_cpu_time = dvode_cpu_end - dvode_cpu_start
        dvode_istate_out = ISTATE
        CALL GET_STATS(dvode_rstats, dvode_istats)

        ! Write solver statistics immediately after EVERY DVODE call (including failures)
        IF (PRESENT(statsarray) .AND. PRESENT(statsarray_size) .AND. PRESENT(dtime)) THEN
            solver_stats_counter = solver_stats_counter + 1

            ! Check for array overflow
            IF (solver_stats_counter > statsarray_size) THEN
                write(*,*) "ERROR: Solver stats array overflow at counter", solver_stats_counter
                write(*,*) "       Allocated size:", statsarray_size
                write(*,*) "       Consider increasing statsarray allocation or reducing finalTime"
                successFlag = SOLVER_STATS_OVERFLOW_ERROR
                RETURN
            END IF

            ! Write stats: column 1 = trajectory index, rest shifted by 1
            statsarray(solver_stats_counter, dstep, 1) = DBLE(dtime)
            statsarray(solver_stats_counter, dstep, 2) = DBLE(dvode_istate_out)
            statsarray(solver_stats_counter, dstep, 3:6) = dvode_rstats(11:14)
            statsarray(solver_stats_counter, dstep, 7:18) = DBLE(dvode_istats(11:22))
            statsarray(solver_stats_counter, dstep, 19) = dvode_cpu_time
        END IF

        IF (f_callback_error .lt. 0) THEN
            successFlag = f_callback_error
            RETURN
        END IF

        SELECT CASE(ISTATE)
            CASE(-1)
                !ISTATE -1 means the integrator can't break the problem into small enough steps
                !We could increase MXSTEP but better to reduce targetTime and get to physics update
                !physical conditions may be easier to solve as time goes by so better to get to that update
                write(*,*) "ISTATE -1: Reducing time step"
                !More steps required for this problem
                !MXSTEP=MXSTEP*2   
                targetTime=currentTime+(targetTime-currentTime)*0.1
            CASE(-2)
                !ISTATE -2 just needs an absol change so let's do that and try again
                write(*,*) "ISTATE -2: Tolerances too small"
                !Tolerances are too small for machine but succesful to current currentTime
                abstol_factor=abstol_factor*10.0
                abstol_ice_factor=abstol_ice_factor*10.0
                reltol_phys=MIN(reltol_phys*10.0, 1.0d-1)
            CASE(-3)
                !ISTATE -3 is unrecoverable so just bail on intergration
                write(*,*) "DVODE found invalid inputs"
                write(*,*) "abstol:"
                write(*,*) abstol
                successFlag=INT_UNRECOVERABLE_ERROR
                RETURN
            CASE(-4)
                !Successful as far as currentTime but many errors.
                !Make targetTime smaller and just go again
                write(*,*) "ISTATE -4 - shortening step"
                targetTime=currentTime+(targetTime-currentTime)*0.1
            CASE(-5)
                write(*,*) "ISTATE -5 - shortening step at time", timeInYears,"years"
                targetTime=currentTime+(targetTime-currentTime)*0.1
            CASE default
                MXSTEP=10000    
        END SELECT
    if (enforceChargeConservation) then
        ! REALLY ensure charge is always conserved (also after integrating)
        abund(nelec,dstep) = sum(abund(ionlist(1:nion),dstep))
    end if
    END SUBROUTINE integrateODESystem

    SUBROUTINE F (NEQUATIONS, T, Y, YDOT)
        USE ODES
        INTEGER, PARAMETER :: WP = KIND(1.0D0)
        INTEGER NEQUATIONS
        REAL(WP) T
        REAL(WP), DIMENSION(NEQUATIONS) :: Y, YDOT
        INTENT(IN)  :: NEQUATIONS, T, Y
        INTENT(OUT) :: YDOT
        REAL(dp) :: D,loss,prod
        REAL(dp) :: surfaceCoverage
        REAL(dp) :: phi,cgr(6),grec,denom
        REAL(dp) :: h2heatfac, h2_denom  ! H&M79 eq. 6.45 thermalization efficiency factor
        INTEGER :: ii
        !Set D to the gas density for use in the ODEs
        D=y(nspec+2)     !Gas density
        ydot=0.0

        ! Column densities are fixed for postprocessing data, so don't do this bit
        if (.not. lusecoldens) then
        !changing abundances of H2 and CO can causes oscillation since their rates depend on their abundances
        !recalculating rates as abundances are updated prevents that.
        !thus these are the only rates calculated each time the ODE system is called.
        cocol=coColToCell+0.5*Y(nco)*D*(cloudSize/real(points))
        h2col=h2ColToCell+0.5*Y(nh2)*D*(cloudSize/real(points))
        rate(nR_H2_hv)=H2PhotoDissRate(h2Col,radField,av(dstep),turbVel) !H2 photodissociation
        rate(nR_CO_hv)=COPhotoDissRate(h2Col,coCol,radField,av(dstep)) !CO photodissociation
        end if

        !recalculate coefficients for ice processes
        safeMantle=MAX(1d-30,Y(nSurface))
        safeBulk=MAX(1d-30,Y(nBulk))
        bulkLayersReciprocal=MIN(1.0,NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
        surfaceCoverage=bulkGainFromMantleBuildUp()

        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        CALL GETYDOT(RATE, Y, bulkLayersReciprocal, ratioSurfaceToBulk, surfaceCoverage, safeMantle,safeBulk, D, YDOT)
        ! get density change from physics module to send to DLSODE
        if (enforceChargeConservation) then 
            ydot(nelec) = sum(ydot(ionlist(1:nion)))
        end if 

        ydot(nspec+2) = densdot(Y(nspec+2))     !Gas density ODE

        IF (heatingFlag) THEN
            ! Check for diverging negative abundances (beyond tolerance).
            ! Values in (-negative_abundance_tol, 0) are solver noise and handled downstream.
            IF (ANY(Y(1:nspec) < -negative_abundance_tol)) THEN
                WRITE(*,'(A,ES12.4,A,I4)') "ERROR: negative abundance(s) entering heating at t=", &
                    timeInYears, " yr, dstep=", dstep
                DO ii = 1, nspec
                    IF (Y(ii) < -negative_abundance_tol) WRITE(*,'(4X,A,A,ES12.4)') &
                        TRIM(specname(ii)), ": ", Y(ii)
                END DO
                f_callback_error = NEGATIVE_ABUNDANCE_ERROR
                RETURN
            END IF
            ! Write(*,*) "Updating heating and cooling rates"
            IF (ABS(y(nspec+1)-oldTemp).gt.MIN(heating_temp_abstol, heating_temp_reltol*oldTemp)) THEN
                gasTemp(dstep)=y(nspec+1)
                IF (gasTemp(dstep) .lt. lower_limit_gastemp) gasTemp(dstep)=lower_limit_gastemp
                IF (gasTemp(dstep) .gt. upper_limit_gastemp) gasTemp(dstep)=upper_limit_gastemp
                ! H&M79 eq. 6.45: critical density for H2 thermalization
                ! (18100 coefficient for consistency with h2FUVPumpHeating in heating.f90)
                h2_denom = 1.6d0*Y(nh)*EXP(-((400.0d0/Y(nspec+1))**2)) &
                         &+ 1.4d0*Y(nh2)*EXP(-(18100.0d0/(Y(nspec+1)+1200.0d0)))
                IF (h2_denom > 0.0d0) THEN
                    h2heatfac = 1.0d0 / (1.0d0 + 1.0d6/(SQRT(Y(nspec+1))*h2_denom*Y(nspec+2)))
                ELSE
                    h2heatfac = 0.0d0
                END IF
                ! Only desorbing products contribute to gas heating; LH/ER remain on grain.
                ! H&M79 eq. 6.43: LHDes gives 0.1 eV kinetic + 4.2 eV vibrational (fraction h2heatfac to gas)
                ! ERDes: 0.6 eV (Bourlot et al. 2012), thermalization-corrected
                ! CT: 1.5 eV (Hollenbach & Tielens 1999), no thermalization correction
                h2form = eV * ( &
                    &  1.5d0 * rate(nR_H2Form_CT) * Y(nspec+2)**2 * Y(nh) &
                    &+ (0.1d0 + 4.2d0*h2heatfac) * rate(nR_H2Form_LHDes) * Y(ngh)**2 * Y(nspec+2) &
                    &+ 0.6d0 * h2heatfac * rate(nR_H2Form_ERDes) * Y(nspec+2)**2 * Y(nh) * Y(ngh) &
                    &  / max(safeMantle, MIN_SURFACE_ABUND) )
                tempDot=getTempDot( &
                            &    timeInYears, &                         ! time
                            &    Y(nspec+1), &                          ! gas temperature
                            &    Y(nspec+2), &                          ! gas density
                            &    colDens(dstep), &                      ! gas column density
                            &    radfield*EXP(-UV_FAC*av(dstep)), &     ! attenuated radiation field
                            &    Y, &                                   ! all number densities
                            &    rate(nR_H2_hv), &                      ! H2 dissociation rate computed in rates.f90
                            &    h2form, &                              ! mechanism-weighted H2 formation heating [erg cm^-3 s^-1]
                            &    zeta, &                                ! cosmic ray ionization rate
                            &    rate(nR_C_hv), &                       ! C-photo rate
                            &    1.0/GAS_DUST_DENSITY_RATIO, &          ! dust-to-gas ratio
                            &    grain_Radius, &                        ! grain radius
                            &    metallicity, &                         ! metallicity
                            &    dusttemp(dstep), &                     ! dust temperature
                            &    turbVel &                              ! turbulence velocity
                        )
                oldTemp=y(nspec+1)
            END IF
            ydot(nspec+1)=tempDot
        ELSE
            ydot(nspec+1)=0.0
        END IF

    END SUBROUTINE F

    ! SUBROUTINE JAC(NEQ, T, Y, ML, MU, J, NROWPD)
    !     INTEGER NEQ,ML,MU,NROWPD
    !     DOUBLE PRECISION T, Y(NEQ), J(NROWPD,NEQ)
    !     REAL(DP) :: D
    !     INTENT(IN)  :: NEQ, T, Y,ML,MU,NROWPD
    !     INTENT(INOUT) :: J
    !     D=y(NEQ)

    !     J=0.0d0
    !     INCLUDE 'jacobian.f90'
    !     J(nh,nh2)=J(nh,nh2)+2.0*h2dis
    !     J(nh2,nh2)=J(nh,nh2)-h2dis
    ! END SUBROUTINE JAC

END MODULE chemistry
