! Chemistry module of UCLCHEM.                                                                          !
! Contains all the core machinery of the code, not really intended to be altered in standard            !
! use. Use a (custom) physics module to alter temperature/density behavior etc.                         !
!                                                                                                       !
! chemistry module contains rates.f90, a series of subroutines to calculate all reaction rate constants !
! when updateChemistry is called from main, these rate constants are calculated, the ODEs are solved    !
! from currentTime to targetTime to get abundances at targetTime and then all abundances are            !
! written to the fullOutput file.                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module chemistry
    use constants, only: dp, SECONDS_PER_YEAR, CONSERVATION_ERROR, SOLVER_STATS_OVERFLOW_ERROR, &
        INT_TOO_MANY_FAILS_ERROR, INT_UNRECOVERABLE_ERROR, NOT_ENOUGH_TIMEPOINTS_ERROR, MIN_ABUND, eV
    use DEFAULTPARAMETERS
    use DVODE_F90_M, only: SET_OPTS, DVODE_F90, VODE_OPTS, GET_STATS  !dvode_f90_m
    !f2py INTEGER, parameter :: dp
    use physicscore, only: points, dstep, cloudsize, radfield, radfield_internal, h2crprateconstant, improvedH2CRPDissociation, &
    & zeta, currentTime, targetTime, timeinyears, freefall, density, ion, get_densdot, gasTemp, dustTemp, av, av_internal, colDens
    use f2py_constants, only: nSpec, nReac, REAC_NOT_PRESENT
    use heating, only: calculateDustTemp, getTempDot, initializeHeating
    ! allow(use-all)
    use network
    use numerics, only: is_equal
    use odes, only: GETYDOT
    use photoreactions, only: getH2PhotoDissRateConstant, getCOPhotoDissRateConstant, UV_FAC
    use postprocess_mod, only: lusecoldens,usepostprocess,tstep,lnh,lnh2,lnco,lnc
    use rates, only: calculateReactionRateConstants, turbVel, lastTemp, rate_constants_lh_unsplit, rate_constants_er_unsplit
    use surfacereactions, only: updateVdiffAndVdes, getNumberMonolayers, getDesorptionFractionBare, &
        getDesorptionFractionFullCoverage, NUM_SITES_PER_GRAIN, GAS_DUST_DENSITY_RATIO, vdes, vdiff, &
        desorptionFractionsBare, getDesorptionFractionIncludingIce, desorptionFractionsFullCoverage, &
        safeMantle, safeBulk, bulkLayersReciprocal, GRAIN_RADIUS, MIN_SURFACE_ABUND, &
        bulkGainFromMantleBuildUp, surfaceToBulkSwappingRateConstants

    implicit none

    public

    !f2py integer, intent(aux) :: points
    !These integers store the array index of important species and reactions, x is for ions
    !loop counters
    integer :: i,j,l,writeCounter=0,loopCounter,failedIntegrationCounter
    integer, parameter :: maxLoops=10,maxConsecutiveFailures=10

    !Array to store reaction rate constants
    real(dp) :: rate_constants(nReac)

    !DLSODE variables
    integer :: ITASK,ISTATE,NEQ
    real(dp), allocatable :: abstol(:)
    real(dp), allocatable :: reltol_vec(:)

    !initial fractional elemental abundances and arrays to store abundances
    real(dp) :: h2col,cocol,ccol,h2colToCell,cocolToCell,ccolToCell
    real(dp), allocatable :: abund(:,:)
    real(dp) :: numMonolayers,ratioSurfaceToBulk, surfGrowthUncorrected

    integer :: nion
    integer, dimension(nSpec) :: ionlist

    real(dp) :: tempDot, oldTemp=0.0_dp, prevIntegrationTemp=0.0_dp
    real(dp) :: h2form

    real(dp) :: lastGasTemp,lastDustTemp

    !DVODE solver statistics (populated by integrateODESystem, read by output)
    !Dimensions are minimum necessary dimensions for each (see dvode.f90)
    ! allow(magic-number-in-array-size)
    real(dp) :: dvode_rstats(22)
    ! allow(magic-number-in-array-size)
    integer :: dvode_istats(31)
    integer :: dvode_istate_out
    real(dp) :: dvode_cpu_start, dvode_cpu_end, dvode_cpu_time

    !Solver statistics counter - tracks all DVODE calls including retries
    integer :: solver_stats_counter

    ! Error code set inside the F callback; cannot use successFlag directly due to fixed DVODE signature.
    integer :: f_callback_error = 0

    ! Initial elemental abundances per parcel, used for runtime conservation check.
    ! Shape: (n_elem_tracked, points) - allocated in initializeChemistry.
    real(dp), allocatable :: initial_elem_abund(:,:)

contains
    pure subroutine get_ion_indices(ionlist, nion)
        ! get list of positive-charged species to conserve charge later
        integer, intent(out) :: ionlist(nSpec)
        integer, intent(out) :: nion

        integer :: i

        nion = 0
        do i=1,nSpec
           if (index(specname(i),"+") /= 0) then
              nion = nion + 1
              ionlist(nion) = i
           end if
        end do
    end subroutine get_ion_indices

    subroutine initializeChemistry(readAbunds, successFlag)
        logical, intent(in) :: readAbunds
        integer, intent(inout) :: successFlag
        !f2py integer, intent(aux) :: points

        ! Sets variables at the start of every run.
        ! Since python module persists, it's not enough to set initial
        ! values in module definitions above. Reset here.
        NEQ=nSpec+2
        if (ALLOCATED(abund)) deallocate(abund,vdiff,vdes)
        if (ALLOCATED(reltol_vec)) deallocate(reltol_vec)
        allocate(abund(NEQ,points),vdiff(N_ICE_SPECIES),vdes(N_ICE_SPECIES))
        !Set abundances to initial elemental if not reading them in.
        if (.NOT. readAbunds) then
            !ensure abund is initially zero
            ! abund= MIN_ABUND
            abund(1:nSpec,:)=MIN_ABUND

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
            !Decide how much carbon is initially ionized using parameters.f90
            select case (ion)
                case(0)
                    abund(nc,:)=fc
                    abund(ncx,:)=1.0e-10_dp
                case(1)
                    abund(nc,:)=fc*0.5_dp
                    abund(ncx,:)=fc*0.5_dp
                case(2)
                    abund(nc,:)=1.0e-10_dp
                    abund(ncx,:)=fc
                case default
                    write(*, *) "ERROR: Invalid value for 'ion':", ion
                    write(*, *) "Valid values are: 0, 1, or 2."
                    stop 1
            end select

            !isotopes
            abund(n18o,:) = f18o
            abund(n15n,:) = f15n
            abund(n13c,:) = f13c

            abund(nelec,:)=abund(ncx,:)+abund(nsix,:)+abund(nsx,:)+abund(nclx,:)+abund(nmgx,:)

            abund=abund*metallicity

            !Total H nuclei is always 1 so put fh into H and whatever is left over in H2
            abund(nh,:) = fh
            abund(nh2,:) = 0.5_dp*(1.0_dp-fh)
            abund(nd,:)=fd

            abund(nhe,:) = fhe
        end if
        abund(nSpec+2,:)=density      !Gas density
        abund(nSpec+1,:)=gasTemp    !Gas temperature
        !Initial calculations of diffusion and desorption frequencies
        !Uses updateVdiffAndVdes which supports both HH1992 and TST treatments
        call updateVdiffAndVdes(dustTemp(1), N_ICE_SPECIES, vdiff, vdes)

        ! Allocate conservation baseline array; values set by setConservationBaseline
        ! after starting abundances are loaded, so the baseline reflects actual initial state.
        if (ALLOCATED(initial_elem_abund)) deallocate(initial_elem_abund)
        allocate(initial_elem_abund(n_elem_tracked, points))
        initial_elem_abund = 0.0_dp

        call get_ion_indices(ionList, nIon)

        !DVODE SETTINGS
        ISTATE=1
        ITASK=1

        !set integration counts
        loopCounter=0
        failedIntegrationCounter=0
        solver_stats_counter=0  !Reset solver statistics counter

        if (.NOT. ALLOCATED(abstol)) then
            allocate(abstol(NEQ))
        end if
        if (.NOT. ALLOCATED(reltol_vec)) then
            allocate(reltol_vec(NEQ))
        end if

        if (heatingFlag) then
            !Initializing heating.f90 --> get coolants
            call initializeHeating(density(dstep),colDens(dstep),cloudSize,successFlag)
            if (successFlag < 0) return
        end if

        !Set rate constants to zero to ensure they don't hold previous values or random ones if we don't set them in calculateReactionRateConstants
        rate_constants=0.0

        !We typically don't recalculate rate constants that only depend on temperature if the temp hasn't changed
        !use arbitrarily high value to make sure they are calculated at least once.
        lastGasTemp=99.0e99_dp
        lastDustTemp=99.0e99_dp
        lastTemp=99.0e99_dp  ! Reset rate constants module lastTemp to force recalculation

        ! If the hydrogen diffusion energy is still its default value (-1.0 in default_parameters.f90),
        ! i.e. no custom value was set in the input dictionary, set it to the correct value
        ! according to the ratio of the ratio of diffusion energy to binding energy.
        if (is_equal(HdiffusionBarrier, -1.0_dp)) then
            species: do i = LBOUND(iceList, 1), UBOUND(iceList, 1)
                if (iceList(i) == ngh) then
                    HdiffusionBarrier = diffToBindRatio*bindingEnergy(i)
                    exit species
                end if
            end do species
        end if

        ! Pre-calculate desorption fractions for LHDES and ERDES reactions
        desorptionFractionsBare = 0.0_dp
        desorptionFractionsFullCoverage = 0.0_dp
        do j = lhdesReacs(1), lhdesReacs(2)
            desorptionFractionsBare(j) = getDesorptionFractionBare(j, j-lhdesReacs(1)+1)
            desorptionFractionsFullCoverage(j) = getDesorptionFractionFullCoverage(j, j-lhdesReacs(1)+1)
        end do
        do j = erdesReacs(1), erdesReacs(2)
            desorptionFractionsBare(j) = getDesorptionFractionBare(j, j-erdesReacs(1)+1)
            desorptionFractionsFullCoverage(j) = getDesorptionFractionFullCoverage(j, j-erdesReacs(1)+1)
        end do

    end subroutine initializeChemistry

    subroutine resetDVODEForNewPoint()
        ! Called at the start of each spatial point's time loop in the (points,time)
        ! loop order. Forces a fresh DVODE BDF restart and resets per-point counters.
        ! Setting ISTATE=1 is sufficient: the prevAbund abundance-change guard
        ! (integrateODESystem lines 427-436) is only evaluated when ISTATE=2.
        ISTATE = 1
        prevIntegrationTemp = 0.0_dp
        failedIntegrationCounter = 0
        solver_stats_counter = 0
    end subroutine resetDVODEForNewPoint

    subroutine updateChemistry(successFlag, statsarray, statsarray_size, dtime)
    !Updates the abundances for the next time step, first updating chemical variables and reaction rate constants,
    !then by solving the ODE system to obtain new abundances.
    !Solving ODEs is complex so we have two checks to try to automatically overcome difficulties and end stalled models
    !Firstly, the integration subroutine is called up to maxLoops times whilst adjusting variables to help integration converge.
    !If it succeeds before maxLoops, we continue as normal, otherwise we'll call it a fail.
    !Secondly, we check for stalls caused by the the solver loop reducing the targetTime to overcome difficulties.
    !That reduction the possibility of the code "succeeding" by integrating tiny target times. We have a counter that resets each time
    !the code integrates to the planned targetTime rather than a reduced one. If the counter reaches maxConsecutiveFailures, we end the code.
        !f2py integer, intent(aux) :: points
        integer, intent(out) :: successFlag
        real(dp), intent(inout), optional, dimension(:,:,:) :: statsarray
        integer, intent(in), optional :: statsarray_size
        integer, intent(in), optional :: dtime
        real(dp) :: originalTargetTime  !targetTime can be altered by integrator but we'd like to know if it was changed
        integer :: ie
        real(dp) :: rel_err
        real(dp) :: total_elem(n_elem_tracked)
        real(dp) :: surfaceCoverage
        real(dp) :: h2form_CT_vol, h2form_LH_vol, h2form_ER_vol  ! per-mechanism volumetric H2 formation rate constant [cm^-3 s^-1]
        real(dp) :: h2form_heat  ! mechanism-weighted H2 formation heating [erg cm^-3 s^-1]
        real(dp) :: h2heatfac, h2_denom  ! H&M79 eq. 6.45 thermalization efficiency factor

        ! write(*,*) "update Chemistry ..."
        !Integration can fail in a way that we can manage. Allow maxLoops tries before giving up.
        loopCounter=0
        successFlag=0
        originalTargetTime=targetTime
        do while((currentTime < targetTime) .and. (loopCounter < maxLoops))
            !allow option for dens to have been changed elsewhere.
            if (.not. freefall) abund(nSpec+2,dstep)=density(dstep)

            !First sum the total column density over all points further towards edge of cloud
            if (dstep>1) then
                h2ColToCell=(sum(abund(nh2,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                coColToCell=(sum(abund(nco,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                cColToCell=(sum(abund(nc,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
            else
                h2ColToCell=0.0_dp
                coColToCell=0.0_dp
                cColToCell=0.0_dp
            end if
            !then add half the column density of the current point to get average in this "cell"
            h2Col=h2ColToCell+0.5_dp*abund(nh2,dstep)*density(dstep)*(cloudSize/real(points))
            coCol=coColToCell+0.5_dp*abund(nco,dstep)*density(dstep)*(cloudSize/real(points))
            cCol=cColToCell+0.5_dp*abund(nc,dstep)*density(dstep)*(cloudSize/real(points))

            ! Postprocessed tracers have column densities provided
            if (lusecoldens) then
               h2col = lnh2(tstep)
               cocol = lnco(tstep)
               ! ccol = lnc(dstep, tstep) ! TODO enable C column density support
               ccol = lnh(tstep) * abund(nc,dstep)  ! No C column densities yet...
            end if

            !Reset surface and bulk values in case of integration error or sputtering
            ! Skip when continuing DVODE (ISTATE=2): recomputing from clamped individual
            ! species (each at MIN_ABUND) gives SUM >> DVODE's nBulk, breaking continuation.
            if (ISTATE /= 2) then
                abund(nBulk,dstep)=sum(abund(bulkList,dstep))
                abund(nSurface,dstep)=sum(abund(surfaceList,dstep))
            end if
            !recalculate coefficients for ice processes
            safeMantle=MAX(MIN_ABUND, abund(nSurface,dstep))
            safeBulk=MAX(MIN_ABUND, abund(nBulk,dstep))

            if (refractoryList(1) > 0) safeBulk=safeBulk-SUM(abund(refractoryList,dstep))

            ratioSurfaceToBulk=MIN(1.0_dp, safeMantle/safeBulk)
            bulkLayersReciprocal=MIN(1.0_dp,NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
            surfaceCoverage=bulkGainFromMantleBuildUp()

            if ((.not. is_equal(dustTemp(dstep), lastDustTemp)) .OR. &
                (.not. is_equal(gasTemp(dstep), lastGasTemp))) then
                call updateVdiffAndVdes(dustTemp(dstep), N_ICE_SPECIES, vdiff, vdes)
            end if

            call calculateReactionRateConstants(abund,safeMantle, h2col, cocol, ccol, rate_constants)
            if (heatingFlag) then
                dustTemp(dstep)=calculateDustTemp( &
                    radfield*EXP(-UV_FAC*av(dstep)) + radfield_internal(dstep)*EXP(-UV_FAC*av_internal(dstep)), &
                    radfield + radfield_internal(dstep), &
                    av(dstep), zeta)


                ! Per-mechanism volumetric H2 formation rate constants [cm^-3 s^-1]
                ! Only accounting for H2 ending up in the gas phase.
                h2form_CT_vol = rate_constants(nR_H2Form_CT) * abund(nSpec+2,dstep)**2 * abund(nh,dstep)
                h2form_LH_vol = (rate_constants(nR_H2Form_LHDes)) &
                              &  * abund(ngh,dstep)**2 * abund(nSpec+2,dstep)
                h2form_ER_vol = (rate_constants(nR_H2Form_ERDes)) &
                              &  * abund(nSpec+2,dstep)**2 * abund(nh,dstep) * abund(ngh,dstep) / safeMantle
                ! H&M79 eq. 6.45: critical density for H2 thermalization
                ! (18100 coefficient for consistency with h2FUVPumpHeating in heating.f90)
                h2_denom = 1.6_dp*abund(nh,dstep)*EXP(-((400.0_dp/gasTemp(dstep))**2)) &
                         &+ 1.4_dp*abund(nh2,dstep)*EXP(-(18100.0_dp/(gasTemp(dstep)+1200.0_dp)))
                if (h2_denom > 0.0_dp) then
                    h2heatfac = 1.0_dp / (1.0_dp + 1.0e6_dp/(SQRT(gasTemp(dstep))*h2_denom*abund(nSpec+2,dstep)))
                else
                    h2heatfac = 0.0_dp
                end if
                ! H&M79 eq. 6.43: LH gives 0.1 eV kinetic + 4.2 eV vibrational (fraction h2heatfac goes to gas)
                ! ER: 0.6 eV (Bourlot et al. 2012), thermalization-corrected
                ! CT: 1.5 eV (Hollenbach & Tielens 1999), no thermalization correction
                h2form_heat = eV * (1.5_dp*h2form_CT_vol &
                            &+ (0.1_dp + 4.2_dp*h2heatfac)*h2form_LH_vol &
                            &+ 0.6_dp*h2heatfac*h2form_ER_vol)
                tempDot= getTempDot(&
                                    abund(nSpec+1,dstep), &              ! gas temperature
                                    abund(nSpec+2,dstep), &              ! gas density
                                    colDens(dstep), &                    ! gas column density
                                    radfield*EXP(-UV_FAC*av(dstep)) + radfield_internal(dstep)*EXP(-UV_FAC*av_internal(dstep)), &   ! attenuated radiation field
                                    abund(:,dstep), &                    ! full abundance vector
                                    rate_constants(nR_H2_hv), &          ! H2 dissociation rate constant
                                    h2form_heat, &                       ! mechanism-weighted H2 formation heating [erg cm^-3 s^-1]
                                    zeta, &                              ! cosmic ray ionization rate
                                    rate_constants(nR_C_hv), &           ! C-photo rate
                                    1.0_dp/GAS_DUST_DENSITY_RATIO, &     ! dust-to-gas ratio
                                    GRAIN_RADIUS, &                      ! grain radius
                                    metallicity, &                       ! metallicity
                                !    heatWriteFlag, &                    ! write flag
                                    dusttemp(dstep), &                   ! dust temperature
                                    turbVel)                             ! turbulence velocity
            end if

            !Integrate_constants chemistry, and return fail if unrecoverable error was reached
            if (PRESENT(statsarray) .AND. PRESENT(statsarray_size) .AND. PRESENT(dtime)) then
                call integrateODESystem(successFlag, statsarray, statsarray_size, dtime)
            else
                call integrateODESystem(successFlag)
            end if
            if (successFlag < 0) then
                write(*,*) "Integration failed, exiting"
                return
            end if

            !1.e-30_dp stops numbers getting too small for fortran.
            ! WHERE(abund<MIN_ABUND) abund=MIN_ABUND
            where(abund(1:nSpec,:)<MIN_ABUND) abund(1:nSpec,:)=MIN_ABUND
            gasTemp(dstep)=abund(nSpec+1,dstep)
            density(dstep)=abund(nSpec+2,dstep)
            ! IF (gasTemp(dstep) .lt. 10) gasTemp(dstep)=10.0
            if (gasTemp(dstep) < lower_limit_gastemp) gasTemp(dstep)=lower_limit_gastemp
            loopCounter=loopCounter+1

            ! For postprocessing, force solver to try and reach original target time
            if (usepostprocess) targettime = originaltargettime
        end do

        ! Postprocessing needs to reach next timestep whatever the cost
        if (.not. usepostprocess) then
            if (loopCounter == maxLoops) successFlag=INT_TOO_MANY_FAILS_ERROR

            !Since targetTime can be altered, eventually leading to "successful" integration we want to
            !check if integrator ever just reaches the planned target time. If it doesn't for many attempts,
            !we will call the run a failure. This stops the target being constantly reduced to tiny increments
            !so that the code all but stalls as the time is increased by seconds each integraiton.
            if (ABS(originalTargetTime- targetTime) < 0.001_dp*originalTargetTime) then
                failedIntegrationCounter=0
            else
                failedIntegrationCounter=failedIntegrationCounter+1
            end if
            if (failedIntegrationCounter > maxConsecutiveFailures) then
                successFlag=INT_TOO_MANY_FAILS_ERROR
            end if
        end if

        ! Runtime element conservation check (every iteration, not inside F)
        if (runtime_conservation_tolerance >= 0.0_dp .AND. successFlag == 0) then
            total_elem = calculate_elemental_abundances(abund(:, dstep))
            do ie = 1, n_elem_tracked
                if (initial_elem_abund(ie, dstep) > 0.0_dp) then
                    rel_err = ABS(total_elem(ie) - initial_elem_abund(ie, dstep)) &
                            & / initial_elem_abund(ie, dstep)
                    if (rel_err > runtime_conservation_tolerance) then
                        write(*,"(A,A2,A,ES10.3,A,ES12.4,A)") &
                            "CONSERVATION ERROR: element ", TRIM(elem_names(ie)), &
                            " changed by ", rel_err*100.0_dp, "% at t=", &
                            currentTime/SECONDS_PER_YEAR, " yr"
                        successFlag = CONSERVATION_ERROR
                        return
                    end if
                end if
            end do
        end if

    end subroutine updateChemistry

    subroutine integrateODESystem(successFlag, statsarray, statsarray_size, dtime)
        integer, intent(out) :: successFlag
        real(dp), intent(inout), optional, dimension(:,:,:) :: statsarray
        integer, intent(in), optional :: statsarray_size
        integer, intent(in), optional :: dtime
        type(VODE_OPTS), save :: OPTIONS  ! SAVE: persists across ISTATE=2 continuation calls
        real(dp), save :: prevAbund(nSpec)  ! end-state snapshot of last ISTATE=1 step; safe: only read when ISTATE=2, which requires a prior successful call that sets this array
        real(dp) :: maxLogChange
        logical :: was_fresh_restart       ! whether this call entered DVODE with ISTATE=1
        ! integer :: ii
        ! species_check_mask was used by the negative-abundance error block (now commented out)
        !LOGICAL :: species_check_mask(nSpec)
        successFlag=0
        f_callback_error=0

    !This subroutine calls DVODE (3rd party ODE solver) until it can reach targetTime with acceptable errors (reltol/abstol)
        !reset parameters for DVODE
        ITASK=1  !try to integrate to targetTime
        if (solverMode == 0) then
            ISTATE = 1                   ! mode 0: always restart fresh
        else
            if (ISTATE < 0) ISTATE = 1  ! reset on solver error; ISTATE=2 carries BDF history forward
            ! Temperature guard: restart if temperature changed significantly at output step boundary
            if (ISTATE == 2) then
                if (ABS(gasTemp(dstep) - prevIntegrationTemp) > 1.0_dp) ISTATE = 1
            end if
            ! Abundance-change guard (mode 2 only): restart if chemistry evolved rapidly since last call.
            ! Large per-step changes mean the frozen abstol and BDF Jacobian are stale.
            if (solverMode == 2 .AND. ISTATE == 2) then
                maxLogChange = MAXVAL( &
                    ABS(LOG10(MAX(abund(1:nSpec,dstep), MIN_ABUND)) - &
                        LOG10(MAX(prevAbund,            MIN_ABUND))), &
                    MASK = abund(1:nSpec,dstep) > MIN_ABUND .AND. &
                           prevAbund            > MIN_ABUND)
                if (maxLogChange > logChangeThreshold) then
                    ISTATE = 1
                end if
            end if
        end if
        ! Setup tolerances and options only on fresh start or error recovery
        if (ISTATE <= 1) then
            !Gas-phase species: absolute tolerances scaled by abundance
            abstol=abstol_factor*abund(:,dstep)
            where(abstol<abstol_min) abstol=abstol_min
            !Ice species (surface + bulk): separate_constants, looser absolute tolerances
            abstol(iceList) = abstol_ice_factor*abund(iceList,dstep)
            where(abstol(iceList)<abstol_ice_min) abstol(iceList)=abstol_ice_min
            !Physical variables: separate_constants tolerance heuristic (T and nH need looser tolerances)
            abstol(nSpec+1) = MAX(abstol_phys_factor * ABS(abund(nSpec+1,dstep)), abstol_T_min)
            abstol(nSpec+2) = MAX(abstol_phys_factor * ABS(abund(nSpec+2,dstep)), abstol_nH_min)
            !Per-component relative tolerances: tight for chemistry, relaxed for physics
            reltol_vec(1:nSpec) = reltol
            reltol_vec(nSpec+1) = reltol_phys
            reltol_vec(nSpec+2) = reltol_phys
            !Call the integrator with ITOL=4 (vector reltol + vector abstol).
            OPTIONS = SET_OPTS(METHOD_FLAG=22, ABSERR_VECTOR=abstol, RELERR_VECTOR=reltol_vec, &
                               USER_SUPPLIED_JACOBIAN=.false.,MXSTEP=MXSTEP)
        end if
        ! Track whether this call enters DVODE as a fresh restart (ISTATE=1).
        ! prevAbund is saved AFTER DVODE succeeds (see success block below), so the guard
        ! compares against the END-state of the last ISTATE=1 step, not the start-state.
        ! This guarantees the first ISTATE=2 step after any ISTATE=1 always sees
        ! maxLogChange=0 (free pass), preventing the immediate re-fire that occurred when
        ! prevAbund was saved before DVODE (change DURING the restart was measured instead
        ! of cumulative drift SINCE the restart).
        was_fresh_restart = (ISTATE == 1)
        call CPU_TIME(dvode_cpu_start)
        call DVODE_F90(F,NEQ,abund(:,dstep),currentTime,targetTime,ITASK,ISTATE,OPTIONS)
        call CPU_TIME(dvode_cpu_end)
        dvode_cpu_time = dvode_cpu_end - dvode_cpu_start
        dvode_istate_out = ISTATE
        call GET_STATS(dvode_rstats, dvode_istats)

        ! Write solver statistics immediately after EVERY DVODE call (including failures)
        if (PRESENT(statsarray) .AND. PRESENT(statsarray_size) .AND. PRESENT(dtime)) then
            solver_stats_counter = solver_stats_counter + 1

            ! Check for array overflow
            if (solver_stats_counter > statsarray_size) then
                write(*,*) "ERROR: Solver stats array overflow at counter", solver_stats_counter
                write(*,*) "       Allocated size:", statsarray_size
                write(*,*) "       Consider increasing statsarray allocation or reducing finalTime"
                successFlag = SOLVER_STATS_OVERFLOW_ERROR
                return
            end if

            ! Write stats: column 1 = trajectory index, rest shifted by 1
            statsarray(solver_stats_counter, dstep, 1) = DBLE(dtime)
            statsarray(solver_stats_counter, dstep, 2) = DBLE(dvode_istate_out)
            statsarray(solver_stats_counter, dstep, 3:6) = dvode_rstats(11:14)
            statsarray(solver_stats_counter, dstep, 7:18) = DBLE(dvode_istats(11:22))
            statsarray(solver_stats_counter, dstep, 19) = dvode_cpu_time
        end if

        if (f_callback_error < 0) then
            successFlag = f_callback_error
            return
        end if

        ! Between-step physical sanity check: after DVODE returns, verify that no
        ! species has a genuinely diverging negative abundance. Tiny negatives from
        ! solver numerics are expected and will be clamped by the WHERE below; we
        ! only abort if a species is significantly negative (beyond negative_abundance_tol).
        ! nBulk and nSurface are pseudo-species (aggregates of @xxx and #xxx) whose DVODE
        ! derivative is set to the sum of component derivatives BEFORE the mantle-retreat
        ! block updates those components. Their values therefore drift from the true component
        ! sum and can go slightly negative. Exclude them from the divergence check and
        ! recompute them from their components after clamping the real species.
        if (ISTATE >= 2) then
            where(abund(1:nSpec,dstep) < MIN_ABUND) abund(1:nSpec,dstep) = MIN_ABUND
            ! Do NOT recompute nBulk/nSurface from sum(clamped bulkList):
            ! at early times that jumps nBulk by ~N_species orders of magnitude,
            ! breaking DVODE's BDF history. Direct clamp is within abstol_ice_min=1e-20.
            ! Save end-state of fresh restarts for the abundance-change guard.
            ! Using the end-state (not start-state) means the first ISTATE=2 step after
            ! any ISTATE=1 always sees maxLogChange=0, measuring real cumulative drift
            ! from this point onward.
            if (was_fresh_restart) prevAbund = abund(1:nSpec, dstep)
            prevIntegrationTemp = gasTemp(dstep)
        end if

        select case(ISTATE)
            case(-1)
                !ISTATE -1 means the integrator can't break the problem into small enough steps
                !We could increase MXSTEP but better to reduce targetTime and get to physics update
                !physical conditions may be easier to solve as time goes by so better to get to that update
                write(*,*) "ISTATE -1: Reducing time step"
                !More steps required for this problem
                !MXSTEP=MXSTEP*2
                targetTime=currentTime+(targetTime-currentTime)*0.1_dp
            case(-2)
                !ISTATE -2 just needs an absol change so let's do that and try again
                write(*,*) "ISTATE -2: Tolerances too small"
                !Tolerances are too small for machine but successful to current currentTime
                abstol_factor=abstol_factor*10.0_dp
                abstol_ice_factor=abstol_ice_factor*10.0_dp
                reltol_phys=MIN(reltol_phys*10.0_dp, 1e-1_dp)
            case(-3)
                !ISTATE -3 is unrecoverable so just bail on integration
                write(*,*) "DVODE found invalid inputs"
                write(*,*) "abstol:"
                write(*,*) abstol
                successFlag=INT_UNRECOVERABLE_ERROR
                return
            case(-4)
                !Successful as far as currentTime but many errors.
                !Make targetTime smaller and just go again
                write(*,*) "ISTATE -4 - shortening step"
                targetTime=currentTime+(targetTime-currentTime)*0.1_dp
            case(-5)
                write(*,*) "ISTATE -5 - shortening step at time", timeInYears,"years"
                targetTime=currentTime+(targetTime-currentTime)*0.1_dp
            case default
                ! Success: MXSTEP stays at whatever param_dict set (do not reset to hardcoded 10000)
        end select

    if (enforceChargeConservation) then
        ! REALLY ensure charge is always conserved (also after integrating)
        abund(nelec,dstep) = sum(abund(ionlist(1:nion),dstep))
    end if
    end subroutine integrateODESystem

    subroutine F (NEQUATIONS, T, Y, YDOT)
        use ODES, only: GETYDOT
        integer, parameter :: WP = KIND(1.0D0)

        integer, intent(in)  :: NEQUATIONS
        real(WP), intent(in) :: T
        real(WP), dimension(NEQUATIONS), intent(in) :: Y
        real(WP), dimension(NEQUATIONS), intent(out) :: YDOT

        real(dp) :: D, gasTemperature
        real(dp) :: surfaceCoverage
        real(dp) :: h2heatfac, h2_denom  ! H&M79 eq. 6.45 thermalization efficiency factor

        integer :: k
        ! Y_safe clamps species abundances to MIN_ABUND during ODE evaluation.
        ! DVODE predictor steps can drive species to small negatives; feeding those
        ! negative values back into destruction terms compounds the overshoot.
        ! Clamping here keeps the RHS physical without altering the accepted step.
        real(WP), dimension(NEQUATIONS) :: Y_safe
        !Set D to the gas density for use in the ODEs
        D=y(nSpec+2)     !Gas density
        gasTemperature=y(nSpec+1)     !Gas temperature
        ydot=0.0_dp

        Y_safe = Y
        ! where(Y_safe(1:nSpec) < MIN_ABUND) Y_safe(1:nSpec) = MIN_ABUND

        ! Column densities are fixed for postprocessing data, so don't do this bit
        if (.not. lusecoldens) then
            !changing abundances of H2 and CO can causes oscillation since their rate_constants depend on their abundances
            !recalculating the rate constants as abundances are updated prevents that.
            !thus these are the only rate constants calculated each time the ODE system is called.
            cocol=coColToCell + 0.5_dp*Y_safe(nco)*D*(cloudSize/real(points))
            h2col=h2ColToCell + 0.5_dp*Y_safe(nh2)*D*(cloudSize/real(points))
            rate_constants(nR_H2_hv)=getH2PhotoDissRateConstant(h2Col,radField,av(dstep),turbVel)  !H2 photodissociation
            rate_constants(nR_CO_hv)=getCOPhotoDissRateConstant(h2Col,coCol,radField,av(dstep))  !CO photodissociation
        end if

        !recalculate coefficients for ice processes
        safeMantle=MAX(MIN_ABUND, Y_safe(nSurface))
        safeBulk=MAX(MIN_ABUND, Y_safe(nBulk))
        bulkLayersReciprocal=MIN(1.0_dp, NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
        surfaceCoverage=bulkGainFromMantleBuildUp()
        ratioSurfaceToBulk=MIN(1.0_dp, safeMantle/safeBulk)

        ! Fix 3: refresh surface-to-bulk swap rate_constants from current safeMantle
        ! (safeMantle was just updated from Y_safe above, but rate_constants(surfSwapReacs) is still
        ! set from the start-of-step call to calculateReactionRateConstants)
        if (THREE_PHASE) rate_constants(surfSwapReacs(1):surfSwapReacs(2)) = surfaceToBulkSwappingRateConstants(dustTemp(dstep))

        ! Fix 1: re-split LH/LHDES and ER/ERDES using the current ice thickness.
        ! desorptionFractionIncludingIce depends on numMonolayers which changes as ice builds up
        ! during DVODE integration. Without this, the fraction is frozen at the start-of-step value.
        numMonolayers = getNumberMonolayers(safeMantle + safeBulk)
        if (lhdesReacs(1) /= REAC_NOT_PRESENT .AND. desorb .AND. chemdesorb &
            .AND. dustTemp(dstep) < maxGrainTemp                            &
            .AND. safeMantle > MIN_SURFACE_ABUND) then
            k = 0
            do i = lhdesReacs(1), lhdesReacs(2)
                k = k + 1
                rate_constants(i) = getDesorptionFractionIncludingIce(i, numMonolayers) &
                          * rate_constants_lh_unsplit(LHDEScorrespondingLHreacs(k))
                if (ANY(bulkList==re1(i))) rate_constants(i) = 0.0
            end do
            k = 0
            do i = lhdesReacs(1), lhdesReacs(2)
                k = k + 1
                rate_constants(LHDEScorrespondingLHreacs(k)) = rate_constants_lh_unsplit(LHDEScorrespondingLHreacs(k)) &
                    - rate_constants(i)
            end do
        end if
        if (erdesReacs(1) /= REAC_NOT_PRESENT .AND. desorb .AND. chemdesorb &
            .AND. dustTemp(dstep) < maxGrainTemp                            &
            .AND. safeMantle > MIN_SURFACE_ABUND) then
            k = 0
            do i = erdesReacs(1), erdesReacs(2)
                k = k + 1
                rate_constants(i) = getDesorptionFractionIncludingIce(i, numMonolayers) &
                          * rate_constants_er_unsplit(ERDEScorrespondingERreacs(k))
                if (ANY(bulkList==re1(i))) rate_constants(i) = 0.0
            end do
            k = 0
            do i = erdesReacs(1), erdesReacs(2)
                k = k + 1
                rate_constants(ERDEScorrespondingERreacs(k)) = rate_constants_er_unsplit(ERDEScorrespondingERreacs(k)) -&
                    rate_constants(i)
            end do
        end if

        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens.
        !Each species ODE is made up of the reactions between it and every other species it reacts with.
        call GETYDOT(rate_constants, Y_safe, surfaceCoverage, D, YDOT, surfGrowthUncorrected)

        ydot(nSpec+2) = get_densdot(Y(nSpec+2))     !Gas density ODE

        if (enforceChargeConservation) then
            ydot(nelec) = sum(ydot(ionlist(1:nion)))
        end if


        if (heatingFlag) then
            ! Species abundances in Y_safe are already clamped; Y used below only
            ! for temperature (nSpec+1) and density (nSpec+2), which are never negative.
            ! Write(*,*) "Updating heating and cooling rate_constantss"
            if (ABS(gasTemperature-oldTemp)>MIN(heating_temp_abstol, heating_temp_reltol*oldTemp)) then
                gasTemp(dstep)=gasTemperature
                if (gasTemp(dstep) < lower_limit_gastemp) gasTemp(dstep)=lower_limit_gastemp
                if (gasTemp(dstep) > upper_limit_gastemp) gasTemp(dstep)=upper_limit_gastemp
                ! Fix 2: update gas-phase two-body rate_constantss for the new temperature.
                ! These rate_constantss are frozen at start-of-step in calculateReactionRateConstants.
                rate_constants(twobodyReacs(1):twobodyReacs(2)) = &
                    alpha(twobodyReacs(1):twobodyReacs(2)) * &
                    ((gasTemp(dstep)/300.0_dp)**beta(twobodyReacs(1):twobodyReacs(2))) * &
                    exp(-gama(twobodyReacs(1):twobodyReacs(2))/gasTemp(dstep))
                ! H&M79 eq. 6.45: critical density for H2 thermalization
                ! (18100 coefficient for consistency with h2FUVPumpHeating in heating.f90)
                h2_denom = 1.6_dp*Y(nh)*EXP(-((400.0_dp/gasTemperature)**2)) &
                         &+ 1.4_dp*Y(nh2)*EXP(-(18100.0_dp/(gasTemperature+1200.0_dp)))
                if (h2_denom > 0.0_dp) then
                    h2heatfac = 1.0_dp / (1.0_dp + 1.0e6_dp/(SQRT(gasTemperature)*h2_denom*D))
                else
                    h2heatfac = 0.0_dp
                end if
                ! Only desorbing products contribute to gas heating; LH/ER remain on grain.
                ! H&M79 eq. 6.43: LHDes gives 0.1 eV kinetic + 4.2 eV vibrational (fraction h2heatfac to gas)
                ! ERDes: 0.6 eV (Bourlot et al. 2012), thermalization-corrected
                ! CT: 1.5 eV (Hollenbach & Tielens 1999), no thermalization correction
                h2form = eV * ( &
                    &  1.5_dp * rate_constants(nR_H2Form_CT) * D**2 * Y(nh) &
                    &+ (0.1_dp + 4.2_dp*h2heatfac) * rate_constants(nR_H2Form_LHDes) * Y(ngh)**2 * D &
                    &+ 0.6_dp * h2heatfac * rate_constants(nR_H2Form_ERDes) * D**2 * Y(nh) * Y(ngh) &
                    &  / max(safeMantle, MIN_SURFACE_ABUND))
                tempDot=getTempDot( &
                               gasTemperature, &                          ! gas temperature
                               D, &                          ! gas density
                               colDens(dstep), &                      ! gas column density
                               radfield*EXP(-UV_FAC*av(dstep)) + radfield_internal(dstep)*EXP(-UV_FAC*av_internal(dstep)), &     ! attenuated radiation field
                               Y, &                                   ! all number densities
                               rate_constants(nR_H2_hv), &            ! H2 dissociation rate_constants computed in rate_constantss.f90
                               h2form, &                              ! mechanism-weighted H2 formation heating [erg cm^-3 s^-1]
                               zeta, &                                ! cosmic ray ionization rate_constants
                               rate_constants(nR_C_hv), &             ! C-photo rate_constants
                               1.0_dp/GAS_DUST_DENSITY_RATIO, &       ! dust-to-gas ratio
                               GRAIN_RADIUS, &                        ! grain radius
                               metallicity, &                         ! metallicity
                               dusttemp(dstep), &                     ! dust temperature
                               turbVel)                               ! turbulence velocity
                oldTemp=gasTemperature
            end if
            ydot(nSpec+1)=tempDot
        else
            ydot(nSpec+1)=0.0
        end if

    end subroutine F

    ! SUBROUTINE JAC(NEQ, T, Y, ML, MU, J, NROWPD)
    !     INTEGER NEQ,ML,MU,NROWPD
    !     DOUBLE PRECISION T, Y(NEQ), J(NROWPD,NEQ)
    !     REAL(DP) :: D
    !     INTENT(IN)  :: NEQ, T, Y,ML,MU,NROWPD
    !     INTENT(INOUT) :: J
    !     D=y(NEQ)

    !     J=0.0_dp
    !     INCLUDE 'jacobian.f90'
    !     J(nh,nh2)=J(nh,nh2)+2.0*h2dis
    !     J(nh2,nh2)=J(nh,nh2)-h2dis
    ! END SUBROUTINE JAC

    pure function calculate_elemental_abundances(abundances) result(elemental_abundances)
        real(dp), intent(in) :: abundances(nSpec)

        real(dp) :: elemental_abundances(n_elem_tracked)

        integer :: element_idx

        do element_idx = 1, n_elem_tracked
            elemental_abundances(element_idx) = &
                sum(real(elem_count(:, element_idx), dp) * abundances)
        end do
    end function calculate_elemental_abundances

    ! Set the conservation baseline for all parcels from the current abund array.
    ! Must be called after all starting abundances have been loaded into abund.
    subroutine setConservationBaseline()
        integer :: point_idx

        do point_idx = 1, points
            call resetConservationBaselineForPoint(point_idx)
        end do
    end subroutine setConservationBaseline

    ! Reset the conservation baseline for a single parcel (used in 1D per-parcel loop).
    ! Must be called after the parcel's starting abundances have been loaded into abund.
    subroutine resetConservationBaselineForPoint(point_idx)
        integer, intent(in) :: point_idx

        initial_elem_abund(:, point_idx) = calculate_elemental_abundances(abund(:, point_idx))
    end subroutine resetConservationBaselineForPoint

end module chemistry
