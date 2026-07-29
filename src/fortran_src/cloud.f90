! Cloud.f90 is our most basic physics module
! It simulates a static or collapsing cloud of isothermal gas.
! It is useful for static models or for producing initial abundances for the other modules.
module cloud_mod
    use constants, only: dp, SECONDS_PER_YEAR, PC, ZERO_INNER_RADIUS_ERROR
    use DEFAULTPARAMETERS
    use f2py_constants, only: nSpec
    !f2py INTEGER, parameter :: dp
    use physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
        zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, &
        gastemp, dusttemp, av, coldens, density_max, get_ngas_r, get_initialDens_r, findcoldens_edge2core, &
        get_coldens_external, initialDens_array, parcel_radius, outer_coldens_for_current_step, get_av

    implicit none

    private
    public :: updateTargetTime, initializePhysics, updatePhysics, sublimation

    real(dp), allocatable :: parcelRadius(:)
    real(dp), allocatable :: coldens_obs(:)

    ! Time sampling control parameters
    real(dp) :: timestep_resolution_factor_early = 0.5_dp   ! For 0 < t < 10 yr: samples per decade (snapped to k*10^n grid)
    real(dp) :: timestep_resolution_factor_mid = 1.0_dp      ! For 10 yr < t < 1 Myr: dt = 10^floor(log_10(t)) / factor
    real(dp) :: timestep_fixed_late_years = 1.0e5_dp            ! For t > 1 Myr: dt = fixed timestep in years

    ! Radial grid spacing: .false. = linear (default), .true. = logarithmic
    logical :: log_radius_sampling = .false.

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called at start of UCLCHEM run
    ! Uses values in defaultparamters.f90 and any inputs to set initial values        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializePhysics(successFlag)
        integer, intent(out) :: successFlag
        successFlag=0

        !Set up basic physics variables
        cloudSize=(rout-rin)*PC

        if (enable_radiative_transfer .AND. points>1 .AND. rin <= 0.0_dp) then
            write(*,*) "ERROR: rin must be > 0 when enable_radiative_transfer=True (innermost parcel would be at r=0)"
            successFlag=ZERO_INNER_RADIUS_ERROR
            return
        end if

        ! Allocate 1D arrays if 1D radiative transfer is enabled
        if (enable_radiative_transfer .AND. points>1) then
            if (ALLOCATED(parcelRadius)) deallocate(parcelRadius)
            allocate(parcelRadius(points))

            if (ALLOCATED(coldens_obs)) deallocate(coldens_obs)
            allocate(coldens_obs(points))

            do dstep=1,points
                if (points > 1) then
                    if (log_radius_sampling) then
                        parcelRadius(dstep)=rin*(rout/rin)**((float(dstep)-1.0e0_dp)/float(points-1))
                    else
                        parcelRadius(dstep)=rin+(dstep-1)*(rout-rin)/float(points-1)
                    end if
                else
                    parcelRadius(dstep)=rout
                end if
                parcel_radius(dstep)=parcelRadius(dstep)
                density_max(dstep)=get_ngas_r(parcelRadius(dstep),finalDens,density_scale_radius,density_power_index)
            end do
        end if

        ! Set up densities (use radial profile if 1D radiative transfer enabled)
        do dstep=1,points
            if (enable_radiative_transfer .AND. points>1) then
                ! Store the raw profile density (without 1.001 bump) so densdot
                ! sees density(dstep) > initialDens_array(dstep) and freefall fires.
                initialDens_array(dstep)=get_initialDens_r(parcelRadius(dstep)*PC,density_power_index)
                density(dstep)=1.001_dp*initialDens_array(dstep)
            else
                density(dstep)=1.001_dp*initialDens
            end if
        end do

        do dstep=1,points
            if (enable_radiative_transfer .AND. points>1) then
                coldens(dstep) = get_coldens_external(parcelRadius(dstep), initialDens)
            else
                coldens(dstep) = real(points-dstep+1)*cloudSize/real(points)*initialDens
            end if
            av(dstep) = get_av(baseAv, coldens(dstep))
        end do

    end subroutine initializePhysics

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine updateTargetTime
        real(dp) :: orderMagnitude, stepSize
        if (timeInYears >= 1.0e6_dp) then
            ! Beyond 1 Myr: use fixed time step
            targetTime = (timeInYears + timestep_fixed_late_years) * SECONDS_PER_YEAR
        else if (timeInYears > 10.0_dp) then
            ! Between 10 years and 1 Myr: linear steps within each decade
            ! Step size is one order of magnitude smaller, divided by steps_per_decade
            orderMagnitude = 10.0_dp**(floor(log10(timeInYears)))
            stepSize = orderMagnitude / timestep_resolution_factor_mid
            targetTime = (timeInYears + stepSize) * SECONDS_PER_YEAR
        else if (timeInYears > 0.0_dp) then
            ! Below 10 years: logarithmic sampling snapped to a k*10^n grid.
            ! orderMagnitude is the decade floor (e.g. 1e-3 when t is in [1e-3, 1e-2)).
            ! stepSize = orderMagnitude / factor gives exactly `factor` steps per decade.
            ! Snapping with floor(t/stepSize)+1 ensures targets land on exact multiples
            ! of stepSize, so decade boundaries (1, 10, 100 ...) are always hit cleanly.
            ! The 1e-10 epsilon guards against floating-point when t is already on a grid point.
            orderMagnitude = 10.0_dp**(floor(log10(timeInYears)))
            stepSize = orderMagnitude / timestep_resolution_factor_early
            targetTime = (floor(timeInYears / stepSize + 1.0e-10_dp) + 1.0_dp) * stepSize * SECONDS_PER_YEAR
        else
            ! Initial timestep: start at the bottom of the first sampled decade
            targetTime = SECONDS_PER_YEAR * 1.0e-7_dp
        end if
    end subroutine updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine updatePhysics
        ! Calculate column densities and update physical parameters for 1D models
        if (enable_radiative_transfer .AND. points>1) then
            call findcoldens_edge2core(coldens_obs(dstep),finalDens,density_scale_radius,density_power_index,parcelRadius(dstep))

            ! Calculate density profiles
            density_max(dstep)=get_ngas_r(parcelRadius(dstep),finalDens,density_scale_radius,density_power_index)

            if (density(dstep)>=density_max(dstep)) density(dstep)=density_max(dstep)

            ! Recompute the coldens from physics-core.f90
            ! coldens should be amount of gas from edge to parcel
            coldens(dstep)=cloudSize/real(points)*density(dstep)

            ! Add previous column densities to current as we move into cloud to get total.
            ! outer_coldens_for_current_step is set by wrap.f90 from coldens_history(dtime)
            ! before this subroutine is called, so the value is exact (no approximation).
            if (dstep < points) then
                coldens(dstep)=coldens(dstep)+outer_coldens_for_current_step
            end if

            av(dstep) = get_av(baseAv, coldens(dstep))

            ! Diagnostic output
            ! print '(A,1PE12.3,A,0PF8.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3)', &
            ! 't=',timeInYears, '  r=',parcelRadius(dstep), '  nH(t)=',density(dstep), '  nH(r)_max=',density_max(dstep), &
            ! '  NH(t)(edge->core)=',coldens(dstep), '  NH(r)(edge->core)=', coldens_obs(dstep), '  av(t)=',av(dstep)
        end if
    end subroutine updatePhysics

    !This is a dummy subroutine.
    subroutine sublimation(abund, lpoints)
        integer, intent(in) :: lpoints
        real(dp), intent(inout) :: abund(nspec+1,lpoints)
    end subroutine sublimation
end module cloud_mod


