! Cloud.f90 is our most basic physics module
! It simulates a static or collapsing cloud of isothermal gas.
! It is useful for static models or for producing initial abundances for the other modules.
MODULE cloud_mod
    USE constants
    USE DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp
    USE physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
    & zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, densdot, gastemp, dusttemp, av,&
    &coldens, density_max, ngas_r, initialDens_r, findcoldens_edge2core
    USE network
    use f2py_constants
    IMPLICIT NONE
    REAL(dp), allocatable :: parcelRadius(:)
    REAL(dp), allocatable :: coldens_obs(:)
CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called at start of UCLCHEM run
    ! Uses values in defaultparamters.f90 and any inputs to set initial values        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        successFlag=0

        !Set up basic physics variables
        cloudSize=(rout-rin)*pc

        ! Allocate 1D arrays if 1D radiative transfer is enabled
        IF (enable_radiative_transfer .AND. points.gt.1) THEN
            IF (ALLOCATED(parcelRadius)) DEALLOCATE(parcelRadius)
            ALLOCATE(parcelRadius(points))

            IF (ALLOCATED(coldens_obs)) DEALLOCATE(coldens_obs)
            ALLOCATE(coldens_obs(points))

            DO dstep=1,points
                parcelRadius(dstep)=dstep*rout/float(points) !unit of parsec -- Note: from core to edge
            END DO
            
            density_max=ngas_r(rin,finalDens,density_scale_radius,density_power_index)
        END IF

        ! Set up densities (use radial profile if 1D radiative transfer enabled)
        DO dstep=1,points
            IF (enable_radiative_transfer .AND. points.gt.1) THEN
                density(dstep)=1.001*initialDens_r(parcelRadius(dstep)*pc,density_power_index)
            ELSE
                density(dstep)=1.001*initialDens
            END IF
        END DO
        
        DO dstep=1,points
            IF (enable_radiative_transfer .AND. points.gt.1) THEN
                coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens_r(parcelRadius(dstep)*pc,density_power_index)
            ELSE
                coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
            END IF
        END DO

    END SUBROUTINE initializePhysics

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE updateTargetTime
        real(dp) :: orderMagnitude, currentValue
        IF (timeInYears .ge. 1.0d6) THEN
            targetTime=(timeInYears+1.0d5)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 10.0) THEN
            orderMagnitude = 10.0_dp**(FLOOR(LOG10(timeInYears)))
            targetTime = ((FLOOR(timeInYears/orderMagnitude) + 1.0_dp) * orderMagnitude)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 0.0) THEN
            targetTime = 10 * timeInYears * SECONDS_PER_YEAR
        ELSE
            targetTime=SECONDS_PER_YEAR*1.0d-7
        ENDIF
    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics
        ! Calculate column densities and update physical parameters for 1D models
        IF (enable_radiative_transfer .AND. points.gt.1) THEN
            call findcoldens_edge2core(coldens_obs(dstep),finalDens,density_scale_radius,density_power_index,parcelRadius(dstep))

            ! Calculate density profiles
            density_max(dstep)=ngas_r(parcelRadius(dstep),finalDens,density_scale_radius,density_power_index)

            IF (density(dstep).ge.density_max(dstep)) density(dstep)=density_max(dstep)

            ! Recompute the coldens from physics-core.f90
            ! coldens should be amount of gas from edge to parcel
            coldens(dstep)=cloudSize/real(points)*density(dstep)
            
            ! Add previous column densities to current as we move into cloud to get total
            IF (dstep .lt. points) THEN 
                coldens(dstep)=coldens(dstep)+coldens(dstep+1)
            END IF

            ! Calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
            av(dstep)= baseAv +coldens(dstep)/1.6d21

            ! Diagnostic output
            ! print '(A,1PE12.3,A,0PF8.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3)', &
            ! 't=',timeInYears, '  r=',parcelRadius(dstep), '  nH(t)=',density(dstep), '  nH(r)_max=',density_max(dstep), &
            ! '  NH(t)(edge->core)=',coldens(dstep), '  NH(r)(edge->core)=', coldens_obs(dstep), '  av(t)=',av(dstep)
        END IF
    END SUBROUTINE updatePhysics

    !This is a dummy subroutine.
    SUBROUTINE sublimation(abund, lpoints)
        INTEGER, INTENT(IN) :: lpoints
        REAL(dp), INTENT(INOUT) :: abund(nspec+1,lpoints)
    END SUBROUTINE sublimation
END MODULE cloud_mod


