! Cloud.f90 is our most basic physics module
! It simulates a static or collapsing cloud of isothermal gas.
! It is useful for static models or for producing initial abundances for the other modules.
MODULE cloud_mod
    USE constants
    USE DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp
    USE physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
    & zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, densdot, gastemp, dusttemp, av,&
    &coldens
    USE network
    use f2py_constants
    IMPLICIT NONE
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


        !Set up freefall.
        IF(freefall) density=1.001*initialDens

    END SUBROUTINE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE updateTargetTime
        IF (timeInYears .gt. 1.0d6) THEN !code in years for readability, targetTime in s
            targetTime=(timeInYears+1.0d5)*SECONDS_PER_YEAR
        ELSE  IF (timeInYears .gt. 1.0d5) THEN
            targetTime=(timeInYears+1.0d4)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1.0d4) THEN
            targetTime=(timeInYears+1000.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1000) THEN
            targetTime=(timeInYears+100.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 0.0) THEN
            targetTime=(timeInYears*10.0)*SECONDS_PER_YEAR
        ELSE
            targetTime=SECONDS_PER_YEAR*1.0d-7
        ENDIF
    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics
        !Nothing to do here :)
    END SUBROUTINE updatePhysics

    !This is a dummy subroutine.
    SUBROUTINE sublimation(abund, lpoints)
        INTEGER, INTENT(IN) :: lpoints
        REAL(dp), INTENT(INOUT) :: abund(nspec+1,lpoints)
    END SUBROUTINE sublimation
END MODULE cloud_mod


