! Cloud.f90 is our most basic physics module
! It simulates a static or collapsing cloud of isothermal gas.
! It is useful for static models or for producing initial abundances for the other modules.
MODULE cloud_mod
    USE physicscore
    USE network
    USE constants
    IMPLICIT NONE
CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called at start of UCLCHEM run
    ! Uses values in defaultparamters.f90 and any inputs to set initial values        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        successFlag=0
        ! Modules not restarted in python wraps so best to reset everything manually.
        IF (ALLOCATED(av)) DEALLOCATE(av,coldens,gasTemp,dustTemp,density)
        ALLOCATE(av(points),coldens(points),gasTemp(points),dustTemp(points),density(points))


        currentTime=0.0
        targetTime=0.0
        timeInYears=0.0

        !Set up basic physics variables
        cloudSize=(rout-rin)*pc
        gasTemp=initialTemp
        dustTemp=gasTemp

        !Set up freefall.
        SELECT CASE(freefall)
            !freefall Rawlings 1992
            CASE(0)
                density=initialDens
            CASE(1)
                density=1.001*initialDens
            CASE DEFAULT
                write(*,*) "freefall must be 0 or 1 for clouds"
        END SELECT

        !calculate initial column density as distance from core edge to current point * density
        DO dstep=1,points
            coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
        END DO
          !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av= baseAv +coldens/1.6d21
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
        !calculate column density. Remember dstep counts from core centre to edge
        !and coldens should be amount of gas from edge to parcel.
        IF (dstep .lt. points) THEN
            !column density of current point + column density of all points further out
            coldens(dstep)=(cloudSize/real(points))*density(dstep)
            coldens(dstep)=coldens(dstep)+sum(coldens(dstep:points))
        ELSE
            coldens(dstep)=cloudSize/real(points)*density(dstep)
        END IF

        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21
        dustTemp=gasTemp
    END SUBROUTINE updatePhysics

    ! This subroutine must be in every physics module so we dummy it here.
    SUBROUTINE sublimation(abund)
        DOUBLE PRECISION :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
    END SUBROUTINE sublimation    
END MODULE cloud_mod


