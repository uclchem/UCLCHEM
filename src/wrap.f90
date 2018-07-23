!WORK IN PROGRESS
!THIS IS A SERIES OF SUBROUTINES THAT CAN BE COMPILED WITH F2PY TO PRODUCE A PYTHON MODULE
!EACH SUBROUTINE BECOMES A PYTHON FUNCTION

SUBROUTINE CLOUD(initDens,finDens,finTime,temp,outFile)
    USE physics
    USE chemistry
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: initDens,finDens,finTime,temp
    CHARACTER(LEN=*), INTENT(IN) :: outFile
    CHARACTER (LEN=100):: abundFile,outputFile,columnFile

    !f2py intent(in) initDens,finDens,finTime,temp,outFile
    include 'defaultparameters.f90'
    close(10)
    close(11)
    close(7)

    open(10,file=outFile,status='unknown') 
    columnFlag=.False.
    open(7,file=outFile//"-start.dat",status='unknown') 
    initialDens=initDens
    IF (ABS(initDens-finDens) .GT. 0.01) THEN
        collapse=1
        finalDens=finDens
    ELSE
        collapse=0
    END IF
    initialTemp=temp
    finalTime=finTime
    switch=0

    CALL initializePhysics
    CALL initializeChemistry

    dstep=1
    currentTime=0.0
    timeInYears=0.0
    !loop until the end condition of the model is reached 
    DO WHILE ((switch .eq. 1 .and. dens(1) < finalDens) .or. (switch .eq. 0 .and. timeInYears < finalTime))

        !store current time as starting point for each depth step
        IF (points .gt. 1) THEN
            currentTimeold=targetTime
            currentTime=currentTimeold
        END IF
        !Each physics module has a subroutine to set the target time from the current time
        CALL updateTargetTime

        !loop over parcels, counting from centre out to edge of cloud
        DO dstep=1,points

            dens=abund(nspec+1,dstep)
            !update physics
            CALL updatePhysics
            !update chemistry
            CALL updateChemistry
            
            !set time to the final time of integrator rather than target     
            targetTime=currentTime
            !reset target for next depth point
            if (points .gt. 1)currentTime=currentTimeold
            !get time in years for output
            timeInYears= currentTime*year
            !write this depth step
            CALL output
        END DO
    END DO 
END SUBROUTINE CLOUD