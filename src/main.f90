! 2015 UCLCHEM by Serena Viti update by Jon Holdship
! Rewritten for FORTRAN 95 and modulised
PROGRAM uclchem
!Everything to do with physics should be stored in a physics module based on physics-template.f90
!UCLCHEM uses density and temperature in chemical calculations so these MUST be provided, everything else is
!user dependent
USE physics
USE chemistry
IMPLICIT NONE
    INTEGER :: ios=0,line=0,fileNum=12,pos,readIndx
    CHARACTER (LEN=100):: paramFile, buffer,label
    CHARACTER (LEN=100):: abundFile,outputFile,columnFile
    LOGICAL :: columnFileRead=.False.
    
    !All model parameters are given a default value in paramters.f90
    !Full explanations of those parameters in the comments of that file
    INCLUDE 'defaultparameters.f90'
    !Any subset of the parameters can be passed in a file on program start
    !see example.inp
    INCLUDE 'readparameters.f90'


    dstep=1
    currentTime=0.0
    timeInYears=0.0

    !Set up with initial values. For chemistry this is setting initial abundances and assigning memory for ODE solver
    CALL initializePhysics
    CALL initializeChemistry
    
    !update physics
    DO dstep=1,points
        CALL updatePhysics
    END DO

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
            !update chemistry
            CALL updateChemistry
            
            !set time to the final time of integrator rather than target     
            targetTime=currentTime
            !reset target for next depth point
            if (points .gt. 1)currentTime=currentTimeold
            !get time in years for output
            timeInYears= currentTime*year
            !write this depth step
            CALL updatePhysics
            CALL output
        END DO
    END DO 
END PROGRAM uclchem