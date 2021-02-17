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

    heatingFlag=.False.

    dstep=1
    currentTime=0.0
    timeInYears=0.0

    !Set up with initial values. For chemistry this is setting initial abundances and assigning memory for ODE solver
    CALL initializePhysics
    write(*,*) "p"
    gastemp=600.0
    ccol= 18058765988.12
    h2col=9205561518560.0
    coldens=632456000000000.0
    av=coldens*6.289E-22
    cloudSize=2000000000.0
    close(10)
    CALL initializeChemistry
    write(*,*) "initialized"
    !loop until the end condition of the model is reached 
    DO WHILE ((switch .eq. 1 .and. density(1) < finalDens) .or. (switch .eq. 0 .and. timeInYears < finalTime))
        !store current time as starting point for each depth step
        currentTimeold=currentTime

        !Each physics module has a subroutine to set the target time from the current time
        CALL updateTargetTime

        !loop over parcels, counting from centre out to edge of cloud
        DO dstep=1,points
            !update chemistry from currentTime to targetTime
            if (timeInYears .gt. 5d5) heatingFlag=.True.
            !write(*,*) heatingFlag
            CALL updateChemistry
            currentTime=targetTime
            !get time in years for output, currentTime is now equal to targetTime
            timeInYears= currentTime/SECONDS_PER_YEAR
            !Update physics so it's correct for new currentTime and start of next time step
            CALL updatePhysics
            !Sublimation checks if Sublimation should happen this time step and does it
            CALL sublimation(abund)
            !write this depth step now time, chemistry and physics are consistent
            CALL output

            !reset time for next depth point
            if (points .gt. 1)currentTime=currentTimeold
        END DO
    END DO 
END PROGRAM uclchem