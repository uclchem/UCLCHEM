! 2015 UCLCHEM by Serena Viti update by Jon Holdship
! Rewritten for FORTRAN 95 and modulised
PROGRAM uclchem
!Everything to do with physics should be stored in a physics module based on physics-template.f90
!UCLCHEM uses density and temperature in chemical calculations so these MUST be provided, everything else is
!user dependent

USE physics
USE chemistry
IMPLICIT NONE

!All variables the user is likely to want to change are in parameters.f90 
include 'parameters.f90'

!This commented out line is a simple read in. This is useful when running large numbers of models, varying only a few paramters
!use it by replacing the variables as needed and running code as "./main < value1 value2 value3"
!read(*,*) variable1,variable2,variable3

!Set up with initial values. For chemistry this is setting initial abundances and assigning memory for ODE solver
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
END PROGRAM uclchem

