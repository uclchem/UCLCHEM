!This module contains all the must have physics variables. It is imported by both
!the other physics modules (so they can modify the physics) and the chemistry
!module (so it can use the physics for reaction rates).

MODULE physicscore
    USE constants
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    INTEGER :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    INTEGER :: freefall,switch

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    INTEGER :: ion,instantSublimation

    !variables either controlled by physics or that user may wish to change
    DOUBLE PRECISION :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime,initialTemp
    DOUBLE PRECISION ::  bc,cloudSize,rout,rin,baseAv
    DOUBLE PRECISION, allocatable :: av(:),coldens(:),gasTemp(:),dustTemp(:),density(:)


CONTAINS

    pure FUNCTION densdot(density)
    ! Returns the time derivative of the density.                                     
    ! Analytical function taken from Rawlings et al. 1992                             
    ! Called from chemistry.f90, density integrated with abundances so this gives ydot
    DOUBLE PRECISION, INTENT(IN) :: density
    DOUBLE PRECISION :: densdot
    !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
    IF ((density .lt. finalDens) .and. (freefall .eq. 1)) THEN
        densdot=bc*(density**4./initialDens)**0.33*&
        &(8.4d-30*initialDens*((density/initialDens)**0.33-1.))**0.5
    ELSE
        densdot=0.0
    ENDIF
    END FUNCTION densdot


    pure FUNCTION dByDnDensdot(density)
    !Defunct function which provides the necessary derivative d(dn/dt)/dn
    !in the case one uses a Jacobian.
    DOUBLE PRECISION, INTENT(IN) :: density
    DOUBLE PRECISION :: dByDnDensdot
    !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
    IF (density .lt. finalDens) THEN
        dByDndensdot=bc*8.4d-30*(density**3)*((9.0d0*((density/initialDens)**0.33))-8.0d0)
        dByDnDensdot=dByDnDensdot/(6.0d0*(((density**4.0)/initialDens)**0.66))
        dByDnDensdot=dByDnDensdot/dsqrt(initialDens*8.4d-30*(((density/initialDens)**0.33))-1.0d0)
    ELSE
        dByDnDensdot=0.0
    ENDIF
    END FUNCTION dByDnDensdot

END MODULE physicscore
