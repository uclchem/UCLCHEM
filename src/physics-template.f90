!Simple physics module. Models points along a 1d line from the centre to edge of a cloud. Assuming the cloud is spherical you can average
!over the points to get a 1d average and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE
    integer :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,volcflag,coflag,tempindx
    
    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,dens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime,grainRadius,initialTemp
    double precision :: cloudSize,rout,rin,baseAv,bc,olddens,maxTemp
    double precision :: tempa(5),tempb(5),codestemp(5),volctemp(5),solidtemp(5)
    double precision, allocatable :: av(:),coldens(:),temp(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18


CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    
    SUBROUTINE initializePhysics
    !Any initialisation logic steps go here
    !cloudSize is important as is allocating space for depth arrays
        allocate(av(points),coldens(points))
        cloudSize=(rout-rin)*pc
        if (collapse .eq. 1) THEN
            dens=1.001*initialDens
        ELSE
            dens=initialDens
        ENDIF 
    END SUBROUTINE

    SUBROUTINE updateTargetTime
    !At each timestep, the time at the end of the step is calculated by calling this function
    !You need to set targetTime in seconds, its initial value will be the time at start of current step
    IF (timeInYears .gt. 1.0d6) THEN
        targetTime=(timeInYears+1.0d5)/year
    ELSE IF (timeInYears .gt. 10000) THEN
        targetTime=(timeInYears+1000.0)/year
    ELSE IF (timeInYears .gt. 1000) THEN
        targetTime=(timeInYears+100.0)/year
    ELSE IF (timeInYears .gt. 0.0) THEN
        targetTime=(timeInYears*10)/year
    ELSE
        targetTime=3.16d7*10.d-8
    ENDIF

    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrator.
    ! DLSODE sorts that out based on chosen error tolerances (RTOL/ATOL) and is simply called repeatedly
    !until it outputs a time >= targetTime.
    END SUBROUTINE updateTargetTime

    !This routine is formed for every parcel at every time step.
    !update any physics here. For example, set density
    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from core to edge
        !and coldens should be amount of gas from edge to parcel.
        coldens(dstep)= cloudSize*((real(points-dstep))/real(points))*dens
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21
    END SUBROUTINE updatePhysics

    pure FUNCTION densdot()
    !Required for collapse=1, works out the time derivative of the density, allowing DVODE
    !to update density with the rest of our ODEs
    !It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DVODE
        double precision :: densdot
        !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (dens .lt. finalDens) THEN
             densdot=bc*(dens**4./initialDens)**0.33*&
             &(8.4d-30*initialDens*((dens/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=1.0d-30       
        ENDIF    
    END FUNCTION densdot
END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
