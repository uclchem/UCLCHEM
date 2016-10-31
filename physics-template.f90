!Simple physics module. Models points along a 1d line from the centre to edge of a cloud. Assuming the cloud is spherical you can average
!over the points to get a 1d average and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: tstep,dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see chem_initialise)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,volcflag,coflag
   
    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,dens,temp,tage,tout,t0,t0old,finalDens,finalTime,av(points),coldens(points)
    double precision :: size,rout,rin,oldtemp,avic,bc,tempa,tempb,olddens,maxTemp,radg

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    
    SUBROUTINE phys_initialise
    !Any initialisation logic steps go here
    !size is important as is allocating space for depth arrays
        allocate(av(points),coldens(points))
        size=(rout-rin)*pc
        if (collapse .eq. 1) THEN
            dens=1.001*initialDens
        ELSE
            dens=initialDens
        ENDIF 
    END SUBROUTINE

    SUBROUTINE timestep
    !At each timestep, the time at the end of the step is calculated by calling this function
    !You need to set tout in seconds, its initial value will be the time at start of current step


    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrator.
    ! DLSODE sorts that out based on chosen error tolerances (RTOL/ATOL) and is simply called repeatedly
    !until it outputs a time >= tout.
    END SUBROUTINE timestep

    SUBROUTINE phys_update
        !Called once for each time and depth step.
        !increase temperature, recalulcate column density (colddens) and AV
        !Can change density if not using a differential equation. ie collapse=0 and set dens another way

        !calculate column density by dens x depth of point 
        coldens(dstep)= size*((real(dstep))/real(points))*dens
        !calculate the Av using an assumed extinction outside of core (avic), depth of point and density

        av(dstep)= avic +coldens(dstep)/1.6d21

        IF (phase .eq. 2) THEN
            !increase temp here.
            temp=10. + (7.8470d-3*tage**0.8395)
            write(*,*) temp
        END IF
    END SUBROUTINE phys_update

    pure FUNCTION densdot()
    !Required for collapse=1, works out the time derivative of the density, allowing DLSODE
    !to update density with the rest of our ODEs
    !It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DLSODE
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
