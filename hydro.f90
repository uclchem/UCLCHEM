!Module to read physical paramters from an input file for every time step and solve chemistry for those paramters
!Used for post-processing other model outputs.
!To Do: Add interpolating routines to allow inputs with time steps too sparse for good chemical modelling

!Set phase=2 to use module to post-process. Set phase=1 to run a cloud model for intiial conditions
MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: tstep,dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see chem_initialise)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,volcflag,coflag,tempindx,io
   
    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,tage,tout,t0,t0old,finalDens,finalTime,grainRadius,initialTemp
    double precision :: size,rout,rin,oldtemp,baseAv,bc,olddens,maxTemp
    double precision :: tempa(6),tempb(6),codestemp(6),volctemp(6),solidtemp(6)
    double precision, allocatable :: av(:),coldens(:),temp(:),dens(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18


CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    
    SUBROUTINE phys_initialise
    !Any initialisation logic steps go here
    !size is important as is allocating space for depth arrays
        allocate(av(points),coldens(points),temp(points),dens(points))
        size=(rout-rin)*pc

        if (collapse .eq. 1) THEN
            dens=1.001*initialDens
        ELSE
            dens=initialDens
        ENDIF 

        !Open a file with time,dens,temp,Av in columns
        !Multiple rows with same time value if multipoint model
        open(64,file='physics_inputs.dat',status='old')

    END SUBROUTINE

    SUBROUTINE timestep
        IF (phase .eq. 1) THEN
            IF (tstep .gt. 2000) THEN
                tout=(tage+20000.0)/year
            ELSE IF (tstep .gt. 1000) THEN
                tout=(tage+10000.0)/year
            ELSE IF (tstep .gt. 1) THEN
                tout=1.58e11*(tstep-0)
            ELSE  IF (tstep .eq. 1) THEN  
                tout=3.16d7*1.0d3
            ELSE
                tout=3.16d7*10.d-8
            ENDIF
        END IF 
   END SUBROUTINE timestep

    SUBROUTINE phys_update

        !Only do post-processing on phase 2, so phase 1can be cloud model for initial conditions
        IF (phase .eq. 2) THEN
            !read in physical inputs from file
            !Tage in years, dens in hydrogen nuclei per cubic cm, temp in kelvin and Av in magnitudes
            !just multiply density by 2.0 if you have density in H2 molecules /cm^3
            read(64,*)tage,dens(dstep),temp(dstep),av(dstep)
            tout=(tage)/year

            coldens(dstep)=1.6d21*av(dstep)
        ELSE
            IF (dstep .lt. points) THEN
                coldens(dstep)= size*((real(points+0.5-dstep))/real(points))*dens(dstep)
            ELSE
                coldens(dstep)= 0.5*(size/real(points))*dens(dstep)
            END IF
            av(dstep)= baseAv +coldens(dstep)/1.6d21
        END IF

    
    END SUBROUTINE phys_update

    pure FUNCTION densdot()
    !Required for collapse=1, works out the time derivative of the density, allowing DLSODE
    !to update density with the rest of our ODEs
    !It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DLSODE
        double precision :: densdot
        !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (dens(dstep) .lt. finalDens .and. phase .eq. 1) THEN
             densdot=bc*(dens(dstep)**4./initialDens)**0.33*&
             &(8.4d-30*initialDens*((dens(dstep)/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=0.00    
        ENDIF    
    END FUNCTION densdot
END MODULE physics

