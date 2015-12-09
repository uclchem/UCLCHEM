!Izaskun Cshock in progress
MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: tstep,dstep
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see chem_initialise)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,monoflag,volcflag,coflag
    
    !Number of depth points included in model
    integer, parameter :: points=1  

    !variables either controlled by physics or that user may wish to change    
    double precision :: d0,dens,temp,tage,tout,t0,dfin,tfin,av(points)
    double precision :: size,oldtemp,avic,bc,tempa,tempb,tstart

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

!This is the time step for outputs from UCL_CHEM NOT the timestep for the integrater. DLSODE sorts that out based on chosen error
!tolerances (RTOL/ATOL) and is simply called repeatedly until it outputs a time >= tout. tout in seconds for DLSODE, tage in
!years for output.
    SUBROUTINE timestep
        IF (phase .eq. 1) THEN
            IF (tstep .gt. 2000) THEN
                tout=(tage+20000.0)/year
            ELSE IF (tstep .gt. 1000) THEN
                tout=(tage+10000.0)/year
            ELSE IF (tstep .gt. 1) THEN
                tout=1.58e11*(tstep-0)
            ELSE
                tout=3.16d7*10.0**(tstep+2)
            ENDIF
        ELSE
            tout=(tage+10000.0)/year
        END IF
        !This is to match Serena's timesteps for testing code.
    END SUBROUTINE timestep

!This is called by main so it should either do any physics the user wants or call the subroutines that do.
!The exception is densdot() as density is integrated with chemistry ODEs.    
    SUBROUTINE phys_update
        !calculate the Av using an assumed extinction outside of core (avic), depth of point and density
        av(dstep)= avic +((size*(real(dstep)/real(points)))*dens)/1.6d21

        IF (phase .eq. 2) THEN
            !temperature increase borrowed from sv for comparison 288.000
            !will add general profile later
            temp=10. + (7.8470d-3*tage**0.8395)
            write(*,*) temp
        END IF

    END SUBROUTINE phys_update

!This FUNCTION works out the time derivative of the density, allowing DLSODE to update density with the rest of our ODEs
!It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DLSODE
!Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot()
        double precision :: densdot

        IF (dens .lt. dfin) THEN
 !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (dens .lt. dfin) THEN
             densdot=bc*(dens**4./d0)**0.33*&
             &(8.4d-30*d0*((dens/d0)**0.33-1.))**0.5
        ELSE
            densdot=1.0d-30       
        ENDIF    
    END FUNCTION densdot
END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
