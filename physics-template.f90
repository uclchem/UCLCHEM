MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: tstep,dstep
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb
    !More complicated controls have multiple options
    !vap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see chem_initialise)
    integer :: evap,ion
    
    !Number of depth points included in model
    integer, parameter :: points=1  

    !variables either controlled by physics or that user may wish to change    
    double precision :: d0,dens,temp,tage,tout,t0,dfin,tfin,av(points)
    double precision :: size,oldtemp,avic,bc

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

!This is the time step for outputs from UCL_CHEM NOT the timestep for the integrater. DLSODE sorts that out based on chosen error
!tolerances (RTOL/ATOL) and is simply called repeatedly until it outputs a time >= tout. tout in seconds for DLSODE, tage in
!years for output.
    SUBROUTINE timestep
        IF (tstep .gt. 2000) THEN
            tout=(tage+20000.0)/year
        ELSE IF (tstep .gt. 1000) THEN
            tout=(tage+10000.0)/year
        ELSE
            tout=(tage+10.0)/year
        ENDIF
        !This is to match Serena's timesteps for testing code.
        IF (tstep .eq. 0) tout=1d-7/year     
    END SUBROUTINE timestep

!This is called by main so it should either do any physics the user wants or call the subroutines that do.
!The exception is densdot() as density is integrated with chemistry ODEs.    
    SUBROUTINE phys_update
        av(dstep)= avic +((size*(real(dstep)/real(points)))*dens)/1.6d21
    END SUBROUTINE phys_update

!This FUNCTION works out the time derivative of the density, allowing DLSODE to update density with the rest of our ODEs
!It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DLSODE
!Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot()
        double precision :: densdot

        IF (dens .lt. dfin) THEN
            densdot=(dens**4)/d0
            densdot=densdot**0.3333
            densdot=densdot*((8.4d-30*d0*(((dens/d0)**0.333)-1.0))**0.5)
            densdot=densdot*bc
        ELSE
            densdot=0.0        
        ENDIF    
    END FUNCTION densdot
END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
