!Physics module based on Holdship & Viti 2015. Models points along a 1d line from the centre
!to edge of a core of gas. models core as it is created by isothermal shock.
!Assuming the cloud is spherical you can average over the points to get a 1d average 
! and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: tstep,dstep,ishock,points
    
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see chem_initialise)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,monoflag,volcflag,coflag

    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,dens,temp,tage,tout,t0,t0old,finalDens,finalTime,grainRadius
    double precision :: size,rout,rin,oldtemp,baseAv,bc,tempa,tempb,olddens,oldt0,maxt

    !old bd model variables
    double precision ::  dshock,rshock, tageold,ffc,mbd,mach
    double precision, allocatable :: tshock(:),av(:),coldens(:)

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18,cs=2.00d4

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

!This is the time step for outputs from UCL_CHEM NOT the timestep for the integrater. DLSODE sorts that out based on chosen error
!tolerances (RTOL/ATOL) and is simply called repeatedly until it outputs a time >= tout. tout in seconds for DLSODE, tage in
!years for output.
    
    SUBROUTINE phys_initialise
        allocate(av(points),coldens(points),tshock(points))
        size=(rout-rin)*pc
        dens=initialDens
        mbd=0.40
        call bdboundaries
    END SUBROUTINE

    SUBROUTINE timestep
            IF (tstep .gt. 2000) THEN
                tout=(tage+20000.0)/year
            ELSE IF (tstep .gt. 1000) THEN
                tout=(tage+10000.0)/year
            ELSE IF (tstep .gt. 1) THEN
                tout=1.58e11*(tstep-0)
            ELSE
                tout=3.16d7*10.0**(tstep+2)
            ENDIF

        !IF (tstep .eq. 0) tout=3.16d7*10.d-8
        !This is to match Serena's timesteps for testing code.
    END SUBROUTINE timestep

!This is called by main so it should either do any physics the user wants or call the subroutines that do.
!The exception is densdot() as density is integrated with chemistry ODEs.    
    SUBROUTINE phys_update
        IF (dens .lt. finalDens .and. phase .eq. 2) THEN
            CALL simpleshock
        END IF
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        !for bd, need to have dstep=1 as core.
        IF (dstep .lt. points) THEN
            coldens(dstep)= size*((real(points+0.5-dstep))/real(points))*dens
        ELSE
            coldens(dstep)= 0.5*(size/real(points))*dens
        END IF
        av(dstep)= baseAv +coldens(dstep)/1.6d21
    END SUBROUTINE phys_update


    pure FUNCTION densdot()
    !This FUNCTION works out the time derivative of the density, allowing DLSODE to update
    !density with the rest of our ODEs. It get's called by F, the SUBROUTINE in chem.f90
    !that sets up the ODEs for DLSODE.
    !Currently set to Rawlings 1992 freefall.
        double precision :: densdot
        IF (dens .lt. finalDens .and. phase .eq. 1) THEN
             densdot=bc*(dens**4./initialDens)**0.33*&
            &(8.4d-30*initialDens*((dens/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=0.00    
        ENDIF    
    END FUNCTION densdot

    subroutine bdboundaries
    !sets the post shock density, shock veloticy and shock times for each point
    !based on final density required for core to be critical at chosen mass.

        ! dshock from padoan nordlund 2004
        dshock=1.0d3*((mbd/3.3)**(-2.0))*((temp/10)**(3.0))
        
        ! total cloud rshock in cm
        rshock=((3.0/(4.0*pi))*((2.0d33*mbd)/(mh*dshock)))
        rshock=rshock**(1.0/3.0)

        ! mach calculated by assuming M ~p/p0
        mach=dsqrt(dshock/initialDens) 
        !time for each point then based on radial position and mach number.
        !tshock in seconds STILL NEEDS N/Ntot
        DO ishock=1,points
            tshock(ishock)=((ishock-1)*rshock)/((points-1.0)*mach*cs)   
        END DO
        !below added for a modified freefall equation in poly, big number is G in cm3/g.year2
        ffc=24.0*pi*mh*initialDens*(6.64d+7)
    END SUBROUTINE bdboundaries


    SUBROUTINE simpleshock    
    !'shock' for turbulent fragmentation, instantenous density increase as shock
    ! front reaches each depth point
        double precision store
        !instantly increases density to final post-shock density at tshock
        IF (tout .lt. tshock(dstep)) THEN
            dens=initialDens
        ENDIF

        IF (tout .gt. tshock(dstep) .AND. dens .lt. dshock) THEN
            dens=dshock+1.0
            tageold=tout
        ENDIF

        !Stops shock collapse at post-shock density and reverts to freefall
        !store is used to store bits of the calculation to break up the code a little, it's convenient not necessary.
        IF (dens .gt. dshock) THEN
            store=((dens/dshock)**(1.0/3.0))-1.0
            store=ffc*store
            store=dsqrt(store)
            store=store*(((dens**4.0)/dshock)**(1.0/3.0))
            !tout and tageold are in seconds, the units ffc are years so 3.171d-8 converts
            dens=dens+((tout-tageold)*(3.171d-8)*0.8*store)
            write(10,*) dens
        ENDIF

        tageold=tout
    END SUBROUTINE simpleshock

END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
