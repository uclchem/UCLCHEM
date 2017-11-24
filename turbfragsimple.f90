!Physics module based on Holdship & Viti 2015. Models points along a 1d line from the centre
!to edge of a core of gas. models core as it is created by isothermal shock.
!Assuming the cloud is spherical you can average over the points to get a 1d average 
! and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE

    !Use main loop counters in calculations so they're kept here
    integer :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,volcflag,coflag,tempindx

    !variables either controlled by physics or that user may wish to change
    double precision :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime,grainRadius,initialTemp
    double precision :: cloudSize,rout,rin,baseAv,bc,olddens,maxTemp
    double precision :: tempa(6),tempb(6),codestemp(6),volctemp(6),solidtemp(6)
    double precision, allocatable :: av(:),coldens(:),temp(:),dens(:)
    
    !variables either controlled by physics or that user may wish to change    
    !BD model specific variables
    double precision ::  dshock,rshock, timeInYearsold,ffc,mbd,mach
    double precision, allocatable ::  tshock(:)


    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18,cs=2.00d4

    character(1) :: modeln
    character(4) :: bdstring
    
CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

!This is the time step for outputs from UCL_CHEM NOT the timestep for the integrater. DLSODE sorts that out based on chosen error
!tolerances (RTOL/ATOL) and is simply called repeatedly until it outputs a time >= targetTime. targetTime in seconds for DLSODE, timeInYears in
!years for output.
    
    SUBROUTINE initializePhysics
        allocate(av(points),coldens(points),tshock(points),temp(points),dens(points))
        cloudSize=(rout-rin)*pc
        dens=1.01*initialDens
        temp=initialTemp
        mbd=0.2
        if (phase .eq. 2) THEN
            collapse=1
            write(*,*) "Collapse flag ignored for phase 2"
            call bdboundaries
        END IF

    END SUBROUTINE

    SUBROUTINE updateTargetTime
        IF (phase .eq. 1) THEN
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
        ELSE
            IF (timeInYears .gt. 1.0d6) THEN
                targetTime=(timeInYears+20000.0)/year
            ELSE IF (timeInYears/year .gt. tshock(points)) THEN
                targetTime=(timeInYears+2500.0)/year
            ELSE
                targetTime=(targetTime+(tshock(points)/50))
            END IF
        END IF
    END SUBROUTINE updateTargetTime

!This is called by main so it should either do any physics the user wants or call the subroutines that do.
!The exception is densdot() as density is integrated with chemistry ODEs.    
    SUBROUTINE updatePhysics
        IF (phase .eq. 2) THEN
            CALL simpleshock
        END IF
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        !for bd, need to have dstep=1 as core.
        IF (dstep .lt. points) THEN
            coldens(dstep)= cloudSize*((real(points+0.5-dstep))/real(points))*dens(dstep)
        ELSE
            coldens(dstep)= 0.5*(cloudSize/real(points))*dens(dstep)
        END IF
        av(dstep)= baseAv +coldens(dstep)/1.6d21
    END SUBROUTINE updatePhysics


    pure FUNCTION densdot()
    !This FUNCTION works out the time derivative of the density, allowing DLSODE to update
    !density with the rest of our ODEs. It get's called by F, the SUBROUTINE in chem.f90
    !that sets up the ODEs for DLSODE.
    !Currently set to Rawlings 1992 freefall.
        double precision :: densdot
        IF (phase .eq. 1) THEN
            IF (dens(dstep) .lt. finalDens) THEN
                densdot=bc*(dens(dstep)**4./initialDens)**0.33*&
                &(8.4d-30*initialDens*((dens(dstep)/initialDens)**0.33-1.))**0.5
            ELSE
               densdot=0.00    
            ENDIF 
        ELSE
            IF (dens(dstep)*1.001 .gt. dshock .and.  dens(dstep) .lt. finalDens) THEN
                densdot=bc*(dens(dstep)**4./dshock)**0.33*&
                &(8.4d-30*dshock*((dens(dstep)/dshock)**0.33-1.))**0.5
            ELSE
                densdot=0.00    
            END IF
        ENDIF    
    END FUNCTION densdot

    subroutine bdboundaries
    !sets the post shock density, shock veloticy and shock times for each point
    !based on final density required for core to be critical at chosen mass.

        ! dshock from padoan nordlund 2004
        dshock=1.0d3*((mbd/3.3)**(-2.0))*((initialTemp/10)**(3.0))
        
        ! total cloud rshock in cm
        rshock=((3.0/(4.0*pi))*((2.0d33*mbd)/(mh*dshock)))
        rshock=rshock**(1.0/3.0)
        cloudSize=rshock

        ! mach calculated by assuming M ~p/p0
        mach=dsqrt(dshock/initialDens) 
        !time for each point then based on radial position and mach number.
        !tshock in seconds
        DO dstep=1,points
            tshock(dstep)=((dstep-1)*rshock)/((points-1.0)*mach*cs)   
        END DO
        !below added for a modified freefall equation in poly, big number is G in cm3/g.year2
        ffc=24.0*pi*mh*initialDens*(6.64d+7)
    END SUBROUTINE bdboundaries


    SUBROUTINE simpleshock    
    !'shock' for turbulent fragmentation, instantenous density increase as shock
    ! front reaches each depth point
        double precision store
        !instantly increases density to final post-shock density at tshock
        IF (targetTime .lt. tshock(dstep)) THEN
            dens(dstep)=initialDens
        ENDIF

        IF (targetTime .gt. tshock(dstep) .and. dens(dstep) .lt. dshock) THEN
            dens(dstep)=dshock*1.01
        END IF

    END SUBROUTINE simpleshock

END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
