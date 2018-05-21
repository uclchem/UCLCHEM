!Module to read physical paramters from an input file for every time step and solve chemistry for those paramters
!Used for post-processing other model outputs.
!To Do: Add interpolating routines to allow inputs with time steps too sparse for good chemical modelling

!Set phase=2 to use module to post-process. Set phase=1 to run a cloud model for intiial conditions
MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,volcflag,coflag,tempindx,io
   
    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime,grainRadius,initialTemp
    double precision :: cloudSize,rout,rin,oldtemp,baseAv,bc,olddens,maxTemp
    double precision, allocatable :: av(:),coldens(:),temp(:),dens(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18

    !Interpolation variables
    integer, parameter :: nInterpPoints=201
    double precision :: times(nInterpPoints),temperatures(nInterpPoints),densities(nInterpPoints)
    double precision :: tempA(nInterpPoints),tempB(nInterpPoints),tempC(nInterpPoints)
    double precision :: densA(nInterpPoints),densB(nInterpPoints),densC(nInterpPoints)
    double precision :: junk ! put this in read loop for rows you don't need
    integer :: iRead


CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    
    SUBROUTINE initializePhysics
    !Any initialisation logic steps go here
    !cloudSize is important as is allocating space for depth arrays
        allocate(av(points),coldens(points),temp(points),dens(points))
        cloudSize=(rout-rin)*pc

        if (collapse .eq. 1) THEN
            dens=1.001*initialDens
        ELSE
            dens=initialDens
        ENDIF 

        !It's not possible to prepare for every use case but the comments in this if statement
        !give the general procedure for using this module
        IF (phase .eq. 2) THEN
            !Possbile to set up code so it reads from command line . 
            !in this way run "uclchem filename.dat" to run with that file
            !CALL getarg(1, name_user_file)
            !open(64,file=name_user_file,status='old')
            
            !Open a file with time,dens,temp in columns
            open(64,file='output/column.dat',status='old')

            !read in those values, use the junk variable defined above to store unused columns
            DO iRead=1,nInterpPoints
                read(64,*) times(iRead),densities(iRead),temperatures(iRead),junk,junk
            END DO

            !Temperatures in Kelvin, densities in hydrogen nuclei / cm^3. Convert if necessary
            !Time in seconds
            times=times/year

            !This sets up the interpolator for temperature and density
            CALL splineSetup(times,densities,densA,densB,densC,nInterpPoints)
            CALL splineSetup(times,temperatures,tempA,tempB,tempC,nInterpPoints)
            initialDens=densities(1)
            dens=initialDens
            initialTemp=temperatures(1)
            temp=initialTemp
            
            !Set model to finish once time exceeds a value
            !Set that maximum time to slightly less than last element of time array
            switch=0
            finalTime=times(UBOUND(times))-1.0

        END IF
    END SUBROUTINE

    !This works fine for phase 1 where user will set initial conditions of cloud.
    !For phase 2, the actual post-processing, the best time step is problem dependent
    SUBROUTINE updateTargetTime
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
   END SUBROUTINE updateTargetTime

    SUBROUTINE updatePhysics
        !Only do post-processing on phase 2, so phase 1can be cloud model for initial conditions
        IF (phase .eq. 2) THEN
            dens(dstep)=getInterp(targetTime,times,densities,densA,densB,densC,nInterpPoints)
            temp(dstep)=getInterp(targetTime,times,temperatures,tempA,tempB,tempC,nInterpPoints)
        END IF


        IF (dstep .lt. points) THEN
            coldens(dstep)= cloudSize*((real(points+0.5-dstep))/real(points))*dens(dstep)
        ELSE
            coldens(dstep)= 0.5*(cloudSize/real(points))*dens(dstep)
        END IF
        av(dstep)= baseAv +coldens(dstep)/1.6d21

    
    END SUBROUTINE updatePhysics

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


   SUBROUTINE splineSetup (x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alex G: January 2010
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    integer n
    double precision x(n), y(n), b(n), c(n), d(n)
    integer i, j, gap
    double precision h

    gap = n-1
    ! check input
    if ( n < 2 ) return
    if ( n < 3 ) then
      b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
    end if
    !
    ! step 1: preparation
    !
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2.0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do
    !
    ! step 2: end conditions 
    !
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0
    c(n) = 0.0
    if(n /= 3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !
    ! step 3: forward elimination 
    !
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do
    !
    ! step 4: back substitution
    !
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !
    ! step 5: compute spline coefficients
    !
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3.*c(i)
    end do
    c(n) = 3.0*c(n)
    d(n) = d(n-1)
    END SUBROUTINE splineSetup

    function getInterp(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    double precision getInterp
    integer n
    double precision  u, x(n), y(n), b(n), c(n), d(n)
    integer i, j, k
    double precision dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u <= x(1)) then
      getInterp = y(1)
      return
    end if
    if(u >= x(n)) then
      getInterp = y(n)
      return
    end if

    !*
    !  binary search for for i, such that x(i) <= u <= x(i+1)
    !*
    i = 1
    j = n+1
    do while (j > i+1)
      k = (i+j)/2
      if(u < x(k)) then
        j=k
        else
        i=k
       end if
    end do
    !*
    !  evaluate spline interpolation
    !*
    dx = u - x(i)
    getInterp = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    end function getInterp

END MODULE physics

