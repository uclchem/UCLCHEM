!Simple physics module. Models points along a 1d line from the centre to edge of a cloud. Assuming the cloud is spherical you can average
!over the points to get a 1d average and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: tstep,dstep,points,NF,ifile
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see chem_initialise)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,volcflag,coflag,tempindx
   
    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,dens,tage,tout,t0,t0old,finalDens,finalTime,radg,initialTemp
    double precision :: size,rout,rin,oldtemp,baseAv,bc,olddens,maxTemp
    double precision :: tempa(5),tempb(5),codestemp(5),volctemp(5),solidtemp(5)
    double precision, allocatable :: av(:),coldens(:),temp(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18

    !Hydro specific stuff
    double precision :: tinc,tincmax,junkdist,junkvel
    !Interpolation variables
    double precision :: temp2(17784),dens2(17784),density(17784),time(17784),temps(17784)


CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    
    SUBROUTINE phys_initialise
    !Any initialisation logic steps go here
    !size is important as is allocating space for depth arrays
        allocate(av(points),coldens(points),temp(points))
        size=(rout-rin)*pc

        if (collapse .eq. 1) THEN
            dens=1.001*initialDens
        ELSE
            dens=initialDens
        ENDIF 

        NF=0
        open(21,file='fp_traj.000263',status='old')
        DO ifile=1,17784
            read(21,*)junkdist,junkdist,junkdist,density(ifile),temps(ifile),junkdist,junkvel,junkvel,junkvel,time(ifile)
            NF=NF+1
            time(ifile)=time(ifile)/year
!            temps(ifile) = temps(ifile)
        END  DO

        !convert density from gcm^-3 to H nuclei cm^-3
        density = 1.1*0.7057*6.02e23 * density

!	write(*,*) time(1), density(1), temps(1)
        t0=time(1)
        finalTime=time(NF)*year
        CALL splinehydro(time, density, NF, 1d30, 1d30, dens2)
        CALL splinehydro(time, temps, NF, 1d30, 1d30, temp2)
        tinc=0.1
        tincmax=(0.1/year)/tinc
        write(*,*) tincmax*year
    END SUBROUTINE

    SUBROUTINE timestep
    !At each timestep, the time at the end of the step is calculated by calling this function
    !You need to set tout in seconds, its initial value will be the time at start of current step
    IF (tout .lt. tincmax) THEN
       !tout=(1.0+tinc)*t0
       tout=t0+(tincmax*tinc)*10
    ELSE
        tout=t0+(tincmax*tinc)*10
    END IF
    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrator.
    ! DLSODE sorts that out based on chosen error tolerances (RTOL/ATOL) and is simply called repeatedly
    !until it outputs a time >= tout.
    END SUBROUTINE timestep

    SUBROUTINE phys_update
        !splint gives temp/dens at time t0 given outputs from spline
        CALL splinthydro(time,temps,temp2,nf,t0,temp(dstep))
        CALL splinthydro(time,density,dens2,nf,t0,dens)

!	write(*,*) 'A', t0, temp, dens
        
        !calculate column density by dens x depth of point 
        coldens(dstep)= size*((real(dstep))/real(points))*dens
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv !+coldens(dstep)/1.6d21
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

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.

   SUBROUTINE splinehydro(x, y, n, yp1, ypn, y2)
!   use nrtype
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for a
! natural spline with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipiated value of n
! (adopted from Numerical Recipes in FORTRAN 77)
!
   INTEGER, PARAMETER :: DP=KIND(1.0D0)
   INTEGER:: n
   INTEGER, PARAMETER:: nmax=20000
   REAL(DP):: yp1, ypn, x(n), y(n), y2(n)
   INTEGER:: i, k
   REAL(DP):: p, qn, sig, un, u(nmax)

     if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
     else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     endif

     !this loop is broken
     do i=2, n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
             & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     enddo

     if (ypn.gt..99e30) then
        qn=0.
        un=0.
     else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
     endif

     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

     do k=n-1, 1, -1
        y2(k)=y2(k)*y2(k+1)+u(k)
     enddo

     return
     END SUBROUTINE splinehydro

   SUBROUTINE splinthydro(xa, ya, y2a, n, x, y)
!   USE nrtype
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! (with the xa(i) in order), and given the array y2a(1:n), which is the output
! from the subroutine spline, and given a value of x, this routine returns a
! cubic spline interpolated value y.
! (adopted from Numerical Recipes in FORTRAN 77)
!
   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   INTEGER:: n
   REAL(DP) :: x, y, xa(n), y2a(n), ya(n)
   INTEGER:: k, khi, klo
   REAL(DP):: a, b, h

     klo=1
     khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     if (h.eq.0.) pause 'bad xa input in splint'

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

     return
     END SUBROUTINE splinthydro


END MODULE physics

