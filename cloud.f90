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
    integer :: evap,ion,solidflag,volcflag,coflag,tempindx
    
    !variables either controlled by physics or that user may wish to change    
    double precision :: initdens,dens,tage,tout,t0,t0old,dfin,tfin,radg,inittemp
    double precision :: size,rout,rin,avic,bc,olddens,maxtemp
    double precision :: tempa(5),tempb(5),codestemp(5),volctemp(5),solidtemp(5)
    double precision, allocatable :: av(:),coldens(:),temp(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN.
!YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    SUBROUTINE phys_initialise
        allocate(av(points),coldens(points),temp(points))
        size=(rout-rin)*pc
        if (collapse .eq. 1) THEN
            dens=1.001*initdens
        ELSE
            dens=initdens
        ENDIF 
        temp=inittemp
    END SUBROUTINE

    SUBROUTINE timestep
    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrator.
    !tout in seconds for DLSODE, tage in years for output.
    
        IF (phase .eq. 1) THEN
            IF (tstep .gt. 2000) THEN
                tout=(tage+20000.0)/year
            ELSE IF (tstep .gt. 1000) THEN
                tout=(tage+10000.0)/year
            ELSE IF (tstep .gt. 1) THEN
                tout=1.58e11*(tstep-0)
            ELSE  IF (tstep .eq. 1) THEN  
                tout=3.16d7*10.0**3
            ELSE
                tout=3.16d7*10.d-8
            ENDIF
        ELSE
            tout=(tage+10000.0)/year
        END IF
    END SUBROUTINE timestep

  
    SUBROUTINE phys_update
        !calculate column density. Remember dstep counts from core to edge
        coldens(dstep)= size*((real(dstep))/real(points))*dens
        !calculate the Av using an assumed extinction outside of core (avic), depth of point and density
        av(dstep)= avic +coldens(dstep)/1.6d21

        IF (phase .eq. 2 .and. temp(dstep) .lt. maxtemp) THEN

        !Below we include a temperature profile for hot cores
        !This is a profile taken from Viti et al. 2004 with an additional distance dependence from Nomura and Millar 2004.
        !It takes the form T=A(t^B)*[(d/R)^-0.5], where A and B are given below for various stellar masses

        ! 28,000    60 Msol
                !temp(depth)=10. + (4.74d-7*tage**1.98)
        ! 70,000    25 Msol
                !td(depth)=10. + (1.706d-4*tage**1.289)
        ! 115,000   15 Msol     
                !td(depth)=10. + (9.6966d-4*tage**1.085)
        ! 288,000   10 Msol
                !td(depth)=10. + (7.8470d-3*tage**0.8395)
        !1.15e6     5 Msol
                !td(depth)=10. + (4.8560d-2*tage**0.6255)

            !temperature increase borrowed from sv for comparison 288.000
            !will add general profile later, this works well for inittemp=10 K
            temp(dstep)=(size/rout)*(real(dstep)/real(points))
            temp(dstep)=temp(dstep)**(-0.5)
            temp(dstep)=inittemp + (tempa(tempindx)*tage**tempb(tempindx))*temp(dstep)
            if (temp(dstep) .gt. solidtemp(tempindx) .and. solidflag .ne. 2) solidflag=1
            if (temp(dstep) .gt. volctemp(tempindx) .and. volcflag .ne. 2) volcflag=1
            if (temp(dstep) .gt. codestemp(tempindx) .and. coflag .ne. 2) coflag=1
        END IF

    END SUBROUTINE phys_update

!This FUNCTION works out the time derivative of the density, allowing DLSODE to update density with the rest of our ODEs
!It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DLSODE
!Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot()
        double precision :: densdot
        !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (dens .lt. dfin) THEN
             densdot=bc*(dens**4./initdens)**0.33*&
             &(8.4d-30*initdens*((dens/initdens)**0.33-1.))**0.5
        ELSE
            densdot=1.0d-30       
        ENDIF    
    END FUNCTION densdot
END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
