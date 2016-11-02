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
    double precision :: initialDens,dens,tage,tout,t0,t0old,finalDens,finalTime,radg,initialTemp
    double precision :: size,rout,rin,baseAv,bc,olddens,maxTemp
    double precision :: tempa(5),tempb(5),codestemp(5),volctemp(5),solidtemp(5)
    double precision, allocatable :: av(:),coldens(:),temp(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18

    !variables for collapse modes
    double precision :: unitrho,unitr,unitt,c_s
    double precision :: dimrho,dimr,dimt,maxdimt
    double precision :: rho0,r0
    double precision,parameter :: G_N = 6.674d-8

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN.
!YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    SUBROUTINE phys_initialise
        allocate(av(points),coldens(points),temp(points))
        size=(rout-rin)*pc
        dens=initialDens
        temp=initialTemp

        !collapse stuff
        !calculate sound speed from temperature and mean molecular mass
        !cs=(1d7*kbolt*initialTemp)/(2*mh) assuming molecular mass from H2
        c_s = sqrt(5d6*kbolt*initialTemp/mh)
        
        !Set up collapse modes.
        !work in dimensionless units so find conversions for density, time and radius.
        !Also find maximum time value at which model must stop.

        
        SELECT CASE(collapse)
            !freefall Rawlings 1992
            CASE(0)
                dens=initialDens
            CASE(1)
                dens=1.001*initialDens
            !foster & chevalier 1993
            CASE(2)
                unitrho = initialDens*mh
                unitt = 1./sqrt(4*pi*G_N*unitrho)
                unitr = c_s*unitt
                maxdimt = 5.75 - (15.1/(14 + log10(finalDens/initialDens)))**(1/0.04)
            !ogino, tomisaka & nakamura 1999
            CASE(3)
                unitrho = initialDens*mh
                unitt = 1./sqrt(G_N*unitrho)
                unitr = c_s*unitt
                maxdimt = 0.33 - (2.5/(2.5 + log10(finalDens/initialDens)))**(1/0.18)
            CASE(4)
            
                unitrho = initialDens*mh
                unitt = 1./sqrt(2*pi*G_N*unitrho)
                unitr = c_s*unitt
                maxdimt = 5.5 - (2.1/(1.35 + log10(finalDens/initialDens)))**(1/0.28)
        END SELECT

        IF (collapse .gt. 1) THEN
             !Enforce maximum time value for collapse modes
             finalTime=maxdimt*year*unitt
        END  IF
       
        IF (switch .eq. 1 .and. collapse .gt. 1) THEN
            write(*,*) "Switch must be 0 for BE collapse, changing to stop at finalTime"
            write(*,*) "Tfin = ",finalTime/year, " years"
            switch=0
        END IF


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
        !and coldens should be amount of gas from edge to parcel.
        coldens(dstep)= size*((real(points-dstep))/real(points))*dens
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21

        IF (phase .eq. 2 .and. temp(dstep) .lt. maxTemp) THEN

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
            !will add general profile later, this works well for initialTemp=10 K
            temp(dstep)=(size/(rout*pc))*(real(dstep)/real(points))
            temp(dstep)=temp(dstep)**(-0.5)
            temp(dstep)=initialTemp + ((tempa(tempindx)*(t0/year)**tempb(tempindx))*temp(dstep))
            
            if (temp(dstep) .gt. solidtemp(tempindx) .and. solidflag .ne. 2) solidflag=1
            if (temp(dstep) .gt. volctemp(tempindx) .and. volcflag .ne. 2) volcflag=1
            if (temp(dstep) .gt. codestemp(tempindx) .and. coflag .ne. 2) coflag=1
        END IF

        !Density update for BE collapse modes
        !first calculate current radius and time in dimensionless units
        dimr = (real(dstep)/real(points))*size/unitr
        dimt = t0/unitt
        if (dimt .gt. maxdimt) dimt = maxdimt
        
        !Then calculate central density and radius of central region
        !use this to find density at larger radii. Three different cases possible.

        !foster & chevalier 1993
        if (collapse .eq. 2) then
            rho0 = 10**(15*(5.75-dimt)**(-0.04) - 14)
            r0 = 10**(-6.7*(5.75-dimt)**(-0.04) + 6.6)
            dimrho = rho0/(1 + (dimr/r0)**2.5)
            if (dimt .lt. 5.75) dens = dimrho*initialDens
        
        !ogino, tomisaka & nakamura 1999
        else if (collapse .eq. 3) then
            rho0 = 10**(2.5*(0.33-dimt)**(-0.18) - 2.5)
            r0 = 10**(-1.2*(0.33-dimt)**(-0.18) + 1.3)
            dimrho = rho0/(1 + (dimr/r0)**2.5)
            if (dimt .lt. 0.33) dens = dimrho*initialDens
        
        !nakamura, hanawa & takano 1995
        else if (collapse .eq. 4) then
            rho0 = 10**(2.1*(5.5-dimt)**(-0.28) - 1.35)
            r0 = 10**(-0.7*(5.5-dimt)**(-0.28) + 0.75)
            dimrho = rho0/(1 + (dimr/r0)**2)**1.5
            if (dimt .lt. 5.5) dens = dimrho*initialDens
        endif

    END SUBROUTINE phys_update

!This FUNCTION works out the time derivative of the density, allowing DVODE to update density with the rest of our ODEs
!It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DVODE
!Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot()
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
