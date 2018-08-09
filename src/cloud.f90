!Simple physics module. Models points along a 1d line from the centre to edge of a cloud. Assuming the cloud is spherical you can average
!over the points to get a 1d average and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    INTEGER :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    INTEGER :: collapse,switch,phase
    INTEGER :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    INTEGER :: evap,ion,solidflag,volcflag,coflag,tempindx

    !variables either controlled by physics or that user may wish to change
    DOUBLE PRECISION :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime,grainRadius,initialTemp
    DOUBLE PRECISION :: cloudSize,rout,rin,baseAv,bc,olddens,maxTemp,vs
    
    !Arrays for phase 2 temp profiles. parameters for equation chosen by index
    !arrays go [1Msun,5, 10, 15, 25,60]
   
    DOUBLE PRECISION,PARAMETER :: tempa(6)=(/1.927d-1,4.8560d-2,7.8470d-3,9.6966d-4,1.706d-4,4.74d-7/)
    DOUBLE PRECISION,PARAMETER :: tempb(6)=(/0.5339,0.6255,0.8395,1.085,1.289,1.98/)
    DOUBLE PRECISION,PARAMETER :: solidtemp(6)=(/20.0,19.6,19.45,19.3,19.5,20.35/)
    DOUBLE PRECISION,PARAMETER :: volctemp(6)=(/84.0,86.3,88.2,89.5,90.4,92.2/)
    DOUBLE PRECISION,PARAMETER :: codestemp(6)=(/95.0,97.5,99.4,100.8,101.6,103.4/)

    DOUBLE PRECISION, allocatable :: av(:),coldens(:),temp(:),dens(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    DOUBLE PRECISION,PARAMETER ::PI=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    DOUBLE PRECISION, PARAMETER :: year=3.16455d-08,pc=3.086d18

    !variables for collapse modes
    DOUBLE PRECISION :: unitrho,unitr,unitt,c_s
    DOUBLE PRECISION :: dimrho,dimr,dimt,maxdimt
    DOUBLE PRECISION :: rho0,r0
    DOUBLE PRECISION,PARAMETER :: G_N = 6.674d-8

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN.
!YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    SUBROUTINE initializePhysics
        IF (ALLOCATED(av)) DEALLOCATE(av,coldens,temp,dens)
        ALLOCATE(av(points),coldens(points),temp(points),dens(points))
        cloudSize=(rout-rin)*pc
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
                unitt = 1.0/sqrt(4.0*pi*G_N*unitrho)
                unitr = c_s*unitt
                maxdimt = 5.75 - (15.1/(14.0 + log10(finalDens/initialDens)))**(1.0/0.04)
            !ogino, tomisaka & nakamura 1999
            CASE(3)
                unitrho = initialDens*mh
                unitt = 1.0/sqrt(G_N*unitrho)
                unitr = c_s*unitt
                maxdimt = 0.33 - (2.5/(2.5 + log10(finalDens/initialDens)))**(1.0/0.18)
            CASE(4)

                unitrho = initialDens*mh
                unitt = 1.0/sqrt(2.0*pi*G_N*unitrho)
                unitr = c_s*unitt
                maxdimt = 5.5 - (2.1/(1.35 + log10(finalDens/initialDens)))**(1.0/0.28)
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

        !calculate initial column density as distance from core edge to current point * density
        DO dstep=1,points
            coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
        END DO

    END SUBROUTINE

    SUBROUTINE updateTargetTime
    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrator.
    !targetTime in seconds for DLSODE, timeInYears in years for output.
        IF (timeInYears .gt. 1.0d6) THEN 
            targetTime=(timeInYears+1.0d5)/year
        ELSE  IF (timeInYears .gt. 1.0d5) THEN
            targetTime=(timeInYears+1.0d4)/year
        ELSE IF (timeInYears .gt. 1.0d4) THEN
            targetTime=(timeInYears+1000.0)/year
        ELSE IF (timeInYears .gt. 1000) THEN
            targetTime=(timeInYears+100.0)/year
        ELSE IF (timeInYears .gt. 0.0) THEN
            targetTime=(timeInYears*10.0)/year
        ELSE
            targetTime=3.16d7*10.d-8
        ENDIF
    END SUBROUTINE updateTargetTime

    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from core to edge
        !and coldens should be amount of gas from edge to parcel.
        IF (dstep .lt. points) THEN
            !column density of current point + column density of all points further out
            coldens(dstep)=(cloudSize/real(points))*dens(dstep)
            coldens(dstep)=coldens(dstep)+sum(coldens(dstep:points))
        ELSE
            coldens(dstep)=cloudSize/real(points)*dens(dstep)
        END IF

        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21

        IF (phase .eq. 2 .and. temp(dstep) .lt. maxTemp) THEN

        !Below we include temperature profiles for hot cores, selected using tempindx in PARAMETERs.f90
        !They are taken from Viti et al. 2004 with an additional distance dependence from Nomura and Millar 2004.
        !It takes the form T=A(t^B)*[(d/R)^-0.5], where A and B are given below for various stellar masses
            temp(dstep)=(cloudSize/(rout*pc))*(real(dstep)/real(points))
            temp(dstep)=temp(dstep)**(-0.5)
            temp(dstep)=initialTemp + ((tempa(tempindx)*(currentTime*year)**tempb(tempindx))*temp(dstep))

            if (temp(dstep) .gt. solidtemp(tempindx) .and. solidflag .ne. 2) solidflag=1
            if (temp(dstep) .gt. volctemp(tempindx) .and. volcflag .ne. 2) volcflag=1
            if (temp(dstep) .gt. codestemp(tempindx) .and. coflag .ne. 2) coflag=1
        END IF

        !Density update for BE collapse modes
        !first calculate current radius and time in dimensionless units
        dimr = (real(dstep)/real(points))*cloudSize/unitr
        dimt = currentTime/unitt
        if (dimt .gt. maxdimt) dimt = maxdimt

        !Then calculate central density and radius of central region
        !use this to find density at larger radii. Three different cases possible.

        !foster & chevalier 1993
        if (collapse .eq. 2) then
            rho0 = 10**(15*(5.75-dimt)**(-0.04) - 14)
            r0 = 10**(-6.7*(5.75-dimt)**(-0.04) + 6.6)
            dimrho = rho0/(1 + (dimr/r0)**2.5)
            if (dimt .lt. 5.75) dens(dstep) = dimrho*initialDens

        !ogino, tomisaka & nakamura 1999
        else if (collapse .eq. 3) then
            rho0 = 10**(2.5*(0.33-dimt)**(-0.18) - 2.5)
            r0 = 10**(-1.2*(0.33-dimt)**(-0.18) + 1.3)
            dimrho = rho0/(1 + (dimr/r0)**2.5)
            if (dimt .lt. 0.33) dens(dstep) = dimrho*initialDens

        !nakamura, hanawa & takano 1995
        else if (collapse .eq. 4) then
            rho0 = 10**(2.1*(5.5-dimt)**(-0.28) - 1.35)
            r0 = 10**(-0.7*(5.5-dimt)**(-0.28) + 0.75)
            dimrho = rho0/(1 + (dimr/r0)**2)**1.5
            if (dimt .lt. 5.5) dens(dstep) = dimrho*initialDens
        endif

    END SUBROUTINE updatePhysics

!This function works out the time derivative of the density, allowing DVODE to update density with the rest of our ODEs
!It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DVODE
!Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot(density)
        DOUBLE PRECISION, INTENT(IN) :: density
        DOUBLE PRECISION :: densdot
        !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (density .lt. finalDens) THEN
             densdot=bc*(density**4./initialDens)**0.33*&
             &(8.4d-30*initialDens*((density/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=0.0
        ENDIF
    END FUNCTION densdot
END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
