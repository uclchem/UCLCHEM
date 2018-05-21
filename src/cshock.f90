!Izaskun Cshock
!Requires slightly different parameters to regular models so include cshock_parameters.f90 instead of parameters.f90 in  main.f90

!Is the same as cloud.f90 if you set phase=1, use this to produce self-consistent starting grain abundances for the cshock
!Set phase=2 to run c-shock
MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,monoflag,volcflag,coflag,tempindx

    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime
    double precision :: cloudSize,rout,rin,baseAv,bc,tstart,maxTemp
    double precision, allocatable :: av(:),coldens(:),temp(:),dens(:)

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18,km=1.d5

    character(2) ::filename
    character(1)  ::densint
    !Cshock specific parameters
    !*******************************************************************
    double precision :: initialTemp, z2,vs,v0,zn,vn,at,z3,targetTime0,tsat
    double precision :: ucm,z1,dv,vi,tempi,vn0,zn0,vA,dlength
    double precision :: grainRadius5,grainRadius,dens6
    double precision, allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)
    !variables for the collisional and radiative heating of grains
    double precision :: mun,tgc0,Frs,tgr0,tgr1,tgr2,tau100,trs0,G0
    double precision :: coshinv1,coshinv2,zmax,a1

    integer :: inrad
    double precision, parameter::nu0=3.0d15,kb2=1.38d-16,bm0=1.e-6,bt=6.
    !*******************************************************************

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    !Set up, calculate cloudSize, give dens a kickstart if collapsing 
    SUBROUTINE initializePhysics
        allocate(av(points),coldens(points),temp(points),dens(points))  

        cloudSize=(rout-rin)*pc
        if (collapse .eq. 1) THEN
            if (phase .eq. 2) THEN
                write(*,*) "Cannot have collapse on during cshock (phase=2)"
                Write(*,*) "setting collapse=0 and continuing"
                collapse=0
                dens=initialDens
            ELSE
                dens=1.001*initialDens
            END IF
        ELSE
            dens=initialDens
        ENDIF 

        !calculate initial column density as distance from core edge to current point * density
        DO dstep=1,points
            coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
        END DO


        IF (phase .eq. 2) THEN
            allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))
            mun=2*mh
            grainRadius5=grainRadius/4.e-5
            dens6=dens(dstep)/1.e6
            targetTime0=0
        END IF

        !maxxtemp set by vs and pre-shock density, polynomial fits to values taken from Draine et al. 1983
        !have been made and coefficients placed here. Tested with log(dens)>3 <6
        IF (initialDens .gt. 10**4.5) THEN
            maxTemp=(2.91731*vs*vs)-(23.78974*vs)+225.204167337
        ELSE
            maxTemp=(0.47258*vs*vs)+(40.44161*vs)-128.635455216
        END IF    
        temp=initialTemp

        !tsat proportional to 1/pre-shock density. Fit to tsats from Jimenez-Serra 2008.
        tsat=(-15.38729*vs*vs*vs)+(2069.56962*vs*vs)-(90272.826991*vs)+1686858.54278
        tsat=tsat/initialDens

        ! The initial parameters that define the C-shock structure
        ! Length of the dissipation region, dlength:
        dlength=12.0*pc*vs/initialDens
        ! Parameters that describe the decoupling between the ion and the neutral
        ! fluids. z2 is obtained by assuming that at z=dlength, the velocity of
        ! the neutrals is 99% (vs-v0). See v0 below and more details in
        ! Jimenez-Serra et al. (2008).
        coshinv1=log((1/0.01)+sqrt((1/0.01)**2-1))
        z2=dlength/coshinv1
        !We assume that z2/z1=4.5 (Jimenez-Serra et al. 2008).
        z1=z2/4.5

        ! zmax is the distance at which Tn reaches its maximum. This happens when
        ! the neutral fluid reaches velocities that are almost 0.85% (vs-v0)
        coshinv2=log((1/0.15)+sqrt((1/0.15)**2-1))
        zmax=dlength/coshinv2

        ! z3 has to be 1/6 zmax
        z3=zmax/6

        ! maxTemp is taken from Fig.9b in Draine et al. (1983) and the at constant is
        ! derived as:
        a1=6.0
        at=(1/zmax)*((maxTemp-initialTemp)*(dexp(a1)-1.))**(1./6.)

        !Second, we calculate v0 that depends on the alfven and the shock velocities
        !Magnetic field in microGauss. We assume strong magnetic field, i.e., bm0=1.microgauss.
        !(Draine, Roberge & Dalgarno 1983)
        !For the general case, the Alfven velocity is calculated as vA=B0/sqrt(4*pi*2*initialDens). If we
        !substitute the expression of B0 on this equation, we obtain that vA=bm0/sqrt(4*pi*mH).
        !B0=bm0*sqrt(2*initialDens)
        vA=bm0/sqrt(4*pi*mh)
        vA=vA/km
    END SUBROUTINE initializePhysics


!This is the time step for outputs from UCL_CHEM NOT the timestep for the integrater. DLSODE sorts that out based on chosen error
!tolerances (RTOL/ATOL) and is simply called repeatedly until it outputs a time >= targetTime. targetTime in seconds for DLSODE, timeInYears in
!years for output.

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
            IF (timeInYears .gt. 1.0d5) THEN
                targetTime=(timeInYears+1.0d4)/year
            ELSE IF (timeInYears.gt. 1.0d4) THEN
                targetTime=(timeInYears+1000.)/year                
            ELSE IF (timeInYears.gt. 1.0d3) THEN
                targetTime=(timeInYears+100.)/year
            ELSE IF (timeInYears .gt. 1000) THEN
                targetTime=(timeInYears+50.)/year
            ELSE IF (timeInYears .gt. 10) THEN
                targetTime=(timeInYears+10.)/year
            ELSE IF  (timeInYears.gt.0.0) THEN
                targetTime=(timeInYears+1.)/year
            ELSE
                targetTime=3.16d6
            ENDIF
        END IF
    END SUBROUTINE updateTargetTime

!This is called by main so it should either do any physics the user wants or call the subroutines that do.
!The exception is densdot() as density is integrated with chemistry ODEs.    

    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from edge of core in to centre
        IF (dstep .lt. points) THEN
            !column density of current point + column density of all points further out
            coldens(dstep)=(cloudSize/real(points))*dens(dstep)
            coldens(dstep)=coldens(dstep)+sum(coldens(dstep:points))
        ELSE
            coldens(dstep)=cloudSize/real(points)*dens(dstep)
        END IF
                !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21

        !emulate cloud.f90 for phase=1, do cshock for phase=2
        IF (phase .eq. 2) THEN
            !First call shst subroutine 
            call shst
            !Below needed for Izaskun's treatment of the shock
            !We introduce the gas temperature curve along the dissipation region of the
            !C-shock. We also take into account that the gas and dust are decoupled. We
            !use the equations for the collisional and radiative heating of grains of
            !Draine, Roberge & Dalgarno (1983) and Hollenbach, Takahashi & Tielens (1991).
            tn(dstep)=initialTemp+((at*zn)**bt)/(dexp(zn/z3)-1)
            ti(dstep)=tn(dstep)+(mun*(dv*km)**2/(3*kb2))

            !grain collisional heating
            tgc(dstep)=15*(dens6/grainRadius5)**(0.1818)*(tn(dstep)/1000.0)**(0.2727)
            !grain radiative heating
            !        Frs=0.25*dens*mun*(vn*km)**3
            !        G0=Frs/Hab
            !        trs0=12.2*G0**0.2
            !        tau100=2.7d2*G0/trs0**5
            !        tgr1=8.9d-11*nu0*G0*dexp(1.8*av(dstep))+2.7**5
            !        tgr2=3.4d-2*(0.42-log(3.5d-2*tau100*trs0))*tau100*trs0**6
            !        tgr(dstep)=(tgr1+tgr2)**0.2
            !If we don't include the radiative heating that is characteristic
            !of J-type shocks
            tgr(dstep)=0.0
            !total grain heating
            tg(dstep)=tgc(dstep)+tgr(dstep)

            !We introduce the variation of the density (nn) as the C-shock evolves
            IF (timeInYears .gt. 0.0) THEN
                dens=initialDens*vs/(vs-vn)
                !dens = 1.d5/((1+5.0d-6*timeInYears)**2)
                !write(6,*)dens
            END IF
            IF (timeInYears .gt. 0.0) THEN
                tn(dstep)=initialTemp+((at*zn)**bt)/(dexp(zn/z3)-1)
                temp(dstep)=tn(dstep)
                ti(dstep)=tn(dstep)+(mun*(dv*km)**2/(3*kb2))
                tempi=ti(dstep)
            ENDIF

            !At tsat, all mantle species evaporated. These flags make chem module aware of it.
            IF (timeInYears .gt. tsat .and. coflag .eq. 0) THEN
                evap=2
                coflag=1
            ENDIF
        ENDIF
    END SUBROUTINE updatePhysics


!This function works out the time derivative of the density, allowing DVODE to update density with the rest of our ODEs
!It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DVODE
!Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot(density)
        double precision, INTENT(IN) :: density
        double precision :: densdot
        !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (density .lt. finalDens) THEN
             densdot=bc*(density**4./initialDens)**0.33*&
             &(8.4d-30*initialDens*((density/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=0.0
        ENDIF
    END FUNCTION densdot

!the subroutine below has been written by Izaskun Jimenez-Serra.
! subroutine that calculates the distance along the dissipation region
!(zn) and the velocity of the gas as the shock evolves with time.
    SUBROUTINE shst
        double precision :: vn1,f1,f0,xcos,acosh,v01,g1,g2
        !Calculation of v0. We initially assume that this velocity is of 2kms-1. v0 is
        !numerically derived as follows:
        v0=2.
        v01=0
        DO WHILE (abs(v0-v01) .ge. 1e-6)
            v01=v0
            g1=-(vA**2*vs**2)/2
            g2=v01**2-v01*vs-vA**2/2

            v0=sqrt(g1/g2)
        END DO

        !We calculate the physical structure of the shock
        !set vn1 arbitrarily high to ensure while loop is done at least once
        vn1=1d30
        vn=vn0

        DO WHILE (abs(vn-vn1).ge.1.e-14)
            vn1=vn
            f1=vs-vn1
            f0=vs-vn0
            zn=zn0+(targetTime-targetTime0)*km*(f1+f0)/2
            xcos=zn/z2
            acosh=0.5*(dexp(xcos)+dexp(-xcos))
            vn=(vs-v0)-((vs-v0)/acosh)
        END  DO

        xcos=zn/z1
        acosh=0.5*(dexp(xcos)+dexp(-xcos))
        vi=(vs-v0)-((vs-v0)/acosh)

        !Store all variables as initial values for next iteration
        dv=vi-vn
        zn0=zn
        vn0=vn
        targetTime0=targetTime
    END SUBROUTINE SHST
END MODULE physics