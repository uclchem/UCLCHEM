!Izaskun Cshock
!Requires slightly different parameters to regular models so include cshock_parameters.f90 instead of parameters.f90 in  main.f90

!Is the same as cloud.f90 if you set phase=1, use this to produce self-consistent starting grain abundances for the cshock
!Set phase=2 to run c-shock
MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: tstep,dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see chem_initialise)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,monoflag,volcflag,coflag,tempindx

    !variables either controlled by physics or that user may wish to change    
    double precision :: initdens,dens,tage,tout,t0,t0old,dfin,tfin
    double precision :: size,rout,rin,avic,bc,tstart,maxtemp
    double precision :: tempa(5),tempb(5),codestemp(5),volctemp(5),solidtemp(5)
    double precision, allocatable :: av(:),coldens(:),temp(:)

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    double precision, parameter :: year=3.16455d-08,pc=3.086d18,km=1.d5

    character(2) ::filename
    character(1)  ::densint
    !Cshock specific parameters
    !*******************************************************************
    double precision :: inittemp, z2,vs,v0,zn,vn,at,z3,tout0,tsat
    double precision :: ucm,z1,dv,vi,tempi,vn0,zn0,vA,dlength
    double precision :: radg5,radg,dens6
    double precision, allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)
    !variables for the collisional and radiative heating of grains
    double precision :: mun,tgc0,Frs,tgr0,tgr1,tgr2,tau100,trs0,G0
    double precision :: coshinv1,coshinv2,zmax,a1

    integer :: inrad
    double precision, parameter::nu0=3.0d15,kb2=1.38d-16,bm0=1.e-6,bt=6.
    !*******************************************************************

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    !Set up, calculate size, give dens a kickstart if collapsing 
    SUBROUTINE phys_initialise
        allocate(av(points),coldens(points),temp(points))  

        size=(rout-rin)*pc
        if (collapse .eq. 1) THEN
            if (phase .eq. 2) THEN
                write(*,*) "Cannot have collapse on during cshock (phase=2)"
                Write(*,*) "setting collapse=0 and continuing"
                collapse=0
                dens=initdens
            ELSE
                dens=1.001*initdens
            END IF
        ELSE
            dens=initdens
        ENDIF 

        IF (phase .eq. 2) THEN
            allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))
            mun=2*mh
            radg5=radg/4.e-5
            dens6=dens/1.e6
            tout0=0
        END IF

        !maxxtemp set by vs and pre-shock density, polynomial fits to values taken from Draine et al. 1983
        !have been made and coefficients placed here. Tested with log(dens)>3 <6
        IF (initdens .gt. 10**4.5) THEN
            maxtemp=(2.91731*vs*vs)-(23.78974*vs)+225.204167337
        ELSE
            maxtemp=(0.47258*vs*vs)+(40.44161*vs)-128.635455216
        END IF    
        temp=inittemp

        !tsat proportional to 1/pre-shock density. Fit to tsats from Jimenez-Serra 2008.
        tsat=(-15.38729*vs*vs*vs)+(2069.56962*vs*vs)-(90272.826991*vs)+1686858.54278
        tsat=tsat/initdens

        write(*,*)tsat,maxtemp

        ! The initial parameters that define the C-shock structure
        ! Length of the dissipation region, dlength:
        dlength=12.0*pc*vs/initdens
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

        ! maxtemp is taken from Fig.9b in Draine et al. (1983) and the at constant is
        ! derived as:
        a1=6.0
        at=(1/zmax)*((maxtemp-inittemp)*(dexp(a1)-1.))**(1./6.)

        !write(92,*) 'L=',dlength,'; zn=',z2,'; zi=',z1,'; zT=',z3,'; at=',at

        !Second, we calculate v0 that depends on the alfven and the shock velocities
        !Magnetic field in microGauss. We assume strong magnetic field, i.e., bm0=1.microgauss.
        !(Draine, Roberge & Dalgarno 1983)
        !For the general case, the Alfven velocity is calculated as vA=B0/sqrt(4*pi*2*initdens). If we
        !substitute the expression of B0 on this equation, we obtain that vA=bm0/sqrt(4*pi*mH).
        !B0=bm0*sqrt(2*initdens)
        vA=bm0/sqrt(4*pi*mh)
        vA=vA/km
    END SUBROUTINE phys_initialise


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
            ELSE  IF (tstep .eq. 1) THEN  
                tout=3.16d7*10.0**3
            ELSE
                tout=3.16d7*10.d-8
            ENDIF
        ELSE
            IF (tstep .gt. 160) THEN
                tout=(tage+500)/3.16d-08
            ELSE IF (tstep .gt. 120) THEN
                tout=(tage+100.)/3.16d-08
            ELSE IF (tstep .gt. 80) THEN
                tout=(tage+50.)/3.16d-08
            ELSE IF (tstep .gt. 40) THEN
                tout=(tage+10.)/3.16d-08
            ELSE IF  (tstep.gt.1) THEN
                tout=(tage+1.)/3.16d-08
            ELSE
                tout=3.16d5*10.**(tstep+1)
            ENDIF
        END IF
    END SUBROUTINE timestep

!This is called by main so it should either do any physics the user wants or call the subroutines that do.
!The exception is densdot() as density is integrated with chemistry ODEs.    

    SUBROUTINE phys_update
        !calculate column density. Remember dstep counts from edge of core in to centre
        coldens(dstep)= size*((real(dstep))/real(points))*dens
        !calculate the Av using an assumed extinction outside of core (avic), depth of point and density
        av(dstep)= avic +coldens(dstep)/1.6d21

        !emulate cloud.f90 for phase=1, do cshock for phase=2
        IF (phase .eq. 2) THEN
            !First call shst subroutine 
            call shst
            !Below needed for Izaskun's treatment of the shock
            !We introduce the gas temperature curve along the dissipation region of the
            !C-shock. We also take into account that the gas and dust are decoupled. We
            !use the equations for the collisional and radiative heating of grains of
            !Draine, Roberge & Dalgarno (1983) and Hollenbach, Takahashi & Tielens (1991).
            tn(dstep)=inittemp+((at*zn)**bt)/(dexp(zn/z3)-1)
            ti(dstep)=tn(dstep)+(mun*(dv*km)**2/(3*kb2))
            write(91,*)tn(dstep), ti(dstep)

            !grain collisional heating
            tgc(dstep)=15*(dens6/radg5)**(0.1818)*(tn(dstep)/1000.0)**(0.2727)
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

            write(91,*) av(dstep),vn,tn(dstep),tg(dstep)


            !We introduce the variation of the density (nn) as the C-shock evolves
            IF (tstep .gt. 1) THEN
                dens=initdens*vs/(vs-vn)
                !dens = 1.d5/((1+5.0d-6*tage)**2)
                !write(6,*)dens
            END IF
            IF (tstep.gt.1.) THEN
                tn(dstep)=inittemp+((at*zn)**bt)/(dexp(zn/z3)-1)
                temp=tn(dstep)
                ti(dstep)=tn(dstep)+(mun*(dv*km)**2/(3*kb2))
                tempi=ti(dstep)

                write(93,1234) zn/pc,vn,vi,(vn-vi),tage
                write(92,1234) zn/pc,tn(dstep),ti(dstep),dens,tage
                1234 format(5(1x,ES12.4e2))
            ENDIF

            !At tsat, all mantle species evaporated. These flags make chem module aware of it.
            IF (tage .gt. tsat .and. coflag .eq. 0) THEN
                evap=2
                coflag=1
            ENDIF
        ENDIF
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

!the subroutine below has been written by Izaskun Jimenez-Serra.
! subroutine that calculates the distance along the dissipation region
!(zn) and the velocity of the gas as the shock evolves with time.
    SUBROUTINE shst
        double precision :: vn1,f1,f0,xcos,acosh,v01,g1,g2
        !Calculation of v0. We initially assume that this velocity is of 2kms-1. v0 is
        !numerically derived as follows:
        v0=2.
        v01=0
        write(*,*)"in loop1"
        DO WHILE (abs(v0-v01) .ge. 1e-6)
            v01=v0
            g1=-(vA**2*vs**2)/2
            g2=v01**2-v01*vs-vA**2/2

            v0=sqrt(g1/g2)
        END DO
        write(*,*)"out loop1"

        !We calculate the physical structure of the shock
        !set vn1 arbitrarily high to ensure while loop is done at least once
        vn1=1d30
        vn=vn0
        write(*,*)"in loop2"

        DO WHILE (abs(vn-vn1).ge.1.e-14)
            vn1=vn
            f1=vs-vn1
            f0=vs-vn0
            zn=zn0+(tout-tout0)*km*(f1+f0)/2
            xcos=zn/z2
            acosh=0.5*(dexp(xcos)+dexp(-xcos))
            vn=(vs-v0)-((vs-v0)/acosh)
            write(*,*) "looping"
        END  DO
        write(*,*)"out loop2"


        xcos=zn/z1
        acosh=0.5*(dexp(xcos)+dexp(-xcos))
        vi=(vs-v0)-((vs-v0)/acosh)

        !Store all variables as initial values for next iteration
        dv=vi-vn
        zn0=zn
        vn0=vn
        tout0=tout
    END SUBROUTINE SHST
END MODULE physics