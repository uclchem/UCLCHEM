!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!J-shock paramterization from James et al. 2019                                         !
!                                                                                              !
!Set phase=2 for jshock. Phase=1 will be a standard cloud model with no change in physical     !
!variables unless collapse is set to 1.                                                        !
!                                                                                              !
!See UCLCHEM manual for reasons to run a phase 1.                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE physics
    USE network
    USE constants
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,phase

    !evap is dummy for defaultparameters.f90, ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when sublimation takes place.
    integer :: instantSublimation,ion,coflag,tempindx

    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime
    double precision :: cloudSize,rout,rin,baseAv,bc,tstart,maxTemp,vMin,mfp,tCool,tShock,d,dMax,maxDens
    double precision :: t_lambda, n_lambda
    double precision, allocatable :: av(:),coldens(:),gasTemp(:),dustTemp(:),density(:)

    character(2) ::filename
    character(1)  ::densint
    !Cshock specific parameters
    !*******************************************************************
    double precision :: initialTemp, z2,vs,v0,zn,vn,at,z3,tsat
    double precision :: ucm,z1,driftVel,vi,tempi,vn0,zn0,vA,dlength
    double precision :: grainRadius5,dens6,grainNumberDensity,dzv,start_vel
    double precision, allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)
    !variables for the collisional and radiative heating of grains
    double precision :: mun,tgc0,Frs,tgr0,tgr1,tgr2,tau100,trs0,G0
    double precision :: coshinv1,coshinv2,zmax,a1,eta,eps,epso,sConst

    integer :: inrad,projectiles(6)
    DOUBLE PRECISION, PARAMETER ::nu0=3.0d15,K_BOLTZ_CGS=1.38d-16,bm0=1.e-6,bt=6.
    DOUBLE PRECISION, PARAMETER :: GAS_DUST_NUMBER_RATIO=1.14d-12,CODES_TEMP=130.0
    DOUBLE PRECISION, PARAMETER :: grainRadius=1.0d-5
    !*******************************************************************

CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks inputs make sense and then calculates a few constants and!
    ! sets up variables for the shock paramterization that follows    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics
        INTEGER :: iLoop

        !Reset variables for python wrap.
        IF (ALLOCATED(av)) deallocate(av,coldens,gasTemp,dustTemp,density)
        allocate(av(points),coldens(points),gasTemp(points),dustTemp(points),density(points))  
        coflag=0 !should reset sputtering
        
        cloudSize=(rout-rin)*pc
        write(*,*) "rout=",rout
        write(*,*) "cloudSize=",cloudSize

        if (collapse .eq. 1) THEN
            write(*,*) "Initialising physics."
            if (phase .eq. 2) THEN
                write(*,*) "Cannot have collapse on during jshock (phase=2)"
                Write(*,*) "setting collapse=0 and continuing"
                collapse=0
                density=initialDens
            ELSE
                density=1.001*initialDens
            END IF
        ELSE
            density=initialDens
        ENDIF

        write(*,*) "Setting temperature and cooling time."
        gasTemp = initialTemp

        ! Determine the maximum temperature
        maxTemp = (5e3)*(vs/10)**2
        write(*,*) "vs=",vs
        write(*,*) "maxTemp=",maxTemp

        ! Determine minimum velocity
        vMin = ((-2.058e-07*(vs**4) + 3.844e-05*(vs**3) - 0.002478*(vs**2) + 0.06183*(vs) - 0.4254)**2)**0.5
        write(*,*) 'vMin=',vMin

        !calculate initial column density as distance from core edge to current point * density
        DO dstep=1,points
            write(*,*) "cloudSize=",cloudSize
            write(*,*) "initialDens=",initialDens
            coldens(dstep)=(real(points-dstep+1)*cloudSize/real(points))*initialDens
            write(*,*) "Determined column density as ", coldens(dstep)
        END DO

        IF (phase .eq. 2) THEN
            if (allocated(tn)) deallocate(tn,ti,tgc,tgr,tg)
            allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))
            mun=2*mh
            grainRadius5=grainRadius/4.e-5
            dens6=density(dstep)/1.e6
            currentTimeOld=0.0

            !Need to find the location of the sputtering projectiles in species arrays
            DO iLoop=1,SIZE(specName)
                IF (specName(iLoop).eq."H2") projectiles(1)=iLoop
                IF (specName(iLoop).eq."HE") projectiles(2)=iLoop
                IF (specName(iLoop).eq."C") projectiles(3)=iLoop
                IF (specName(iLoop).eq."O") projectiles(4)=iLoop
                IF (specName(iLoop).eq."SI") projectiles(5)=iLoop
                IF (specName(iLoop).eq."CO") projectiles(6)=iLoop
            END DO
        END IF

        !maxxtemp set by vs and pre-shock density, polynomial fits to values taken from Draine et al. 1983
        !have been made and coefficients placed here. Tested with log(dens)>3 <6
        ! IF (initialDens .gt. 10**4.5) THEN
        !     maxTemp=(2.91731*vs*vs)-(23.78974*vs)+225.204167337
        ! ELSE
        !     maxTemp=(0.47258*vs*vs)+(40.44161*vs)-128.635455216
        ! END IF
        gasTemp=initialTemp
        dustTemp=gasTemp
        !tsat proportional to 1/pre-shock density. Fit to tsats from Jimenez-Serra 2008.
        tsat=(-15.38729*vs*vs*vs)+(2069.56962*vs*vs)-(90272.826991*vs)+1686858.54278
        tsat=tsat/initialDens        
    END SUBROUTINE initializePhysics

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updateTargetTime
        IF (phase .eq. 1) THEN
            IF (timeInYears .gt. 1.0d6) THEN
                targetTime=(timeInYears+1.0d4)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 10000) THEN
                targetTime=(timeInYears+1000.0)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 1000) THEN
                targetTime=(timeInYears+100.0)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 0.0) THEN
                targetTime=(timeInYears*10)*SECONDS_PER_YEAR
            ELSE
                targetTime=3.16d7*10.d-8
            ENDIF
        ELSE
            IF (timeInYears .gt. 1e6) THEN
                targetTime=(timeInYears+1e5)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 1.0d4) THEN
                targetTime=(timeInYears+1000)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 1.0d3) THEN
                targetTime=(timeInYears+100.)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 1.) THEN
                targetTime=(timeInYears+1.)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 0.001) THEN
                targetTime=(timeInYears+0.1)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 0.00001) THEN
                targetTime=(timeInYears+0.0001)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 0.000001) THEN
               targetTime=(timeInYears+0.0000001)*SECONDS_PER_YEAR
            ELSE IF  (timeInYears.gt.0.0) THEN
                targetTime=(timeInYears+0.0000000001)*SECONDS_PER_YEAR
            ELSE
                targetTime=3.16d-03
            ENDIF
        END IF
    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate shock properties for current time and set density, temperature and Av  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from edge of core in to centre
        !calculate column density. Remember dstep counts from edge of core in to centre
        IF (dstep .lt. points) THEN
            !column density of current point + column density of all points further out
            coldens(dstep)=(cloudSize/real(points))*density(dstep)
            coldens(dstep)=coldens(dstep)+sum(coldens(dstep:points))
        ELSE
            coldens(dstep)=cloudSize/real(points)*density(dstep)
        END IF
      
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21

        ! phase=2 for the J-shock specific calculations
        IF (phase .eq. 2) THEN
            ! Determine the shock velocity at the current time
            v0 = vs*(exp(LOG(vMin/vs)*(currentTime/(finalTime*60*60*24*365))))
            IF (v0 .lt. vMin) THEN
                v0 = vMin
            END IF
            ! Determine the shock width (of the order of the mean free path)
            mfp = ((SQRT(2.0)*(1e3)*(pi*(2.4e-8)**2))**(-1))/1d4
            tShock = mfp/(vs*1d5)
            ! Determine shock width
            tCool = (1/initialDens)*1d6*(60*60*24*365)
            ! Determine the maximum density attained
            maxDens = vs*initialDens*(1d2)
            ! Determine the rate constants
            t_lambda = LOG(maxTemp/initialTemp)
            n_lambda = LOG(maxDens/initialDens)
            ! Determine whether shock is still increasing the temperature
            ! Or whether it is in the post-shock cooling phase
            ! Or whether the temperature is now constant
            IF (currentTime .le. tShock) THEN
                !write(*,*) "currentTime (",currentTime,") < tShock (",tShock,")"
                tn(dstep) = ((currentTime/tShock)**2)*(maxTemp) + initialTemp
                density = (((currentTime/tShock)**3)*(4*initialDens))

                IF (density(1) .lt. initialDens) THEN
                    density = initialDens
                END IF

            ELSE IF (currentTime .gt. tShock .AND. currentTime .le. tCool) THEN
                ! Otherwise we're in the cooling phase
                tn(dstep) = maxTemp*EXP(-t_lambda*(currentTime/(tCool)))
                density = (4*initialDens)*EXP(n_lambda*(currentTime/(tCool)))

                ! Ensure the gas does not cool below around 10 K
                IF (tn(dstep) .le. 10) THEN
                    tn(dstep) = 10
                END IF

                IF (density(dstep) .gt. maxDens) THEN
                    density = maxDens
                END IF

            ELSE
                tn(dstep) = 10
                density = maxDens
            END IF
            gasTemp(dstep)=tn(dstep)
            dustTemp(dstep)=gasTemp(dstep)
            IF (timeInYears .gt. 0) THEN
                write(92,1234) tn(dstep),density(dstep),timeInYears
                1234 format(3(e16.9))
            ENDIF
        ENDIF
    END SUBROUTINE updatePhysics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine must be in every physics module.                                !
    ! It receives the abundance array and performs any sublimation related activity   !
    ! In hot core that means following thermalEvaporation subroutine.                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE sublimation(abund)
        DOUBLE PRECISION :: abund(nspec+1,points)
        INTENT(INOUT) :: abund

        IF (coflag .ne. 1 .and. phase .eq. 2) THEN
            IF (gasTemp(dstep) .gt. CODES_TEMP) THEN
                coflag=1
                abund(gasIceList,dstep)=abund(gasIceList,dstep)+abund(iceList,dstep)
                abund(iceList,dstep)=1d-30
            ELSE
                IF ((sum(abund(iceList,dstep)) .gt. 1d-25) .AND. (driftVel .gt. 0)) CALL sputtering(abund)
            END IF
        END IF
        WHERE(abund<1.0d-30) abund=1.0d-30
    END SUBROUTINE sublimation


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine that will sputter the ices based in their entirety
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE sputtering(abund)
        DOUBLE PRECISION :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
        abund(gasIcelist,dstep)=abund(gasIcelist,dstep)+abund(iceList,dstep)
        abund(iceList,dstep)=abund(iceList,dstep)-abund(iceList,dstep)
    END SUBROUTINE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Returns the time derivative of the density.                                     !
    ! Analytical function taken from Rawlings et al. 1992                             !
    ! Called from chemistry.f90, density integrated with abundances so this gives ydot!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
