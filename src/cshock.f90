!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Cshock paramterization from Jimenez-Serra et al. 2008                                         !
!                                                                                              !
!Set phase=2 for cshock. Phase=1 will be a standard cloud model with no change in physical     !
!variables unless collapse is set to 1.                                                        !
!                                                                                              !
!See UCLCHEM manual for reasons to run a phase 1.                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE physics
    USE network
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,phase
    integer :: h2desorb,crdesorb,uvcr,desorb

    !evap is dummy for defaultparameters.f90, ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when sublimation takes place.
    integer :: instantSublimation,ion,coflag,tempindx

    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime
    double precision :: cloudSize,rout,rin,baseAv,bc,tstart,maxTemp
    double precision, allocatable :: av(:),coldens(:),temp(:),density(:)

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::PI=3.141592654,MH=1.674d-24
    double precision, parameter :: YEAR=3.16455d-08,PC=3.086d18,KM=1.d5,SECONDS_PER_YEAR=3.16d7

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
        DOUBLE PRECISION :: v01,g1,g2
        !Reset variables for python wrap.
        IF (ALLOCATED(av)) deallocate(av,coldens,temp,density)
        allocate(av(points),coldens(points),temp(points),density(points))  

        !check input sanity and set inital values
        cloudSize=(rout-rin)*pc
        if (collapse .eq. 1) THEN
            if (phase .eq. 2) THEN
                write(*,*) "Cannot have collapse on during cshock (phase=2)"
                Write(*,*) "setting collapse=0 and continuing"
                collapse=0
                density=initialDens
            ELSE
                density=1.001*initialDens
            END IF
        ELSE
            density=initialDens
        ENDIF 
        temp=initialTemp


        !calculate initial column density as distance from core edge to current point * density
        DO dstep=1,points
            coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
        END DO

        !cshock initialization
        IF (phase .eq. 2) THEN
            !No freeze out during shock.
            evap=1

            IF (ALLOCATED(tn)) deallocate(tn,ti,tgc,tgr,tg)
            allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))
            mun=2*mh
            grainRadius5=grainRadius/4.e-5
            dens6=density(dstep)/1.e6
            currentTimeOld=0.0
            driftVel=0.0
            zn0=0.0
            vn0=0.0

            !maxtemp set by vs and pre-shock density, polynomial fits to values taken from Draine et al. 1983
            !have been made and coefficients placed here. Tested with log(dens)>3 <6
            !Fits only available for density of 1e4 and 1e6 so in between we average
            IF (initialDens .gt. 10**5.5) THEN
                maxTemp=(2.91731*vs*vs)-(23.78974*vs)+225.204167337
            ELSE IF (initialDens .gt. 10.0**4.5) THEN
                maxTemp=(3.38989*vs*vs)+(16.6519*vs)+96.569
                maxTemp=0.5*maxTemp
            ELSE
                maxTemp=(0.47258*vs*vs)+(40.44161*vs)-128.635455216
            END IF    
        
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

            !Calculation of v0, final velocity of ions/neutrals in the shock frame
            v0=2.
            v01=0
            DO WHILE (abs(v0-v01) .ge. 1e-6)
                v01=v0
                g1=-(vA**2*vs**2)/2
                g2=v01**2-v01*vs-vA**2/2
                v0=sqrt(g1/g2)
            END DO

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
    END SUBROUTINE initializePhysics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updateTargetTime
        IF (phase .eq. 1) THEN
            IF (timeInYears .gt. 1.0d6) THEN
                targetTime=(timeInYears+1.0d5)*SECONDS_PER_YEAR
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
            IF (timeInYears .gt. 1.0d5) THEN
                targetTime=(timeInYears+1.0d4)*SECONDS_PER_YEAR
            ELSE IF (timeInYears.gt. 1.0d4) THEN
                targetTime=(timeInYears+1000.)*SECONDS_PER_YEAR                
            ELSE IF (timeInYears .gt. 1000) THEN
                targetTime=(timeInYears+50.)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 100) THEN
                targetTime=(timeInYears+1.)*SECONDS_PER_YEAR
            ELSE IF (timeInYears .gt. 10) THEN
                targetTime=(timeInYears+.1)*SECONDS_PER_YEAR
            ELSE IF  (timeInYears.lt.tsat) THEN
                targetTime=(timeInYears+0.5)*SECONDS_PER_YEAR
            ELSE
                targetTime=(timeInYears+.1)*SECONDS_PER_YEAR
            ENDIF
        END IF
    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate shock properties for current time and set density, temperature and Av  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics
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

        !All cshock logic only happens if user sets phase=2 otherwise: basic cloud
        IF (phase .eq. 2) THEN
            !First calculate shock velocity and position of shock front at currentTime
            call shst

            dzv=sqrt(2*vn*(vs-v0)-vn**2)
            dzv=1./(dzv*km)
            dzv=z2*((vs-v0)/(vs-v0-vn))*dzv

            !number density of dust grains, rearrange following:
            !100.0*dust grain mass*dustgrainnumberdensity=hydrogen mass * h nuclei number density
            grainNumberDensity=density(dstep)*GAS_DUST_NUMBER_RATIO
            !We introduce the gas temperature curve along the dissipation region of the
            !C-shock. We also take into account that the gas and dust are decoupled. We
            !use the equations for the collisional and radiative heating of grains of
            !Draine, Roberge & Dalgarno (1983) and Hollenbach, Takahashi & Tielens (1991).
            tn(dstep)=initialTemp+((at*zn)**bt)/(dexp(zn/z3)-1)
            ti(dstep)=tn(dstep)+(mun*(driftVel*km)**2/(3*K_BOLTZ_CGS))

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

            !Density change as shock evolves
            IF (timeInYears .gt. 0.0) THEN
                density=initialDens*vs/(vs-vn)
            END IF

            !temperature change as shock evolves
            IF (timeInYears .gt. 0.0) THEN
                tn(dstep)=initialTemp+((at*zn)**bt)/(dexp(zn/z3)-1)
                temp(dstep)=tn(dstep)
                ti(dstep)=tn(dstep)+(mun*(driftVel*km)**2/(3*K_BOLTZ_CGS))
                tempi=ti(dstep)
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
            IF (temp(dstep) .gt. CODES_TEMP) THEN
                coflag=1
                abund(gasGrainList,dstep)=abund(gasGrainList,dstep)+abund(grainList,dstep)
                abund(grainList,dstep)=1d-30
            ELSE
                IF (sum(abund(grainList,dstep)) .gt. 1d-25) CALL sputtering(abund)
            END IF
        END IF
        WHERE(abund<1.0d-30) abund=1.0d-30
    END SUBROUTINE sublimation


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


!the subroutine below has been written by Izaskun Jimenez-Serra.
! subroutine that calculates the distance along the dissipation region
!(zn) and the velocity of the gas as the shock evolves with time.
    SUBROUTINE shst
        double precision :: vn1,f1,f0,xcos,acosh
        !We calculate the physical structure of the shock
        !set vn1 arbitrarily high to ensure while loop is done at least once
        vn1=1d30
        vn=vn0
        start_vel=vn0
        DO WHILE (abs(vn-vn1).ge.1.e-14)
            vn1=vn
            f1=vs-vn1
            f0=vs-vn0
            zn=zn0+(currentTime-currentTimeOld)*km*(f1+f0)/2
            xcos=zn/z2
            acosh=0.5*(dexp(xcos)+dexp(-xcos))
            vn=(vs-v0)-((vs-v0)/acosh)
        END  DO

        xcos=zn/z1
        acosh=0.5*(dexp(xcos)+dexp(-xcos))
        vi=(vs-v0)-((vs-v0)/acosh)

        !Store all variables as initial values for next iteration
        driftVel=vi-vn
        zn0=zn
        vn0=vn
    END SUBROUTINE SHST

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Mock subroutine for cshock that will sputter the ices based on Jimenez-Serra!
    ! 2008 paper.                                                                 !
    !
    ! TO DO:          
    !  - actual loop to get masses and abundances for each projectile and do calculation!
    !  - More careful decisions over what should be module variable rather than Function
    !  - Tests against the izaskun paper    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE sputtering(abund)
        DOUBLE PRECISION :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
        DOUBLE PRECISION :: sputterRate,abundChangeFrac,totalMantle
        INTEGER :: iSpec
        !Constant relating mass and speed of projectile to energy
        sConst=(driftVel*driftVel*km*km)/(2.0*temp(dstep)*K_BOLTZ_CGS)
        sConst=sqrt(sConst)

        !loop over projectile species and get rates of change of mantle for each, summing them
        sputterRate=0.0
        DO iSpec=1,SIZE(projectiles) !!!! Make projectiles array in initialize
            sputterRate=sputterRate+iceYieldRate(mass(projectiles(iSpec))*MH,density(dstep)*abund(projectiles(iSpec),dstep))
        END DO

        grainNumberDensity=density(dstep)*GAS_DUST_NUMBER_RATIO
        !Total rate/cm3 (ie released particles /cm3/s) is sputterRate (per grain) multiplied by grain number density
        sputterRate=sputterRate*grainNumberDensity

        !integrate that forward from currentTimeOld to currentTime. to get total number of particles released
        abundChangeFrac=sputterRate*(currentTime-currentTimeOld)!/density(dstep)
        !I think that commented out dens is required for units. However, sputtering doesn't happen if it is uncommented
        !and sputtering matches Jimenez-Serra et al. 2008 curves when it's commented out.


        !if M particles are released and there are N particles on the grain total
        !then a species with X particles on the grain will release M*(X/N)
        !this is M/N and we'll multiply by X below
        totalMantle=sum(abund(grainList,dstep))
        abundChangeFrac=abundChangeFrac/totalMantle
        if (abundChangeFrac .gt. 1) abundChangeFrac=1.0
        !multiply M/N by x and add to gas phase
        abund(gasGrainList,dstep)=abund(gasGrainList,dstep)+abundChangeFrac*abund(grainList,dstep)
        abund(grainList,dstep)=abund(grainList,dstep)-abundChangeFrac*abund(grainList,dstep)
    END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function calculates rate of change of ice mantle abundance of a species!
!due to the impact of molecules of a given mass. actual rate is         !
!proportional to projectile abundance                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION iceYieldRate(projectileMass,projectileDensity)
        DOUBLE PRECISION :: iceYieldRate
        DOUBLE PRECISION projectileMass,projectileDensity
        DOUBLE PRECISION :: lowerLimit,upperLimit,s

        DOUBLE PRECISION, PARAMETER :: iceBindingEnergy=0.53*1.6d-12,targetMass=18.0*MH   
        DOUBLE PRECISION, PARAMETER :: iceYieldEfficiency=0.8 !

        !eta is effectively reduced mass of the collision
        eta=4.*iceYieldEfficiency*projectileMass*targetMass*((projectileMass+targetMass)**(-2.0))
        epso=max(1.,4.*eta)
        s=sConst*sqrt(projectileMass)


        !Lower limit is xth in Jimenez-Serra et al. 2008
        lowerLimit=sqrt(epso*iceBindingEnergy/(eta*K_BOLTZ_CGS*temp(dstep)))

        !Upper limit is just where the integrand goes to zero
        upperLimit=iceYieldIntegralLimit(lowerLimit,projectileMass)

        !calculate eq B.1 from Jimenez-Serra et al. 2008
        IF ((upperlimit-lowerLimit) .gt. 1d-4) THEN
            !first get integral from Eq B.1 including 1/s factor
            iceYieldRate=trapezoidIntegrate(iceYieldIntegrand,lowerLimit,upperLimit,projectileMass)/s
            !multiply through by constants
            iceYieldRate=iceYieldRate*grainRadius*grainRadius*sqrt(8.0*K_BOLTZ_CGS*temp(dstep)*pi/projectileMass)
            !need projectile number density
            iceYieldRate=iceYieldRate*projectileDensity
        ELSE
            iceYieldRate=0.0
        ENDIF
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function calculates integrand from Eq B.1 of Jimenez-Serra et al. 2008 !
!                                                                       !
!Inputs are mass of projectile and x. Returns value of integrand at x   !
!allowing trapezium rule to integrate from xth to infinity              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
FUNCTION iceYieldIntegrand(x,projectileMass)
      DOUBLE PRECISION :: iceYieldIntegrand,x,projectileMass
      DOUBLE PRECISION :: yield,s

      DOUBLE PRECISION, PARAMETER :: yieldConst=8.3d-4
      DOUBLE PRECISION, PARAMETER :: iceBindingEnergy=0.53*1.6d-12

      !this is s from exp(x+s) in eq B.1, varies only with mass so constant precalculated in initialize
      s=sConst*sqrt(projectileMass)

      !epsilon is calculated from inmpact energy (Ep)
      eps=(x**2)*K_BOLTZ_CGS*temp(dstep)
      !and some other factors
      eps=eta*eps/iceBindingEnergy
      !this yield is for ice averaged over all angles. There's a different one for cores (Appendix B Jimenez-Serra 2008)
      !it's 2 times the normal incidence yield, but there's a factor of 0.5 in integrand so we drop both
      yield=yieldConst*((eps-epso)**2)/(1.+((eps/30.)**(1.3333)))
      iceYieldIntegrand=yield*(x**2)*(DEXP(-((x-s)**2))-DEXP(-((x+s)**2)))
END FUNCTION iceYieldIntegrand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function to calculate the upper limit beyond which there's no point   !
!evaluating the ice yield integrand. Ie trapezoids from upper limit to !
!upperlimit+dx will have an area~0                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION iceYieldIntegralLimit(xth,projectileMass)
      DOUBLE PRECISION iceYieldIntegralLimit,xth,projectileMass
      INTEGER :: i
      i=1
      !Take upperlimit to be half way between lowerLimit and 1000.
      iceYieldIntegralLimit=xth+(1d3-xth)*(0.5**i)
      !decrease upper limit for as long as f(upperlimit) is <1.0e-20 and 
      !difference between lower and upper limit is not zero.
      DO WHILE (iceYieldIntegrand(iceYieldIntegralLimit,projectileMass) .lt. 1d-200 .and.&
        & (iceYieldIntegralLimit-xth) .gt. 1.0d-3)
            i=i+1
            iceYieldIntegralLimit=xth+(1d3-xth)*(0.5**i)
      END DO
END FUNCTION iceYieldIntegralLimit


   !Subroutine that calculates an integral using the trapezoidal method.
    Function trapezoidIntegrate(func,lowerLimit,upperlimit,projectileMass)
        DOUBLE PRECISION :: trapezoidIntegrate
        INTEGER JMAX
        DOUBLE PRECISION lowerLimit,upperlimit,func,tolerance,projectileMass
        PARAMETER (tolerance=1.e-3, JMAX=25)
        INTEGER j
        DOUBLE PRECISION olds
        external func
        olds=-1.e30
        DO j=1,JMAX
            call trapzd(func,lowerLimit,upperlimit,trapezoidIntegrate,j,projectileMass)
            if (abs(trapezoidIntegrate-olds).le.tolerance*abs(olds)) RETURN
            olds=trapezoidIntegrate
        END DO
    END

    SUBROUTINE trapzd(func,a,b,s,n,func_arg)
        INTEGER n
        DOUBLE PRECISION a,b,s,func,func_arg
        INTEGER it,j
        DOUBLE PRECISION del,sum,tnm,x
        external func
        IF(n.eq.1) THEN
            s=0.5*(b-a)*(func(a,func_arg)+func(b,func_arg))
        ELSE
            it=2**(n-2)
            tnm=it
            del=(b-a)/tnm
            x=a+0.5*del
            sum=0.
            DO j=1,it
                sum=sum+func(x,func_arg)
                x=x+del
            END DO
            s=0.5*(s+(b-a)*sum/tnm)
        ENDIF
        RETURN
    END
END MODULE physics