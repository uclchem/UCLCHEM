
!Cshock parametrization
!Based on Jimenez-Serra et al. 2008 A&A 482
!http://adsabs.harvard.edu/abs/2008A&A...482..549J
module cshock_mod
    use constants, only: dp, MH, SECONDS_PER_YEAR, K_BOLTZ, KM, PC, PI
    use DEFAULTPARAMETERS
    use f2py_constants, only: nSpec
    !f2py INTEGER, parameter :: dp
    use network, only: iceList
    use physicscore, only: points, dstep, cloudsize, radfield, h2CRPRateConstant, improvedH2CRPDissociation, &
        & zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, &
        & gastemp, dusttemp, av, coldens
    use sputtering, only: sputterIces, sputteringSetup

    implicit none

    private
    public :: initializePhysics, updatePhysics, updateTargetTime, sublimation, dissipationTime, &
        timestepFactor, vs, minimumPostshockTemp

    real(dp) :: tstart,maxTemp
    real(dp) :: timestepFactor=0.01_dp
    real(dp) :: z2,vs,v0,zn,vn,at,z3,tsat
    real(dp) :: ucm,z1,driftVel,vi,vn0,zn0,vA,dlength,dissipationTime
    real(dp) :: grainRadius5,dens6,dzv
    real(dp), allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)
    logical :: postShock
    real(dp) :: minimumPostshockTemp=0.0_dp
    !variables for the collisional and radiative heating of grains
    real(dp) :: mun,tgc0,Frs,tgr0,tgr1,tgr2,tau100,trs0,G0
    real(dp) :: coshinv1,coshinv2,zmax,eps

    integer :: inrad
    real(dp), parameter :: nu0=3.0e15_dp,bt=6.0_dp
    real(dp), parameter :: grainRadius=1.0e-5_dp

contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks inputs make sense and then calculates a few constants and!
    ! sets up variables for the shock parametrization that follows    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializePhysics(successFlag)
        integer, intent(out) :: successFlag
        !f2py integer, intent(aux) :: points
        real(dp) :: v01,g1,g2

        successFlag=0
        driftVel=0.0_dp
        zn0=0.0_dp
        vn0=0.0_dp

        ! Set cooling variables to off by default (change if reqd)
        postShock = .false.

        !check input sanity and set initial values
        cloudSize=(rout-rin)*PC
        if (freefall) then
            write(*,*) "Cannot have freefall on during cshock"
            write(*,*) "setting freefall=0 and continuing"
            freefall=.false.
        end if
        if (points > 1) then
            write(*,*) "Cannot have more than one point in cshock"
            successFlag=-1
            return
        end if
        density=initialDens


        !cshock initialization
        if (ALLOCATED(tn)) deallocate(tn,ti,tgc,tgr,tg)
        allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))
        mun=2*MH
        grainRadius5=grainRadius/4.0e-5_dp
        dens6=density(dstep)/1.0e6_dp
        currentTimeOld=0.0_dp
        driftVel=0.0_dp
        zn0=0.0_dp
        vn0=0.0_dp

        !maxtemp set by vs and pre-shock density, polynomial fits to values taken from Draine et al. 1983
        !have been made and coefficients placed here. Tested with log(dens)>3 <6
        !Fits only available for density of 1e4 and 1e6 so in between we average
        if (initialDens > 10**5.5_dp) then
            maxTemp=(2.91731_dp*vs**2)-(23.78974_dp*vs)+225.204167337_dp
        else if (initialDens > 10**4.5_dp) then
            maxTemp=(3.38989_dp*vs**2)+(16.6519_dp*vs)+96.569_dp
            maxTemp=0.5_dp*maxTemp
        else
            maxTemp=(0.47258_dp*vs**2)+(40.44161_dp*vs)-128.635455216_dp
        end if

        !tsat proportional to 1/pre-shock density. Fit to tsats from Jimenez-Serra 2008.
        tsat=(-15.38729_dp*vs**3)+(2069.56962_dp*vs**2)-(90272.826991_dp*vs)+1686858.54278_dp
        tsat=tsat/initialDens

        ! The initial parameters that define the C-shock structure
        ! Length of the dissipation region, dlength:
        dlength=12.0_dp*PC*vs/initialDens
        dissipationTime=(dlength*1.0e-5_dp/vs)/SECONDS_PER_YEAR

        ! Parameters that describe the decoupling between the ion and the neutral
        ! fluids. z2 is obtained by assuming that at z=dlength, the velocity of
        ! the neutrals is 99% (vs-v0). See v0 below and more details in
        ! Jimenez-Serra et al. (2008).
        coshinv1=log((1/0.01_dp)+sqrt((1/0.01_dp)**2-1))
        z2=dlength/coshinv1
        !We assume that z2/z1=4.5 (Jimenez-Serra et al. 2008).
        z1=z2/4.5_dp

        ! zmax is the distance at which Tn reaches its maximum. This happens when
        ! the neutral fluid reaches velocities that are almost 0.85% (vs-v0)
        coshinv2=log((1/0.15_dp)+sqrt((1/0.15_dp)**2-1))
        zmax=dlength/coshinv2

        ! z3 has to be 1/6 zmax
        z3=zmax/6

        ! maxTemp is taken from Fig.9b in Draine et al. (1983) and the at constant is
        ! derived as:
        ! a1=6.0_dp
        at=(1/zmax)*((maxTemp-initialTemp)*(6.0_dp-1.0_dp))**(1.0_dp/6.0_dp)

        !Second, we calculate v0 that depends on the alfven and the shock velocities
        !Magnetic field in microGauss. We assume strong magnetic field, i.e., bm0=1.microgauss.
        !(Draine, Roberge & Dalgarno 1983)
        !For the general case, the Alfven velocity is calculated as vA=B0/sqrt(4*PI*2*initialDens). If we
        !substitute the expression of B0 on this equation, we obtain that vA=bm0/sqrt(4*PI*mH).
        !B0=bm0*sqrt(2*initialDens)
        bm0 = bm0*1e-6_dp
        vA=bm0/sqrt(4*PI*MH)
        vA=vA/KM

        !Calculation of v0, final velocity of ions/neutrals in the shock frame
        v0=2.0_dp
        v01=0
        do while (abs(v0-v01) >= 1e-6_dp)
            v01=v0
            g1=-(vA**2*vs**2)/2
            g2=v01**2-v01*vs-vA**2/2
            v0=sqrt(g1/g2)
        end do

        call sputteringSetup
    end subroutine initializePhysics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine updateTargetTime
        if (timeInYears < 2.0_dp*dissipationTime) then
            !get a nice sampling along the shock
            targetTime=(timeInYears+timestepFactor*dissipationTime)*SECONDS_PER_YEAR
        else
            targetTime=(1.1_dp*timeInYears)*SECONDS_PER_YEAR
        end if
    end subroutine updateTargetTime

    !Calculate shock properties for current time and set density, temperature and Av
    subroutine updatePhysics
        !First calculate velocity of neutrals and position of shock front at currentTime
        call shst

        dzv=sqrt(2*vn*(vs-v0)-vn**2)
        dzv=1.0_dp/(dzv*KM)
        dzv=z2*((vs-v0)/(vs-v0-vn))*dzv

        !We introduce the gas temperature curve along the dissipation region of the
        !C-shock. We also take into account that the gas and dust are decoupled. We
        !use the equations for the collisional and radiative heating of grains of
        !Draine, Roberge & Dalgarno (1983) and Hollenbach, Takahashi & Tielens (1991).
        tn(dstep)=initialTemp+((at*zn)**bt)/((zn/z3)-1)
        ti(dstep)=tn(dstep)+(mun*(driftVel*KM)**2/(3*K_BOLTZ))

        !grain collisional heating
        tgc(dstep)=15*(dens6/grainRadius5)**(0.1818_dp)*(tn(dstep)/1000.0_dp)**(0.2727_dp)
        !grain radiative heating
        ! Frs=0.25*density(1)*mun*(vn*KM)**3
        ! G0=Frs
        ! trs0=12.2*G0**0.2
        ! tau100=2.7e2_dp*G0/trs0**5
        ! tgr1=8.9e-11_dp*nu0*G0*(1.8*av(dstep))+2.7**5
        ! tgr2=3.4e-2_dp*(0.42-log(3.5e-2_dp*tau100*trs0))*tau100*trs0**6
        ! tgr(dstep)=(tgr1+tgr2)**0.2
        !If we don't include the radiative heating that is characteristic
        !of J-type shocks
        tgr(dstep)=0.0_dp
        !total grain heating
        tg(dstep)=tgc(dstep)+tgr(dstep)

        !Density change as shock evolves
        if (timeInYears > 0.0_dp) then
            density=initialDens*vs/(vs-vn)
        end if

        !temperature change as shock evolves
        if (timeInYears > 0.0_dp) then
            tn(dstep)=initialTemp+((at*zn)**bt)/((zn/z3)-1)
            gasTemp(dstep)=tn(dstep)
            ti(dstep)=tn(dstep)+(mun*(driftVel*KM)**2/(3*K_BOLTZ))

        end if
        postShock = (timeInYears > dissipationTime)

        if ((gasTemp(dstep) < minimumPostshockTemp) .AND. (postShock)) then
            gasTemp(dstep) = minimumPostshockTemp
        end if
        dustTemp=gasTemp
    end subroutine updatePhysics

    !For c-shock, sublimation is simply the sputtering subroutine
    subroutine sublimation(abund,lpoints)
        real(dp), intent(inout) :: abund(nSpec+1,lpoints)
        integer, intent(in) :: lpoints
        real(dp) :: timeDelta
        timeDelta=(currentTime-currentTimeOld)
        if ((sum(abund(iceList,dstep)) > 1e-25_dp) .AND. (driftVel > 0)) then
          call sputterIces(abund(:,dstep),driftVel,gasTemp(dstep),density(dstep),timeDelta)
        end if
        where(abund < 1.0e-50_dp) abund=0.0_dp
    end subroutine sublimation


    !the subroutine below has been written by Izaskun Jimenez-Serra.
    ! subroutine that calculates the distance along the dissipation region
    !(zn) and the velocity of the gas as the shock evolves with time.
    subroutine shst
        real(dp) :: vn1,f1,f0,xcos,acosh
        integer :: loopCount
        !We calculate the physical structure of the shock
        !set vn1 arbitrarily high to ensure while loop is done at least once
        vn1=1e30_dp
        vn=vn0
        loopCount=0
        do while ((abs(vn-vn1)>=1.0e-10_dp) .and. (loopCount < 100))
            vn1=vn
            f1=vs-vn1
            f0=vs-vn0
            zn=zn0+(currentTime-currentTimeOld)*KM*(f1+f0)/2
            xcos=zn/z2
            acosh=0.5_dp*((xcos)+exp(-xcos))
            vn=(vs-v0)-((vs-v0)/acosh)
            loopCount=loopCount+1
        end  do

        xcos=zn/z1
        acosh=0.5_dp*((xcos)+exp(-xcos))
        vi=(vs-v0)-((vs-v0)/acosh)

        !Store all variables as initial values for next iteration
        driftVel=vi-vn
        zn0=zn
        vn0=vn
    end subroutine shst


end module cshock_mod
