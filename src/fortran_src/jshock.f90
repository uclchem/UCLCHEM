
!J-shock parametrization
!Based on James et al. 2019 A&A 634
!https://ui.adsabs.harvard.edu/abs/2020A%26A...634A..17J/abstract
module jshock_mod
    use constants, only: dp, PI, PC, SECONDS_PER_YEAR
    use DEFAULTPARAMETERS
    use f2py_constants, only: nSpec
    use network, only: iceList
    use numerics, only: evaluate_polynomial
    use physicscore, only: points, dstep, cloudsize, radfield, h2CRPRateConstant, improvedH2CRPDissociation, &
        zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, &
        gastemp, dusttemp, av, coldens
    use sputtering, only: sputterIces, sputteringSetup
    implicit none

    private
    public :: initializePhysics, updatePhysics, updateTargetTime, sublimation, vs

    real(dp) :: maxTemp,vMin,mfp,tCool,tShock,maxDens
    real(dp) :: t_lambda, n_lambda

    real(dp) :: vs,v0
    real(dp), allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)

    !*******************************************************************

contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks inputs make sense and then calculates a few constants and!
    ! sets up variables for the shock parametrization that follows    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializePhysics(successFlag)
        !f2py integer, intent(aux) :: points
        integer, intent(out) :: successFlag
        ! allow(C031)
        real(dp), parameter, dimension(5) :: vMinCoeffs = (/0.4254_dp, 0.06183_dp, 0.002478_dp, 3.844e-5_dp, -2.058e-7_dp/)

        successFlag=0
        !Reset variables for python wrap.

        cloudSize=(rout-rin)*PC

        if (freefall) then
            write(*,*) "Cannot have freefall on during jshock"
            write(*,*) "setting freefall=0 and continuing"
            freefall=.false.
        end if
        if (points > 1) then
            write(*,*) "Cannot have more than one point in shock"
            successFlag=-1
            return
        end if

        density=initialDens

        ! Determine the maximum temperature
        maxTemp = (5e3_dp)*(vs/10.0_dp)**2
        currentTimeOld=0.0_dp

        ! Determine minimum velocity
        vMin = abs(evaluate_polynomial(vMinCoeffs, vs))

        ! Determine the shock width (of the order of the mean free path)
        mfp = ((sqrt(2.0_dp)*(1e3_dp)*(PI*(2.4e-8_dp)**2))**(-1))/1e4_dp
        tShock = mfp/(vs*1e5_dp)
        ! Determine shock width
        tCool = (1/initialDens)*1e6_dp*SECONDS_PER_YEAR
        ! Determine the maximum density attained
        maxDens = vs*initialDens*(1e2_dp)
        ! Determine the rate constants
        t_lambda = log(maxTemp/initialTemp)
        n_lambda = log(maxDens/initialDens)


        if (allocated(tn)) deallocate(tn,ti,tgc,tgr,tg)
        allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))

        currentTimeOld=0.0
        call sputteringSetup
    end subroutine initializePhysics

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine updateTargetTime
        if (timeInYears > 1e6_dp) then
            targetTime=(timeInYears+1e5_dp)*SECONDS_PER_YEAR
        else if (timeInYears > 1.0e4_dp) then
            targetTime=(timeInYears+1000.0_dp)*SECONDS_PER_YEAR
        else if (timeInYears > 1.0e3_dp) then
            targetTime=(timeInYears+100.0_dp)*SECONDS_PER_YEAR
        else if (timeInYears*SECONDS_PER_YEAR < tShock) then
            targetTime=currentTime+0.05_dp*tShock
        else
            targetTime=1.1_dp*currentTime
        end if
    end subroutine updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate shock properties for current time and set density, temperature and Av  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine updatePhysics

        ! Determine the shock velocity at the current time
        v0 = max(vMin, vs*(exp(log(vMin/vs)*(currentTime/(finalTime*SECONDS_PER_YEAR)))))

        ! Determine whether shock is still increasing the temperature
        ! Or whether it is in the post-shock cooling phase
        ! Or whether the temperature is now constant
        if (currentTime <= tShock) then
            tn(dstep) = ((currentTime/tShock)**2)*(maxTemp) + initialTemp
            density = (((currentTime/tShock)**3)*(4*initialDens))
            where (density < initialDens) density = initialDens
        else if (currentTime > tShock .AND. currentTime <= tCool) then
            ! Otherwise we're in the cooling phase
            tn(dstep) = maxTemp*exp(-t_lambda*(currentTime/tCool))
            density = (4*initialDens)*exp(n_lambda*(currentTime/tCool))

            ! Ensure the gas does not cool below around 10 K
            tn(dstep) = max(10.0_dp, tn(dstep))

            where(density > maxDens) density = maxDens
        else
            tn(dstep) = 10.0_dp
            density = maxDens
        end if
        gasTemp(dstep)=tn(dstep)
        dustTemp(dstep)=gasTemp(dstep)
    end subroutine updatePhysics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine must be in every physics module.                                !
    ! It receives the abundance array and performs any sublimation related activity   !
    ! In hot core that means following thermalEvaporation subroutine.                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sublimation(abund, lpoints)
        integer, intent(in) :: lpoints
        real(dp), intent(inout) :: abund(nSpec+1,lpoints)
        real(dp) :: timeDelta
        timeDelta=(currentTime-currentTimeOld)

        if ((sum(abund(iceList,dstep)) > 1e-25_dp) .AND. (v0 > 0)) then
          call sputterIces(abund(:,dstep),v0,gasTemp(dstep),density(dstep),timeDelta)
        end if
        where(abund < 1.0e-50_dp) abund=0.0_dp
    end subroutine sublimation

end module jshock_mod
