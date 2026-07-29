! hotcore simulates a hot core/corino by following the increase in temperature as a function of time
! It can also reproduce the episodic thermal sublimation seen in laboratory experiments for two phase models
! It is based on Viti et al. 2004 and Collings et al. 2004
module hotcore
    use constants, only: dp, PI, AMU, PC, SECONDS_PER_YEAR, K_BOLTZ_SI, MIN_ABUND, &
        ZERO_INNER_RADIUS_ERROR, Lsun
    use DEFAULTPARAMETERS
    use f2py_constants, only: nSpec
    use network, only: bindingEnergy, gasIceList, iceList, THREE_PHASE, solidFractions, &
        volcanicFractions, monoFractions, mass
    !f2py INTEGER, parameter :: dp
    use physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
    & zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, densdot, gastemp, dusttemp, av,&
    &coldens, density_max, get_ngas_r, coldens_internal, get_coldens_external, parcel_radius, radiation, &
    &outer_coldens_for_current_step, av_internal, radfield_internal, get_G0_internal_at_r

    implicit none

    private
    public :: initializePhysics, updatePhysics, updateTargetTime, sublimation, maxTemp, tempIndx, &
        tempA, tempB

    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: solidflag,volcflag,coflag

    !Arrays for phase 2 temp profiles. parameters for equation chosen by index
    !arrays go [1Msun,5, 10, 15, 25,60]
    integer, parameter :: nMasses = 6
    integer :: tempIndx
    ! Initialize with dummy values
    real(dp) :: tempa(nMasses) = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
    real(dp) :: tempb(nMasses) = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
    ! Deprecated solid, volc and codesorption values, can be removed at some point.
    real(dp),parameter :: solidtemp(nMasses)=(/20.0_dp,19.6_dp,19.45_dp,19.3_dp,19.5_dp,20.35_dp/)
    real(dp),parameter :: volctemp(nMasses)=(/84.0_dp,86.3_dp,88.2_dp,89.5_dp,90.4_dp,92.2_dp/)
    real(dp),parameter :: codestemp(nMasses)=(/95.0_dp,97.5_dp,99.4_dp,100.8_dp,101.6_dp,103.4_dp/)
    real(dp), allocatable :: monoFracCopy(:)
    real(dp) :: maxTemp

    ! 1D radiative transfer arrays
    real(dp) :: Td_r
    real(dp) :: U_r
    real(dp), allocatable :: parcelRadius(:)
    real(dp), allocatable :: maximum_Temp(:)

    ! Time sampling control parameters
    ! During heating (gasTemp < maxTemp): log-based stepping with finer resolution
    real(dp) :: timestep_resolution_factor_heating = 2.0_dp  ! steps per decade during heating
    ! After max_temp reached: coarser log-based or fixed stepping
    real(dp) :: timestep_resolution_factor_hot = 1.0_dp      ! steps per decade post-heating
    real(dp) :: timestep_fixed_late_years = 1.0e5_dp            ! for t > 1 Myr: fixed step in years

    ! Radial grid spacing: .false. = linear (default), .true. = logarithmic
    logical :: log_radius_sampling = .false.

contains

    subroutine initializePhysics(successFlag)
        integer, intent(out) :: successFlag
        successFlag=0

        if (enable_radiative_transfer .AND. points>1 .AND. rin <= 0.0_dp) then
            write(*,*) "ERROR: rin must be > 0 when enable_radiative_transfer=True (G0 diverges at r=0)"
            successFlag=ZERO_INNER_RADIUS_ERROR
            return
        end if

        ! Modules not restarted in python wraps so best to reset everything manually.
        if (ALLOCATED(monoFracCopy)) deallocate(monoFracCopy)
        allocate(monoFracCopy(size(monoFractions)))
        coFlag=0  !reset sublimation
        solidFlag=0
        volcFlag=0
        monoFracCopy=monoFractions  !reset monofractions

        ! Allocate 1D arrays if 1D radiative transfer is enabled
        if (enable_radiative_transfer .AND. points>1) then
            if (ALLOCATED(parcelRadius)) deallocate(parcelRadius)
            allocate(parcelRadius(points))

            if (ALLOCATED(maximum_Temp)) deallocate(maximum_Temp)
            allocate(maximum_Temp(points))

            do dstep=1,points
                if (points > 1) then
                    if (log_radius_sampling) then
                        parcelRadius(dstep)=rin*(rout/rin)**((float(dstep)-1.0e0_dp)/float(points-1))
                    else
                        parcelRadius(dstep)= rin + (dstep-1)*(rout-rin)/float(points-1)
                    end if
                else
                    parcelRadius(dstep)= rout
                end if
                parcel_radius(dstep)=parcelRadius(dstep)
            end do
            ! Better fit for 1D:
            tempa(:) = (/3.1417e-2_dp,3.5495e-2_dp,4.9653e-4_dp,9.5928e-4_dp,1.4158e-3_dp,2.817e-3_dp/)
            tempb(:) = (/0.5329_dp,0.5324_dp,0.9_dp,0.9_dp,0.9_dp,0.9_dp/)
        else
            ! 0D/single-point: parcelRadius must always be allocated because updatePhysics
            ! uses it unconditionally for the radial temperature factor (r/rout)^-0.5.
            ! For 0D the single parcel sits at rout, giving factor = 1.0 (no radial correction).
            if (ALLOCATED(parcelRadius)) deallocate(parcelRadius)
            allocate(parcelRadius(max(1, points)))
            parcelRadius(:) = rout
            ! Default values for 0D:
            tempa(:) = (/1.927e-1_dp,4.8560e-2_dp,7.8470e-3_dp,9.6966e-4_dp,1.706e-4_dp,4.74e-7_dp/)
            tempb(:) = (/0.5339_dp,0.6255_dp,0.8395_dp,1.085_dp,1.289_dp,1.98_dp/)
        end if

        if (freefall) density=1.001_dp*initialDens

        if (enable_radiative_transfer .AND. points>1) then
            do dstep=1,points
                density_max(dstep)=get_ngas_r(parcelRadius(dstep),finalDens,density_scale_radius,density_power_index)
                density(dstep)=density_max(dstep)
                maximum_Temp(dstep) = maxTemp
                ! Internal shielding: from protostar to parcel (core-to-edge).
                coldens(dstep)           = coldens_internal(parcelRadius(dstep))
                av_internal(dstep)       = coldens(dstep) / 1.6e21_dp
                ! External shielding: from cloud edge to parcel (edge-to-core), includes baseAv.
                av(dstep)                = baseAv + get_coldens_external(parcelRadius(dstep), finalDens) / 1.6e21_dp
                radfield_internal(dstep) = get_G0_internal_at_r(lum_star*Lsun, parcelRadius(dstep)*PC)
            end do
        end if

        if (tempindx > nMasses) then
            write(*,*) "tempindx was ",tempindx
            write(*,*) "tempindx must be less than",nMasses
            write(*,*) "1=1Msol, 2=5M, 3=10M, 4=15M, 5=25M, 6=60M"
            successFlag=-1
            return
        end if
    end subroutine initializePhysics

    !Called every time loop in main.f90. Sets the timestep for the next output from
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90
    !but the integrator itself chooses an integration timestep.
    subroutine updateTargetTime
        real(dp) :: orderMagnitude, stepSize, resolutionFactor
        logical :: heating_complete

        ! Determine whether the innermost point (closest to protostar, hottest) has reached maxTemp.
        ! For 0D (points=1) this reduces to gasTemp(1) >= maxTemp.
        ! Once heating is complete we switch to coarser sampling; until then we use finer
        ! steps to capture the rapid chemistry changes during the warm-up phase.
        heating_complete = (gasTemp(1) >= maxTemp)

        if (timeInYears >= 1.0e6_dp) then
            ! Beyond 1 Myr: fixed step regardless of heating state
            targetTime = (timeInYears + timestep_fixed_late_years) * SECONDS_PER_YEAR
        else if (timeInYears > 0.0_dp) then
            ! Log-based stepping: step = 10^floor(log10(t)) / factor
            ! Pick resolution factor based on whether heating is still ongoing.
            ! Heating phase uses finer steps to resolve the rapid temperature rise;
            ! post-heating uses coarser steps since chemistry evolves more slowly.
            if (heating_complete) then
                resolutionFactor = timestep_resolution_factor_hot
            else
                resolutionFactor = timestep_resolution_factor_heating
            end if

            orderMagnitude = 10.0_dp**(floor(log10(timeInYears)))
            stepSize = orderMagnitude / resolutionFactor
            ! Snap to exact multiples of stepSize so decade boundaries are always hit cleanly.
            targetTime = (floor(timeInYears / stepSize + 1.0e-10_dp) + 1.0_dp) * stepSize * SECONDS_PER_YEAR
        else
            targetTime = SECONDS_PER_YEAR * 1.0e-7_dp
        end if
    end subroutine updateTargetTime

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine updatePhysics
        !f2py integer, intent(aux) :: points

        ! 1D radiative transfer calculations
        if (enable_radiative_transfer .AND. points>1) then
            density_max(dstep)=get_ngas_r(parcelRadius(dstep),finalDens,density_scale_radius,density_power_index)
            density(dstep)=density_max(dstep)

            ! Internal shielding: from protostar to parcel (core-to-edge).
            coldens(dstep)           = coldens_internal(parcelRadius(dstep))
            av_internal(dstep)       = coldens(dstep) / 1.6e21_dp
            ! External shielding: from cloud edge to parcel (edge-to-core), includes baseAv.
            av(dstep)                = baseAv + get_coldens_external(parcelRadius(dstep), finalDens) / 1.6e21_dp
            radfield_internal(dstep) = get_G0_internal_at_r(lum_star*Lsun, parcelRadius(dstep)*pc)

            ! Dust temperature from internal protostellar radiation; use av_internal for attenuation.
            call radiation(parcelRadius(dstep)*pc, lum_star*Lsun, temp_star, av_internal(dstep), Td_r, U_r)
            maximum_Temp(dstep)=Td_r  !get global variable value for wrap.f90
            maxTemp=maximum_Temp(dstep)  !set the local variable in this routine
        end if

         if (gasTemp(dstep) < maxTemp) then
        !Below we include temperature profiles for hot cores, selected using tempindx
        !They are taken from Viti et al. 2004 with an additional distance dependence from Nomura and Millar 2004.
        !It takes the form T=A(t^B)*[(d/R)^-0.5], where A and B are given below for various stellar masses
            gasTemp(dstep)=parcelRadius(dstep)/rout
            gasTemp(dstep)=gasTemp(dstep)**(-0.5_dp)
            gasTemp(dstep)=initialTemp + ((tempa(tempindx)*(currentTime/SECONDS_PER_YEAR)**tempb(tempindx))*gasTemp(dstep))
            if (gasTemp(dstep) > maxTemp) gasTemp(dstep)=maxTemp
        end if
        dustTemp=gasTemp
        ! IF (.not. heatingFlag) THEN
        !     dustTemp=gasTemp
        ! END IF

        ! 1D diagnostic output
        if (enable_radiative_transfer .AND. points>1) then
            print "(A,I4,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3)", &
                " pt=",dstep, "  t=",timeInYears, "  r=",parcelRadius(dstep), "  ngas=",density(dstep), &
                "  Tdust=",dustTemp(dstep), "  Tgas=",gasTemp(dstep), "  Av_ext=",av(dstep), &
                "  Av_int=",av_internal(dstep), "  Td_max=",Td_r, "  U_max=",U_r
        end if

    end subroutine updatePhysics

    subroutine sublimation(abund, lpoints)
        ! This subroutine mimics episodic thermal desorption if the network is two pahse
        real(dp), intent(inout) :: abund(nspec+1,lpoints)
        integer, intent(in) :: lpoints
        if (.not. THREE_PHASE) then
            if (instantSublimation) then
                instantSublimation=.false.
                call totalSublimation(abund, lpoints)
            else if (coflag /= 2) then
                if (dustTemp(dstep) > solidtemp(tempindx) .and. solidflag /= 2) solidflag=1
                if (dustTemp(dstep) > volctemp(tempindx) .and. volcflag /= 2) volcflag=1
                if (dustTemp(dstep) > codestemp(tempindx)) coflag=1
                call thermalEvaporation(abund, lpoints)
            end if
        end if
    end subroutine sublimation

    subroutine thermalEvaporation(abund, lpoints)
        !Evaporation is based on Viti et al. 2004. A proportion of the frozen species is released into the gas phase
        !in specific events. These events are activated by flags (eg solidflag) which can be set in physics module.
        !The species evaporated are in lists, created by Makerates and based on groupings. see the viti 2004 paper.
        !f2py integer, intent(aux) :: points
        real(dp), intent(inout) :: abund(nspec+1,lpoints)
        integer, intent(in) :: lpoints

            if (sum(abund(iceList,dstep)) > MIN_ABUND) then
                !Solid Evap
                if (solidflag == 1) then
                    call partialSublimation(solidFractions,abund, lpoints)
                    solidflag=2
                end if

                !monotonic evaporation at binding energy of species
                call bindingEnergyEvap(abund, lpoints)

                !Volcanic evap
                if (volcflag == 1) then
                    call partialSublimation(volcanicFractions,abund, lpoints)
                    volcflag=2  !Set flag to 2 to stop it being recalled
                end if

                !Co-desorption
                if (coflag == 1) then
                    call totalSublimation(abund, lpoints)
                    coflag=2
                end if
            end if
    end subroutine thermalEvaporation

    subroutine partialSublimation(fractions, abund, lpoints)
        real(dp), intent(inout) :: abund(nspec+1,lpoints)
        integer, intent(in) :: lpoints
        real(dp), intent(in) :: fractions(:)

        abund(gasIceList,dstep)=abund(gasIceList,dstep)+fractions*abund(iceList,dstep)
        abund(iceList,dstep)=(1.0_dp-fractions)*abund(iceList,dstep)

    end subroutine partialSublimation

    subroutine totalSublimation(abund, lpoints)
        real(dp), intent(inout) :: abund(nspec+1,lpoints)
        integer, intent(in) :: lpoints
        abund(gasIceList,dstep)=abund(gasIceList,dstep)+abund(iceList,dstep)
        abund(iceList,dstep)=MIN_ABUND
    end subroutine totalSublimation

    subroutine bindingEnergyEvap(abund, lpoints)
        real(dp), intent(inout) :: abund(nspec+1,lpoints)
        integer, intent(in) :: lpoints
        real(dp), parameter :: SURFACE_SITE_DENSITY = 1.5e15_dp
        integer :: i
        !Subroutine to handle mono-evaporation. See viti 2004
        real(dp) :: en,newm,expdust,freq,kevap
        integer :: speci
        !mono evaporation at the binding energy of each species
        do i=lbound(iceList,1),ubound(iceList,1)
            speci=iceList(i)
            en=bindingEnergy(i)*K_BOLTZ_SI
            expdust=bindingEnergy(i)/dustTemp(dstep)
            newm = mass(speci) * AMU / 1e3_dp
            freq = sqrt((2*(SURFACE_SITE_DENSITY)*en)/((PI**2)*newm))
            kevap=freq*exp(-expdust)
            if (kevap >= 0.99_dp) then
                abund(gasIceList(i),dstep)=abund(gasIceList(i),dstep)+(monoFracCopy(i)*abund(speci,dstep))
                abund(speci,dstep)=(1.0_dp-monoFracCopy(i))*abund(speci,dstep)
                monoFracCopy(i)=0.0
            end if
        end do
    end subroutine bindingEnergyEvap

end module hotcore
