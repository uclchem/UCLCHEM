! hotcore simulates a hot core/corino by following the increase in temperature as a function of time
! It can also reproduce the episodic thermal sublimation seen in laboratory experiments for two phase models
! It is based on Viti et al. 2004 and Collings et al. 2004
MODULE hotcore
    USE constants
    USE DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp    
    USE physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
    & zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, densdot, gastemp, dusttemp, av,&
    &coldens, density_max, ngas_r, findcoldens_core2edge
    USE network
    USE f2py_constants
    USE extinction_module
    IMPLICIT NONE
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    INTEGER :: solidflag,volcflag,coflag
    
    !Arrays for phase 2 temp profiles. parameters for equation chosen by index
    !arrays go [1Msun,5, 10, 15, 25,60]
    INTEGER, PARAMETER :: nMasses= 6 
    INTEGER :: tempIndx
    REAL(dp),PARAMETER :: tempa(nMasses)=(/1.927d-1,4.8560d-2,7.8470d-3,9.6966d-4,1.706d-4,4.74d-7/)
    REAL(dp),PARAMETER :: tempb(nMasses)=(/0.5339,0.6255,0.8395,1.085,1.289,1.98/)
    REAL(dp),PARAMETER :: solidtemp(nMasses)=(/20.0,19.6,19.45,19.3,19.5,20.35/)
    REAL(dp),PARAMETER :: volctemp(nMasses)=(/84.0,86.3,88.2,89.5,90.4,92.2/)
    REAL(dp),PARAMETER :: codestemp(nMasses)=(/95.0,97.5,99.4,100.8,101.6,103.4/)
    REAL(dp), allocatable :: monoFracCopy(:)
    REAL(dp) :: maxTemp
    
    ! 1D radiative transfer arrays
    REAL(dp) :: av_max
    REAL(dp) :: Td_r
    REAL(dp) :: U_r
    REAL(dp), allocatable :: parcelRadius(:)
    REAL(dp), allocatable :: coldens_max(:)
    REAL(dp), allocatable :: maximum_Temp(:)

contains

    SUBROUTINE initializePhysics(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        successFlag=0

        ! Modules not restarted in python wraps so best to reset everything manually.
        IF (ALLOCATED(monoFracCopy)) DEALLOCATE(monoFracCopy)
        ALLOCATE(monoFracCopy(size(monoFractions)))
        coFlag=0 !reset sublimation
        solidFlag=0
        volcFlag=0
        monoFracCopy=monoFractions !reset monofractions
        
        ! Allocate 1D arrays if 1D radiative transfer is enabled
        IF (enable_radiative_transfer .AND. points.gt.1) THEN
            IF (ALLOCATED(parcelRadius)) DEALLOCATE(parcelRadius)
            ALLOCATE(parcelRadius(points))

            IF (ALLOCATED(coldens_max)) DEALLOCATE(coldens_max)
            ALLOCATE(coldens_max(points))
            
            IF (ALLOCATED(maximum_Temp)) DEALLOCATE(maximum_Temp)
            ALLOCATE(maximum_Temp(points))

            DO dstep=1,points
                parcelRadius(dstep)=dstep*rout/float(points) !unit of parsec -- note: parcelRadius is from core to edge
            END DO
        END IF

        IF (freefall) density=1.001*initialDens
        
        IF (enable_radiative_transfer .AND. points.gt.1) THEN
            DO dstep=1,points 
                density_max(dstep)=ngas_r(parcelRadius(dstep),finalDens,density_scale_radius,density_power_index)
                density(dstep)=density_max(dstep)
                maximum_Temp(dstep) = maxTemp
            END DO
        END IF

        IF (tempindx .gt. nMasses) THEN
            write(*,*) "tempindx was ",tempindx
            write(*,*) "tempindx must be less than",nMasses
            write(*,*) "1=1Msol, 2=5M, 3=10M, 4=15M, 5=25M, 6=60M"
            successFlag=-1
            RETURN
        END IF 
    END SUBROUTINE

    !Called every time loop in main.f90. Sets the timestep for the next output from   
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 
    !but the integrator itself chooses an integration timestep.                       
    SUBROUTINE updateTargetTime
        IF (timeInYears .gt. 1.0d6) THEN !code in years for readability, targetTime in s
            targetTime=(timeInYears+1.0d5)*SECONDS_PER_YEAR
        ELSE  IF (timeInYears .gt. 1.0d5) THEN
            targetTime=(timeInYears+1.0d4)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1.0d4) THEN
            targetTime=(timeInYears+1000.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1000) THEN
            targetTime=(timeInYears+100.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 100) THEN
            targetTime=(timeInYears+10.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 0.0) THEN
            targetTime=(timeInYears*10.0)*SECONDS_PER_YEAR
        ELSE
            targetTime=SECONDS_PER_YEAR*1.0d-7
        ENDIF
    END SUBROUTINE updateTargetTime

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics
        !f2py integer, intent(aux) :: points
        
        ! 1D radiative transfer calculations
        IF (enable_radiative_transfer .AND. points.gt.1) THEN
            density_max(dstep)=ngas_r(parcelRadius(dstep),finalDens,density_scale_radius,density_power_index)
            density(dstep)=density_max(dstep)

            call findcoldens_core2edge(coldens_max(dstep),rin,finalDens,density_scale_radius,density_power_index,parcelRadius(dstep))
            av_max = coldens_max(dstep)/1.6d21 + baseAv !note: av_max from core to edge

            coldens(dstep)=cloudSize/real(points)*density(dstep) !Note: Ngas from edge to core
            ! add previous column densities to current as we move into cloud to get total
            IF (dstep .lt. points) coldens(dstep)=coldens(dstep)+coldens(dstep+1)

            !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
            av(dstep)= baseAv +coldens(dstep)/1.6d21

            ! Get maximum temperature at a given position
            call radiation(parcelRadius(dstep)*pc, lum_star*Lsun, temp_star, av_max, Td_r, U_r)
            maximum_Temp(dstep)=Td_r !get global variable value for wrap.f90
            maxTemp=maximum_Temp(dstep) !set the local variable in this routine
        END IF

         IF (gasTemp(dstep) .lt. maxTemp) THEN
        !Below we include temperature profiles for hot cores, selected using tempindx
        !They are taken from Viti et al. 2004 with an additional distance dependence from Nomura and Millar 2004.
        !It takes the form T=A(t^B)*[(d/R)^-0.5], where A and B are given below for various stellar masses
            gasTemp(dstep)=(cloudSize/(rout*pc))*(real(dstep)/real(points))
            gasTemp(dstep)=gasTemp(dstep)**(-0.5)
            gasTemp(dstep)=initialTemp + ((tempa(tempindx)*(currentTime/SECONDS_PER_YEAR)**tempb(tempindx))*gasTemp(dstep))
            if (gasTemp(dstep) .gt. maxTemp) gasTemp(dstep)=maxTemp
        END IF
        dustTemp=gasTemp
        ! IF (.not. heatingFlag) THEN 
        !     dustTemp=gasTemp
        ! END IF
        
        ! 1D diagnostic output
        IF (enable_radiative_transfer .AND. points.gt.1) THEN
            PRINT '(A,1PE12.3,A,0PF8.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3,A,1PE12.3)', &
                't=',timeInYears, '  r=',parcelRadius(dstep), '  ngas(t)=',density(dstep), '  Ngas(t)=',coldens(dstep), &
                '  Td(t)=',gasTemp(dstep), '  Td(r)_max=',Td_r, ' U(r)_max=',U_r, ' av(core2edge)=',av_max, ' av(edge2core)=',av(dstep)
        END IF

    END SUBROUTINE updatePhysics

    SUBROUTINE sublimation(abund, lpoints)
        ! This subroutine mimics episodic thermal desorption if the network is two pahse
        REAL(dp), INTENT(INOUT) :: abund(nspec+1,lpoints)
        INTEGER, INTENT(IN) :: lpoints
        IF (.not. THREE_PHASE) THEN
            IF (instantSublimation) THEN
                instantSublimation=.False.
                CALL totalSublimation(abund, lpoints)
            ELSE IF (coflag .ne. 2) THEN
                IF (dustTemp(dstep) .gt. solidtemp(tempindx) .and. solidflag .ne. 2) solidflag=1
                IF (dustTemp(dstep) .gt. volctemp(tempindx) .and. volcflag .ne. 2) volcflag=1
                IF (dustTemp(dstep) .gt. codestemp(tempindx)) coflag=1
                CALL thermalEvaporation(abund, lpoints)
            END IF
        END IF
    END SUBROUTINE sublimation

    SUBROUTINE thermalEvaporation(abund, lpoints)
        !Evaporation is based on Viti et al. 2004. A proportion of the frozen species is released into the gas phase
        !in specific events. These events are activated by flags (eg solidflag) which can be set in physics module.
        !The species evaporated are in lists, created by Makerates and based on groupings. see the viti 2004 paper.
        !f2py integer, intent(aux) :: points
        REAL(dp), INTENT(INOUT) :: abund(nspec+1,lpoints)
        INTEGER, INTENT(IN) :: lpoints
       
            IF (sum(abund(iceList,dstep)) .gt. 1d-30) THEN
                !Solid Evap
                IF (solidflag .eq. 1) THEN
                    CALL partialSublimation(solidFractions,abund, lpoints)
                    solidflag=2
                ENDIF
    
                !monotonic evaporation at binding energy of species
                CALL bindingEnergyEvap(abund, lpoints)
    
                !Volcanic evap
                IF (volcflag .eq. 1) THEN
                    CALL partialSublimation(volcanicFractions,abund, lpoints)
                    volcflag=2 !Set flag to 2 to stop it being recalled
                ENDIF
    
                !Co-desorption
                IF (coflag .eq. 1) THEN
                    CALL totalSublimation(abund, lpoints)
                    coflag=2
                ENDIF
            ENDIF
    END SUBROUTINE thermalEvaporation

    SUBROUTINE partialSublimation(fractions, abund, lpoints)
        REAL(dp), INTENT(INOUT) :: abund(nspec+1,lpoints)
        INTEGER, INTENT(IN) :: lpoints
        REAL(dp), INTENT(IN) :: fractions(:)

        abund(gasiceList,dstep)=abund(gasiceList,dstep)+fractions*abund(iceList,dstep)
        abund(iceList,dstep)=(1.0-fractions)*abund(iceList,dstep)

    END SUBROUTINE partialSublimation

    SUBROUTINE totalSublimation(abund, lpoints)
        REAL(dp), INTENT(INOUT) :: abund(nspec+1,lpoints)
        INTEGER, INTENT(IN) :: lpoints
        abund(gasiceList,dstep)=abund(gasiceList,dstep)+abund(iceList,dstep)
        abund(iceList,dstep)=1d-30
    END SUBROUTINE totalSublimation

    SUBROUTINE bindingEnergyEvap(abund, lpoints)
        REAL(dp), INTENT(INOUT) :: abund(nspec+1,lpoints)
        INTEGER, INTENT(IN) :: lpoints
        REAL(dp), parameter :: SURFACE_SITE_DENSITY = 1.5d15
        INTEGER :: i
        !Subroutine to handle mono-evaporation. See viti 2004
        REAL(dp) en,newm,expdust,freq,kevap
        integer speci
        !mono evaporation at the binding energy of each species
        DO i=lbound(iceList,1),ubound(iceList,1)
            speci=iceList(i)
            en=bindingEnergy(i)*K_BOLTZ_SI
            expdust=bindingEnergy(i)/dustTemp(dstep)
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(SURFACE_SITE_DENSITY)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                abund(gasiceList(i),dstep)=abund(gasiceList(i),dstep)+(monoFracCopy(i)*abund(speci,dstep))
                abund(speci,dstep)=(1.0-monoFracCopy(i))*abund(speci,dstep)
                monoFracCopy(i)=0.0
            END IF 
        END DO
    END SUBROUTINE bindingEnergyEvap

    ! 1D Radiative transfer subroutines
    SUBROUTINE radiation(r, Lstar, Tstar, Avs, temperatures, U)
        IMPLICIT NONE
        REAL(dp) :: Lstar, Tstar, r
        REAL(dp), INTENT(OUT) :: temperatures, U
        INTEGER :: i
        REAL(dp), DIMENSION(:), ALLOCATABLE :: wave, wave_cm, Istar, uwave_star, uwave_red, tau_wave
        REAL(dp) :: ZZ, wave1, Avs, wave2, rsub
        REAL(dp) :: h, sum_odd, sum_even
        REAL(dp) :: uISRF, NH_EBV, RV, Tdsub
        REAL(8)  ::  urad_red
        INTEGER, PARAMETER :: nw=129
        REAL(dp), DIMENSION(2, nw) :: ext_curves
        CHARACTER(LEN=10) :: model

        uISRF = 8.64d-13 !erg cm-3
        Tdsub=2300.d0 !K
        RV = 4.0d0
        NH_EBV = 5.8d21
        
        ! sublimation distance
        rsub = 155.3d0*(Lstar/1.0d6/Lsun)**(0.5) * (Tdsub/1500.d0)**(-5.6/2.0) * aunit !in cm

        ZZ = HP * C / (K_BOLTZ * Tstar)

        wave1 = HP * C * 1.d4 / (13.6d0 * EV)
        wave2 = 20.0d0

        !logspace for wave in micron
        ALLOCATE(wave(nw))
        CALL logspace(LOG10(wave1), LOG10(wave2), nw, wave)

        ! Call the function from the module
        CALL extcurve_obs(wave, RV, NH_EBV, model, ext_curves)

        wave_cm = wave*1.d-4 !in cm

        Istar   = (2.d0*HP*C**2.0/wave_cm**5.0)*(1.d0/(EXP(ZZ/wave_cm)-1.0d0)) !an array
        uwave_star = (4.d0*PI*wave_cm/C)*(Istar)/wave_cm !an array

        tau_wave = Avs * ext_curves(1,:)/1.086d0 !an array
        uwave_red = uwave_star*EXP(-tau_wave) !an array
        
        ! Initialize the integral value
        urad_red = 0.0d0

        ! Apply trapezoidal rule for numerical integration
        DO i = 1, 129-1
            urad_red = urad_red + 0.5d0 * (wave_cm(i+1) - wave_cm(i)) * (uwave_red(i+1) + uwave_red(i))
        END DO
        ! Outputs
        U = urad_red / uISRF * (r / rsub)**(-2.0)
        temperatures = 16.4d0 * U**(1.0/6.0)  ! For silicate grains
    END SUBROUTINE radiation

    SUBROUTINE logspace(start, stop, num, result)
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: start, stop
        INTEGER, INTENT(IN) :: num
        REAL(dp), DIMENSION(num), INTENT(OUT) :: result
        INTEGER :: i

        DO i = 1, num
            result(i) = 10.0d0**(start + (i-1)*(stop-start)/DBLE(num-1))
        END DO
    END SUBROUTINE logspace
END MODULE hotcore 
