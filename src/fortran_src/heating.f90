!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module that provides heating and cooling rates		  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module heating
    use CONSTANTS, only: dp, eV, K_BOLTZ, GRAV_G, SB_CONST, T_CMB, MH, PI, &
        COOLANT_CONFIG_ERROR
    use COOLANT_MODULE, only: NCOOLANTS, coolantIndices, coolants, &
        CLOUD_COLUMN, CLOUD_DENSITY, CLOUD_SIZE, coolant_populations_initialized, &
        CHECK_CONVERGENCE, MANAGE_COOLANT_POPULATIONS, UPDATE_COOLANT_LINEWIDTHS, &
        UPDATE_COOLANT_ABUNDANCES, CALCULATE_LEVEL_POPULATIONS, CALCULATE_LINE_OPACITIES, &
        CALCULATE_LINE_OPACITIES, CALCULATE_LAMBDA_OPERATOR, READ_COOLANTS
    use DEFAULTPARAMETERS
    use F2PY_CONSTANTS, only: nSpec, nReac, &
        coolantParentNames, MAX_COOLANTS, coolant_active, coolantNames
    use network, only: specName, exothermicities, REACTIONRATE, enableChemicalHeating, &
        nelec, nh2, nhex, nhx, nh, nhe, nc

    implicit none

    private
    public :: initializeHeating, getTempDot, calculateDustTemp, nHeatingTerms, coolingLabels, &
        coolingValues, heatingValues, NCOOLING, NHEATING, lineCoolingArray, se_coolant_iterations, &
        se_coolant_max_rel_change, chemheating, median_line_index, heatingLabels, &
        cooling_modules, heating_modules, dust_gas_coupling_method, LINE_SOLVER_ATTEMPTS, pahAbund

    real(dp) :: pahAbund=6e-7_dp
    real(dp) :: chemheating

    ! Heating mechanisms:
    ! only ONE photoelectric mechanism should be enabled at once
    ! 1 = Photoelectric - Bakes method
    ! 2 = Photoelectric - Weingartner method
    ! 3 = H2 Formation Heating
    ! 4 = H2 Photodissociation Heating
    ! 5 = H2 FUV Pumping Heating
    ! 6 = Carbon Ionization Heating
    ! 7 = Cosmic Ray Heating
    ! 8 = Turbulent Heating
    ! 9 = Gas-Grain Collisional Heating/Cooling
    integer, parameter :: NHEATING = 9
    real(dp) :: heatingValues(NHEATING)
    logical :: heating_modules(NHEATING)=(/ .true.,.false.,.true.,.true.,.true.,.true.,.true.,.true.,.true./)
    character(LEN=30), parameter :: heatingLabels(NHEATING) = (/ &
        "PhotoelectricBakes            ", &
        "PhotoelectricWeingartner      ", &
        "H2Formation                   ", &
        "H2Photodissociation           ", &
        "H2FUVPumping                  ", &
        "CarbonIonization              ", &
        "CosmicRay                     ", &
        "Turbulent                     ", &
        "GasGrainCollisions            " /)

    ! Cooling mechanisms:
    ! 1 = Atomic Line Cooling
    ! 2 = H2 Collisionally Induced Emission
    ! 3 = Compton Cooling
    ! 4 = Continuum Emission
    ! 5 = Molecular Line Cooling
    integer, parameter :: NCOOLING = 5
    real(dp) :: coolingValues(NCOOLING)
    logical :: cooling_modules(NCOOLING)=(/ .true.,.true.,.true.,.true.,.true./)
    character(LEN=30), parameter :: coolingLabels(NCOOLING) = (/ &
        "AtomicLineEmission            ", &
        "H2CollisionallyInduced        ", &
        "Compton                       ", &
        "ContinuumEmission             ", &
        "MolecularLine                 " /)

    integer, parameter :: nHeatingTerms = 2 + NHEATING + NCOOLING + NCOOLANTS  !Total number of heating and cooling terms tracked including time and chemical heating.

    ! Treatment of the dust-gas temperature coupling
    ! 1 = Simple treatment Hocuk et al. 2017
    ! 2 = Detailed balance method Hollenbach 1991
    ! 3 = Ivlev et al. 2019 (Hocuk + CR heating correction)
    integer :: dust_gas_coupling_method = 3

    ! LINE_SOLVER_ATTEMPTS: Number of times to solve line cooling and take median, must be <= MAX_LINE_SOLVE_ATTEMPTS
    integer :: LINE_SOLVER_ATTEMPTS = 1
    ! Maximum number of line solver attempts (fixed size for arrays)
    integer, parameter :: MAX_LINE_SOLVE_ATTEMPTS = 3

    ! Arrays for the line cooling and median calculation
    real(dp) :: lineCoolingArray(MAX_LINE_SOLVE_ATTEMPTS, MAX_COOLANTS)
    integer :: permutationArray(MAX_LINE_SOLVE_ATTEMPTS)
    real(dp) :: lineCoolingSum(MAX_LINE_SOLVE_ATTEMPTS)
    ! Median line index - must be module-level so io.f90 can access it after getCoolingRate returns
    integer :: median_line_index = 1

    ! SE solver statistics (per-coolant tracking for output arrays)
    integer, dimension(MAX_COOLANTS) :: se_coolant_iterations
    real(dp), dimension(MAX_COOLANTS) :: se_coolant_max_rel_change
    real(dp) :: se_cpu_start, se_cpu_end

 contains

    subroutine initializeHeating(gasTemperature, gasDensity,abundances,columnDensity,cloudSize,successFlag)
        real(dp), intent(in) :: gasTemperature,gasDensity,columnDensity,cloudSize
        real(dp), intent(in) :: abundances(:)
        integer, intent(inout) :: successFlag
        integer :: i, j

        ! write(*,*) "Initializing heating.f90 ..."
        call READ_COOLANTS(successFlag)
        if (successFlag < 0) return

        ! Reset population initialization flag for new model run
        coolant_populations_initialized = .false.

        do i=1,NCOOLANTS
            do j=1,nSpec
                if (coolantParentNames(i) == specName(j)) coolantIndices(i)=j
            end do
        end do

        ! Validate that all coolant indices were successfully initialized
        do i=1,NCOOLANTS
            if (coolantIndices(i) == 0) then
                write(*,"(A,I3,A,A,A,A,A)") &
                    "ERROR: Coolant #", i, " ('", TRIM(coolantNames(i)), &
                    "') could not find parent species '", TRIM(coolantParentNames(i)), &
                    "' in the network. Check your coolant configuration."
                successFlag = COOLANT_CONFIG_ERROR
                return
            end if
        end do

        CLOUD_COLUMN=columnDensity
        CLOUD_DENSITY=gasDensity
        cloud_size=cloudSize

        ! WRITE(*,'(A45,*(L1))') "Coolants enabled:", coolant_active
        ! Moved IO handling to io.f90
    end subroutine initializeHeating


    function getTempDot(time,gasTemperature,gasDensity,gasCols,habingField,abundances,h2dis,h2form,zeta,cIonRate, &
                                & dustAbundance,dustRadius,metallicity, &
                                & dustTemp,turbVel) result(tempDot)
        !Habing field is radfield diminished by Av
        real(dp), intent(in) :: time,gasTemperature,gasDensity,gasCols,habingField,h2dis,h2form,metallicity
        real(dp), intent(in) :: zeta,cIonRate,dustAbundance,dustRadius,dustTemp,turbVel
        real(dp), intent(in) :: abundances(:)  !,exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
        real(dp) :: tempDot

        real(dp) :: adiabaticIdx,heating,cooling

        !First calculate adiabatic index - should use number density but that's just an additional common factor
        adiabaticIdx=5.0_dp*(abundances(nh)+abundances(nhe)+abundances(nelec)+abundances(nh2))+2.0_dp*abundances(nh2)
        adiabaticIdx=adiabaticIdx/&
            (3.0_dp*(abundances(nh)+abundances(nhe)+abundances(nelec)+abundances(nh2))+2.0_dp*abundances(nh2))

        !then calculate overall heating/cooling rate
        heating=getHeatingRate(time,gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate, &
                                & dustAbundance,dustRadius,metallicity, dustTemp,turbVel)

        cooling=0.0
        if (gasTemperature > 3.0_dp) then
            cooling=getCoolingRate(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,h2dis,turbVel)  !6.951290e-17_dp!!
        end if

        chemheating=0.0
        if (enableChemicalHeating) then
            ! Here we need the "flux" in cm^-3 s^-1 to get the heating per reaction.
            ! Energy is in kcal so convert to eV per reaction by multiplying by 1.6022e-12
            chemheating= sum(reactionrate(:nReac) * exothermicities(:))
        end if

        tempDot=heating+chemheating-cooling
        !write(*,*) "Temp Dot",getTempDot
         !and convert to dT/dt
        tempDot=((adiabaticIdx-1.0_dp)*tempDot)/(K_BOLTZ*gasDensity)
    end function getTempDot


    function getHeatingRate(time, gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate, &
            dustAbundance,dustRadius,metallicity,dustTemp,turbVel) &
            result(heatingRate)
        real(dp), intent(in) :: time,gasTemperature,gasDensity,habingField,h2dis,metallicity
        real(dp), intent(in) :: h2form,zeta,cIonRate,dustAbundance,dustRadius,dustTemp,turbVel
        real(dp), intent(in) :: abundances(:)

        real(dp) :: heatingRate

        heatingValues = 0.0_dp
        heatingValues(1) = getPhotoelectricHeatingBakes(gasTemperature,gasDensity,habingField,abundances(nelec),metallicity)
        heatingValues(2) = getPhotoelectricHeatingWeingartner(gasTemperature,gasDensity,habingField,abundances(nelec),metallicity)
        heatingValues(3) = getH2FormationHeating(h2form)
        heatingValues(4) = getH2PhotoDissHeating(gasDensity,abundances(nh2),h2dis)
        heatingValues(5) = getH2FUVPumpHeating(abundances(nh),abundances(nh2),gasTemperature,gasDensity,h2dis)
        heatingValues(6) = getCarbonIonizationHeating(cIonRate,abundances(nc),gasDensity)
        heatingValues(7) = getCosmicRayHeating(zeta,gasDensity,abundances(nh2))
        heatingValues(8) = getTurbulentHeatingRate(gasDensity,turbVel)
        heatingValues(9) = getGasGrainCollisionsHeatingRate(gasTemperature,gasDensity,dustAbundance,dustRadius,dustTemp)
        ! Mask values we do not need.
        where(heating_modules .eqv. .false.) heatingValues = 0.0_dp
        heatingRate=SUM(heatingValues)

    end function getHeatingRate

    function getCoolingRate(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,h2dis,turbVel) &
            result(coolingRate)

        real(dp), intent(in) :: time,gasTemperature,gasDensity,gasCols,dustTemp,h2dis,turbVel
        real(dp), intent(in) :: abundances(:)
        real(dp) :: coolingRate
        integer :: ti, num_attempts
        coolingValues = 0.0_dp
        lineCoolingArray = 0.0_dp
        lineCoolingSum = 0.0_dp

        coolingValues(1)=getAtomicCooling(gasTemperature,gasDensity,abundances(nh),abundances(nhe),&
                        &abundances(nelec),abundances(nhx),abundances(nhex))
        coolingValues(2)=getCollisionallyInducedEmission(gasTemperature,gasDensity,abundances(nh2))
        coolingValues(3)=getComptonCooling(gasTemperature,gasDensity,abundances(nelec))
        coolingValues(4)=getContinuumEmission(gasTemperature,gasDensity)

        !getH2VibrationalCooling is handled by heating FUV pumping function

        ! Only compute expensive line cooling if enabled (guard clause for performance)
        if (cooling_modules(5)) then
            ! Use LINE_SOLVER_ATTEMPTS (up to MAX_LINE_SOLVE_ATTEMPTS) for actual number of solves
            num_attempts = MIN(LINE_SOLVER_ATTEMPTS, MAX_LINE_SOLVE_ATTEMPTS)

            !We do the line cooling multiple times and take median value since solver will occasionally do something wild
            do ti=1,num_attempts
                lineCoolingArray(ti, :NCOOLANTS)= &
                    getLineCooling(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,turbVel)
                lineCoolingSum(ti) = sum(lineCoolingArray(ti, :NCOOLANTS))
            end do
            call pair_insertion_sort_with_perm(lineCoolingSum(1:num_attempts), permutationArray(1:num_attempts))
            median_line_index = permutationArray((num_attempts + 1) / 2)
            coolingValues(5) = lineCoolingSum(median_line_index)
        else
            coolingValues(5) = 0.0_dp
        end if

        ! Mask the first 4 cooling values (line cooling already handled above)
        where(cooling_modules(1:4) .eqv. .false.) coolingValues(1:4) = 0.0_dp
        coolingRate = sum(coolingValues)
    end function getCoolingRate


    function getLineCooling(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,turbVel) &
            result(moleculeCooling)
        real(dp), intent(in) :: time,gasTemperature,gasDensity,gasCols,dustTemp,turbVel
        real(dp), intent(in) :: abundances(:)

        real(dp)  :: moleculeCooling(NCOOLANTS)

        integer :: N,I,level  !, collisionalIndices(5)=(/nh,nhx,nh2,nhe,nelec/)
        real(dp) :: rel_change

        moleculeCooling = 0.0_dp

        ! Update CLOUD_DENSITY and CLOUD_COLUMN
        CLOUD_DENSITY = gasDensity
        CLOUD_COLUMN  = gasCols

        call UPDATE_COOLANT_LINEWIDTHS(gasTemperature,turbVel)
        call UPDATE_COOLANT_ABUNDANCES(gasDensity,gasTemperature,abundances)

        call MANAGE_COOLANT_POPULATIONS(gasTemperature)

        call CALCULATE_LINE_OPACITIES()
        call CALCULATE_LAMBDA_OPERATOR()

        ! Initialize SE solver statistics
        se_coolant_iterations(1:NCOOLANTS) = 0
        se_coolant_max_rel_change(1:NCOOLANTS) = 0.0_dp
        call CPU_TIME(se_cpu_start)

         tries: do I=1,500  !while not converged and less than 100 tries:
            coolant_scf: do N=1,NCOOLANTS
                if (.NOT. coolant_active(N)) cycle coolant_scf
                if (coolants(N)%CONVERGED) then
                    if (se_coolant_iterations(N) == 0) then
                        se_coolant_iterations(N) = I-1  ! Track if converged last iteration
                    end if
                else
                    call CALCULATE_LEVEL_POPULATIONS(coolants(N),gasTemperature,gasDensity,&
                        &abundances,dustTemp)
                end if
                ! Calculate max relative change for this coolant (inline)
                se_coolant_max_rel_change(N) = 0.0_dp
                if (ALLOCATED(coolants(N)%POPULATION) .AND. ALLOCATED(coolants(N)%PREVIOUS_POPULATION)) then
                    do level = 1, SIZE(coolants(N)%POPULATION)
                        if (coolants(N)%POPULATION(level) > 0.0_dp) then
                            rel_change = ABS(coolants(N)%POPULATION(level) - coolants(N)%PREVIOUS_POPULATION(level)) / &
                                        coolants(N)%POPULATION(level)
                            se_coolant_max_rel_change(N) = MAX(se_coolant_max_rel_change(N), rel_change)
                        end if
                    end do
                end if
            end do coolant_scf
            call CALCULATE_LINE_OPACITIES()
            call CALCULATE_LAMBDA_OPERATOR()
            if (CHECK_CONVERGENCE()) then
                where(se_coolant_iterations == 0) se_coolant_iterations = I
                exit tries
            end if


            !IF (I .eq. 499) write(*,*) "Failed convergence"
        end do tries
        call CPU_TIME(se_cpu_end)

        !  Calculate the cooling rate due to the Lyman-alpha emission for each particle
        !  using the analytical expression of Spitzer (1978) neglecting photon trapping
        coolant: do N=1,NCOOLANTS
            if (.NOT. coolant_active(N)) cycle coolant
            if (.not. allocated(coolants(N)%EMISSIVITY)) then
                write(*,*) "ERROR: EMISSIVITY not allocated for coolant ", N
            end if

            if(coolants(N)%NAME=="H") then
               coolants(N)%EMISSIVITY(2,1) = 7.3e-19_dp*(abundances(nelec)*gasDensity) &
                                                          & *(abundances(nH)*gasDensity) &
                                                          & *exp(-118400.0_dp/gasTemperature)
            end if
        end do coolant

        !Calculate the cooling rates
        coolant_2: do N=1,NCOOLANTS
            if (.NOT. coolant_active(N)) then
                moleculeCooling(N) = 0.0_dp
                cycle coolant_2
            end if
            where(coolants(N)%EMISSIVITY < -HUGE(1.0_dp)) coolants(N)%EMISSIVITY = 0.0_dp
            moleculeCooling(N)=SUM(coolants(N)%EMISSIVITY,MASK=.NOT.ISNAN(coolants(N)%EMISSIVITY))
            if (moleculeCooling(N) < 0.0_dp) moleculeCooling(N)=0.0_dp
        end do coolant_2
        where(moleculeCooling < 0.0_dp) moleculeCooling=0.0_dp

    end function getLineCooling

    ! !-----------------------------------------------------------------------
    ! !  Atomic and ionic cooling rates
    ! !  from Neal et al. 1995 based on Cen (1992) via Grassi et al. (2014)
    ! !-----------------------------------------------------------------------
    function getAtomicCooling(gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund) &
            result(atomicCooling)
        real(dp), intent(in) :: gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund
        real(dp) :: atomicCooling
        real(dp) :: t5,invT,rootT,collTFactor  !temp/10^5, 1/T and a weird factor from the table
        real(dp) :: hDens,elecDens,heDens,hxDens,hexDens,gauntFactor
        hDens=gasDensity*hAbund
        elecDens=gasDensity*electronAbund
        heDens=gasDensity*heAbund
        hxDens=gasDensity*hxAbund
        hexDens=gasDensity*hexAbund
        t5=1.0e-5_dp*gasT
        invT=1.0_dp/gasT
        rootT=sqrt(gasT)
        collTFactor=1.0_dp/(1.0_dp+sqrt(t5))

        !gauntFactor from Neal et al. 1995
        gauntFactor=1.1_dp+(0.34_dp*exp(-((5.5_dp-log10(gasT))**2.0_dp)/3.0_dp))
        !Neal et al. 1995 lists several fits to cooling each in ergs/cm3/s so we'll just sum them
        !see table 1 of that paper
        !These are just numerical fits so there's loads of magic numbers
        !I've shorted variable names to make it easier to write/read (tn is temperature/10^n)

        !collisional excitation and ionization
        atomicCooling=(7.5e-19_dp*collTFactor*exp(-118348.0_dp*invT)*elecDens*hDens) &
            &+(5.54e-17_dp*(gasT**-0.397_dp)*collTFactor*exp(-473638.0_dp*invT)*elecDens*hexDens)&
            &+(1.27e-21_dp*rootT*exp(-157809.1_dp*invT)*elecDens*hDens*collTFactor)&
            &+(9.38e-22_dp*rootT*exp(-285335.4_dp*invT)*elecDens*heDens*collTFactor)&
            &+(4.95e-22_dp*rootT*exp(-631515.0_dp*invT)*elecDens*hexDens*collTFactor)&
            !dielectric
            &+(1.24e-13_dp*(gasT**-1.5_dp)*exp(-470000.0_dp*invT)*(1.0_dp+0.3_dp*exp(-94000.0_dp*invT))*elecDens*hexDens)
        if (gasT > 1.0e5_dp) then
        !recombination
            atomicCooling=atomicCooling&
            &+(8.7e-27_dp*rootT*((1.0e-3_dp*gasT)**-0.2_dp)*elecDens*hxDens/(1.0_dp+((0.1_dp*t5)**0.7_dp)))&
            &+(1.55e-26_dp*(gasT**0.3647_dp)*elecDens*hexDens)&
            !&+(3.48e-26_dp*rootT*((0.001*gasT)**-0.2)*nelec*nhexI/(1+(0.1*t5)**0.7))&
            !Free-free emission
            &+(1.42e-27_dp*rootT*nelec*(nhex+nhx)*gauntFactor)
        end if

    end function getAtomicCooling

    !!-----------------------------------------------------------------------
    !! Collisionally Induced Emission
    !! Hirano & Yoshida (2013) and Ripamonti & Abel 2004 via Grassi 2012
    !!-----------------------------------------------------------------------
    function getCollisionallyInducedEmission(gasTemperature,gasDensity,h2Abund) &
            result(collisionallyInducedEmission)
        real(dp), intent(in) :: gasTemperature,gasDensity,h2Abund
        real(dp) :: collisionallyInducedEmission
        integer, parameter :: NUM_FIT_CONSTANTS = 6
        real(dp), parameter :: aConsts(NUM_FIT_CONSTANTS)=(/-30.3314216559651_dp,19.0004016698518_dp,-17.1507937874082_dp&
                                                            &,9.49499574218739_dp,-2.54768404538229_dp,0.265382965410969_dp/)
        real(dp), parameter :: bConsts(NUM_FIT_CONSTANTS)=(/-180.992524120965_dp,168.471004362887_dp,-67.499549702687_dp,&
                                                            &13.5075841245848_dp,-1.31983368963974_dp,0.0500087685129987_dp/)
        real(dp),parameter :: c=3.0_dp,d=21.2968837223113_dp
        real(dp) :: tau,logt
        integer :: i

        logt=log10(gasTemperature)

        tau=(gasDensity*h2Abund/7.0e15_dp)**2.8_dp
        !if (tau.lt.0.2) THEN
            !avoid numerical problems, tau tends to 1 for low tau but fortran can't do it
            !Taylor series is fine for tau<0.2 and that's well above the area we get an issue
         !   tau=1.0-(0.5*tau)+((tau**2.0)/6.0)-((tau**3.0)/24.0)+((tau**4)/120.0)
        !ELSE
            tau=(1.0_dp-exp(-tau))/tau
        !END IF
        tau=min(1.0_dp,tau)


        collisionallyInducedEmission=0.0_dp

        if (gasTemperature >= 1.0e5_dp) then
            collisionallyInducedEmission=(c*logt)-d
        else if (gasTemperature >= 891.0_dp) then
            do i=1,NUM_FIT_CONSTANTS
                collisionallyInducedEmission=collisionallyInducedEmission+(bConsts(i)*(logt**(i-1)))
            end do
        !technically fit below is ok down to 100 K but bad fit seems better than no cooling at 70 K
        else if (gasTemperature >= 100.0_dp) then
            do i=1,NUM_FIT_CONSTANTS
                collisionallyInducedEmission=collisionallyInducedEmission+(aConsts(i)*(logt**(i-1)))
            end do
        end if

        if (gasTemperature >= 100.0_dp) then
            collisionallyInducedEmission=(10.0_dp**collisionallyInducedEmission)*tau
        end if
        !collisionallyInducedEmission=(10.0**collisionallyInducedEmission)*tau
    end function getCollisionallyInducedEmission


    !!-----------------------------------------------------------------------
    !! Continuum Emission
    !! Hirano & Yoshida (2013) and Ripamonti & Abel 2004 via Grassi 2012
    !!-----------------------------------------------------------------------
    function getContinuumEmission(gasTemperature,gasDensity) result(continuumEmission)
        real(dp), intent(in) :: gasTemperature,gasDensity
        real(dp) :: continuumEmission

        real(dp) :: massDensity,opacity,opticalDepth

        massDensity=min(0.5_dp,gasDensity*MH*1.22_dp)  !assume mean molecular weight 1.22 and give max
        opacity=10.0_dp**(1.000042_dp*log(massDensity)+2.14989_dp)  !Lenzuni opacity fit

        opticalDepth=sqrt(3.14159_dp*K_BOLTZ*gasTemperature/(massDensity*MH*1.22_dp*GRAV_G))
        opticalDepth=opticalDepth*opacity*massDensity+1.0e-40_dp  ! stop it going to zero

        continuumEmission=4.0_dp*SB_CONST*(gasTemperature**4)*opacity*massDensity*min((opticalDepth**(-2)),1.0_dp)
    end function getContinuumEmission

    !!-----------------------------------------------------------------------
    !! Compton cooling
    !! Cen 1992 via Grassi 2012
    !! Cooling due to compton scattering of CMB photons.
    !! Shouldn't be important in near universe but include for completion
    !!-----------------------------------------------------------------------
    function getComptonCooling(gasTemperature,gasDensity,elecAbund) result(comptonCooling)
        real(dp), intent(in) :: gasTemperature,gasDensity,elecAbund
        real(dp) :: comptonCooling

        comptonCooling=1.017e-37_dp*(T_CMB**4.0_dp)*(gasTemperature-T_CMB)*elecAbund*gasDensity
    end function getComptonCooling

    ! !-----------------------------------------------------------------------
    ! !  Grain + PAH photoelectric heating (with graphitic and silicate grains)
    !
    !  Use the treatment of Bakes & Tielens (1994, ApJ, 427, 822) with the
    !  modifications suggested by Wolfire et al. (2003, ApJ, 587, 278) to
    !  account for the revised PAH abundance estimate from Spitzer data.
    !
    ! !-----------------------------------------------------------------------
    function getPhotoelectricHeatingBakes(gasTemperature,gasDensity,habingField,electronAbund,metallicity) &
            result(photoelectricHeatingBakes)
        real(dp), intent(in) :: gasTemperature,gasDensity,habingField,electronAbund,metallicity
        real(dp) :: photoelectricHeatingBakes
        real(dp), parameter :: PHI_PAH=0.4_dp,ALPHA=0.944_dp
        real(dp) :: beta,delta,epsilon,nElec,PAH_HEATING_RATE,PAH_COOLING_RATE

        !Bakes & Tielens 1994 with updates from Wolfire 2008
        !  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680, 384)
        !  Setting this factor to 1.0 gives the standard Bakes & Tielens expression

        !Skip photoelectric heating for now if electron abundance is zero because cooling is infinite
        !The only way to have no E- is bad initial conditions so this will resolve itself within a time step
        if (electronAbund > 1.0e-20_dp) then
            nElec=electronAbund*gasDensity
            BETA=0.735_dp/gasTemperature**0.068_dp
            DELTA=habingField*sqrt(gasTemperature)/(nElec*PHI_PAH)
            EPSILON=4.87e-2_dp/(1.0_dp+4.0e-3_dp*DELTA**0.73_dp) + &
                3.65e-2_dp*(gasTemperature/1.0e4_dp)**0.7_dp/(1.0_dp+2.0e-4_dp*DELTA)

            PAH_HEATING_RATE=1.30e-24_dp*EPSILON*habingField*gasDensity
            PAH_COOLING_RATE=4.65e-30_dp*gasTemperature**ALPHA*(DELTA**BETA)*nElec*PHI_PAH*gasDensity
            photoelectricHeatingBakes=PAH_HEATING_RATE-PAH_COOLING_RATE

            !Assume the PE heating rate scales linearly with PAH abundance
            photoelectricHeatingBakes=photoelectricHeatingBakes*metallicity  !*(pahAbund/6.0e-7_dp)
        else
            photoelectricHeatingBakes=0.0
        end if
    end function getPhotoelectricHeatingBakes

    function getPhotoelectricHeatingWeingartner(gasTemperature,gasDensity,habingField,electronAbund,metallicity) &
            result(photoelectricHeatingWeingartner)
        real(dp), intent(in) :: gasTemperature,gasDensity,habingField,electronAbund,metallicity
        real(dp) :: photoelectricHeatingWeingartner
        real(dp), parameter :: C0=5.72e+0_dp,C1=3.45e-2_dp,C2=7.08e-3_dp
        real(dp), parameter :: C3=1.98e-2_dp, C4=4.95e-1_dp,C5=6.92e-1_dp
        real(dp), parameter :: C6=5.20e-1_dp
        !Weingartner & Draine2001
        photoelectricHeatingWeingartner=1.0e-26_dp*(habingField*gasDensity)*(C0+C1*gasTemperature**C4) &
                & /(1.0_dp+C2*(habingField*sqrt(gasTemperature)/(gasDensity*electronAbund))**C5  &
                & *(1.0_dp+C3*(habingField*sqrt(gasTemperature)/(gasDensity*electronAbund))**C6))
    end function getPhotoelectricHeatingWeingartner

    ! !-----------------------------------------------------------------------
    ! !  H2 formation heating
    ! !
    ! !  Energy released to gas depends on formation mechanism and thermalization efficiency.
    ! !  Only H2 that desorbs to the gas phase contributes; LH/ER products remain on grain.
    ! !  Pre-computed in chemistry.f90 using Hollenbach & McKee (1979) eq. 6.43:
    ! !
    ! !    CT    = 1.5 eV (Hollenbach & Tielens 1999), no thermalization correction
    ! !    LHDes = (0.1 + 4.2 * h2heatfac) eV  (H&M79 eq. 6.43)
    ! !              0.1 eV kinetic, 4.2 eV vibrational thermalized at fraction h2heatfac
    ! !    ERDes = 0.6 eV * h2heatfac  (Bourlot et al. 2012)
    ! !
    ! !  where h2heatfac = 1/(1 + n_cr/n), n_cr from H&M79 eq. 6.45:
    ! !    n_cr = 10^6 T^{-1/2} / (1.6*x_H*exp(-(400/T)^2) + 1.4*x_H2*exp(-18100/(T+1200)))
    ! !
    ! !  This function receives the already-computed heating rate [erg cm^-3 s^-1].
    ! !-----------------------------------------------------------------------
    function getH2FormationHeating(h2formHeat) result(H2FormationHeating)
        real(dp), intent(in) :: h2formHeat  ! mechanism-weighted H2 formation heating [erg cm^-3 s^-1]
        real(dp) :: H2FormationHeating
        H2FormationHeating = h2formHeat
    end function getH2FormationHeating


    !-----------------------------------------------------------------------
    !  H2 vibrational heating/cooling
    !
    !  Treat the vibrationally excited levels of H2 as a single pseudo level
    !  with effective rates of spontaneous emission, collisional excitation,
    !  FUV pumping and photodissociation that describe the behavior of all
    !  the vibrational levels combined.
    !
    !  Use the treatment of Rollig et al. (2006, A&A, 451, 917)
    !-----------------------------------------------------------------------
    function getH2VibrationalCooling(gasTemperature,gasDensity,h2Abund,h2dis) &
            result(H2VibrationalCooling)
        real(dp), intent(in) :: gasTemperature,gasDensity,h2Abund,h2dis
        real(dp) :: H2VibrationalCooling
        real(dp) :: photoDisRate,DELTA_E_10,A_COEFF_10,C_COEFF_10
        real(dp) :: DELTA_E_EFF,A_COEFF_EFF,R_PUMP_EFF,R_PHOTO_EFF

        DELTA_E_10=6587.0_dp  ! Energy gap (K) between the v=1 and v=0 levels of H2
        A_COEFF_10=8.6e-7_dp  ! Einstein A-coefficient (s^-1) for emission from the v=1 to v=0 level of H2
        C_COEFF_10=5.4e-13_dp*sqrt(gasTemperature)  ! Collisional rate coefficient (cm^3 s^-1) for v=0 to v=1
        photoDisRate=h2dis  ! Photodissociation rate (s^-1) for the v=1 level of H2

        DELTA_E_EFF=23500.0_dp  ! Characteristic vibrational level energy (K)
        A_COEFF_EFF=1.9e-6_dp  ! Effective Einstein A-coefficient (s^-1)

        R_PUMP_EFF=11.2_dp*h2dis  ! Effective vibrational pumping rate (s^-1)
        R_PHOTO_EFF=18.0_dp*h2dis  ! Effective photodissociation rate (s^-1)

        H2VibrationalCooling=K_BOLTZ*DELTA_E_10*C_COEFF_10*gasDensity*exp(-DELTA_E_10/gasTemperature)* &
            h2Abund*gasDensity*(A_COEFF_10+h2dis)/(C_COEFF_10*gasDensity+A_COEFF_10+h2dis)

        !Some heating from H2 vibrational interactions so subtract from cooling rate
        H2VibrationalCooling=H2VibrationalCooling-(h2Abund*gasDensity*(R_PUMP_EFF*K_BOLTZ*DELTA_E_EFF) &
                           & /(1.0_dp+(A_COEFF_10+R_PHOTO_EFF)/(C_COEFF_10*gasDensity)))
    end function getH2VibrationalCooling

    ! !-----------------------------------------------------------------------
    ! !  H2 photodissociation heating
    ! !
    ! !  On average, 0.4 eV of kinetic energy per photodissociated molecule
    ! !
    ! !  Use the H2 photodissociation rate determined by the subroutine
    ! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRH2) (s^-1)
    ! !  JH: again, grabbed h2dis from chemistry.f90 for consistency.
    ! !-----------------------------------------------------------------------
    function getH2PhotoDissHeating(gasDensity,h2Abund,h2dis) result(H2PhotoDissHeating)
        real(dp), intent(in) :: gasDensity,h2Abund,h2dis
        real(dp) :: H2PhotoDissHeating
        H2PhotoDissHeating=(0.4_dp*eV)*h2dis*gasDensity*h2Abund
    end function getH2PhotoDissHeating

    ! !-----------------------------------------------------------------------
    ! !  Cosmic-ray ionization heating
    ! !
    ! !  20.0 eV of kinetic energy deposited per H2 ionization,
    ! !  based on the estimate of Goldsmith (2001, ApJ, 557,736)
    ! !
    ! !  See also:
    ! !  Clavel et al. (1978, A&A, 65, 435)
    ! !  Tielens & Hollenbach (1985, ApJ, 291, 722)
    ! !  Shull & van Steenberg (1985, ApJ, 298, 268)
    ! !  Kamp & van Zadelhoff (2001)
    ! !-----------------------------------------------------------------------
    function getCosmicRayHeating(zeta,gasDensity,h2Abund) result(cosmicRayHeating)
        real(dp), intent(in) :: zeta,gasDensity,h2Abund
        real(dp) :: cosmicRayHeating
        ! cosmicRayHeating=(20.0*eV)*(1.3e-17_dp*zeta)*h2Abund*gasDensity
        ! According to Ivlev et al. 2019
        cosmicRayHeating=(16.0_dp*eV)*(1.3e-17_dp*zeta)*h2Abund*gasDensity
    end function getCosmicRayHeating



    ! !-----------------------------------------------------------------------
    ! !  H2 FUV pumping heating
    ! !
    ! !  On average, 2.2 eV released per vibrationally excited H2* molecule
    ! !
    ! !  Use the treatment of Hollenbach & McKee (1979, ApJ)

    ! !  Use the H2 photodissociation rate determined by the subroutine
    ! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRH2) (s^-1)
    ! !
    ! !  Use the H2 critical density expression from Hollenbach & McKee (1979)
    ! !  NOTE: The equation for the collisional de-excitation rate coefficient
    ! !  for the v=2-1 and v=1-0 transitions by collisions with H2 was wrongly
    ! !  stated in the Hollenbach & McKee (1979) paper (equation 6.29), but is
    ! !  corrected in Hollenbach & McKee (1989, equation 2.8) and used below.

    ! ! JH: h2dis instead of Rate(nRH2) again
    ! !-----------------------------------------------------------------------
    function getH2FUVPumpHeating(hAbund,h2Abund,gasTemperature,gasDensity,h2dis) &
            result(H2FUVPumpHeating)
        real(dp), intent(in) :: hAbund,h2Abund,gasTemperature,gasDensity,h2dis
        real(dp) :: H2FUVPumpHeating

        real(dp) :: NCRIT_H2
        NCRIT_H2=1.0e6_dp/sqrt(gasTemperature)/(1.6_dp*hAbund*exp(-((400.0_dp/gasTemperature)**2)) &
                                      & + 1.4_dp*h2Abund*exp(-(18100.0_dp/(gasTemperature+1200.0_dp))))

        H2FUVPumpHeating=(2.2_dp*eV)*9.0_dp*h2dis*gasDensity*h2Abund/(1.0_dp+NCRIT_H2/gasDensity)
        ! !  If vibrationally excited H2 (H2*) is included in the chemical network,
        ! !  then use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722)
        !    IF(nH2v.NE.0) THEN
        !       H2_FUV_PUMPING_HEATING_RATE=(DENSITY(nH)*1.0e-12_dp*sqrt(gasTemperature)*exp(-1000.0_dp/gasTemperature) &
        !                                & +DENSITY(nH2)*1.4e-12_dp*sqrt(gasTemperature)*exp(-18100.0_dp/(gasTemperature+1200.0_dp))) &
        !                                  *(2.6*eV)*DENSITY(nH2v)
        !    END IF

    end function getH2FUVPumpHeating



! !-----------------------------------------------------------------------
! !  Carbon photoionization heating
! !
! !  On average, 1 eV of kinetic energy released per carbon ionization
! !  Use the carbon photoionization rate determined by the subroutine
! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRCI) (s^-1)
! !-----------------------------------------------------------------------

function getCarbonIonizationHeating(cIonizationRate,carbonAbund,gasDensity) &
        result(carbonIonizationHeating)
    real(dp), intent(in) :: cIonizationRate,carbonAbund,gasDensity
    real(dp) :: carbonIonizationHeating
    carbonIonizationHeating=(1.0_dp*eV)*cIonizationRate*carbonAbund*gasDensity
end function getCarbonIonizationHeating

! !-----------------------------------------------------------------------
! !  Exothermic chemical reaction heating
! !
! !  See:
! !  Clavel et al. (1978, A&A,  65, 435)
! !  Meijerink & Spaans (2005, A&A, 436, 397)
! !  Glassgold & Langer (1973, ApJ, 179, L147)
! !
! !  Recombination reactions:
! !     H2+ (10.9 eV); H3+ (9.23+4.76 eV); H3O+ (1.16+5.63+6.27 eV); HCO+ (7.51 eV)
! !
! !  Ion-neutral reactions:
! !     H2+ + H (0.94 eV); He+ + H2 (6.51 eV); He+ + CO (2.22 eV)
! !
! !  For each reaction, the heating rate is given by: n(1) * n(2) * K * E
! !  where n(1) and n(2) are the number densities, K the rate coefficient
! !  (cm^3 s^-1), and E the energy released (erg).
! !-----------------------------------------------------------------------

function getChemicalHeatingRate(gasDensity,exoReactants1,exoReactants2,exoRates,exothermicities) &
        result(chemicalHeatingRate)
    real(dp), intent(in) :: gasDensity
    real(dp), intent(in), dimension(:) :: exoReactants1,exoReactants2,exoRates,exothermicities
    real(dp) :: chemicalHeatingRate

    chemicalHeatingRate=SUM(exoReactants1*exoReactants2*exoRates*exothermicities)
    chemicalHeatingRate=chemicalHeatingRate*gasDensity*gasDensity*eV  !each abundance should be a number density to multiply through
end function getChemicalHeatingRate
! !-----------------------------------------------------------------------
! !  Gas-grain collisional heating
! !
! !  Use the treatment of Burke & Hollenbach (1983, ApJ, 265, 223)
! !
! !  Other relevant references:
! !  Hollenbach & McKee (1979, ApJS, 41, 555)
! !  Tielens & Hollenbach (1985, ApJ, 291, 722)
! !  Goldsmith (2001, ApJ, 557, 736)
! !  Young et al. (2004, ApJ, 614, 252)
! !
! !  This process is insignificant for the energy balance of the dust,
! !  but can influence the gas temperature. If the dust temperature is
! !  lower than the gas temperature, this becomes a cooling mechanism.
! !-----------------------------------------------------------------------

function getGasGrainCollisionsHeatingRate(gasTemperature,gasDensity,dustAbundance,dustRadius,dustTemp) &
        result(gasGrainCollisionsHeatingRate)
    real(dp), intent(in) :: gasTemperature,gasDensity,dustAbundance,dustRadius,dustTemp
    real(dp) :: gasGrainCollisionsHeatingRate
    real(dp) :: nGrain,accommodation,C_GRAIN
    nGrain=dustAbundance*gasDensity

    !nGrain=2.0e-12_dp*gasDensity
    C_GRAIN=PI*dustRadius**2

    !!$!  Accommodation fitting formula of Groenewegen (1994, A&A, 290, 531)
    !!$   ACCOMMODATION=0.35_dp*exp(-sqrt((dustTemperature+gasTemperature)/5.0e2_dp))+0.1_dp

    !  Accommodation coefficient of Burke & Hollenbach (1983, ApJ, 265, 223)
    accommodation=0.37_dp*(1.0_dp-0.8_dp*exp(-75.0_dp/gasTemperature))

    gasGrainCollisionsHeatingRate=nGrain*C_GRAIN*gasDensity*sqrt(8.0_dp*K_BOLTZ*gasTemperature/(PI*MH)) &
                       & *accommodation*(2.0_dp*K_BOLTZ*dustTemp-2.0_dp*K_BOLTZ*gasTemperature)
end function getGasGrainCollisionsHeatingRate

function calculateDustTemp(localField,surfaceField,Av,zeta) result(dustTemperature)
    real(dp), intent(in) :: localField,surfaceField,Av
    real(dp), intent(in), optional :: zeta
    real(dp) :: dustTemperature

    !Choose which dust temperature calculation to use
    select case (dust_gas_coupling_method)
    case (1)
        dustTemperature=calculateDustTempHollenbach(localField,surfaceField)
    case (2)
        dustTemperature=calculateDustTempHocuk(surfaceField,Av)
    case (3)
        ! Ivlev et al. 2019: CR heating correction on top of Hocuk base temperature
        dustTemperature=calculateDustTempHocuk(surfaceField,Av)
        if (present(zeta)) then
            dustTemperature=calculateDustTempIvlev(dustTemperature,zeta)
        end if
    case default
        write(*,*) "Unimplemented dust temperature calculation method choose 1, 2, or 3. Exiting."
        stop
    end select
end function calculateDustTemp

!=======================================================================
!
!  Calculate the dust temperature for each particle using the treatment
!  of Hollenbach, Takahashi & Tielens (1991, ApJ, 377, 192, eqns 5 & 6)
!  for the heating due to the incident FUV photons and the treatment of
!  Meijerink & Spaans (2005, A&A, 436, 397, eqn B.6) for heating due to
!  the incident flux of X-ray photons.
!
!  Among other things, the dust temperature can influence:
!
!     1) Cooling budget by emitting FIR photons that
!        interact with the line radiative transfer;
!     2) Gas-grain collisional heating or cooling rate;
!     3) H2 formation by changing the sticking probability;
!     4) Evaporation and condensation of molecules on grains.
!
!  The formula derived by Hollenbach, Takahashi & Tielens (1991) has
!  been modified to include the attenuation of the IR radiation. The
!  incident FUV radiation is absorbed and re-emitted in the infrared
!  by dust at the surface of the cloud (up to Av ~ 1mag). In the HTT
!  derivation, this IR radiation then serves as a second heat source
!  for dust deeper into the cloud. However, in their treatment, this
!  second re-radiated component is not attenuated with distance into
!  the cloud so it is *undiluted* with depth, leading to higher dust
!  temperatures deep within the cloud which in turn heat the gas via
!  collisions to unrealistically high temperatures. Models with high
!  gas densities and high incident FUV fluxes (e.g. n_H = 10^5 cm-3,
!  X_0 = 10^8 Draine) can produce T_gas ~ 100 K at Av ~ 50 mag!
!
!  Attenuation of the FIR radiation has therefore been introduced by
!  using an approximation for the infrared-only dust temperature from
!  Rowan-Robinson (1980, eqn 30b):
!
!  T_dust = T_0*(r/r_0)^(-0.4)
!
!  where r_0 is the cloud depth at which T_dust = T_0, corresponding
!  to an A_V of ~ 1 mag, the assumed size of the outer region of the
!  cloud that processes the incident FUV radiation and then re-emits
!  it in the FIR (see the original HTT 1991 paper for details). This
!  should prevent the dust temperature from dropping off too rapidly
!  with distance and maintain a larger warm dust region (~50-100 K).
!
!-----------------------------------------------------------------------

function calculateDustTempHollenbach(localField,surfaceField) result(dustTemperature)
    !UV field in Habing at this depth and at surface required for this calculation
    !Both in Habing as required for this treatment
    real(dp), intent(in) :: localField,surfaceField
    real(dp) :: dustTemperature

    real(KIND=DP) :: NU_0,R_0,T_0,TAU_100

    !Parameters used in the HHT equations (see their paper for details)
    NU_0=2.65e15_dp
    TAU_100=1.0e-3_dp
    R_0=1.0_dp/1.6e-21_dp  !avFactor

    !Calculate the contribution to the dust temperature from the local FUV flux and the CMB background
    !UCLPDR had afactor of 1.7 which I assume was Habing conversion so removed
    dustTemperature=8.9e-11_dp*NU_0*localField+T_CMB**5


    !The minimum dust temperature is related to the incident FUV flux along each ray
    T_0=12.2_dp*surfaceField**0.2_dp

    !!Attenuate the FIR radiation produced in the surface layer
    !!JH why is this commented?
    !IF(PARTICLE(P)%TOTAL_COLUMN(J).GT.R_0) THEN
    !     T_0=T_0*(PARTICLE(P)%TOTAL_COLUMN(J)/R_0)**(-0.4)
    !END IF

    !        Add the contribution to the dust temperature from the FUV flux incident along this ray
    if(T_0>0) dustTemperature=dustTemperature &
          & + (0.42_dp-log(3.45e-2_dp*TAU_100*T_0))*(3.45e-2_dp*TAU_100*T_0)*T_0**5

    !Convert from total dust emission intensity to dust temperature
    dustTemperature=dustTemperature**0.2_dp

    !Calculate the contribution to the dust temperature from the local X-ray flux (assuming a fixed grain abundance of 1.6E-8)
    !JH We have no xrays sthis
    !dustTemperature=dustTemperature+1.5e4_dp*(XRAY_ENERGY_DEPOSITION_RATE/1.6e-8_dp)**0.2

    !Impose a lower limit on the dust temperature, since values below 10 K can dramatically
    !limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
    if(dustTemperature<lower_limit_dusttemp) then
        dustTemperature=lower_limit_dusttemp
    end if

    !     Check that the dust temperature is physical
    if(dustTemperature>upper_limit_dusttemp) then
        write(*,*) localField, surfaceField  !WRITE(6,*) 'ERROR! Calculated dust temperature exceeds upper limit'
        dustTemperature=upper_limit_dusttemp
    end if
end function calculateDustTempHollenbach

! Using the new parametric formulation for dust temperature in Hocuk et al. 2017
! See equation 8 in DOI: 10.1051/0004-6361/201629944
function calculateDustTempHocuk(surfaceField,Av) result(dustTemperature)
    real(dp), intent(in) :: surfaceField, Av
    real(dp) :: tanh_term
    real(dp) :: dustTemperature
    tanh_term = 0.61_dp - log10(Av)
    dustTemperature = (11.0_dp + 5.7_dp*tanh(tanh_term)) * (1.7_dp*surfaceField)**(1.0_dp/5.9_dp)

    !Impose a lower limit on the dust temperature, since values below 10 K can dramatically
    !limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
    if(dustTemperature<lower_limit_dusttemp) then
        dustTemperature=lower_limit_dusttemp
    end if

    !     Check that the dust temperature is physical
    if(dustTemperature>upper_limit_dusttemp) then
        write(*,*) Av, surfaceField  !WRITE(6,*) 'ERROR! Calculated dust temperature exceeds upper limit'
        dustTemperature=upper_limit_dusttemp
    end if
end function calculateDustTempHocuk

! Ivlev et al. 2019 cosmic-ray dust heating correction
! Equation 30: T_d,eff = T__dp * [1 + 0.202 * (zeta_ion/1e-16)*(T__dp/6)^(-6)]^(1/6)
! where T__dp is the base dust temperature (e.g. from Hocuk 2017)
! and zeta_ion is the cosmic ray ionization rate in s^-1 (= 1.3e-17_dp * zeta)
function calculateDustTempIvlev(T_dp, zeta) result(dustTemperature)
    real(dp), intent(in) :: T_dp   ! Base dust temperature from Hocuk
    real(dp), intent(in) :: zeta  ! Dimensionless CR ionization rate scaling factor
    real(dp) :: dustTemperature
    real(dp) :: zeta_ion  ! CR ionization rate in s^-1

    zeta_ion = 1.3e-17_dp * zeta
    dustTemperature = T_dp * (1.0_dp + 0.202_dp * (zeta_ion / 1.0e-16_dp) * (T_dp / 6.0_dp)**(-6))**(1.0_dp/6.0_dp)

    ! Apply dust temperature limits
    if (dustTemperature<lower_limit_dusttemp) then
        dustTemperature = lower_limit_dusttemp
    end if
    if (dustTemperature>upper_limit_dusttemp) then
        dustTemperature = upper_limit_dusttemp
    end if
end function calculateDustTempIvlev
!=======================================================================



! !-----------------------------------------------------------------------
! !  Coulomb heating
! !
! !  Use the treatment of Meijerink & Spaans (2005, A&A, 436, 397)
! !
! !  Other relevant references:
! !  Dalgarno et al. (1999, ApJS, 125, 237)
! !
! !  This is an X-ray heating mechanism. When X-rays are absorbed, fast
! !  electrons are produced. These fast electrons lose part of their
! !  energy through Coulomb interactions with thermal electrons.
! !-----------------------------------------------------------------------

!    R=DENSITY(nH2)/DENSITY(nH) ! n(H2)/n(H) ratio
!    X_PRIME=1.83_dp*ABUNDANCE(nelect)/(1.0_dp+0.83_dp*ABUNDANCE(nelect)) ! Correction to the electron abundance for a pure H2-He mixture

!    ETA_H2_He=1.0_dp+(0.055_dp-1.0_dp)/(1.0_dp+2.17_dp*X_PRIME**0.366) ! Heating efficiency for a pure H2-He mixture
!    ETA_H_He =1.0_dp+(0.117_dp-1.0_dp)/(1.0_dp+7.95_dp*ABUNDANCE(nelect)**0.678) ! Heating efficiency for a pure H-He mixture
!    ETA=(10.0_dp*R*ETA_H2_He+ETA_H_He)/(10.0_dp*R+1.0_dp) ! Total heating efficiency for mixed atomic and molecular gas
!    H_X=XRAY_ENERGY_DEPOSITION_RATE ! X-ray energy deposition rate per hydrogen nucleus (erg s^-1)

!    COULOMB_HEATING_RATE=ETA*GAS_DENSITY*H_X



! !-----------------------------------------------------------------------
! !  Supersonic turbulent decay heating
! !
! !  Most relevant for the inner parsecs of galaxies (Black)
! !  Black, in Interstellar Processes, 1987, p731
! !  See also: Rodriguez-Fernandez et al., 2001, A&A, 365, 174
! !
! !  V_TURB = turbulent velocity (km/s); Galactic center ~ 15 km/s
! !  L_TURB = turbulent scale length (pc); typically 5 pc
! !-----------------------------------------------------------------------
function getTurbulentHeatingRate(gasDensity,V_TURB) result(turbulentHeatingRate)
    real(dp), intent(in) :: gasDensity,V_TURB
    real(dp) :: turbulentHeatingRate
    real(dp) :: L_TURB

   L_TURB=5.0_dp
   turbulentHeatingRate=3.5e-28_dp*((V_TURB/1.0e5_dp)**3)*(1.0_dp/L_TURB)*gasDensity

end function getTurbulentHeatingRate



subroutine pair_insertion_sort(array)
    real(dp), intent(inout) :: array(:)
    integer :: i,j,last
    real(dp) :: t1,t2

    last=size(array)
    do i=2,last-1,2
       t1=min(array(i),array(i+1))
       t2=max(array(i),array(i+1))
       j=i-1
       do while((j>=1).and.(array(j)>t2))
          array(j+2)=array(j)
          j=j-1
       end do
       array(j+2)=t2
       do while((j>=1).and.(array(j)>t1))
          array(j+1)=array(j)
          j=j-1
       end do
       array(j+1)=t1
    end do

    if(mod(last,2)==0)then
       t1=array(last)
       loop: do j=last-1,1,-1
          if (array(j)<=t1) exit loop
          array(j+1)=array(j)
       end do loop
       array(j+1)=t1
    end if
end subroutine pair_insertion_sort

!-----------------------------------------------------------------------
! pair_insertion_sort_with_perm - Sort a 1D array and track permutation
!
! Sorts a 1D array using pair insertion sort while maintaining a
! permutation array that tracks where each element originally came from.
! This allows you to reorder other arrays using the same permutation.
!
! Arguments:
!   array(:)    - 1D array to sort (modified in place)
!-----------------------------------------------------------------------
subroutine pair_insertion_sort_with_perm(array, perm)
    real(dp), intent(inout) :: array(:)
    integer, intent(out) :: perm(:)
    integer :: i, j, last
    real(dp) :: t1, t2
    integer :: p1, p2

    last = SIZE(array)

    ! Initialize permutation array (caller must allocate with correct size)
    perm = [(i, i=1, last)]

    ! Pair insertion sort - process elements two at a time
    do i = 2, last-1, 2
       ! Get the pair and their permutation indices
       if (array(i) <= array(i+1)) then
          t1 = array(i)
          t2 = array(i+1)
          p1 = perm(i)
          p2 = perm(i+1)
       else
          t1 = array(i+1)
          t2 = array(i)
          p1 = perm(i+1)
          p2 = perm(i)
       end if

       ! Find position for larger element (t2)
       j = i - 1
       do while ((j >= 1) .AND. (array(j) > t2))
          array(j+2) = array(j)
          perm(j+2) = perm(j)
          j = j - 1
       end do
       array(j+2) = t2
       perm(j+2) = p2

       ! Find position for smaller element (t1)
       do while ((j >= 1) .AND. (array(j) > t1))
          array(j+1) = array(j)
          perm(j+1) = perm(j)
          j = j - 1
       end do
       array(j+1) = t1
       perm(j+1) = p1
    end do

    ! Handle last element if array has even number of elements
    if (MOD(last, 2) == 0) then
       t1 = array(last)
       p1 = perm(last)
       loop: do j = last-1, 1, -1
          if (array(j) <= t1) exit loop
          array(j+1) = array(j)
          perm(j+1) = perm(j)
       end do loop
       array(j+1) = t1
       perm(j+1) = p1
    end if
end subroutine pair_insertion_sort_with_perm

end module heating

!Abandoned heating and cooling mechanisms

! ! !-----------------------------------------------------------------------
! ! !  Molecular Hydrogen cooling
! ! !  from Galli & Palla (1998) via Grassi et al. (2014)
! ! !  Acceptable up to 10^5 K can easily use  Glover & Abel (2008) fits instead
! ! !  but they're more complex.
! ! ! JH: Dropped in favor of UCLPDR treatment of H2 rotational cooling + h2 Vibrational
! ! !-----------------------------------------------------------------------
! REAL(dp) FUNCTION h2Cooling(gasDensity,gasTemperature,hAbund,h2Abund)
!     REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,hAbund,h2Abund
!     REAL(dp) :: highDensLimit,lowDensLimit,T1000,logT
!     ! Temperature in kilokelvin
!     T1000=gasTemperature*0.001

!     !high density limit is same in all models: sum of vibrational and rotational cooling
!     highDensLimit=((9.5e-22_dp*T1000**3.76)/(1+0.12*(T1000**2.1)))*exp(-(0.13/T1000)**3.0)
!     highDensLimit=highDensLimit+(3.0e-24_dp*exp(-0.51/T1000))

!     highDensLimit=highDensLimit+(6.7e-19_dp*exp(-5.86/T1000))+(1.6e-18_dp*exp(-11.7/T1000))

!     !I'm using Galli & Palli limit here which is ok up to 10^5.
!     !I think Glover and Abel is more accurate but is many fits so hard to imlpement
!     logT=log10(gasTemperature)
!     lowDensLimit=-103.0+(97.59*logT)-(48.05*logT*logT)+(10.8*logT**3.0)-(0.9032*logT**4.0)
!     lowDensLimit=10.0**(lowDensLimit*gasDensity*hAbund)

!     !Combine together for final cooling rate
!     h2Cooling=gasDensity*h2Abund*highDensLimit
!     h2Cooling=h2Cooling/(1.0+(highDensLimit/lowDensLimit))
! END FUNCTION h2Cooling




!=========================================
!     JH: Alternative photoelectric heating rates. The correct one depends on dust distribution
!     The one coded as photoelectricHeating() is the current one in ucl in UCLPDR
!-----------------------------------------------------------------------
!  Grain photoelectric heating (large grains only; r approx 100 Angstrom)
!
!  Use the treatment of

!  which in turn follows de Jong (1977, 1980)
!
!  The charge of a dust grain can be found by equating the rate of
!  photo-ejection of electrons from the dust grain to the rate of
!  recombination of electrons with the dust grain (Spitzer)
!
!  The various parameter values are taken from Table 2 of the paper
!-----------------------------------------------------------------------
!     REAL(dp) FUNCTION photoelectricHeating(gasTemperature,gasDensity,habingField,electronAbund)
!         REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,habingField,electronAbund
!         REAL(dp), PARAMETER :: DELTA_D=1.0_dp
!         REAL(dp), PARAMETER :: DELTA_UV=1.8_dp
!         REAL(dp), PARAMETER :: Y=0.1_dp
!         REAL(dp), PARAMETER :: HNU_D=6.0_dp
!         REAL(dp), PARAMETER :: HNU_H=13.6_dp
!         REAL(dp) :: delta,gamma,XK,XD,X,XX
!         INTEGER :: ITERATION


!         XK=K_BOLTZ*gasTemperature/(HNU_H*eV)
!         XD=HNU_D/HNU_H
!         gamma=2.9e-4_dp*Y*sqrt(gasTemperature)*habingField/(gasDensity*electronAbund)
!         delta=XK-XD+gamma

!         !  Iterate to determine X by finding the zero of the function F
!         X=0.5_dp
!         DO ITERATION=1,100
!           XX=X-(grainChargeFunc(X,DELTA,GAMMA)/deltaGrainChargeFunc(X,DELTA))
!           IF(ABS(XX-X).LT.1.0e-2_dp) EXIT
!           X=XX
!         END DO
!         X=XX

!         IF(ITERATION.GE.100) THEN
!           WRITE(10,*)'WARNING! Grain parameter X not found in PE heating'
!           WRITE(10,*)'Using final value from iteration loop: X =',X
!         END IF

!         photoelectricHeating=2.7e-25_dp*DELTA_UV*DELTA_D*gasDensity*Y*habingField &
!                                      & *(((1.0_dp-X)**2)/X + XK*((X**2)-1.0_dp)/(X**2))

!         !  Assume the PE heating rate scales linearly with metallicity
!         !TH85_PHOTOELECTRIC_HEATING_RATE=TH85_PHOTOELECTRIC_HEATING_RATE*METALLICITY
!     END FUNCTION  photoelectricHeating

! !=======================================================================
! !  X is the grain charge parameter and is the solution to F(X)=0
! !-----------------------------------------------------------------------
!    FUNCTION grainChargeFunc(X,DELTA,GAMMA)
!       IMPLICIT NONE
!       REAL(dp) :: grainChargeFunc
!       REAL(dp), INTENT(IN) :: X,DELTA,GAMMA
!       grainChargeFunc=(X**3)+DELTA*(X**2)-GAMMA
!    END FUNCTION grainChargeFunc

! !=======================================================================
! !  FF(X) is the derivative of F(X) with respect to X
! !-----------------------------------------------------------------------
!    FUNCTION deltaGrainChargeFunc(X,DELTA)
!       IMPLICIT NONE
!       REAL(dp) :: deltaGrainChargeFunc
!       REAL(dp), INTENT(IN) :: X,DELTA
!       deltaGrainChargeFunc=3*(X**2)+DELTA*(2*X)
!    END FUNCTION deltaGrainChargeFunc

! !-----------------------------------------------------------------------
! !  Grain + PAH photoelectric heating (MRN size distribution; r = 3-100 Angstrom)
! !
! !  Use the treatment of Bakes & Tielens (1994, ApJ, 427, 822) with the
! !  modifications suggested by Wolfire et al. (2003, ApJ, 587, 278) to
! !  account for the revised PAH abundance estimate from Spitzer data.
! !
! !  See also:
! !  Wolfire et al. (1995, ApJ, 443, 152)
! !  Le Page, Snow & Bierbaum (2001, ApJS, 132, 233)
! !-----------------------------------------------------------------------

! !  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680, 384)
! !  Setting this factor to 1.0 gives the standard Bakes & Tielens expression
!    PHI_PAH=0.4_dp

!    ALPHA=0.944_dp
!    BETA=0.735_dp/gasTemperature**0.068
!    DELTA=HABING_FIELD*sqrt(gasTemperature)/(DENSITY(nelect)*PHI_PAH)
!    EPSILON=4.87e-2_dp/(1.0_dp+4.0e-3_dp*DELTA**0.73) + 3.65e-2_dp*(gasTemperature/1.0e4_dp)**0.7/(1.0_dp+2.0e-4_dp*DELTA)

!    PAH_HEATING_RATE=1.30e-24_dp*EPSILON*HABING_FIELD*GAS_DENSITY
!    PAH_COOLING_RATE=4.65e-30_dp*gasTemperature**ALPHA*(DELTA**BETA)*DENSITY(nelect)*PHI_PAH*GAS_DENSITY

!    BT94_PHOTOELECTRIC_HEATING_RATE=PAH_HEATING_RATE - PAH_COOLING_RATE

! !  Assume the PE heating rate scales linearly with metallicity
!    BT94_PHOTOELECTRIC_HEATING_RATE=BT94_PHOTOELECTRIC_HEATING_RATE*METALLICITY

 ! REAL(dp) FUNCTION getEquilibriumTemp(gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate,dustAbundance&
 !                                &,exoReactants1,exoReactants2,exoRates,exothermicities,writeFlag,dustTemp,turbVel,fixedCooling,coolingFlag&
 !                                ,fixedHeating,heatingFixFlag)
 !        !Habing field is radfield diminished by Av
 !        REAL(dp), INTENT(in) :: gasTemperature,gasDensity,habingField,h2dis,h2form,zeta,cIonRate,dustAbundance,dustTemp,turbVel,fixedCooling,fixedHeating
 !        REAL(dp), INTENT(in) :: abundances(:),exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
 !        LOGICAL, INTENT(IN) :: writeFlag,coolingFlag,heatingFixFlag
 !        REAL(dp) :: previousTemp,previousDifference,thigh,tlow
 !        LOGICAL :: binaryChopSearch,BRACKET_expANDED,TEMPERATURE_CONVERGED
 !        REAL(dp) :: heating,cooling,temperatureDiff,difference,outTemp,relative_difference
 !        REAL(dp),parameter :: TDIFF=0.01, FCRIT=0.1,TMIN=10.0, TMAX=1.0e5_dp
 !        INTEGER :: tempLoops = 0.0

 !        previousTemp=0.0
 !        previousDifference=0.0
 !        thigh=TMAX
 !        tlow=TMIN

 !        binaryChopSearch=.False.
 !        BRACKET_expANDED=.False.
 !        TEMPERATURE_CONVERGED=.False.
 !        tempLoops=0
 !        getEquilibriumTemp=gasTemperature
 !        DO WHILE (.NOT. TEMPERATURE_CONVERGED .AND. tempLoops .lt. 100)
 !            tempLoops=tempLoops+1

 !            !absolute difference between current temperature and previous
 !            temperatureDiff=ABS(getEquilibriumTemp-previousTemp)

 !            !then calculate overall heating/cooling rate
 !            heating=fixedHeating
 !            if (heatingFixFlag) heating=getHeatingRate(getEquilibriumTemp,gasDensity,habingField,abundances,h2dis,h2form,zeta,&
 !              &cIonRate,dustAbundance,exoReactants1,exoReactants2,exoRates,exothermicities,dustTemp,turbVel)
 !            !write(*,*) "Total Heating", heating

 !            !set cooling rate to fixed cooling rate then overwrite if we want real cooling
 !            cooling=fixedCooling
 !            if (coolingFlag) cooling=getCoolingRate(getEquilibriumTemp,gasDensity,dustTemp,abundances,h2dis,turbVel,writeFlag)

 !            !Calculate the difference between the total heating and total cooling rates
 !            !and the absolute value of the relative difference between the two rates
 !            difference=heating-cooling
 !            relative_difference=2.0_dp*ABS(difference)/ABS(heating+cooling)

 !            ! !Quick fix to get fixed T whilst calculating cooling
 !            ! TEMPERATURE_CONVERGED=.TRUE.
 !            ! EXIT
 !            !Check if we've converged heating/cooling balanace
 !            IF (relative_difference .lt. FCRIT) THEN
 !                previousTemp=getEquilibriumTemp
 !                TEMPERATURE_CONVERGED=.True.
 !                CYCLE
 !            END IF


 !        !  Determine the temperature bracket to begin searching within by first increasing
 !        !  or decreasing the temperature by 30% according to the heating-cooling imbalance
 !            IF(.NOT. binaryChopSearch) THEN
 !                !If the heating continues to outweigh the cooling, increase the temperature by 30%
 !                IF(DIFFERENCE.GT.0 .AND. previousDifference.GE.0) THEN
 !                    TLOW=getEquilibriumTemp ! Update the value of T_low
 !                    getEquilibriumTemp=1.3_dp*getEquilibriumTemp
 !                    previousDifference=difference
 !                    THIGH=getEquilibriumTemp ! Update the value of T_high
 !        !     If the cooling continues to outweigh the heating, decrease the temperature by 30%
 !              ELSE IF(DIFFERENCE.LT.0 .AND. previousDifference.LE.0) THEN
 !                 THIGH=getEquilibriumTemp ! Update the value of T_high
 !                 getEquilibriumTemp=0.7_dp*getEquilibriumTemp
 !                 previousDifference=DIFFERENCE
 !                 TLOW=getEquilibriumTemp ! Update the value of T_low


 !        !     If the heating-cooling balance has reversed (either from net heating to net cooling or
 !        !     vice-versa) then switch to the binary chop search method to determine the temperature
 !              ELSE
 !                 getEquilibriumTemp=(THIGH+TLOW)/2.0_dp
 !                 previousDifference=DIFFERENCE
 !                 binaryChopSearch=.TRUE. ! From now on
 !              END IF

 !        !  Perform a binary chop search (the min-max range was found by the 30% increase/decrease method)
 !           ELSE

 !              IF(DIFFERENCE.GT.0) THEN
 !                TLOW=getEquilibriumTemp ! Update the value of T_low
 !                 getEquilibriumTemp=(getEquilibriumTemp+THIGH)/2.0_dp
 !                 previousDifference=DIFFERENCE

 !              END IF
 !              IF(DIFFERENCE.LT.0) THEN
 !                 THIGH=getEquilibriumTemp !Update the value of T_high
 !                 getEquilibriumTemp=(getEquilibriumTemp+TLOW)/2.0_dp
 !                 previousDifference=DIFFERENCE
 !              END IF

 !           END IF

 !        !  If the search routine is unable to converge on a temperature that satisfies the thermal balance
 !        !  criterion, expand the min-max search bracket asymmetrically and begin to narrow the search again
 !        !  If the repeated search fails to converge once more, force convergence at the current temperature
 !           IF(temperatureDiff.LE.TDIFF) THEN
 !              IF(.NOT.BRACKET_expANDED) THEN
 !                 THIGH=THIGH+sqrt(PI)
 !                 TLOW=TLOW-sqrt(2.0)
 !                 BRACKET_expANDED=.TRUE.
 !              ELSE
 !                 previousTemp=getEquilibriumTemp
 !                 TEMPERATURE_CONVERGED=.TRUE.
 !                 CYCLE
 !              END IF
 !           END IF

 !        !  Check if the temperature falls outside of the allowed limits and force convergence if so
 !           IF(getEquilibriumTemp.LE.TMIN .AND. DIFFERENCE.LT.0) THEN
 !              getEquilibriumTemp=TMIN
 !              TEMPERATURE_CONVERGED=.TRUE.
 !           END IF
 !           IF(getEquilibriumTemp.GE.TMAX .AND. DIFFERENCE.GT.0) THEN
 !              getEquilibriumTemp=TMAX
 !              TEMPERATURE_CONVERGED=.TRUE.
 !           END IF

 !        !  Replace the previous temperature with the current value
 !           previousTemp=getEquilibriumTemp
 !        END DO

 !    END FUNCTION getEquilibriumTemp
