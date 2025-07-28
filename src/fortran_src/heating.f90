!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module that provides heating and cooling rates		  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE heating
    USE CONSTANTS
    USE F2PY_CONSTANTS, only : nspec
    USE DEFAULTPARAMETERS
    USE NETWORK
    USE COOLANT_MODULE
IMPLICIT NONE

    REAL(dp) :: pahAbund=6e-7

 CONTAINS

    SUBROUTINE initializeHeating(gasTemperature, gasDensity,abundances,columnDensity,cloudSize)
        REAL(dp), INTENT(in) :: gasTemperature,gasDensity,columnDensity,cloudSize
        REAL(dp), INTENT(in) :: abundances(:)
        CHARACTER(5) :: coolName
        INTEGER ::i,j

        ! write(*,*) "Initializing heating.f90 ..."
        CALL READ_COOLANTS()
        DO i=1,ncool
            coolName=coolantNames(i)
            if (coolName .eq. "p-H2") coolName="H2"
            if (coolName .eq. "o-H2") coolName="H2"
            if (coolName .eq. "p-H2O") coolName="H2O"
            if (coolName .eq. "o-H2O") coolName="H2O"
            DO j=1,nspec
                if (coolName .eq. specName(j)) coolantIndices(i)=j
            END DO
        END DO

        CLOUD_COLUMN=columnDensity
        CLOUD_DENSITY=gasDensity
        cloud_size=cloudSize
        ! Moved IO handling to io.f90
    END SUBROUTINE initializeHeating


    ! REAL(dp) FUNCTION getTempDot(gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate,dustAbundance,dustRadius,metallicity&
    !                             &,exoReactants1,exoReactants2,exoRates,exothermicities,,dustTemp,turbVel)
    REAL(dp) FUNCTION getTempDot(time,gasTemperature,gasDensity,gasCols,habingField,abundances,h2dis,h2form,zeta,cIonRate, &
                                & dustAbundance,dustRadius,metallicity, &
                                & dustTemp,turbVel)
        !Habing field is radfield diminished by Av
        REAL(dp), INTENT(in) :: time,gasTemperature,gasDensity,gasCols,habingField,h2dis,h2form,metallicity
        REAL(dp), INTENT(in) :: zeta,cIonRate,dustAbundance,dustRadius,dustTemp,turbVel
        REAL(dp), INTENT(in) :: abundances(:)!,exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)

        REAL(dp) adiabaticIdx,heating,cooling

        !First calculate adiabatic index - should use number density but that's just an additional common factor
        adiabaticIdx=5.0*(abundances(nh)+abundances(nhe)+abundances(nelec)+abundances(nh2))+2.0*abundances(nh2)
        adiabaticIdx=adiabaticIdx/(3.0*(abundances(nh)+abundances(nhe)+abundances(nelec)+abundances(nh2))+2.0*abundances(nh2))

        !then calculate overall heating/cooling rate
        ! heating=getHeatingRate(gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate,dustAbundance,dustRadius,metallicity&
        !     &,exoReactants1,exoReactants2,exoRates,exothermicities,dustTemp,turbVel,writeFlag)
        heating=getHeatingRate(time,gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate, &
                                & dustAbundance,dustRadius,metallicity, dustTemp,turbVel)

        IF (gasTemperature .gt. 3.0) THEN
            cooling=getCoolingRate(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,h2dis,turbVel)!6.951290d-17!!
        ELSE
            Cooling=0.0
        END IF

        getTempDot=heating-cooling
        !write(*,*) "Temp Dot",getTempDot
         !and convert to dT/dt
        getTempDot=((adiabaticIdx-1.0)*getTempDot)/(K_BOLTZ*gasDensity)
    END FUNCTION getTempDot


    ! REAL(dp) FUNCTION getHeatingRate(gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate,dustAbundance,dustRadius,metallicity&
    !                                 &,exoReactants1,exoReactants2,exoRates,exothermicities,dustTemp,turbVel,writeFlag)
    REAL(dp) FUNCTION getHeatingRate(time, gasTemperature,gasDensity,habingField,abundances,h2dis,h2form,zeta,cIonRate, & 
                                      & dustAbundance,dustRadius,metallicity,dustTemp,turbVel)
        REAL(dp), INTENT(in) :: time,gasTemperature,gasDensity,habingField,h2dis,metallicity
        REAL(dp), INTENT(in) :: h2form,zeta,cIonRate,dustAbundance,dustRadius,dustTemp,turbVel
        REAL(dp), INTENT(IN):: abundances(:)!,exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
        REAL(dp) :: turbHeating,L_TURB=5.0d0
        REAL(dp) :: photoelec,h2forming,fuvpumping,photodis,cionizing,crheating,chemheating,gasgraincolls

        !Photoelectric heating due to PAHs and large dust grains
        photoelec=photoelectricHeating(gasTemperature,gasDensity,habingField,abundances(nelec),metallicity)
        ! write(*,*) "    [heating.f90] photoelectric-heating=", photoelec

        !heating due to H2 formation
        h2forming=H2FormationHeating(gasTemperature,gasDensity,abundances(nh),h2form)
        ! write(*,*) "    [heating.f90] H2formation-heating=", h2forming

        !heating due to photodissociation of H2
        photodis=H2PhotodisHeating(gasDensity,abundances(nh2),h2dis)
        ! write(*,*) "    [heating.f90] H2dissocitation-heating=", photodis
        
        !heating due to FUV pumping of H2
        fuvpumping=h2FUVPumpHeating(abundances(nh),abundances(nh2),gasTemperature,gasDensity,h2dis)
        ! write(*,*) "    [heating.f90] FUVpumping-heating=", fuvpumping

        !heating due to carbon photo ionization          
        cionizing=CarbonIonizationHeating(cIonRate,abundances(nc),gasDensity)
        ! write(*,*) "    [heating.f90] carbon-photo-ionization-heating=", cionizing

        !CR dissociation of H2
        crheating=cosmicRayHeating(zeta,gasDensity,abundances(nh2))
        ! write(*,*) "    [heating.f90] CR-heating=", crheating

        !Exothermic reactions -- keep it for the moment because need to make the indices in newtork.f90 consistent
        ! chemheating=chemicalHeating(gasDensity,exoReactants1,exoReactants2,exoRates,exothermicities)
        !turbulent heating
        ! turbHeating=3.5D-28*((turbVel/1.0D5)**3)*(1.0D0/L_TURB)*gasDensity
        turbHeating=3.5D-28*(turbVel**3)*(1.0D0/L_TURB)*gasDensity  !check the unit of turbVel

        ! write(*,*) "    [heating.f90] turbulence-heating=", turbHeating
                
        !gas-grain collisional heaing/cooling
        gasgraincolls=gasGrainCollisions(gasTemperature,gasDensity,dustAbundance,dustRadius,dustTemp)
        ! write(*,*) "    [heating.f90] gas-grain-collisional-heating=", gasgraincolls

        !sum all heating types
        ! getHeatingRate=photoelec+h2forming+fuvpumping+photodis+cionizing+crheating+chemheating+turbHeating+gasgraincolls
        getHeatingRate=photoelec+h2forming+fuvpumping+photodis+cionizing+crheating+turbHeating+gasgraincolls
        !write to file

    END FUNCTION getHeatingRate

    REAL(dp) FUNCTION getCoolingRate(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,h2dis,turbVel)
        REAL(dp), INTENT(IN) :: time,gasTemperature,gasDensity,gasCols,dustTemp,h2dis,turbVel
        REAL(dp), INTENT(IN) :: abundances(:)
        real(dp) :: coolingMode, coolings(5)
        INTEGER :: ti


        coolingMode=atomicCooling(gasTemperature,gasDensity,abundances(nh),abundances(nhe),&
                        &abundances(nelec),abundances(nhx),abundances(nhex))
        getCoolingRate=coolingMode
        ! write(*,*) "    [heating.f90] atomic-cooling=", coolingMode

        coolingMode=collionallyInducedEmission(gasTemperature,gasDensity,abundances(nh2))
        getCoolingRate=getCoolingRate+coolingMode
        ! write(*,*) "    [heating.f90] collisional-induced-emisison-cooling=", coolingMode

        coolingMode=comptonCooling(gasTemperature,gasDensity,abundances(nelec))
        getCoolingRate=getCoolingRate+coolingMode
        ! write(*,*) "    [heating.f90] compton-cooling=", coolingMode

        coolingMode=continuumEmission(gasTemperature,gasDensity)
        getCoolingRate=getCoolingRate+coolingMode
        ! write(*,*) "    [heating.f90] continuum-emission-cooling=", coolingMode

        !This is already handed by heating FUV pumping function
        !coolingMode=H2VibrationalCooling(gasTemperature,gasDensity,abundances(nh2),h2dis)
        !getCoolingRate=getCoolingRate+coolingMode
        
        !We do the line cooling 5 times and take median value since solver will occasionally do something wild
        DO ti=1,5
            coolings(ti)=lineCooling(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,turbVel)
        END DO
        call pair_insertion_sort(coolings)

        coolingMode=coolings(3)
        getCoolingRate=coolingMode+getCoolingRate
        ! write(*,*) "    [heating.f90] line-cooling=", coolingMode
    END FUNCTION getCoolingRate


    REAL(dp) FUNCTION lineCooling(time,gasTemperature,gasDensity,gasCols,dustTemp,abundances,turbVel)
        REAL(dp), INTENT(IN) :: time,gasTemperature,gasDensity,gasCols,dustTemp,abundances(:),turbVel

        INTEGER :: N,I!, collisionalIndices(5)=(/nh,nhx,nh2,nhe,nelec/)
        real(dp) :: moleculeCooling(NCOOL)=0.0

        ! Update CLOUD_DENSITY and CLOUD_COLUMN
        CLOUD_DENSITY = gasDensity
        CLOUD_COLUMN  = gasCols

        CALL UPDATE_COOLANT_LINEWIDTHS(gasTemperature,turbVel)
        CALL UPDATE_COOLANT_ABUNDANCES(gasDensity,gasTemperature,abundances)

        DO N=1,NCOOL
            ! IF (TRIM(coolants(N)%NAME) == "C+") THEN
            !     write(*,*) "    DEBUG C+: DENSITY=", coolants(N)%DENSITY
            !     write(*,*) "    DEBUG C+: ENERGY(1:5)=", coolants(N)%ENERGY(1:MIN(5,coolants(N)%NLEVEL))
            !     write(*,*) "    DEBUG C+: WEIGHT(1:5)=", coolants(N)%WEIGHT(1:MIN(5,coolants(N)%NLEVEL))
            !     write(*,*) "    DEBUG C+: gasTemperature=", gasTemperature
            ! END IF

            CALL CALCULATE_LTE_POPULATIONS(coolants(N)%NLEVEL,coolants(N)%ENERGY,coolants(N)%WEIGHT, &
                                          & coolants(N)%DENSITY,gasTemperature, &
                                          & coolants(N)%POPULATION)
            ! IF (TRIM(coolants(N)%NAME) == "C+") THEN
            !     write(*,*) "    DEBUG C+: LTE population=", coolants(N)%POPULATION

            ! END IF

        END DO
        ! After LTE populations

        CALL CALCULATE_LINE_OPACITIES()
        CALL CALCULATE_LAMBDA_OPERATOR()

        !!write(*,*)  "lte done"
        !I should then do LVG interactions
         DO I=1,500!while not converged and less than 100 tries:
            DO N=1,NCOOL
                ! write(*,*) "-----------------------------------------------------------------"
                ! write(*,*) "  cloud_column: ", CLOUD_COLUMN
                ! write(*,*) "  cloud_density: ", CLOUD_DENSITY
                ! write(*,*) "  Coolant: ", N, TRIM(coolants(N)%NAME)
                ! write(*,*) "  DENSITY: ", coolants(N)%DENSITY
                ! write(*,*) "  ENERGY(1:5): ", coolants(N)%ENERGY(1:MIN(5,coolants(N)%NLEVEL))
                ! write(*,*) "  WEIGHT(1:5): ", coolants(N)%WEIGHT(1:MIN(5,coolants(N)%NLEVEL))
                ! write(*,*) "  gasTemperature: ", gasTemperature
                ! write(*,*) "  gasDensity: ", gasDensity
                ! write(*,*) "  abundance: ", abundances(coolantIndices(N))

                CALL CALCULATE_LEVEL_POPULATIONS(coolants(N),gasTemperature,gasDensity,&
                    &abundances,dustTemp)
                ! write(*,*) "LVG iter", I, "' Coolant: '", TRIM(coolants(N)%name), "' POPULATION(1:10): '", &
                !             coolants(N)%POPULATION(1:MIN(10, coolants(N)%NLEVEL))
            END DO
         !!write(*,*) "after Lvg",abundances(nh)

            !!write(*,*) I

            CALL CALCULATE_LINE_OPACITIES()
            CALL CALCULATE_LAMBDA_OPERATOR()
            IF (CHECK_CONVERGENCE()) EXIT
            !IF (I .eq. 499) write(*,*) "Failed convergence"
        END DO 

        !  Calculate the cooling rate due to the Lyman-alpha emission for each particle
        !  using the analytical expression of Spitzer (1978) neglecting photon trapping
        DO N=1,NCOOL
            if (.not. allocated(coolants(N)%EMISSIVITY)) then
                write(*,*) "ERROR: EMISSIVITY not allocated for coolant ", N
            end if

            IF(coolants(N)%NAME.EQ."H") THEN
               coolants(N)%EMISSIVITY(2,1) = 7.3D-19*(abundances(nelec)*gasDensity) &
                                                          & *(abundances(nH)*gasDensity) &
                                                          & *EXP(-118400.0D0/gasTemperature)
            END IF
        END DO

        !Calculate the cooling rates
        lineCooling=0.0
        DO N=1,NCOOL
            WHERE(coolants(N)%EMISSIVITY .lt. -HUGE(1.0)) coolants(N)%EMISSIVITY = 0.0
            moleculeCooling(N)=SUM(coolants(N)%EMISSIVITY,MASK=.NOT.ISNAN(coolants(N)%EMISSIVITY))
            ! !we get these wild changes in cooling rate so let's force it not to change too much in a timestep
            ! IF ( (abs(moleculeCooling-coolants(N)%previousCooling)/coolants(N)%previousCooling) .gt. 2.0d0) THEN
            !     !unless it's step 1, use old cooling for this time step
            !     IF (coolants(N)%previousCooling .ne. 0.0d0) moleculeCooling=coolants(N)%previousCooling
            ! ELSE
            !     coolants(N)%previousCooling=moleculeCooling
            ! END IF
            IF (moleculeCooling(N) .lt. 0.0) moleculeCooling(N)=0.0d0
            !IF (moleculeCooling(N) .gt. -1.0d-30 .and. abundances(coolantIndices(N)) .gt. 1.0d-20) 
            lineCooling= lineCooling+moleculeCooling(N)
        END DO


    END FUNCTION lineCooling

    ! !-----------------------------------------------------------------------
    ! !  Atomic and ionic cooling rates
    ! !  from Neal et al. 1995 based on Cen (1992) via Grassi et al. (2014)
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION atomicCooling(gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund)
        REAL(dp), INTENT(IN) :: gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund
        REAL(dp) :: t5,invT,rootT,collTFactor !temp/10^5, 1/T and a weird factor from the table
        REAL(dp) :: hDens,elecDens,heDens,hxDens,hexDens,gauntFactor
        hDens=gasDensity*hAbund
        elecDens=gasDensity*electronAbund
        heDens=gasDensity*heAbund
        hxDens=gasDensity*hxAbund
        hexDens=gasDensity*hexAbund
        t5=1.0d-5*gasT
        invT=1.0/gasT
        rootT=SQRT(gasT)
        collTFactor=1.0/(1.0+SQRT(t5))

        !gauntFactor from Neal et al. 1995
        gauntFactor=1.1+(0.34*EXP(-((5.5-LOG10(gasT))**2.0)/3.0))
        !Neal et al. 1995 lists several fits to cooling each in ergs/cm3/s so we'll just sum them
        !see table 1 of that paper
        !These are just numerical fits so there's loads of magic numbers
        !I've shorted variable names to make it easier to write/read (tn is temperature/10^n)

        !collisional excitation and ionization
        atomicCooling=(7.5d-19*collTFactor*EXP(-118348.0*invT)*elecDens*hDens) &
            &+(5.54d-17*(gasT**-0.397)*collTFactor*EXP(-473638.0*invT)*elecDens*hexDens)&
            &+(1.27d-21*rootT*EXP(-157809.1*invT)*elecDens*hDens*collTFactor)&
            &+(9.38d-22*rootT*EXP(-285335.4*invT)*elecDens*heDens*collTFactor)&
            &+(4.95d-22*rootT*EXP(-631515.0*invT)*elecDens*hexDens*collTFactor)&
            !dielectric
            &+(1.24d-13*(gasT**-1.5)*EXP(-470000.0*invT)*(1.0+0.3*EXP(-94000.0*invT))*elecDens*hexDens)
        IF (gasT .gt. 1.0d5) THEN
        !recombination
            atomicCooling=atomicCooling&
            &+(8.7d-27*rootT*((1.0d-3*gasT)**-0.2)*elecDens*hxDens/(1.0+((0.1*t5)**0.7)))&
            &+(1.55d-26*(gasT**0.3647)*elecDens*hexDens)&
            !&+(3.48d-26*rootT*((0.001*gasT)**-0.2)*nelec*nhexI/(1+(0.1*t5)**0.7))&
            !Free-free emission
            &+(1.42d-27*rootT*nelec*(nhex+nhx)*gauntFactor)
        END IF

    END FUNCTION atomicCooling

    !!-----------------------------------------------------------------------
    !! Collisionally Induced Emission
    !! Hirano & Yoshida (2013) and Ripamonti & Abel 2004 via Grassi 2012
    !!-----------------------------------------------------------------------
    REAL(dp) FUNCTION collionallyInducedEmission(gasTemperature,gasDensity,h2Abund)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,h2Abund
        REAL(dp), PARAMETER :: aConsts(6)=(/-30.3314216559651,19.0004016698518,-17.1507937874082&
                                                            &,9.49499574218739,-2.54768404538229,0.265382965410969/)
        REAL(dp), PARAMETER :: bConsts(6)=(/-180.992524120965,168.471004362887,-67.499549702687,&
                                                            &13.5075841245848,-1.31983368963974,0.0500087685129987/)
        REAL(dp),PARAMETER :: c=3.0,d=21.2968837223113
        REAL(dp) :: tau,logt
        INTEGER :: i        

        logt=LOG10(gasTemperature)

        tau=(gasDensity*h2Abund/7.0d15)**2.8
        !if (tau.lt.0.2) THEN
            !avoid numerical problems, tau tends to 1 for low tau but fortran can't do it
            !Taylor series is fine for tau<0.2 and that's well above the area we get an issue
         !   tau=1.0-(0.5*tau)+((tau**2.0)/6.0)-((tau**3.0)/24.0)+((tau**4)/120.0)
        !ELSE
            tau=(1.0-dexp(-tau))/tau
        !END IF
        tau=min(1.0,tau)


        collionallyInducedEmission=0.0

        IF (gasTemperature .ge. 1.0d5) THEN
            collionallyInducedEmission=(c*logt)-d            
        ELSE IF (gasTemperature .ge. 891.0) THEN
            DO i=1,SIZE(bConsts)
                collionallyInducedEmission=collionallyInducedEmission+(bConsts(i)*(logt**(i-1)))
            END DO
        !technically fit below is ok down to 100 K but bad fit seems better than no cooling at 70 K
        ELSE IF (gasTemperature .ge. 100.0) THEN
            DO i=1,SIZE(aConsts)
                collionallyInducedEmission=collionallyInducedEmission+(aConsts(i)*(logt**(i-1)))
            END DO
        END IF

        if (gasTemperature .ge. 100.0) collionallyInducedEmission=(10.0**collionallyInducedEmission)*tau
        !collionallyInducedEmission=(10.0**collionallyInducedEmission)*tau
    END FUNCTION collionallyInducedEmission


    !!-----------------------------------------------------------------------
    !! Continuum Emission
    !! Hirano & Yoshida (2013) and Ripamonti & Abel 2004 via Grassi 2012
    !!-----------------------------------------------------------------------
    REAL(dp) FUNCTION continuumEmission(gasTemperature,gasDensity)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity
        REAL(dp) :: massDensity,opacity,opticalDepth
        massDensity=min(0.5,gasDensity*MH*1.22) !assume mean molecular weight 1.22 and give max
        opacity=10.0**(1.000042*log(massDensity)+2.14989) !Lenzuni opacity fit

        opticalDepth=SQRT(3.14159*K_BOLTZ*gasTemperature/(massDensity*MH*1.22*GRAV_G))
        opticalDepth=opticalDepth*opacity*massDensity+1.0d-40 ! stop it going to zero

        continuumEmission=4.0*SB_CONST*(gasTemperature**4.0)*opacity*massDensity*min((opticalDepth**(-2.0)),1.0)
    END FUNCTION continuumEmission

    !!-----------------------------------------------------------------------
    !! Compton cooling
    !! Cen 1992 via Grassi 2012
    !! Cooling due to compton scattering of CMB photons.
    !! Shouldn't be important in near universe but include for completion
    !!-----------------------------------------------------------------------
    REAL(dp) FUNCTION comptonCooling(gasTemperature,gasDensity,elecAbund)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,elecAbund
        REAL(dp), PARAMETER :: cmbTemp=2.73
            comptonCooling=1.017d-37*(cmbTemp**4.0)*(gasTemperature-cmbTemp)*elecAbund*gasDensity
    END FUNCTION comptonCooling

    ! !-----------------------------------------------------------------------
    ! !  Grain + PAH photoelectric heating (with graphitic and silicate grains)
    !
    !  Use the treatment of Bakes & Tielens (1994, ApJ, 427, 822) with the
    !  modifications suggested by Wolfire et al. (2003, ApJ, 587, 278) to
    !  account for the revised PAH abundance estimate from Spitzer data.
    ! 
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION photoelectricHeating(gasTemperature,gasDensity,habingField,electronAbund,metallicity)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,habingField,electronAbund,metallicity
        REAL(dp), PARAMETER :: PHI_PAH=0.4d0,ALPHA=0.944D0
        REAL(dp) :: beta,delta,epsilon,nElec,PAH_HEATING_RATE,PAH_COOLING_RATE
        REAL(dp), PARAMETER ::C0=5.72D+0,C1=3.45D-2,C2=7.08D-3
        REAL(dp), PARAMETER ::C3=1.98D-2, C4=4.95D-1,C5=6.92D-1
        REAL(dp), PARAMETER ::C6=5.20D-1

        ! !Weingartner & Draine2001
        ! photoelectricHeating=1.0D-26*(habingField*gasDensity)*(C0+C1*gasTemperature**C4) &
        !         & /(1.0D0+C2*(habingField*SQRT(gasTemperature)/(gasDensity*electronAbund))**C5  &
        !         & *(1.0D0+C3*(habingField*SQRT(gasTemperature)/(gasDensity*electronAbund))**C6))

        !Bakes & Tielens 1994 with updates from Wolfire 2008
        !  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680, 384)
        !  Setting this factor to 1.0 gives the standard Bakes & Tielens expression

        !Skip photoelectric heating for now if electron abundance is zero because cooling is infinite
        !The only way to have no E- is bad initial conditions so this will resolve itself within a time step
        IF (electronAbund .gt. 1.0d-20) THEN
            nElec=electronAbund*gasDensity
            BETA=0.735D0/gasTemperature**0.068
            DELTA=habingField*SQRT(gasTemperature)/(nElec*PHI_PAH)
            EPSILON=4.87D-2/(1.0D0+4.0D-3*DELTA**0.73) + 3.65D-2*(gasTemperature/1.0D4)**0.7/(1.0D0+2.0D-4*DELTA)

            PAH_HEATING_RATE=1.30D-24*EPSILON*habingField*gasDensity
            PAH_COOLING_RATE=4.65D-30*gasTemperature**ALPHA*(DELTA**BETA)*nElec*PHI_PAH*gasDensity
            photoelectricHeating=PAH_HEATING_RATE-PAH_COOLING_RATE

            !Assume the PE heating rate scales linearly with PAH abundance
            photoelectricHeating=photoelectricHeating*metallicity!*(pahAbund/6.0d-7)
        ELSE
            photoelectricHeating=0.0
        END IF
    END FUNCTION photoelectricHeating

    ! !-----------------------------------------------------------------------
    ! !  H2 formation heating
    ! !
    ! !  Assume that 1.5 eV is liberated as heat during H2 formation
    ! !
    ! !  See: Hollenbach & Tielens (Review of Modern Physics, 1999, 71, 173)
    ! !  Use the H2 formation rate determined by subroutine H2_FORMATION_RATE
    ! !  and stored as REACTION_RATE(nRGR) (cm^3 s^-1)
    !! JH: I replaced REACTION_RATE(nRGR) with chemistry.f90's h2form=1.0d-17*dsqrt(Y(NEQ-1))

    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION H2FormationHeating(gasTemperature,gasDensity,hAbund,h2form)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,hAbund,h2form
        H2FormationHeating=(1.5*eV)*h2form*hAbund*gasDensity*gasDensity
    END FUNCTION H2FormationHeating


    !-----------------------------------------------------------------------
    !  H2 vibrational heating/cooling
    !
    !  Treat the vibrationally excited levels of H2 as a single pseudo level
    !  with effective rates of spontaneous emission, collisional excitation,
    !  FUV pumping and photodissociation that describe the behaviour of all
    !  the vibrational levels combined.
    !
    !  Use the treatment of Rollig et al. (2006, A&A, 451, 917)
    !-----------------------------------------------------------------------
    REAL(dp) FUNCTION H2VibrationalCooling(gasTemperature,gasDensity,h2Abund,h2dis)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,h2Abund,h2dis
        REAL(dp) :: photoDisRate,DELTA_E_10,A_COEFF_10,C_COEFF_10
        REAL(dp) ::DELTA_E_EFF,A_COEFF_EFF,R_PUMP_EFF,R_PHOTO_EFF
        DELTA_E_10=6587.0 ! Energy gap (K) between the v=1 and v=0 levels of H2
        A_COEFF_10=8.6D-7 ! Einstein A-coefficient (s^-1) for emission from the v=1 to v=0 level of H2
        C_COEFF_10=5.4D-13*SQRT(gasTemperature) ! Collisional rate coefficient (cm^3 s^-1) for v=0 to v=1
        photoDisRate=h2dis ! Photodissociation rate (s^-1) for the v=1 level of H2

        DELTA_E_EFF=23500.0 ! Characteristic vibrational level energy (K)
        A_COEFF_EFF=1.9D-6  ! Effective Einstein A-coefficient (s^-1)

        R_PUMP_EFF=11.2*h2dis ! Effective vibrational pumping rate (s^-1)
        R_PHOTO_EFF=18.0*h2dis ! Effective photodissociation rate (s^-1)

        H2VibrationalCooling=K_BOLTZ*DELTA_E_10*C_COEFF_10*gasDensity*EXP(-DELTA_E_10/gasTemperature)*h2Abund*gasDensity&
                           & *(A_COEFF_10+h2dis)/(C_COEFF_10*gasDensity+A_COEFF_10+h2dis)

        !Some heating from H2 vibrational interactions so subtract from cooling rate
        H2VibrationalCooling=H2VibrationalCooling-(h2Abund*gasDensity*(R_PUMP_EFF*K_BOLTZ*DELTA_E_EFF) &
                           & /(1.0D0+(A_COEFF_10+R_PHOTO_EFF)/(C_COEFF_10*gasDensity)))
    END FUNCTION H2VibrationalCooling

    ! !-----------------------------------------------------------------------
    ! !  H2 photodissociation heating
    ! !
    ! !  On average, 0.4 eV of kinetic energy per photodissociated molecule
    ! !
    ! !  Use the H2 photodissociation rate determined by the subroutine
    ! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRH2) (s^-1)
    ! !  JH: again, grabbed h2dis from chemistry.f90 for consistency.
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION H2PhotodisHeating(gasDensity,h2Abund,h2dis)
        REAL(dp), INTENT(IN) :: gasDensity,h2Abund,h2dis
        H2PhotodisHeating=(0.4*eV)*h2dis*gasDensity*h2Abund
    END FUNCTION H2PhotodisHeating

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
    REAL(dp) FUNCTION cosmicRayHeating(zeta,gasDensity,h2Abund)
        REAL(dp), INTENT(IN) :: zeta,gasDensity,h2Abund
        cosmicRayHeating=(20.0*eV)*(1.3D-17*zeta)*h2Abund*gasDensity
    END FUNCTION cosmicRayHeating


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
    REAL(dp) FUNCTION h2FUVPumpHeating(hAbund,h2Abund,gasTemperature,gasDensity,h2dis)
        REAL(dp), INTENT(IN) :: hAbund,h2Abund,gasTemperature,gasDensity,h2dis
        REAL(dp) :: NCRIT_H2
        NCRIT_H2=1.0D6/SQRT(gasTemperature)/(1.6D0*hAbund*EXP(-((400.0D0/gasTemperature)**2)) &
                                      & + 1.4D0*h2Abund*EXP(-(18100.0D0/(gasTemperature+1200.0D0))))

        h2FUVPumpHeating=(2.2*eV)*9.0D0*h2dis*gasDensity*h2Abund/(1.0D0+NCRIT_H2/gasDensity)
        ! !  If vibrationally excited H2 (H2*) is included in the chemical network,
        ! !  then use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722)
        !    IF(nH2v.NE.0) THEN
        !       H2_FUV_PUMPING_HEATING_RATE=(DENSITY(nH)*1.0D-12*SQRT(gasTemperature)*EXP(-1000.0D0/gasTemperature) &
        !                                & +DENSITY(nH2)*1.4D-12*SQRT(gasTemperature)*EXP(-18100.0D0/(gasTemperature+1200.0D0))) &
        !                                  *(2.6*eV)*DENSITY(nH2v)
        !    END IF

    END FUNCTION h2FUVPumpHeating   



! !-----------------------------------------------------------------------
! !  Carbon photoionization heating
! !
! !  On average, 1 eV of kinetic energy released per carbon ionization
! !  Use the carbon photoionization rate determined by the subroutine
! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRCI) (s^-1)
! !-----------------------------------------------------------------------

FUNCTION CarbonIonizationHeating(cIonizationRate,carbonAbund,gasDensity)
    real(dp), intent(in) :: cIonizationRate,carbonAbund,gasDensity
    real(dp) :: CarbonIonizationHeating
    CarbonIonizationHeating=(1.0*eV)*cIonizationRate*carbonAbund*gasDensity
END FUNCTION

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

Function chemicalHeating(gasDensity,exoReactants1,exoReactants2,exoRates,exothermicities)
REAL(dp), INTENT(IN) :: gasDensity,exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
REAL(dp) :: chemicalHeating

    chemicalHeating=SUM(exoReactants1*exoReactants2*exoRates*exothermicities)
    chemicalHeating=chemicalHeating*gasDensity*gasDensity*EV !each abundance should be a number dnesity to multiply through
  END FUNCTION chemicalHeating
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

FUNCTION gasGrainCollisions(gasTemperature,gasDensity,dustAbundance,dustRadius,dustTemp)
    real(dp), intent(in) :: gasTemperature,gasDensity,dustAbundance,dustRadius,dustTemp
    REAL(dp) :: gasGrainCollisions
    REAL(dp) :: nGrain,accommodation,C_GRAIN
    nGrain=dustAbundance*gasDensity

    !nGrain=2.0d-12*gasDensity
    C_GRAIN=PI*dustRadius**2

    !!$!  Accommodation fitting formula of Groenewegen (1994, A&A, 290, 531)
    !!$   ACCOMMODATION=0.35D0*EXP(-SQRT((dustTemperature+gasTemperature)/5.0D2))+0.1D0

    !  Accommodation coefficient of Burke & Hollenbach (1983, ApJ, 265, 223)
    accommodation=0.37D0*(1.0D0-0.8D0*EXP(-75.0D0/gasTemperature))

    gasGrainCollisions=nGrain*C_GRAIN*gasDensity*SQRT(8.0*K_BOLTZ*gasTemperature/(PI*MH)) &
                       & *accommodation*(2.0*K_BOLTZ*dustTemp-2.0*K_BOLTZ*gasTemperature)
END FUNCTION gasGrainCollisions

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

! FUNCTION calculateDustTemp(localField,surfaceField) result(dustTemperature)
!     USE constants
!     IMPLICIT NONE
!     !UV field in Habing at this depth and at surface required for this calculation
!     !Both in Habing as required for this treatment
!     REAL(dp), INTENT(IN) :: localField,surfaceField 
!     REAL(dp) :: dustTemperature

!     REAL(KIND=DP) :: NU_0,R_0,T_0,TAU_100

!     !Parameters used in the HHT equations (see their paper for details)
!     NU_0=2.65D15
!     TAU_100=1.0D-3
!     R_0=1.0D0/1.6d-21!avFactor

!     !Calculate the contribution to the dust temperature from the local FUV flux and the CMB background
!     !UCLPDR had afactor of 1.7 which I assume was Habing conversion so removed
!     dustTemperature=8.9D-11*NU_0*localField+T_CMB**5


!     !The minimum dust temperature is related to the incident FUV flux along each ray
!     T_0=12.2*surfaceField**0.2

!     !!Attenuate the FIR radiation produced in the surface layer
!     !!JH why is this commented?
!     !IF(PARTICLE(P)%TOTAL_COLUMN(J).GT.R_0) THEN
!     !     T_0=T_0*(PARTICLE(P)%TOTAL_COLUMN(J)/R_0)**(-0.4)
!     !END IF

!     !        Add the contribution to the dust temperature from the FUV flux incident along this ray
!     IF(T_0.GT.0) dustTemperature=dustTemperature &
!           & + (0.42-LOG(3.45D-2*TAU_100*T_0))*(3.45D-2*TAU_100*T_0)*T_0**5

!     !Convert from total dust emission intensity to dust temperature
!     dustTemperature=dustTemperature**0.2

!     !Calculate the contribution to the dust temperature from the local X-ray flux (assuming a fixed grain abundance of 1.6E-8)
!     !JH We have no xrays sthis
!     !dustTemperature=dustTemperature+1.5D4*(XRAY_ENERGY_DEPOSITION_RATE/1.6D-8)**0.2

!     !Impose a lower limit on the dust temperature, since values below 10 K can dramatically
!     !limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
!     IF(dustTemperature.LT.10) THEN
!         dustTemperature=10.0D0
!     END IF

!     !     Check that the dust temperature is physical
!     IF(dustTemperature.GT.1000) THEN
!         write(*,*) localField, surfaceField!WRITE(6,*) 'ERROR! Calculated dust temperature exceeds 1000 K'
!         dustTemperature=1000.0
!     END IF
!     RETURN
! END FUNCTION calculateDustTemp

! Using the new parametric formulation for dust temperature in Hocuk et al. 2017
! See equation 8 in DOI: 10.1051/0004-6361/201629944
FUNCTION calculateDustTemp(surfaceField,Av) result(dustTemperature)
    IMPLICIT NONE
    REAL(dp):: surfaceField, Av
    REAL(dp):: tanh_term
    REAL(dp):: dustTemperature
    tanh_term = 0.61D0 - LOG10(Av)
    dustTemperature = (11.0D0 + 5.7D0*TANH(tanh_term)) * (1.7D0*surfaceField)**(1.D0/5.9D0)

    !Impose a lower limit on the dust temperature, since values below 10 K can dramatically
    !limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
    IF(dustTemperature.LT.2.73) THEN
        dustTemperature=2.73D0
    END IF

    !     Check that the dust temperature is physical
    IF(dustTemperature.GT.1000) THEN
        write(*,*) Av, surfaceField!WRITE(6,*) 'ERROR! Calculated dust temperature exceeds 1000 K'
        dustTemperature=1000.0
    END IF
    RETURN
END FUNCTION calculateDustTemp
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
!    X_PRIME=1.83D0*ABUNDANCE(nelect)/(1.0D0+0.83D0*ABUNDANCE(nelect)) ! Correction to the electron abundance for a pure H2-He mixture

!    ETA_H2_He=1.0D0+(0.055D0-1.0D0)/(1.0D0+2.17D0*X_PRIME**0.366) ! Heating efficiency for a pure H2-He mixture
!    ETA_H_He =1.0D0+(0.117D0-1.0D0)/(1.0D0+7.95D0*ABUNDANCE(nelect)**0.678) ! Heating efficiency for a pure H-He mixture
!    ETA=(10.0D0*R*ETA_H2_He+ETA_H_He)/(10.0D0*R+1.0D0) ! Total heating efficiency for mixed atomic and molecular gas
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

!    L_TURB=5.0D0
!    TURBULENT_HEATING_RATE=3.5D-28*((V_TURB/1.0D5)**3)*(1.0D0/L_TURB)*GAS_DENSITY



SUBROUTINE pair_insertion_sort(array)
    REAL(dp), INTENT(inout) :: array(:)
    INTEGER :: i,j,last
    REAL(dp) :: t1,t2

    last=size(array)
    DO i=2,last-1,2
       t1=min(array(i),array(i+1))
       t2=max(array(i),array(i+1))
       j=i-1
       DO while((j.ge.1).and.(array(j).gt.t2))
          array(j+2)=array(j)
          j=j-1
       ENDDO
       array(j+2)=t2
       DO while((j.ge.1).and.(array(j).gt.t1))
          array(j+1)=array(j)
          j=j-1
       ENDDO
       array(j+1)=t1
    END DO

    IF(mod(last,2).eq.0)then
       t1=array(last)
       DO j=last-1,1,-1
          IF (array(j).le.t1) exit
          array(j+1)=array(j)
       END DO
       array(j+1)=t1
    ENDIF
END SUBROUTINE pair_insertion_sort

END MODULE heating

!Abandoned heating and cooling mechanisms

! ! !-----------------------------------------------------------------------
! ! !  Molecular Hydrogen cooling 
! ! !  from Galli & Palla (1998) via Grassi et al. (2014)
! ! !  Acceptable up to 10^5 K can easily use  Glover & Abel (2008) fits instead
! ! !  but they're more complex.
! ! ! JH: Dropped in favour of UCLPDR treatment of H2 rotational cooling + h2 Vibrational
! ! !-----------------------------------------------------------------------
! REAL(dp) FUNCTION h2Cooling(gasDensity,gasTemperature,hAbund,h2Abund)
!     REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,hAbund,h2Abund
!     REAL(dp) :: highDensLimit,lowDensLimit,T1000,logT
!     ! Temperature in kilokelvin
!     T1000=gasTemperature*0.001

!     !high density limit is same in all models: sum of vibrational and rotational cooling
!     highDensLimit=((9.5d-22*T1000**3.76)/(1+0.12*(T1000**2.1)))*dexp(-(0.13/T1000)**3.0)
!     highDensLimit=highDensLimit+(3.0d-24*dexp(-0.51/T1000))
    
!     highDensLimit=highDensLimit+(6.7d-19*dexp(-5.86/T1000))+(1.6d-18*dexp(-11.7/T1000))

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
!  Grain photoelectric heating (large grains only; r ~ 100 Å)
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
!         REAL(dp), PARAMETER :: DELTA_D=1.0D0
!         REAL(dp), PARAMETER :: DELTA_UV=1.8D0
!         REAL(dp), PARAMETER :: Y=0.1D0
!         REAL(dp), PARAMETER :: HNU_D=6.0D0
!         REAL(dp), PARAMETER :: HNU_H=13.6D0
!         REAL(dp) :: delta,gamma,XK,XD,X,XX
!         INTEGER :: ITERATION


!         XK=K_BOLTZ*gasTemperature/(HNU_H*eV)
!         XD=HNU_D/HNU_H
!         gamma=2.9D-4*Y*DSQRT(gasTemperature)*habingField/(gasDensity*electronAbund)
!         delta=XK-XD+gamma

!         !  Iterate to determine X by finding the zero of the function F
!         X=0.5D0
!         DO ITERATION=1,100
!           XX=X-(grainChargeFunc(X,DELTA,GAMMA)/deltaGrainChargeFunc(X,DELTA))
!           IF(ABS(XX-X).LT.1.0D-2) EXIT
!           X=XX
!         END DO
!         X=XX

!         IF(ITERATION.GE.100) THEN
!           WRITE(10,*)'WARNING! Grain parameter X not found in PE heating'
!           WRITE(10,*)'Using final value from interation loop: X =',X
!         END IF

!         photoelectricHeating=2.7D-25*DELTA_UV*DELTA_D*gasDensity*Y*habingField &
!                                      & *(((1.0D0-X)**2)/X + XK*((X**2)-1.0D0)/(X**2))

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
! !  Grain + PAH photoelectric heating (MRN size distribution; r = 3-100 Å)
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
!    PHI_PAH=0.4D0

!    ALPHA=0.944D0
!    BETA=0.735D0/gasTemperature**0.068
!    DELTA=HABING_FIELD*SQRT(gasTemperature)/(DENSITY(nelect)*PHI_PAH)
!    EPSILON=4.87D-2/(1.0D0+4.0D-3*DELTA**0.73) + 3.65D-2*(gasTemperature/1.0D4)**0.7/(1.0D0+2.0D-4*DELTA)

!    PAH_HEATING_RATE=1.30D-24*EPSILON*HABING_FIELD*GAS_DENSITY
!    PAH_COOLING_RATE=4.65D-30*gasTemperature**ALPHA*(DELTA**BETA)*DENSITY(nelect)*PHI_PAH*GAS_DENSITY

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
 !        LOGICAL :: binaryChopSearch,BRACKET_EXPANDED,TEMPERATURE_CONVERGED
 !        REAL(dp) :: heating,cooling,temperatureDiff,difference,outTemp,relative_difference
 !        REAL(dp),parameter :: TDIFF=0.01, FCRIT=0.1,TMIN=10.0, TMAX=1.0d5
 !        INTEGER :: tempLoops = 0.0

 !        previousTemp=0.0
 !        previousDifference=0.0
 !        thigh=TMAX
 !        tlow=TMIN

 !        binaryChopSearch=.False.
 !        BRACKET_EXPANDED=.False.
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
 !            relative_difference=2.0D0*ABS(difference)/ABS(heating+cooling)

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
 !                    getEquilibriumTemp=1.3D0*getEquilibriumTemp
 !                    previousDifference=difference
 !                    THIGH=getEquilibriumTemp ! Update the value of T_high
 !        !     If the cooling continues to outweigh the heating, decrease the temperature by 30%
 !              ELSE IF(DIFFERENCE.LT.0 .AND. previousDifference.LE.0) THEN
 !                 THIGH=getEquilibriumTemp ! Update the value of T_high
 !                 getEquilibriumTemp=0.7D0*getEquilibriumTemp
 !                 previousDifference=DIFFERENCE
 !                 TLOW=getEquilibriumTemp ! Update the value of T_low
                

 !        !     If the heating-cooling balance has reversed (either from net heating to net cooling or
 !        !     vice-versa) then switch to the binary chop search method to determine the temperature
 !              ELSE
 !                 getEquilibriumTemp=(THIGH+TLOW)/2.0D0
 !                 previousDifference=DIFFERENCE
 !                 binaryChopSearch=.TRUE. ! From now on
 !              END IF

 !        !  Perform a binary chop search (the min-max range was found by the 30% increase/decrease method)
 !           ELSE

 !              IF(DIFFERENCE.GT.0) THEN
 !                TLOW=getEquilibriumTemp ! Update the value of T_low
 !                 getEquilibriumTemp=(getEquilibriumTemp+THIGH)/2.0D0
 !                 previousDifference=DIFFERENCE
                 
 !              END IF
 !              IF(DIFFERENCE.LT.0) THEN
 !                 THIGH=getEquilibriumTemp !Update the value of T_high
 !                 getEquilibriumTemp=(getEquilibriumTemp+TLOW)/2.0D0
 !                 previousDifference=DIFFERENCE
 !              END IF

 !           END IF

 !        !  If the search routine is unable to converge on a temperature that satisfies the thermal balance
 !        !  criterion, expand the min-max search bracket asymmetrically and begin to narrow the search again
 !        !  If the repeated search fails to converge once more, force convergence at the current temperature
 !           IF(temperatureDiff.LE.TDIFF) THEN
 !              IF(.NOT.BRACKET_EXPANDED) THEN
 !                 THIGH=THIGH+SQRT(PI)
 !                 TLOW=TLOW-SQRT(2.0)
 !                 BRACKET_EXPANDED=.TRUE.
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