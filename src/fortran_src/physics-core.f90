!This module contains all the must have physics variables. It is imported by both
!the other physics modules (so they can modify the physics) and the chemistry
!module (so it can use the physics for reaction rates).

module physicscore
    use constants, only: dp, C, SECONDS_PER_YEAR, PC, PI, HP, K_BOLTZ, aunit, eV, &
        uISRF, uISRF_UV, Lsun, SB_CONST
    use DEFAULTPARAMETERS
    use extinction_module, only: extcurve_obs
    !f2py INTEGER, parameter :: dp

    implicit none

    public

    !Use main loop counters in calculations so they're kept here
    integer :: dstep

    !Optional CR attenuation with column density and better H2 dissociation rates
    real(dp) :: h2CRPRate,zetaScale

    !variables either controlled by physics or that user may wish to change
    real(dp) :: timeInYears,targetTime,currentTimeOld
    real(dp) ::  cloudSize
    real(dp), allocatable :: av(:),coldens(:),gasTemp(:),dustTemp(:),density(:),density_max(:)
    ! Per-parcel internal Av shielding (get_coldens_internal/1.6e21); 0 for models without a protostar.
    real(dp), allocatable :: av_internal(:)
    ! Per-parcel internal radiation field [Habing units]; 0 for models without a protostar.
    real(dp), allocatable :: radfield_internal(:)
    ! Per-point initial density used by get_densdot for multi-point radial profile models.
    ! Defaults to global initialDens; overridden by cloud.f90 when enable_radiative_transfer=T.
    real(dp), allocatable :: initialDens_array(:)
    ! Radial position of each parcel (pc). Set by collapse/cloud/hotcore for 1D models; 0 otherwise.
    real(dp), allocatable :: parcel_radius(:)
    ! coldens time series of the last-processed outer shell; used in (points,time) loop order
    ! so inner shells can look up coldens(dstep+1, dtime) exactly. Allocated in wrap.f90.
    real(dp), allocatable :: coldens_history(:)
    ! Set by wrap.f90 before each modelUpdatePhysics call. Holds coldens(dstep+1) at the
    ! current timestep (from coldens_history), so the exact edge-to-core accumulation is
    ! preserved without accessing a live coldens(dstep+1) that may not yet be computed.
    real(dp) :: outer_coldens_for_current_step = 0.0_dp

    !Arrays for calculating rates
    !if ionModel = L use the L model coefficients, if = H use the H model

    ! REAL(dp),PARAMETER :: ckLIon(10)=(/1.545456645800e7_dp, -6.307708626617e6_dp, 1.142680666041e6_dp, -1.205932302621e5_dp,&
    ! 8.170913352693e3_dp, -3.686121296079e2_dp,1.107203722057e1_dp, -2.135293914267e-1_dp,&
    ! 2.399219033781e-3_dp, -1.196664901916e-5_dp/)

    ! REAL(dp),PARAMETER :: ckHIon(10)=(/1.223529865309e7_dp, -5.013766644305e6_dp, 9.120125566763e5_dp, -9.665446168847e4_dp,&
    !     6.576930812109e3_dp, -2.979875686226e2_dp,8.989721355058_dp, -1.741300519598e-1_dp,&
    !     1.965098116126e-3_dp, -9.844203439473e-6_dp/)

    !ckLIon and ckHIon are the coefficients for the cosmic ray ionization rate
    ! Values from Padovani et al. 2018
    integer, parameter :: k_max = 10
    real(dp), parameter :: ckLIon(k_max) = (/-3.331056497233e6_dp,&
                                        1.207744586503e6_dp,&
                                        -1.913914106234e5_dp,&
                                        1.731822350618e4_dp,&
                                        -9.790557206178e2_dp,&
                                        3.543830893824e1_dp,&
                                        -8.034869454520e-1_dp,&
                                        1.04880859308e-2_dp,&
                                        -6.188760100997e-5_dp,&
                                        3.122820990797e-8_dp/)

    real(dp), parameter :: ckHion(k_max) = (/1.001098610761e7_dp,&
                                        -4.231294690194e6_dp,&
                                        7.921914432011e5_dp,&
                                        -8.623677095423e4_dp,&
                                        6.015889127529e3_dp,&
                                        -2.789238383353e2_dp,&
                                        8.595814402406_dp,&
                                        -1.698029737474e-1_dp,&
                                        1.951179287567e-3_dp,&
                                        -9.937499546711e-6_dp/)

    !ckLDiss and ckHDiss are the coefficients for the H2 dissociation rate
    real(dp), parameter :: ckLDiss(k_max)=(/1.582911005330e7_dp,-6.465722684896e6_dp, 1.172189025424e6_dp, -1.237950798073e5_dp, &
        8.393404654312e3_dp, -3.788811358130e2_dp, 1.138688455029e1_dp, -2.197136304567e-1_dp, &
        2.469841278950e-3_dp, -1.232393620924e-5_dp/)
    real(dp), parameter :: ckHDiss(k_max)=(/1.217227462831e7_dp,-4.989649250304e6_dp, 9.079152156645e5_dp, -9.624890825395e4_dp, &
        6.551161486120e3_dp, -2.968976216187e2_dp, 8.959037875226_dp, -1.735757324445e-1_dp, &
        1.959267277734e-3_dp, -9.816996707980e-6_dp/)

    ! parameters for 1D radiation field calculation
    real(dp), parameter :: Tdsub = 1500.0_dp  !sublimation/melted temperature for dust grains
    ! ! range of wavelength for integration
    real(dp), parameter :: wave1 = HP * C * 1.0e4_dp / (13.6_dp * eV)  !in micron
    real(dp), parameter :: wave2 = 20.0_dp  ! in micron

    interface get_av
        module procedure get_av_array, get_av_scalar
    end interface get_av

contains
    !basic initialization of physics. All physics modules should call this and then
    !do their own initialization.
    subroutine coreInitializePhysics(successFlag)
        integer, intent(out) :: successFlag
        timeInYears=currentTime/SECONDS_PER_YEAR

        ! Modules not restarted in python wraps so best to reset everything manually.
        if (ALLOCATED(av)) deallocate(av,coldens,gasTemp,dustTemp,density,density_max)
        allocate(av(points),coldens(points),gasTemp(points),dustTemp(points),density(points),density_max(points))
        if (ALLOCATED(av_internal))       deallocate(av_internal)
        if (ALLOCATED(radfield_internal)) deallocate(radfield_internal)
        allocate(av_internal(points), radfield_internal(points))
        av_internal       = 0.0_dp
        radfield_internal = 0.0_dp
        if (ALLOCATED(initialDens_array)) deallocate(initialDens_array)
        allocate(initialDens_array(points))
        initialDens_array = initialDens   ! default: same for all points (single-point / no radial profile)
        if (ALLOCATED(parcel_radius)) deallocate(parcel_radius)
        allocate(parcel_radius(points))
        if (ALLOCATED(coldens_history)) deallocate(coldens_history)
        ! coldens_history is allocated in wrap.f90 where timePoints is known
        parcel_radius = 0.0_dp

        cloudSize = (rout-rin)*pc
        gasTemp=initialTemp
        dustTemp=gasTemp
        density=initialDens
        currentTimeOld=0.0_dp
        if (.not. ((ionModel == "L") .or. (ionModel == "H"))) then
            successFlag=-1
            write(*,*) "Error: ionModel must be either L or H"
            return
        end if
        if ((improvedH2CRPDissociation) .and. (.not. cosmicRayAttenuation)) then
            successFlag=-1
            write(*,*) "Error: improvedH2CRPDissociation requires cosmicRayAttenuation to also be True"
            return
        end if
        if ((freefall) .and. (finalDens < initialDens)) then
            successFlag=-1
            write(*,*) "Error: freefall finalDens (", finalDens, ") is less than initialDens (", initialDens, ")"
            return
        end if

        !calculate initial column density as distance from core edge to current point * density
        do dstep=1,points
            coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
        end do
        av = get_av(baseAv, coldens)
        zetaScale=zeta
    end subroutine coreInitializePhysics

    subroutine coreUpdatePhysics
        ! In the 1D RT (points,time) loop, modelUpdatePhysics (cloud.f90) owns
        ! coldens and av using the edge-to-core accumulation with coldens_history.
        ! Skip the center-to-edge accumulation here to avoid clobbering those values.
        if (.NOT. (enable_radiative_transfer .AND. points>1)) then
            !calculate column density. Remember dstep counts from core center to edge
            !and coldens should be amount of gas from edge to parcel.
            coldens(dstep)=cloudSize/real(points)*density(dstep)

            ! add previous column densities to current as we move into cloud to get total
            if (dstep < points) coldens(dstep)=coldens(dstep)+coldens(dstep-1)

            av(dstep)= get_av(baseAv, coldens(dstep))
        end if
        if (.not. heatingFlag) then
            dustTemp(dstep)=gasTemp(dstep)
        end if

        if (cosmicRayAttenuation) call ionizationDependency
    end subroutine coreUpdatePhysics

    pure function get_densdot(density) result(densdot)
    ! Returns the time derivative of the density.
    ! Analytical function taken from Rawlings et al. 1992
    ! Called from chemistry.f90, density integrated with abundances so this gives ydot
    ! Uses initialDens_array(dstep) as the per-point reference density so that
    ! multi-point radial-profile models (enable_radiative_transfer=T) undergo
    ! freefall correctly from each parcel's own starting density.
    real(dp), intent(in) :: density
    real(dp) :: densdot
    real(dp) :: n0_pt   ! per-point initial density for this dstep
    n0_pt = initialDens_array(dstep)
    !Rawlings et al. 1992 freefall collapse. With freefallFactor for B-field etc
    if ((density < finalDens) .and. (freefall) .and. (density > n0_pt)) then
        densdot=freefallFactor*(density**4.0_dp/n0_pt)**0.33_dp*&
        &(8.4e-30_dp*n0_pt*((density/n0_pt)**0.33_dp-1.0_dp))**0.5_dp
    else
        densdot=0.0
    end if
    end function get_densdot

    pure function get_av_array(baseAv, coldens) result(av)
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        real(dp), intent(in) :: baseAv
        real(dp), intent(in), dimension(:) :: coldens

        real(dp), dimension(size(coldens)) :: av

        av = baseAv + coldens/1.6e21_dp
    end function get_av_array

    pure function get_av_scalar(baseAv, coldens) result(av)
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        real(dp), intent(in) :: baseAv
        real(dp), intent(in) :: coldens

        real(dp) :: av

        av = baseAv + coldens/1.6e21_dp
    end function get_av_scalar


    pure function getDByDnDensdot(density) result(dByDnDensdot)
    !Defunct function which provides the necessary derivative d(dn/dt)/dn
    !in the case one uses a Jacobian.
    real(dp), intent(in) :: density
    real(dp) :: dByDnDensdot
    !Rawlings et al. 1992 freefall collapse. With freefallFactor for B-field etc
    if (density < finalDens) then
        dByDnDensdot=freefallFactor*8.4e-30_dp*(density**3)*((9.0_dp*((density/initialDens)**0.33_dp))-8.0_dp)
        dByDnDensdot=dByDnDensdot/(6.0_dp*(((density**4)/initialDens)**0.66_dp))
        dByDnDensdot=dByDnDensdot/sqrt(initialDens*8.4e-30_dp*(((density/initialDens)**0.33_dp))-1.0_dp)
    else
        dByDnDensdot=0.0_dp
    end if
    end function getDByDnDensdot

    subroutine ionizationDependency
        real(dp) :: dissSum,dRate,zSum,ionRate
        integer :: k
        !Attenuate CR by column density
        zeta = 1.0_dp
        zSum = 0.0_dp
        do k=0,9,1
            if (ionModel == "L") then
                ionRate=ckLIon(k+1)*log10(coldens(dstep))**k
            else if (ionModel == "H") then
                ionRate=ckHIon(k+1)*log10(coldens(dstep))**k
            else
                write(*,*) "WARNING: ionModel switch must be 0 or 1"
            end if
            zSum=zSum+ionRate
        end do

        ! update/overwrite zeta with attenuated value
        zeta = ((10**zSum)/1.3e-17_dp)* zetaScale

        !rate calculation for H2 dissociation
        if (improvedH2CRPDissociation) then
            dissSum = 0.0_dp
            do k=0,9,1
                if (ionModel == "L") then
                    dRate=ckLDiss(k+1)*log10(coldens(dstep))**k
                else if (ionModel == "H") then
                    dRate=ckHDiss(k+1)*log10(coldens(dstep))**k
                else
                    write(*,*) "WARNING: ionModel switch must be L or H"
                end if
                dissSum=dissSum+dRate
            end do
            h2CRPRate=(10**dissSum)*zetaScale
        end if
    end subroutine ionizationDependency

    ! Analytical column density from center (r=0) to r [cm^-2].
    pure subroutine findcoldens_core2edge(coldens,rho0,density_scale_radius,density_power_index,r)
      real(dp),intent(in) :: r,rho0,density_scale_radius,density_power_index
      real(dp),intent(out) :: coldens

      if (r <= density_scale_radius) then
          coldens = rho0 * r * pc
      else
          coldens = rho0*density_scale_radius*pc * (1.0_dp + (1.0_dp/(density_power_index-1.0_dp)) * &
              & (1.0_dp - (r/density_scale_radius)**(1.0_dp-density_power_index)))
      end if

    end subroutine findcoldens_core2edge

    pure subroutine findcoldens_edge2core(coldens,rho0,density_scale_radius,density_power_index,r)
        real(dp),intent(in) :: rho0,density_scale_radius,density_power_index,r
        real(dp),intent(out) :: coldens
        if (r>density_scale_radius) then
            coldens = rho0*density_scale_radius*pc/(density_power_index-1.0_dp) &
                & * (r/density_scale_radius)**(1.0_dp-density_power_index)
        else
            coldens = rho0*density_scale_radius*pc*(density_power_index/(density_power_index-1.0_dp)-r/density_scale_radius)
        end if
    end subroutine findcoldens_edge2core

    ! Column density shielding from external UV (stage 1 / cloud): edge-to-parcel integral.
    pure function get_coldens_external(r, rho0) result(coldens_external)
        real(dp), intent(in) :: r    ! parcel radius [pc]
        real(dp), intent(in) :: rho0  ! reference density [cm-3]
        real(dp) :: coldens_external
        call findcoldens_edge2core(coldens_external, rho0, density_scale_radius, &
                                   density_power_index, r)
    end function get_coldens_external

    ! Column density shielding from central protostar (stage 2 / hotcore): integral from center to parcel.
    pure function get_coldens_internal(r) result(coldens_internal)
        real(dp), intent(in) :: r    ! parcel radius [pc]
        real(dp) :: coldens_internal
        call findcoldens_core2edge(coldens_internal, finalDens, &
                                   density_scale_radius, density_power_index, r)
    end function get_coldens_internal

    ! The profile of the gas volume density
    ! REAL(dp) FUNCTION rhofit(r,rho0,r0,a)
    pure function get_ngas_r(r,rho0,density_scale_radius,density_power_index) result(ngas_r)
      real(dp), intent(in) :: r,rho0,density_scale_radius,density_power_index
      real(dp) :: ngas_r
      ! [r] in pc, [density_scale_radius] in pc
      ngas_r = rho0/(1.0_dp + (r/density_scale_radius)**density_power_index)

    end function get_ngas_r

    pure function get_initialDens_r(r,p) result(initialDens_r)
        real(dp), intent(in) :: r,p
        real(dp) :: initialDens_r

        real(dp) :: t, logn0, logr0,n0_init,r0_init

        t = 0.0_dp
        logn0=61.8_dp*(1.175e6_dp-t)**(-0.01_dp) - 49.4_dp
        logr0=-28.5_dp*(1.175e6_dp-t)**(-0.01_dp) + 28.93_dp
        n0_init=10**(logn0)
        r0_init=10**(logr0) * aunit
        initialDens_r = 1.0_dp+(r/r0_init)**p
        initialDens_r = n0_init/initialDens_r
    end function get_initialDens_r

    subroutine radiation(r, Lstar, Tstar, Avs, Temp_dust, U)
        real(dp), intent(in) :: Lstar, Tstar, Avs, r
        real(dp), intent(out) :: Temp_dust, U

        real(dp)  :: rsub, U_star, U_shell
        integer, parameter :: nw=129

        ! sublimation distance
        rsub = get_rsub(Lstar)

        ! radiation from the star
        call radiation_star(r, Lstar, Tstar, Avs, U_star)

        ! radiation from the shell
        if (r<rsub) then
            U_shell = 0.0_dp
        else
            call radiation_shell(r, Lstar, Tstar, Avs, U_shell)
        end if

        ! total radiation field
        U = U_star + U_shell

        ! dust temperature at equilibrium with the radiation field
        Temp_dust=get_temp_average_grain(U)

    end subroutine radiation

    pure subroutine radiation_star(r, Lstar, Tstar, Avs, U)
        real(dp), intent(in) :: Lstar, Tstar, r, Avs
        real(dp), intent(out) :: U
        integer :: i
        real(dp), dimension(:), allocatable :: wave, wave_cm, Istar, uwave_star, tau_wave, uwave_red
        real(dp) :: ZZ, rsub, Rstar
        real(dp) :: NH_EBV, RV
        real(dp)  :: urad_red
        integer, parameter :: nw=129
        real(dp), dimension(2, nw) :: ext_curves
        character(len=10) :: model

        RV = 4.0_dp
        NH_EBV = 5.8e21_dp

        ! sublimation distance
        rsub=get_rsub(Lstar)

        ZZ = HP * C / (K_BOLTZ * Tstar)

        !logspace for wave in micron
        allocate(wave(nw))
        call logspace(log10(wave1), log10(wave2), nw, wave)

        ! convert wave from micron to cm
        wave_cm = wave*1.0e-4_dp  !in cm

        ! Call the function from the module
        call extcurve_obs(wave, RV, NH_EBV, model, ext_curves)

        Istar     = (2.0_dp*HP*C**2.0_dp/wave_cm**5.0_dp)*(1.0_dp/(exp(ZZ/wave_cm)-1.0_dp))
        uwave_star = (4.0_dp*PI*wave_cm/C)*(Istar)/wave_cm

        tau_wave = Avs * ext_curves(1,:)/1.086_dp
        uwave_red = uwave_star*exp(-tau_wave)

        ! Apply trapezoidal rule for numerical integration
        urad_red = 0.0_dp
        do i = 1, nw-1
            urad_red = urad_red + 0.5_dp * (wave_cm(i+1) - wave_cm(i)) * (uwave_red(i+1) + uwave_red(i))
        end do

        ! The stellar radius
        Rstar=get_Rstar_rsub(Tstar,Tdsub)*rsub

        ! The total radiation field (dimensionless)
        U = urad_red / uISRF * (r / Rstar)**(-2.0_dp)
    end subroutine radiation_star

    pure subroutine radiation_shell(r, Lstar, Tstar, Avs, U)
        real(dp), intent(in) :: Lstar, Tstar, r, Avs
        real(dp), intent(out) :: U
        integer :: i
        real(dp), dimension(:), allocatable :: wave, wave_cm, Istar, uwave_star, tau_wave, uwave_red
        real(dp) :: ZZ, rsub
        real(dp) :: NH_EBV, RV
        real(dp)  :: urad_red, Tshell
        integer, parameter :: nw=129
        real(dp), dimension(2, nw) :: ext_curves
        character(len=10) :: model

        RV = 4.0_dp
        NH_EBV = 5.8e21_dp

        ! sublimation distance
        rsub=get_rsub(Lstar)

        ! calculate Tshell
        Tshell = get_Tshell(Tstar)

        ! parameters for Planck function
        ZZ = HP * C / (K_BOLTZ * Tshell)

        !logspace for wave in micron
        allocate(wave(nw))
        call logspace(log10(wave1), log10(wave2), nw, wave)

        ! convert wave from micron to cm
        wave_cm = wave*1.0e-4_dp  !in cm

        ! Call the function from the module
        call extcurve_obs(wave, RV, NH_EBV, model, ext_curves)

        Istar     = (2.0_dp*HP*C**2.0_dp/wave_cm**5.0_dp)*(1.0_dp/(exp(ZZ/wave_cm)-1.0_dp))
        uwave_star = (4.0_dp*PI*wave_cm/C)*(Istar)/wave_cm

        tau_wave = Avs * ext_curves(1,:)/1.086_dp
        uwave_red = uwave_star*exp(-tau_wave)

        ! Apply trapezoidal rule for numerical integration
        urad_red = 0.0_dp
        do i = 1, nw-1
            urad_red = urad_red + 0.5_dp * (wave_cm(i+1) - wave_cm(i)) * (uwave_red(i+1) + uwave_red(i))
        end do
        U = urad_red / uISRF * (r / rsub)**(-2.0_dp)
    end subroutine radiation_shell

    pure subroutine logspace(start, stop, num, result)
        real(dp), intent(in) :: start, stop
        integer, intent(in) :: num
        real(dp), dimension(num), intent(out) :: result
        integer :: i

        do i = 1, num
            result(i) = 10.0_dp**(start + (i-1)*(stop-start)/DBLE(num-1))
        end do
    end subroutine logspace

    pure function get_Rstar_rsub(Tstar,Tdmax) result(Rstar_rsub)
        real(dp), intent(in) :: Tstar,Tdmax
        real(dp) :: Rstar_rsub
        real(dp) :: f1,f2

        f1 = sqrt(1.0e6_dp * Lsun)
        f2 = 155.3_dp * aunit * (Tdmax/1500.0_dp)**(-5.6_dp/2.0_dp) * sqrt(4.0_dp*PI*SB_CONST) * Tstar**2
        Rstar_rsub= f1/f2
    end function get_Rstar_rsub

    pure function get_temp_average_grain(U) result(temp_average)
        real(dp), intent(in) :: U
        real(dp) :: temp_average
        real(dp) :: Td_sil,Td_car
        Td_sil = 16.4_dp * U**(1.0_dp/6.0_dp)  ! For silicate grains
        Td_car = 19.5_dp * U**(1.0_dp/5.6_dp)  ! For carbon grains
        temp_average = Td_sil**(4.0_dp) + Td_car**(4.0_dp)
        temp_average = (0.5_dp*temp_average)**(1.0_dp/4.0_dp)
    end function get_temp_average_grain

    pure function get_Tshell(Tstar) result(Tshell)
        real(dp), intent(in) :: Tstar
        real(dp) :: Tshell

        Tshell = get_Rstar_rsub(Tstar,Tdsub)
        Tshell = Tshell**(0.5_dp) * Tstar
    end function get_Tshell

    pure function get_rsub(Lstar) result(rsub)
        real(dp), intent(in) :: Lstar
        real(dp) :: rsub
        rsub = 155.3_dp*(Lstar/1.0e6_dp/Lsun)**(0.5_dp) * (Tdsub/1500.0_dp)**(-5.6_dp/2.0_dp) * aunit  !in cm
    end function get_rsub

    ! Unattenuated UV radiation field from central protostar at radius r_cm [Habing units].
    ! Uses 45% of the bolometric luminosity as the UV fraction and scales as r^-2.
    pure function get_G0_internal_at_r(Lstar, r_cm) result(G0_internal_at_r)
        real(dp), intent(in) :: Lstar   ! bolometric luminosity [erg s^-1]
        real(dp), intent(in) :: r_cm    ! parcel radius [cm]
        real(dp) :: G0_internal_at_r
        real(dp) :: Luv
        Luv = 0.45_dp * Lstar
        G0_internal_at_r = Luv / (4.0_dp * PI * C * r_cm**2) / uISRF_UV
    end function get_G0_internal_at_r

end module physicscore
