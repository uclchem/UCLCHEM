!This module contains all the must have physics variables. It is imported by both
!the other physics modules (so they can modify the physics) and the chemistry
!module (so it can use the physics for reaction rates).

module physicscore
    use constants
    use DEFAULTPARAMETERS
    use extinction_module
    !f2py INTEGER, parameter :: dp
    implicit none
    !Use main loop counters in calculations so they're kept here
    integer :: dstep


    !Optional CR attenuation with column density and better H2 dissociation rates
    real(dp) :: h2CRPRate,zetaScale

    !variables either controlled by physics or that user may wish to change
    real(dp) :: timeInYears,targetTime,currentTimeold
    real(dp) ::  cloudSize
    real(dp), allocatable :: av(:),coldens(:),gasTemp(:),dustTemp(:),density(:),density_max(:)
    ! Per-parcel internal Av shielding (coldens_internal/1.6e21); 0 for models without a protostar.
    real(dp), allocatable :: av_internal(:)
    ! Per-parcel internal radiation field [Habing units]; 0 for models without a protostar.
    real(dp), allocatable :: radfield_internal(:)
    ! Per-point initial density used by densdot for multi-point radial profile models.
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

    ! REAL(dp),PARAMETER :: ckLIon(10)=(/1.545456645800d7, -6.307708626617d6, 1.142680666041d6, -1.205932302621d5,&
    ! 8.170913352693d3, -3.686121296079d2,1.107203722057d1, -2.135293914267d-1,&
    ! 2.399219033781d-3, -1.196664901916d-5/)

    ! REAL(dp),PARAMETER :: ckHIon(10)=(/1.223529865309d7, -5.013766644305d6, 9.120125566763d5, -9.665446168847d4,&
    !     6.576930812109d3, -2.979875686226d2,8.989721355058d0, -1.741300519598d-1,&
    !     1.965098116126d-3, -9.844203439473d-6/)

    !ckLIon and ckHIon are the coefficients for the cosmic ray ionization rate
    ! VALUES FROM Padovani et al. 2018
    real(dp), parameter :: ckLIon(10) = (/-3.331056497233d6,&
                                        1.207744586503d6,&
                                        -1.913914106234d5,&
                                        1.731822350618d4,&
                                        -9.790557206178d2,&
                                        3.543830893824d1,&
                                        -8.034869454520d-1,&
                                        1.04880859308d-2,&
                                        -6.188760100997d-5,&
                                        3.122820990797d-8/)

    real(dp), parameter :: ckHion(10) = (/1.001098610761d7,&
                                        -4.231294690194d6,&
                                        7.921914432011d5,&
                                        -8.623677095423d4,&
                                        6.015889127529d3,&
                                        -2.789238383353d2,&
                                        8.595814402406d0,&
                                        -1.698029737474d-1,&
                                        1.951179287567d-3,&
                                        -9.937499546711d-6/)

    !ckLDiss and ckHDiss are the coefficients for the H2 dissociation rate
    real(dp), parameter :: ckLDiss(10)=(/1.582911005330d7,-6.465722684896d6, 1.172189025424d6, -1.237950798073d5, &
        8.393404654312d3, -3.788811358130d2, 1.138688455029d1, -2.197136304567d-1, &
        2.469841278950d-3, -1.232393620924d-5/)
    real(dp), parameter :: ckHDiss(10)=(/1.217227462831d7,-4.989649250304d6, 9.079152156645d5, -9.624890825395d4, &
        6.551161486120d3, -2.968976216187d2, 8.959037875226d0, -1.735757324445d-1, &
        1.959267277734d-3, -9.816996707980d-6/)

    ! parameters for 1D radiation field calculation
    real(dp), parameter :: Tdsub = 1500.0d0  !sublimation/melted temperature for dust grains
    ! ! range of wavelength for integration
    real(dp), parameter :: wave1 = HP * C * 1.0d4 / (13.6d0 * EV)  !in micron
    real(dp), parameter :: wave2 = 20.0d0  ! in micron

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
        currentTimeOld=0.0
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
          !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av= baseAv + coldens/1.6d21
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

            !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
            av(dstep)= baseAv + coldens(dstep)/1.6d21
        end if
        if (.not. heatingFlag) then
            dustTemp(dstep)=gasTemp(dstep)
        end if

        if (cosmicRayAttenuation) call ionizationDependency
    end subroutine coreUpdatePhysics

    function densdot(density)
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
        densdot=freefallFactor*(density**4.0/n0_pt)**0.33*&
        &(8.4d-30*n0_pt*((density/n0_pt)**0.33-1.0))**0.5
    else
        densdot=0.0
    end if
    end function densdot


    pure function dByDnDensdot(density)
    !Defunct function which provides the necessary derivative d(dn/dt)/dn
    !in the case one uses a Jacobian.
    real(dp), intent(in) :: density
    real(dp) :: dByDnDensdot
    !Rawlings et al. 1992 freefall collapse. With freefallFactor for B-field etc
    if (density < finalDens) then
        dByDndensdot=freefallFactor*8.4d-30*(density**3)*((9.0d0*((density/initialDens)**0.33))-8.0d0)
        dByDnDensdot=dByDnDensdot/(6.0d0*(((density**4.0)/initialDens)**0.66))
        dByDnDensdot=dByDnDensdot/sqrt(initialDens*8.4d-30*(((density/initialDens)**0.33))-1.0d0)
    else
        dByDnDensdot=0.0
    end if
    end function dByDnDensdot

    subroutine ionizationDependency
        real(dp) :: dissSum,dRate,zSum,ionRate
        integer :: k
        !Attenuate CR by column density
        zeta = 1.0
        zSum = 0
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
        zeta = ((10**zSum)/1.3d-17)* zetaScale

        !rate calculation for H2 dissociation
        if (improvedH2CRPDissociation) then
            dissSum = 0
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
    subroutine findcoldens_core2edge(coldens,rho0,density_scale_radius,density_power_index,r)
      real(dp),intent(in) :: r,rho0,density_scale_radius,density_power_index
      real(dp),intent(out) :: coldens

      if (r <= density_scale_radius) then
          coldens = rho0 * r * pc
      else
          coldens = rho0*density_scale_radius*pc * (1.0d0 + (1.0d0/(density_power_index-1.0d0)) * &
              & (1.0d0 - (r/density_scale_radius)**(1.0d0-density_power_index)))
      end if

    end subroutine findcoldens_core2edge

    subroutine findcoldens_edge2core(coldens,rho0,density_scale_radius,density_power_index,r)
        real(dp),intent(in) :: rho0,density_scale_radius,density_power_index,r
        real(dp),intent(out) :: coldens
        if (r>density_scale_radius) then
            coldens = rho0*density_scale_radius*pc/(density_power_index-1.0d0) &
                & * (r/density_scale_radius)**(1.0d0-density_power_index)
        else
            coldens = rho0*density_scale_radius*pc*(density_power_index/(density_power_index-1.0d0)-r/density_scale_radius)
        end if
    end subroutine findcoldens_edge2core

    ! Column density shielding from external UV (stage 1 / cloud): edge-to-parcel integral.
    real(dp) function coldens_external(r, rho0)
        real(dp), intent(in) :: r    ! parcel radius [pc]
        real(dp), intent(in) :: rho0  ! reference density [cm-3]
        call findcoldens_edge2core(coldens_external, rho0, density_scale_radius, &
                                   density_power_index, r)
    end function coldens_external

    ! Column density shielding from central protostar (stage 2 / hotcore): integral from center to parcel.
    real(dp) function coldens_internal(r)
        real(dp), intent(in) :: r    ! parcel radius [pc]
        call findcoldens_core2edge(coldens_internal, finalDens, &
                                   density_scale_radius, density_power_index, r)
    end function coldens_internal

    ! The profile of the gas volume density
    ! REAL(dp) FUNCTION rhofit(r,rho0,r0,a)
    real(dp) function ngas_r(r,rho0,density_scale_radius,density_power_index)
      real(dp) :: r,rho0,density_scale_radius,density_power_index
      ! [r] in pc, [density_scale_radius] in pc
      ngas_r = rho0/(1.0d0 + (r/density_scale_radius)**density_power_index)

    end function ngas_r

    real(dp) function initialDens_r(r,p)
        real(dp) :: logn0, logr0,n0_init,r0_init
        real(dp) :: r,t,p
        t = 0.0d0
        logn0=61.8d0*(1.175d6-t)**(-0.01) - 49.4d0
        logr0=-28.5d0*(1.175e6-t)**(-0.01) + 28.93d0
        n0_init=10**(logn0)
        r0_init=10**(logr0) * aunit
        initialDens_r=1.0+(r/r0_init)**p
        initialDens_r = n0_init/initialDens_r
    end function initialDens_r

    subroutine radiation(r, Lstar, Tstar, Avs, Temp_dust, U)
        real(dp) :: Lstar, Tstar, Avs, r, U_star, U_shell
        real(dp), intent(out) :: Temp_dust, U
        real(8)  :: rsub
        integer, parameter :: nw=129

        ! sublimation distance
        rsub = get_rsub(Lstar)

        ! radiation from the star
        call radiation_star(r, Lstar, Tstar, Avs, U_star)

        ! radiation from the shell
        if (r<rsub) then
            U_shell = 0.0d0
        else
            call radiation_shell(r, Lstar, Tstar, Avs, U_shell)
        end if

        ! total radiation field
        U = U_star + U_shell

        ! dust temperature at equilibrium with the radiation field
        Temp_dust=Temp_average(U)

    end subroutine radiation

    subroutine radiation_star(r, Lstar, Tstar, Avs, U)
        real(dp) :: Lstar, Tstar, Rstar, r
        real(dp), intent(out) :: U
        integer :: i
        real(dp), dimension(:), allocatable :: wave, wave_cm, Istar, uwave_star, tau_wave, uwave_red
        real(dp) :: ZZ, Avs, rsub
        real(dp) :: NH_EBV, RV
        real(8)  :: urad_red
        integer, parameter :: nw=129
        real(dp), dimension(2, nw) :: ext_curves
        character(len=10) :: model

        RV = 4.0d0
        NH_EBV = 5.8d21

        ! sublimation distance
        rsub=get_rsub(Lstar)

        ZZ = HP * C / (K_BOLTZ * Tstar)

        !logspace for wave in micron
        allocate(wave(nw))
        call logspace(log10(wave1), log10(wave2), nw, wave)

        ! convert wave from micron to cm
        wave_cm = wave*1.0d-4  !in cm

        ! Call the function from the module
        call extcurve_obs(wave, RV, NH_EBV, model, ext_curves)

        Istar     = (2.0d0*HP*C**2.0/wave_cm**5.0)*(1.0d0/(exp(ZZ/wave_cm)-1.0d0))
        uwave_star = (4.0d0*PI*wave_cm/C)*(Istar)/wave_cm

        tau_wave = Avs * ext_curves(1,:)/1.086d0
        uwave_red = uwave_star*exp(-tau_wave)

        ! Apply trapezoidal rule for numerical integration
        urad_red = 0.0d0
        do i = 1, nw-1
            urad_red = urad_red + 0.5d0 * (wave_cm(i+1) - wave_cm(i)) * (uwave_red(i+1) + uwave_red(i))
        end do

        ! The stellar radius
        Rstar=Rstar_rsub(Tstar,Tdsub)*rsub

        ! The total radiation field (dimensionless)
        U = urad_red / uISRF * (r / Rstar)**(-2.0)
    end subroutine radiation_star

    subroutine radiation_shell(r, Lstar, Tstar, Avs, U)
        real(dp) :: Lstar, Tstar, r
        real(dp), intent(out) :: U
        integer :: i
        real(dp), dimension(:), allocatable :: wave, wave_cm, Istar, uwave_star, tau_wave, uwave_red
        real(dp) :: ZZ, Avs, rsub
        real(dp) :: NH_EBV, RV
        real(8)  :: urad_red, Tshell
        integer, parameter :: nw=129
        real(dp), dimension(2, nw) :: ext_curves
        character(len=10) :: model

        RV = 4.0d0
        NH_EBV = 5.8d21

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
        wave_cm = wave*1.0d-4  !in cm

        ! Call the function from the module
        call extcurve_obs(wave, RV, NH_EBV, model, ext_curves)

        Istar     = (2.0d0*HP*C**2.0/wave_cm**5.0)*(1.0d0/(exp(ZZ/wave_cm)-1.0d0))
        uwave_star = (4.0d0*PI*wave_cm/C)*(Istar)/wave_cm

        tau_wave = Avs * ext_curves(1,:)/1.086d0
        uwave_red = uwave_star*exp(-tau_wave)

        ! Apply trapezoidal rule for numerical integration
        urad_red = 0.0d0
        do i = 1, nw-1
            urad_red = urad_red + 0.5d0 * (wave_cm(i+1) - wave_cm(i)) * (uwave_red(i+1) + uwave_red(i))
        end do
        U = urad_red / uISRF * (r / rsub)**(-2.0)
    end subroutine radiation_shell

    subroutine logspace(start, stop, num, result)
        real(dp), intent(in) :: start, stop
        integer, intent(in) :: num
        real(dp), dimension(num), intent(out) :: result
        integer :: i

        do i = 1, num
            result(i) = 10.0d0**(start + (i-1)*(stop-start)/DBLE(num-1))
        end do
    end subroutine logspace

    real(dp) function Rstar_rsub(Tstar,Tdmax)
        real(dp) :: Tstar,Tdmax
        real(dp) :: f1,f2

        f1 = sqrt(1.0d6 * Lsun)
        f2 = 155.3d0 * aunit * (Tdmax/1500.0)**(-5.6/2) * sqrt(4*PI*SB_CONST) * Tstar**(2.0)
        Rstar_rsub= f1/f2
    end function Rstar_rsub

    real(dp) function Temp_average(U)
        real(dp) :: U
        real(dp) :: Td_sil,Td_car
        Td_sil = 16.4d0 * U**(1.0/6.0)  ! For silicate grains
        Td_car = 19.5d0 * U**(1.0/5.6)  ! For carbon grains
        Temp_average = Td_sil**(4.0) + Td_car**(4.0)
        Temp_average = (0.5*Temp_average)**(1.0/4.0)
    end function Temp_average

    real(dp) function get_Tshell(Tstar)
        real(dp) :: Tstar

        get_Tshell = Rstar_rsub(Tstar,Tdsub)
        get_Tshell = get_Tshell**(0.5) * Tstar
    end function get_Tshell

    real(dp) function get_rsub(Lstar)
        real(dp) :: Lstar
        get_rsub = 155.3d0*(Lstar/1.0d6/Lsun)**(0.5) * (Tdsub/1500.0d0)**(-5.6/2.0) * aunit  !in cm
    end function get_rsub

    ! Unattenuated UV radiation field from central protostar at radius r_cm [Habing units].
    ! Uses 45% of the bolometric luminosity as the UV fraction and scales as r^-2.
    real(dp) function G0_internal_at_r(Lstar, r_cm)
        real(dp), intent(in) :: Lstar   ! bolometric luminosity [erg s^-1]
        real(dp), intent(in) :: r_cm    ! parcel radius [cm]
        real(dp) :: Luv
        Luv = 0.45_dp * Lstar
        G0_internal_at_r = Luv / (4.0_dp * PI * C * r_cm**2) / uISRF_UV
    end function G0_internal_at_r

end module physicscore
