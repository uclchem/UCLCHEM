! Models the chemistry in a collapsing prestellar core
! F. D. Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)
! Uses the following parameterizations of MHD models.
! collapse = 2: Bonnor-Ebert sphere, overdensity factor 1.1 (Aikawa+2005)
! collapse = 3: Bonnor-Ebert sphere, overdensity factor 4 (Aikawa+2005)
! collapse = 4: magnetized filament, initially unstable to collapse (Nakamura+1995)
! collapse = 5: magnetized cloud, initially stable, collapse due to ambipolar diffusion (Fiedler+1993)
module collapse_mod
   use constants, only: dp, MH, PI, SECONDS_PER_YEAR, PC, AU
   use DEFAULTPARAMETERS
   use F2PY_CONSTANTS
   use network
   !f2py INTEGER, parameter :: dp
   use physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
   & zeta, currentTime, currentTimeold, targetTime, timeinyears, freefall, density, ion, densdot, gastemp, dusttemp, av, &
   &coldens, parcel_radius
   implicit none

   integer :: collapse_mode
   real(dp) :: maxTime
   real(dp) :: collapseFinalTime
   real(dp), allocatable :: massInRadius(:),parcelRadius(:)
   real(dp) :: dt,drad
contains

    subroutine initializePhysics(successFlag)
        !f2py integer, intent(aux) :: points
        integer, intent(out) :: successFlag

        if (ALLOCATED(parcelRadius)) deallocate(parcelRadius,massInRadius)
        allocate(parcelRadius(points),massInRadius(points))

         select case(collapse_mode)
            case(1)
                maxTime=1.175e6_dp
                collapseFinalTime=1.173387e6_dp
            case(2)
                maxTime=1.855e5_dp
                collapseFinalTime=1.84265e5_dp
            case(3)
                collapseFinalTime=1.393761e6_dp
            case(4)
                collapseFinalTime=1.6132984e7_dp
            case default
                write(*,*) "unacceptable collapse mode"
                successFlag=-1
                return
         end select

         if (points == 1) then
               parcelRadius(1)=rout
               parcel_radius(1)=rout
         else
            do dstep=1,points
                  parcelRadius(dstep)=rin*(rout/rin)**((float(dstep)-1.0_dp)/(float(points)-1.0_dp))
                  parcel_radius(dstep)=parcelRadius(dstep)
            end do
         end if

         density=rhofit(parcelRadius(dstep),rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears))
         if (collapse_mode <= 2) call findmassInRadius
    end subroutine initializePhysics

    subroutine updateTargetTime
        real(dp) :: remaining
        remaining = collapseFinalTime - timeInYears

        ! ===== REGIME 1: POST-COLLAPSE (timeInYears > collapseFinalTime) =====
        ! After the collapse endpoint, density is frozen. Use coarse timesteps to advance to finalTime.
        if (timeInYears > collapseFinalTime) then
            if (timeInYears > 1.0e6_dp) then
                targetTime=(timeInYears+1.0e5_dp)*SECONDS_PER_YEAR    ! 100 kyr steps beyond 1 Myr
            else
                targetTime=(timeInYears+2.5e4_dp)*SECONDS_PER_YEAR    ! 25 kyr steps below 1 Myr
            end if
            return
        end if

        ! ===== REGIME 2: ENDING (remaining time approaching collapseFinalTime) =====
        ! Fine adaptive steps near the singularity to capture raPId density evolution.
        ! The collapse accelerates strongly in the final approach to the singularity,
        ! so step size is reduced as remaining time shrinks.
        if (remaining > 0.0_dp .AND. remaining < 1.0e3_dp) then
            targetTime=(timeInYears+1.0e2_dp)*SECONDS_PER_YEAR    ! 100 yr steps in final 1 kyr
        else if (remaining > 0.0_dp .AND. remaining < 1.0e4_dp) then
            targetTime=(timeInYears+1.0e3_dp)*SECONDS_PER_YEAR    ! 1 kyr steps in final 10 kyr
        else if (remaining > 0.0_dp .AND. remaining < 1.0e5_dp) then
            targetTime=(timeInYears+1.0e4_dp)*SECONDS_PER_YEAR    ! 10 kyr steps in final 100 kyr
        ! ===== REGIME 3: STARTING (early times before ending phase) =====
        else if (timeInYears > 1.0e6_dp) then
            targetTime=(timeInYears+1.0e5_dp)*SECONDS_PER_YEAR    ! 100 kyr steps beyond 1 Myr
        else if (timeInYears > 1.0e5_dp) then
            targetTime=(timeInYears+1.0e4_dp)*SECONDS_PER_YEAR    ! 10 kyr steps beyond 100 kyr
        else if (timeInYears > 10000) then
            targetTime=(timeInYears+1000.0_dp)*SECONDS_PER_YEAR   ! 1 kyr steps beyond 10 kyr
        else if (timeInYears > 1000) then
            targetTime=(timeInYears+100.0_dp)*SECONDS_PER_YEAR    ! 100 yr steps beyond 1 kyr
        else if (timeInYears > 0.0_dp) then
            targetTime=(timeInYears*10)*SECONDS_PER_YEAR       ! *10 exponential early on
        else
            targetTime=3.16e7_dp*10.0e-8_dp
        end if

        ! Cap at collapseFinalTime to avoid overshooting the singularity during collapse phase.
        if (remaining > 0.0_dp .AND. targetTime > collapseFinalTime*SECONDS_PER_YEAR) then
          targetTime=collapseFinalTime*SECONDS_PER_YEAR
        end if
    end subroutine updateTargetTime

    !This routine is formed for every parcel at every time step.
    !update any physics here. For example, set density
    subroutine updatePhysics
         !f2py integer, intent(aux) :: points
        real(dp) :: effectiveTime
        !Freeze density and parcel radius at collapseFinalTime; chemistry continues beyond
        effectiveTime = MIN(timeInYears, collapseFinalTime)
        !calculate column density. Remember dstep counts from core to edge
        !and coldens should be amount of gas from edge to parcel.
        call findcoldens(coldens(dstep),rin,rho0fit(effectiveTime),r0fit(effectiveTime),afit(effectiveTime),rout)
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6e21_dp
        !If collapse is one of the parameterized modes, find new density and radius

        if ((collapse_mode <= 2)) then
            !I changed rin to rout
            call findNewRadius(massInRadius(dstep),rout,rho0fit(effectiveTime),&
                &r0fit(effectiveTime),afit(effectiveTime),parcelRadius(dstep))
        else
            dt = targetTime - currentTime
            if (timeInYears < collapseFinalTime) then
                drad = vrfit(parcelRadius(dstep),rminfit(effectiveTime),vminfit(effectiveTime),avfit(effectiveTime))*dt/PC
                parcelRadius(dstep) = parcelRadius(dstep) + drad
            end if
        end if
        parcel_radius(dstep) = parcelRadius(dstep)
        density(dstep)=rhofit(parcelRadius(dstep),rho0fit(effectiveTime),r0fit(effectiveTime),afit(effectiveTime))
        ! Apply hard density of n_H=1e8 limit to prevent unphysical behavior
        density(dstep) = MIN(density(dstep), 1e8_dp)
    end subroutine updatePhysics

    !This module is isothermal and as such, no sublimation occurs.
    !This is a dummy subroutine.
    subroutine sublimation(abund, lpoints)
      real(dp), intent(inout) :: abund(nspec+1,lpoints)
      integer, intent(in) :: lpoints
    end subroutine sublimation


    ! finds initial mass within starting radius, assuming spherical symmetry
    subroutine findMassInRadius
      !f2py integer, intent(aux) :: points
      real(dp) :: rho0,r0,a
      integer :: i,np,dstep
      real(dp) :: dr,drho

        rho0=rho0fit(timeInYears)
        r0=r0fit(timeInYears)
        a=afit(timeInYears)
      do dstep=1,points
        np = 1000
        dr = parcelRadius(dstep)/np
        massInRadius(dstep) = 0.0_dp

        do i=1,np
           drho = 0.5_dp*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
           massInRadius(dstep) = massInRadius(dstep) + drho*dr*(i*dr)**2
        end do
      end do
    end subroutine findMassInRadius

! finds radius enclosing a mass of massInRadius
    subroutine findNewRadius(massInRadius,r,rho0,r0,a,newRadius)
      real(dp),intent(in) :: massInRadius,r,rho0,r0,a
      real(dp),intent(out) :: newRadius
      integer :: i
      real(dp) :: dr,drho,m1

      i=1
      dr = r/1.0e4_dp
      m1 = 0.0_dp
      do while (m1 < massInRadius)
         drho = 0.5_dp*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
         m1 = m1 + drho*dr*(i*dr)**2
         newRadius = i*dr
         i=i+1
      end do

    end subroutine findNewRadius

! finds column density to edge of cloud based on density profile
    subroutine findcoldens(coldens,rin,rho0,r0,a,rout)
      real(dp),intent(in) :: rin,rout,rho0,r0,a
      real(dp),intent(out) :: coldens
      integer :: i,np
      real(dp) :: dr,drho,size,r1,r2

      np = 10000
      size = rout-rin
      dr = size/np
      coldens = 0.0_dp
      if (size <= 0.0_dp) return

      do i=1,np
         r1 = rin + (i-1)*dr
         r2 = rin + i*dr
         drho = 0.5_dp*(rhofit(r2,rho0,r0,a)+rhofit(r1,rho0,r0,a))
         coldens = coldens + drho*dr*PC
      end do

    end subroutine findcoldens

! fit to density profile of hydrodynamic simulations
    real(dp) function rhofit(r,rho0,r0,a)
      real(dp), intent(in) :: r,rho0,r0,a
      real(dp) :: rau,unitrho,unitr,r75

      if (collapse_mode == 1) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      else if (collapse_mode == 2) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      else if (collapse_mode == 3) then
         unitrho = 2.2e4_dp
         unitr = sqrt(1.38e-16_dp*10/2/MH)*(2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! distance unit equal to c_s * (2 PI G rho0)**-1/2
         unitr = unitr/PC
         rhofit = unitrho*rho0/(1+(r/unitr/r0)**2)**a
      else if (collapse_mode == 4) then
         r75 = r/7.5e-1_dp
         rhofit = rho0/(1.0_dp + (r75/r0)**a)
      end if

    end function rhofit

! fit to time evolution of central density
    real(dp) function rho0fit(t)
      real(dp), intent(in) :: t
      real(dp) :: logrho0, unitt
      if (collapse_mode == 1) then
         logrho0 = 61.8_dp*(maxTime-t)**(-0.01_dp) - 49.4_dp
         rho0fit = 10**logrho0
      else if (collapse_mode == 2) then
         logrho0 = 68.4_dp*(maxTime-t)**(-0.01_dp) - 55.7_dp
         rho0fit = 10**logrho0
      else if (collapse_mode == 3) then
         unitt = (2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! time unit equal to (2 PI G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         logrho0 = 3.54_dp*(5.47_dp-t/unitt)**(-0.15_dp) - 2.73_dp
         rho0fit = 10**logrho0
      else if (collapse_mode == 4) then
         if (t <= 6.0_dp) then
            rho0fit = 2.0e3_dp + 1.7e3_dp*(t/6.0_dp - 1.0_dp)
         else
            logrho0 = 5.3_dp*(16.138_dp-1e-6_dp*t)**(-0.1_dp) - 1.0_dp
            rho0fit = 10**logrho0
         end if
      end if

    end function rho0fit

! fit to time evolution of radius parameter
    real(dp) function r0fit(t)
      real(dp), intent(in) :: t
      real(dp) :: logr0,unitt

      if (collapse_mode == 1) then
         logr0 = -28.5_dp*(maxTime-t)**(-0.01_dp) + 28.93_dp
         r0fit = 10**logr0
      else if (collapse_mode == 2) then
         logr0 = -39.0_dp*(maxTime-t)**(-0.01_dp) + 38.7_dp
         r0fit = 10**logr0
      else if (collapse_mode == 3) then
         unitt = (2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! time unit equal to (2 PI G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         logr0 = -1.34_dp*(5.47_dp-t/unitt)**(-0.15_dp) + 1.47_dp
         r0fit = 10**logr0
      else if (collapse_mode == 4) then
         logr0 = -2.57_dp*(16.138_dp-1e-6_dp*t)**(-0.1_dp) + 1.85_dp
         r0fit = 10**logr0
      end if

    end function r0fit

! fit to time evolution of density slope parameter
    real(dp) function afit(t)
      real(dp), intent(in) :: t
      real(dp) :: unitt

      if (collapse_mode == 1) then
         afit = 2.4_dp
      else if (collapse_mode == 2) then
         afit = 1.9_dp + 0.5_dp*exp(-t/1e5_dp)
      else if (collapse_mode == 3) then
         unitt = (2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! time unit equal to (2 PI G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         afit = 2.0_dp - 0.5_dp*(t/unitt/5.47_dp)**9
      else if (collapse_mode == 4) then
         afit = 2.4_dp - 0.2_dp*(1e-6_dp*t/16.138_dp)**40
      end if

    end function afit

! fit to radial velocity of hydrodynamical simulation
    real(dp) function vrfit(r,rmin,vmin,a)
      real(dp), intent(in) :: r,rmin,vmin,a
      real(dp) :: unitr,newRadius,rmid,r75

      if (collapse_mode == 3) then
         unitr = sqrt(1.38e-16_dp*10/2/MH)*(2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! distance unit equal to c_s * (2 PI G rho0)**-1/2
         unitr = unitr/PC
         newRadius = r/unitr - rmin
         if (newRadius < 0.0_dp) then
            vrfit = vmin*((newRadius/rmin)**2 -1)
         else
            vrfit = vmin*(exp(-2.0_dp*a*newRadius) - 2*exp(-a*newRadius))
         end if
         vrfit = sqrt(1.38e-16_dp*10/2/MH)*vrfit  ! convert to cm s-1 using c_s
      else if (collapse_mode == 4) then
         rmid = 0.5_dp
         r75 = r/7.5e-1_dp
         newRadius = r75 - rmin
         if (r75 < rmin) then
            vrfit = vmin*((newRadius/rmin)**2 - 1.0_dp)
         else if (r75 <= rmid) then
            vrfit = (vmin-a)*(newRadius/(rmid-rmin))**0.3_dp - vmin
         else
            vrfit = a/(1.0_dp-rmid)*(r75-rmid) - a
         end if
         vrfit = 1e3_dp*vrfit  ! convert to cm s-1 from 1e-2 km s-1
      end if

    end function vrfit

! fit to time evolution of radius of minimum velocity
    real(dp) function rminfit(t)
      real(dp), intent(in) :: t
      real(dp) :: unitt,tnew
      real(dp) :: t6

      if (collapse_mode == 3) then
         unitt = (2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! time unit equal to (2 PI G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         if (tnew == 0.0_dp) then
            rminfit = 7.2_dp
         else if (log(tnew) < 1.6_dp) then
            rminfit = -1.149_dp*tnew + 7.2_dp
         else if (log(tnew) < 1.674_dp) then
            rminfit = -9.2_dp*log(tnew) + 16.25_dp
         else
            rminfit = -22.0_dp*log(tnew) + 37.65_dp
         end if
      else if (collapse_mode == 4) then
         t6 = 1e-6_dp*t
         if (t6 <= 10.2_dp) then
            rminfit = -0.0039_dp*t6 + 0.49_dp
         else if (t6 <= 15.1_dp) then
            rminfit = -0.0306_dp*(t6-10.2_dp) + 0.45_dp
         else
            rminfit = -0.282_dp*(t6-15.1_dp) + 0.3_dp
         end if
      end if

    end function rminfit

! fit to time evolution of minimum velocity
    real(dp) function vminfit(t)
      real(dp), intent(in) :: t
      real(dp) :: unitt,tnew
      real(dp) :: t6

      if (collapse_mode == 3) then
         unitt = (2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! time unit equal to (2 PI G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         if (tnew == 0.0_dp) then
            vminfit = 0.0_dp
         else if (log(tnew) < 1.6_dp) then
            vminfit = 0.0891_dp*tnew
         else if (log(tnew) < 1.674_dp) then
            vminfit = 5.5_dp*log(tnew) - 8.37_dp
         else
            vminfit = 18.9_dp*log(tnew) - 30.8_dp
         end if
      else if (collapse_mode == 4) then
         t6 = 1e-6_dp*t
         vminfit = 3.44_dp*(16.138_dp-t6)**(-0.35_dp) - 0.7_dp
      end if

    end function vminfit

! fit to time evolution of velocity a-parameter (collapse 4) or velocity at r=0.5 (collapse 5)
    real(dp) function avfit(t)
      real(dp), intent(in) :: t
      real(dp) :: unitt,tnew
      real(dp) :: t6

      if (collapse_mode == 3) then
         unitt = (2*PI*6.67e-8_dp*2.2e4_dp*MH)**(-0.5_dp)  ! time unit equal to (2 PI G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         if (tnew == 0.0_dp) then
            avfit = 0.4_dp
         else if (log(tnew) < 1.6_dp) then
            avfit = 0.0101_dp*tnew + 0.4_dp
         else if (log(tnew) < 1.674_dp) then
            avfit = 0.695_dp*log(tnew) - 0.663_dp
         else
            avfit = 2.69_dp*log(tnew) - 4.0_dp
         end if
      else if (collapse_mode == 4) then
         t6 = 1e-6_dp*t
         if (t6 <= 10.2_dp) then
            avfit = 0.143_dp*t6
         else
            avfit = 0.217_dp*(t6-10.2_dp) + 1.46_dp
         end if
      end if

    end function avfit
end module collapse_mod

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
