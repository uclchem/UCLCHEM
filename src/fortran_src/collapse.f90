! Models the chemistry in a collapsing prestellar core
! F. D. Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)
! Uses the following parameterizations of MHD models.
! collapse = 2: Bonnor-Ebert sphere, overdensity factor 1.1 (Aikawa+2005)
! collapse = 3: Bonnor-Ebert sphere, overdensity factor 4 (Aikawa+2005)
! collapse = 4: magnetized filament, initially unstable to collapse (Nakamura+1995)
! collapse = 5: magnetized cloud, initially stable, collapse due to ambipolar diffusion (Fiedler+1993)
module collapse_mod
   use constants
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
                maxTime=1.175d6
                collapseFinalTime=1.173387d6
            case(2)
                maxTime=1.855d5
                collapseFinalTime=1.84265d5
            case(3)
                collapseFinalTime=1.393761d6
            case(4)
                collapseFinalTime=1.6132984d7
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
                  parcelRadius(dstep)=rin*(rout/rin)**((float(dstep)-1.0d0)/(float(points)-1.0d0))
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
            if (timeInYears > 1.0d6) then
                targetTime=(timeInYears+1.0d5)*SECONDS_PER_YEAR    ! 100 kyr steps beyond 1 Myr
            else
                targetTime=(timeInYears+2.5d4)*SECONDS_PER_YEAR    ! 25 kyr steps below 1 Myr
            end if
            return
        end if

        ! ===== REGIME 2: ENDING (remaining time approaching collapseFinalTime) =====
        ! Fine adaptive steps near the singularity to capture rapid density evolution.
        ! The collapse accelerates strongly in the final approach to the singularity,
        ! so step size is reduced as remaining time shrinks.
        if (remaining > 0.0d0 .AND. remaining < 1.0d3) then
            targetTime=(timeInYears+1.0d2)*SECONDS_PER_YEAR    ! 100 yr steps in final 1 kyr
        else if (remaining > 0.0d0 .AND. remaining < 1.0d4) then
            targetTime=(timeInYears+1.0d3)*SECONDS_PER_YEAR    ! 1 kyr steps in final 10 kyr
        else if (remaining > 0.0d0 .AND. remaining < 1.0d5) then
            targetTime=(timeInYears+1.0d4)*SECONDS_PER_YEAR    ! 10 kyr steps in final 100 kyr
        ! ===== REGIME 3: STARTING (early times before ending phase) =====
        else if (timeInYears > 1.0d6) then
            targetTime=(timeInYears+1.0d5)*SECONDS_PER_YEAR    ! 100 kyr steps beyond 1 Myr
        else if (timeInYears > 1.0d5) then
            targetTime=(timeInYears+1.0d4)*SECONDS_PER_YEAR    ! 10 kyr steps beyond 100 kyr
        else if (timeInYears > 10000) then
            targetTime=(timeInYears+1000.0)*SECONDS_PER_YEAR   ! 1 kyr steps beyond 10 kyr
        else if (timeInYears > 1000) then
            targetTime=(timeInYears+100.0)*SECONDS_PER_YEAR    ! 100 yr steps beyond 1 kyr
        else if (timeInYears > 0.0) then
            targetTime=(timeInYears*10)*SECONDS_PER_YEAR       ! *10 exponential early on
        else
            targetTime=3.16d7*10.0d-8
        end if

        ! Cap at collapseFinalTime to avoid overshooting the singularity during collapse phase.
        if (remaining > 0.0d0 .AND. targetTime > collapseFinalTime*SECONDS_PER_YEAR) then
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
        av(dstep)= baseAv +coldens(dstep)/1.6d21
        !If collapse is one of the parameterized modes, find new density and radius

        if ((collapse_mode <= 2)) then
            !I changed rin to rout
            call findNewRadius(massInRadius(dstep),rout,rho0fit(effectiveTime),&
                &r0fit(effectiveTime),afit(effectiveTime),parcelRadius(dstep))
        else
            dt = targetTime - currentTime
            if (timeInYears < collapseFinalTime) then
                drad = vrfit(parcelRadius(dstep),rminfit(effectiveTime),vminfit(effectiveTime),avfit(effectiveTime))*dt/pc
                parcelRadius(dstep) = parcelRadius(dstep) + drad
            end if
        end if
        parcel_radius(dstep) = parcelRadius(dstep)
        density(dstep)=rhofit(parcelRadius(dstep),rho0fit(effectiveTime),r0fit(effectiveTime),afit(effectiveTime))
        ! Apply hard density of n_H=1e8 limit to prevent unphysical behavior
        density(dstep) = MIN(density(dstep), 1e8)
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
        massInRadius(dstep) = 0.0d0

        do i=1,np
           drho = 0.5d0*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
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
      dr = r/1.0d4
      m1 = 0.0d0
      do while (m1 < massInRadius)
         drho = 0.5d0*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
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
      coldens = 0.0d0
      if (size <= 0.0d0) return

      do i=1,np
         r1 = rin + (i-1)*dr
         r2 = rin + i*dr
         drho = 0.5d0*(rhofit(r2,rho0,r0,a)+rhofit(r1,rho0,r0,a))
         coldens = coldens + drho*dr*pc
      end do

    end subroutine findcoldens

! fit to density profile of hydrodynamic simulations
    real(dp) function rhofit(r,rho0,r0,a)
      real(dp) :: r,rho0,r0,a
      real(dp) :: rau,unitrho,unitr,r75

      if (collapse_mode == 1) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      else if (collapse_mode == 2) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      else if (collapse_mode == 3) then
         unitrho = 2.2d4
         unitr = sqrt(1.38d-16*10/2/mh)*(2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! distance unit equal to c_s * (2 pi G rho0)**-1/2
         unitr = unitr/pc
         rhofit = unitrho*rho0/(1+(r/unitr/r0)**2)**a
      else if (collapse_mode == 4) then
         r75 = r/7.5d-1
         rhofit = rho0/(1 + (r75/r0)**a)
      end if

    end function rhofit

! fit to time evolution of central density
    real(dp) function rho0fit(t)
      real(dp) :: t,logrho0,unitt
      if (collapse_mode == 1) then
         logrho0 = 61.8*(maxTime-t)**(-0.01) - 49.4
         rho0fit = 10**logrho0
      else if (collapse_mode == 2) then
         logrho0 = 68.4*(maxTime-t)**(-0.01) - 55.7
         rho0fit = 10**logrho0
      else if (collapse_mode == 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         logrho0 = 3.54*(5.47-t/unitt)**(-0.15) - 2.73
         rho0fit = 10**logrho0
      else if (collapse_mode == 4) then
         if (t <= 6.0d0) then
            rho0fit = 2.0d3 + 1.7d3*(t/6.0 - 1.0)
         else
            logrho0 = 5.3*(16.138-1d-6*t)**(-0.1) - 1.0
            rho0fit = 10**logrho0
         end if
      end if

    end function rho0fit

! fit to time evolution of radius parameter
    real(dp) function r0fit(t)
      real(dp) :: t,logr0,unitt

      if (collapse_mode == 1) then
         logr0 = -28.5*(maxTime-t)**(-0.01) + 28.93
         r0fit = 10**logr0
      else if (collapse_mode == 2) then
         logr0 = -39.0*(maxTime-t)**(-0.01) + 38.7
         r0fit = 10**logr0
      else if (collapse_mode == 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         logr0 = -1.34*(5.47-t/unitt)**(-0.15) + 1.47
         r0fit = 10**logr0
      else if (collapse_mode == 4) then
         logr0 = -2.57*(16.138-1d-6*t)**(-0.1) + 1.85
         r0fit = 10**logr0
      end if

    end function r0fit

! fit to time evolution of density slope parameter
    real(dp) function afit(t)
      real(dp) :: t,unitt

      if (collapse_mode == 1) then
         afit = 2.4d0
      else if (collapse_mode == 2) then
         afit = 1.9 + 0.5*exp(-t/1e5)
      else if (collapse_mode == 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         afit = 2.0 - 0.5*(t/unitt/5.47)**9
      else if (collapse_mode == 4) then
         afit = 2.4 - 0.2*(1d-6*t/16.138)**40
      end if

    end function afit

! fit to radial velocity of hydrodynamical simulation
    real(dp) function vrfit(r,rmin,vmin,a)
      real(dp) :: r,rmin,vmin,a
      real(dp) :: unitr,newRadius,rmid,r75

      if (collapse_mode == 3) then
         unitr = sqrt(1.38d-16*10/2/mh)*(2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! distance unit equal to c_s * (2 pi G rho0)**-1/2
         unitr = unitr/pc
         newRadius = r/unitr - rmin
         if (newRadius < 0.0d0) then
            vrfit = vmin*((newRadius/rmin)**2 -1)
         else
            vrfit = vmin*(exp(-2.0d0*a*newRadius) - 2*exp(-a*newRadius))
         end if
         vrfit = sqrt(1.38d-16*10/2/mh)*vrfit  ! convert to cm s-1 using c_s
      else if (collapse_mode == 4) then
         rmid = 0.5
         r75 = r/7.5d-1
         newRadius = r75 - rmin
         if (r75 < rmin) then
            vrfit = vmin*((newRadius/rmin)**2 - 1)
         else if (r75 <= rmid) then
            vrfit = (vmin-a)*(newRadius/(rmid-rmin))**0.3 - vmin
         else
            vrfit = a/(1.0-rmid)*(r75-rmid) - a
         end if
         vrfit = 1d3*vrfit  ! convert to cm s-1 from 1e-2 km s-1
      end if

    end function vrfit

! fit to time evolution of radius of minimum velocity
    real(dp) function rminfit(t)
      real(dp) :: t,unitt,tnew
      real(dp) :: t6

      if (collapse_mode == 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         if (tnew == 0.0d0) then
            rminfit = 7.2d0
         else if (log(tnew) < 1.6d0) then
            rminfit = -1.149*tnew + 7.2
         else if (log(tnew) < 1.674d0) then
            rminfit = -9.2*log(tnew) + 16.25
         else
            rminfit = -22.0*log(tnew) + 37.65
         end if
      else if (collapse_mode == 4) then
         t6 = 1d-6*t
         if (t6 <= 10.2) then
            rminfit = -0.0039*t6 + 0.49
         else if (t6 <= 15.1) then
            rminfit = -0.0306*(t6-10.2) + 0.45
         else
            rminfit = -0.282*(t6-15.1) + 0.3
         end if
      end if

    end function rminfit

! fit to time evolution of minimum velocity
    real(dp) function vminfit(t)
      real(dp) :: t,unitt,tnew
      real(dp) :: t6

      if (collapse_mode == 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         if (tnew == 0.0d0) then
            vminfit = 0.0d0
         else if (log(tnew) < 1.6d0) then
            vminfit = 0.0891*tnew
         else if (log(tnew) < 1.674d0) then
            vminfit = 5.5*log(tnew) - 8.37
         else
            vminfit = 18.9*log(tnew) - 30.8
         end if
      else if (collapse_mode == 4) then
         t6 = 1d-6*t
         vminfit = 3.44*(16.138-t6)**(-0.35) - 0.7
      end if

    end function vminfit

! fit to time evolution of velocity a-parameter (collapse 4) or velocity at r=0.5 (collapse 5)
    real(dp) function avfit(t)
      real(dp) :: t,unitt,tnew
      real(dp) :: t6

      if (collapse_mode == 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5)  ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         if (tnew == 0.0d0) then
            avfit = 0.4d0
         else if (log(tnew) < 1.6d0) then
            avfit = 0.0101*tnew + 0.4
         else if (log(tnew) < 1.674d0) then
            avfit = 0.695*log(tnew) - 0.663
         else
            avfit = 2.69*log(tnew) - 4.0
         end if
      else if (collapse_mode == 4) then
         t6 = 1d-6*t
         if (t6 <= 10.2) then
            avfit = 0.143*t6
         else
            avfit = 0.217*(t6-10.2) + 1.46
         end if
      end if

    end function avfit
end module collapse_mod

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
