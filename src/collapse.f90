!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics module from Priestley et al. 2018, models the chemistry in a collapsing prestellar core using parameterizations
! of MHD models. To run, use collapse=1 to go from a density of 1e2/cm3 to your initial core density and create an abundance file
! (readabunds=0). This is the phase 1 described in most UCLCHEM publications and the manual.
! Then run the code again with (readabunds=1), setting rout to the radius of the core and points to the number of points from the core
! centre to rout that you want to model. Then set collapse to a value from 2-5 depending on the model you wish to use.
!
! collapse modes > 1 are parameterizations:
! collapse = 2: Bonnor-Ebert sphere, overdensity factor 1.1 (Aikawa+2005)
! collapse = 3: Bonnor-Ebert sphere, overdensity factor 4 (Aikawa+2005)
! collapse = 4: magnetised filament, initially unstable to collapse (Nakamura+1995)
! collapse = 5: magnetised cloud, initially stable, collapse due to ambipolar diffusion (Fiedler+1993)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE physics
    IMPLICIT NONE
    INTEGER :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    INTEGER :: collapse,switch,first,phase
    INTEGER :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/DOne
    INTEGER :: evap,ion,solidflag,volcflag,coflag,tempindx
    
    !variables either controlled by physics or that user may wish to change    
    DOUBLE PRECISION :: initialDens,timeInYears,targetTime,currentTime,currentTimeold,finalDens,finalTime,grainRadius,initialTemp
    DOUBLE PRECISION :: cloudSize,rout,rin,baseAv,bc,olddens,maxTemp,vs,maxTime
    DOUBLE PRECISION :: tempa(5),tempb(5),codestemp(5),volctemp(5),solidtemp(5)
    DOUBLE PRECISION, allocatable :: av(:),coldens(:),temp(:),dens(:),massInRadius(:),parcelRadius(:)
    !Everything should be in cgs units. Helpful constants and conversions below
    DOUBLE PRECISION,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23
    DOUBLE PRECISION, parameter :: year=3.16455d-08,pc=3.086d18,au=2.063d5

    DOUBLE PRECISION :: dt,drad


CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.

    
    SUBROUTINE initializePhysics
    !Any initialisation logic steps go here
    !cloudSize is important as is allocating space for depth arrays
      IF (ALLOCATED(av)) DEALLOCATE(av,coldens,temp,dens)
      ALLOCATE(av(points),coldens(points),dens(points),temp(points),parcelRadius(points),massInRadius(points))
      cloudSize=(rout-rin)*pc

      temp = initialTemp

      SELECT CASE(collapse)
        CASE(0) !no collapse
          dens=initialDens
        CASE(1) !standard freefall collapse
          dens=1.001*initialDens
        CASE DEFAULT !collapse modes detailed at top of module
          phase=2
          SELECT CASE(collapse)
            CASE(2)
              maxTime=1.175d6
              finalTime=0.97*maxTime
            CASE(3) 
              maxTime=1.855d5
              finalTime=0.97*maxTime
            CASE DEFAULT
              write(*,*) "unacceptable collapse mode"
              switch=0
              finalTime=0.0
          END SELECT

          DO dstep=1,points
            parcelRadius(dstep)=dstep*rout/float(points)
          END DO
               
          open(unit=66,file='output/radius.dat',status='unknown')
          dens=rhofit(rin,rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears))
          IF ((collapse .eq. 2) .or. (collapse .eq. 3)) then
            CALL findmassInRadius
          END IF
      END SELECT
    END SUBROUTINE initializePhysics

    SUBROUTINE updateTargetTime
    !At each timestep, the time at the END of the step is calculated by calling this function
    !You need to set targetTime in seconds, its initial value will be the time at start of current step
    IF (timeInYears .gt. 10000) THEN
        targetTime=(timeInYears+1000.0)/year
    ELSE IF (timeInYears .gt. 1000) THEN
        targetTime=(timeInYears+100.0)/year
    ELSE IF (timeInYears .gt. 0.0) THEN
        targetTime=(timeInYears*10)/year
    ELSE
        targetTime=3.16d7*10.d-8
    ENDIF

    IF (targetTime .gt. finalTime/year)targetTime=finalTime/year

    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrator.
    ! DLSODE sorts that out based on chosen error tolerances (RTOL/ATOL) and is simply called repeatedly
    !until it outputs a time >= targetTime.
    END SUBROUTINE updateTargetTime

    !This routine is formed for every parcel at every time step.
    !update any physics here. For example, set density
    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from core to edge
        !and coldens should be amount of gas from edge to parcel.
        call findcoldens(coldens(dstep),rin,rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears),rout)
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21

        !If collapse is one of the parameterized modes, find new density and radius
        IF ((collapse .gt. 1).and.(phase .eq.2)) THEN
           IF ((collapse .eq. 2) .or. (collapse .eq. 3)) THEN
              CALL findNewRadius(massInRadius(dstep),rin,rho0fit(timeInYears),&
                &r0fit(timeInYears),afit(timeInYears),parcelRadius(dstep))
           ELSE IF ((collapse .eq. 4) .or. (collapse .eq. 5)) THEN
              dt = targetTime - currentTime
              drad = vrfit(parcelRadius(dstep),rminfit(timeInYears),vminfit(timeInYears),avfit(timeInYears))*dt/pc
              parcelRadius(dstep) = parcelRadius(dstep) + drad
              write(66,*) timeInYears,parcelRadius(dstep),rhofit(parcelRadius(dstep),&
                &rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears)),&
              &vrfit(parcelRadius(dstep),rminfit(timeInYears),vminfit(timeInYears),avfit(timeInYears))
           END IF
           dens(dstep)=rhofit(parcelRadius(dstep),rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears))
        END IF
        
    END SUBROUTINE updatePhysics

    !This function works out the time derivative of the density, allowing DVODE to update density with the rest of our ODEs
    !It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DVODE
    !Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot(density)
        DOUBLE PRECISION, INTENT(IN) :: density
        DOUBLE PRECISION :: densdot
        !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (density .lt. finalDens) THEN
             densdot=bc*(density**4./initialDens)**0.33*&
             &(8.4d-30*initialDens*((density/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=0.0
        ENDIF
    END FUNCTION densdot

    ! finds initial mass within starting radius, assuming spherical symmetry
    SUBROUTINE findMassInRadius
      DOUBLE PRECISION :: rho0,r0,a
      INTEGER :: i,np,dstep
      DOUBLE PRECISION :: dr,drho

        rho0=rho0fit(timeInYears)
        r0=r0fit(timeInYears)
        a=afit(timeInYears)
      DO dstep=1,points
        np = 1000
        dr = parcelRadius(dstep)/np
        massInRadius(dstep) = 0.0d0

        DO i=1,np
           drho = 0.5d0*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
           massInRadius(dstep) = massInRadius(dstep) + drho*dr*(i*dr)**2
        END DO
      END DO
    END SUBROUTINE findMassInRadius

! finds radius enclosing a mass of massInRadius
    SUBROUTINE findNewRadius(massInRadius,r,rho0,r0,a,newRadius)
      DOUBLE PRECISION,intent(in) :: massInRadius,r,rho0,r0,a
      DOUBLE PRECISION,intent(out) :: newRadius
      INTEGER :: i
      DOUBLE PRECISION :: dr,drho,m1

      i=1
      dr = r/1.0d4
      m1 = 0.0d0
      DO WHILE (m1 .lt. massInRadius)
         drho = 0.5d0*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
         m1 = m1 + drho*dr*(i*dr)**2
         newRadius = i*dr
         i=i+1
      END DO
      write(66,*) timeInYears,newRadius,rhofit(newRadius,rho0,r0,a),m1

    END SUBROUTINE findNewRadius

! finds column density to edge of cloud based on density profile
    SUBROUTINE findcoldens(coldens,rin,rho0,r0,a,rout)
      DOUBLE PRECISION,intent(in) :: rin,rout,rho0,r0,a
      DOUBLE PRECISION,intent(out) :: coldens
      INTEGER :: i,np
      DOUBLE PRECISION :: dr,drho,size,r1,r2

      np = 10000
      size = rout-rin
      dr = size/np
      coldens = 0.0d0
      IF (size .le. 0.0d0) return

      DO i=1,np
         r1 = rin + (i-1)*dr
         r2 = rin + i*dr
         drho = 0.5d0*(rhofit(r2,rho0,r0,a)+rhofit(r1,rho0,r0,a))
         coldens = coldens + drho*dr*pc
      END DO

    END SUBROUTINE findcoldens

! fit to density profile of hydrodynamic simulations
    DOUBLE PRECISION function rhofit(r,rho0,r0,a)
      DOUBLE PRECISION :: r,rho0,r0,a
      DOUBLE PRECISION :: rau,unitrho,unitr,r75

      IF (collapse .eq. 2) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      ELSE IF (collapse .eq. 3) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      ELSE IF (collapse .eq. 4) then
         unitrho = 2.2d4
         unitr = sqrt(1.38d-16*10/2/mh)*(2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! distance unit equal to c_s * (2 pi G rho0)**-1/2
         unitr = unitr/pc
         rhofit = unitrho*rho0/(1+(r/unitr/r0)**2)**a
      ELSE IF (collapse .eq. 5) then
         r75 = r/7.5d-1
         rhofit = rho0/(1 + (r75/r0)**a)
      END IF

    END function rhofit

! fit to time evolution of central density
    DOUBLE PRECISION function rho0fit(t)
      DOUBLE PRECISION :: t,logrho0,unitt
      IF (collapse .eq. 2) then
         logrho0 = 61.8*(maxTime-t)**(-0.01) - 49.4
         rho0fit = 10**logrho0
      ELSE IF (collapse .eq. 3) then
         logrho0 = 68.4*(maxTime-t)**(-0.01) - 55.7
         rho0fit = 10**logrho0
      ELSE IF (collapse .eq. 4) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt*year
         logrho0 = 3.54*(5.47-t/unitt)**(-0.15) - 2.73
         rho0fit = 10**logrho0
      ELSE IF (collapse .eq. 5) then
         IF (t .le. 6.0d0) then
            rho0fit = 2.0d3 + 1.7d3*(t/6.0 - 1.0)
         ELSE
            logrho0 = 5.3*(16.138-1d-6*t)**(-0.1) - 1.0
            rho0fit = 10**logrho0
         END IF
      END IF

    END function rho0fit

! fit to time evolution of radius parameter
    DOUBLE PRECISION function r0fit(t)
      DOUBLE PRECISION :: t,logr0,unitt

      IF (collapse .eq. 2) then
         logr0 = -28.5*(maxTime-t)**(-0.01) + 28.93
         r0fit = 10**logr0
      ELSE IF (collapse .eq. 3) then
         logr0 = -39.0*(maxTime-t)**(-0.01) + 38.7
         r0fit = 10**logr0
      ELSE IF (collapse .eq. 4) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt*year
         logr0 = -1.34*(5.47-t/unitt)**(-0.15) + 1.47
         r0fit = 10**logr0
      ELSE IF (collapse .eq. 5) then
         logr0 = -2.57*(16.138-1d-6*t)**(-0.1) + 1.85
         r0fit = 10**logr0
      END IF

    END function r0fit

! fit to time evolution of density slope parameter
    DOUBLE PRECISION function afit(t)
      DOUBLE PRECISION :: t,unitt

      IF (collapse .eq. 2) then
         afit = 2.4d0
      ELSE IF (collapse .eq. 3) then
         afit = 1.9 + 0.5*exp(-t/1e5)
      ELSE IF (collapse .eq. 4) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt*year
         afit = 2.0 - 0.5*(t/unitt/5.47)**9
      ELSE IF (collapse .eq. 5) then
         afit = 2.4 - 0.2*(1d-6*t/16.138)**40
      END IF

    END function afit

! fit to radial velocity of hydrodynamical simulation
    DOUBLE PRECISION function vrfit(r,rmin,vmin,a)
      DOUBLE PRECISION :: r,rmin,vmin,a
      DOUBLE PRECISION :: unitr,newRadius,rmid,r75

      IF (collapse .eq. 4) then
         unitr = sqrt(1.38d-16*10/2/mh)*(2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! distance unit equal to c_s * (2 pi G rho0)**-1/2
         unitr = unitr/pc
         newRadius = r/unitr - rmin
         IF (newRadius .lt. 0.0d0) then
            vrfit = vmin*((newRadius/rmin)**2 -1)
         ELSE
            vrfit = vmin*(exp(-2.0d0*a*newRadius) - 2*exp(-a*newRadius))
         END IF
         vrfit = sqrt(1.38d-16*10/2/mh)*vrfit ! convert to cm s-1 using c_s
      ELSE IF (collapse .eq. 5) then
         rmid = 0.5
         r75 = r/7.5d-1
         newRadius = r75 - rmin
         IF (r75 .lt. rmin) then
            vrfit = vmin*((newRadius/rmin)**2 - 1)
         ELSE IF (r75 .le. rmid) then
            vrfit = (vmin-a)*(newRadius/(rmid-rmin))**0.3 - vmin
         ELSE
            vrfit = a/(1.0-rmid)*(r75-rmid) - a
         END IF
         vrfit = 1d3*vrfit ! convert to cm s-1 from 1e-2 km s-1
      END IF

    END function vrfit

! fit to time evolution of radius of minimum velocity
    DOUBLE PRECISION function rminfit(t)
      DOUBLE PRECISION :: t,unitt,tnew
      DOUBLE PRECISION :: t6

      IF (collapse .eq. 4) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt*year
         tnew = t/unitt
         IF (tnew .eq. 0.0d0) then
            rminfit = 7.2d0
         ELSE IF (log(tnew) .lt. 1.6d0) then
            rminfit = -1.149*tnew + 7.2
         ELSE IF (log(tnew) .lt. 1.674d0) then
            rminfit = -9.2*log(tnew) + 16.25
         ELSE
            rminfit = -22.0*log(tnew) + 37.65
         END IF
      ELSE IF (collapse .eq. 5) then
         t6 = 1d-6*t
         IF (t6 .le. 10.2) then
            rminfit = -0.0039*t6 + 0.49
         ELSE IF (t6 .le. 15.1) then
            rminfit = -0.0306*(t6-10.2) + 0.45
         ELSE
            rminfit = -0.282*(t6-15.1) + 0.3
         END IF
      END IF

    END function rminfit

! fit to time evolution of minimum velocity
    DOUBLE PRECISION function vminfit(t)
      DOUBLE PRECISION :: t,unitt,tnew
      DOUBLE PRECISION :: t6

      IF (collapse .eq. 4) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt*year
         tnew = t/unitt
         IF (tnew .eq. 0.0d0) then
            vminfit = 0.0d0
         ELSE IF (log(tnew) .lt. 1.6d0) then
            vminfit = 0.0891*tnew
         ELSE IF (log(tnew) .lt. 1.674d0) then
            vminfit = 5.5*log(tnew) - 8.37
         ELSE
            vminfit = 18.9*log(tnew) - 30.8
         END IF
      ELSE IF (collapse .eq. 5) then
         t6 = 1d-6*t
         vminfit = 3.44*(16.138-t6)**(-0.35) - 0.7
      END IF

    END function vminfit

! fit to time evolution of velocity a-parameter (collapse 4) or velocity at r=0.5 (collapse 5)
    DOUBLE PRECISION function avfit(t)
      DOUBLE PRECISION :: t,unitt,tnew
      DOUBLE PRECISION :: t6

      IF (collapse .eq. 4) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt*year
         tnew = t/unitt
         IF (tnew .eq. 0.0d0) then
            avfit = 0.4d0
         ELSE IF (log(tnew) .lt. 1.6d0) then
            avfit = 0.0101*tnew + 0.4
         ELSE IF (log(tnew) .lt. 1.674d0) then
            avfit = 0.695*log(tnew) - 0.663
         ELSE
            avfit = 2.69*log(tnew) - 4.0
         END IF
      ELSE IF (collapse .eq. 5) then
         t6 = 1d-6*t
         IF (t6 .le. 10.2) then
            avfit = 0.143*t6
         ELSE
            avfit = 0.217*(t6-10.2) + 1.46
         END IF
      END IF

    END function avfit
END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
