! Cshock paramterization
! Based on Jimenez-Serra et al. 2008 A&A 482
! http://adsabs.harvard.edu/abs/2008A&A...482..549J
! Latest version: 11-2024 by Dutkowska, K.M.

MODULE cshock_mod
    USE network
    USE physicscore
    USE constants
    USE sputtering
    IMPLICIT NONE

    LOGICAL               :: postShock                           ! current time > dissipation time
    !--------------------
    REAL(dp)              :: dlength, z1, z2, z3, zmax, zn       ! variables related to the shock length scale
    REAL(dp)              :: coshinv1, coshinv2                  ! variables describing the decoupling between ions and neutrals
    REAL(dp)              :: vs, v0, vn, vi, vA, driftVel, dzv   ! velocity variables
    REAL(dp)              :: zn0, vn0, v_out                     ! updatable length and velocity variables
    REAL(dp)              :: grainRadius5, dens6                 ! grain radius and medium density fixed for unit
    REAL(dp)              :: tsat, dissipationTime               ! saturation and dissipation time
    REAL(dp)              :: mun                                 ! mean molecular weight of neutrals (2*H nucleus mass in cgs)
    REAL(dp)              :: bm0                                 ! magnetic parameter; bm0 = B0/sqrt(initialDens); in microGauss
    REAL(dp)              :: at, a1, maxTemp                     ! temperature constants (Eq. (4); Jimenez-Serra et al. 2008) and T_n,max
    REAL(dp)              :: minimumPostshockTemp = 0.0          ! defined post-shock conditon 
    REAL(dp)              :: timestepFactor = 0.01               ! time increment 
    !--------------------
    REAL(dp), PARAMETER   :: bt = 6.                             ! b_T parameter (Eq. (4); Jimenez-Serra et al. 2008)
    REAL(dp), PARAMETER   :: grainRadius=1.0d-5                  ! a constant grain radius of 1um 
    !--------------------
    REAL(dp), allocatable :: tn(:), ti(:), tgc(:), tgr(:), tg(:) ! fluids and grains heating variables
    !--------------------
    ! Variables related to the collisional and radiative heating of grains. Approach follows
    ! that of Draine, Roberge & Dalgarno (1983) and Hollenbach, Takahashi & Tielens (1991)
    ! Note: currently not considered!
    ! REAL(dp)              :: Frs, tau100, G0
    ! REAL(dp)              :: tgr1, tgr2, trs0      
    ! REAL(dp), PARAMETER   :: nu0 = 3.0d15             

CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks if inputs make sense, calculates various constants and !
    ! sets up variables for the shock paramterization               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics(successFlag)
        REAL(dp)             :: v01, g1, g2
        INTEGER, INTENT(OUT) :: successFlag
        
        successFlag = 1
        driftVel    = 0.0
        zn0         = 0.0
        vn0         = 0.0

        ! Set cooling variables to off by default (change if required)
        postShock   = .False.

        !Check input sanity and set inital values
        cloudSize   = (rout-rin)*pc

        IF (freefall) THEN
            write(*,*) "Cannot have freefall on during cshock"
            Write(*,*) "setting freefall=0 and continuing"
            freefall = .False.
        ENDIF

        IF (points .gt. 1) THEN
            WRITE(*,*) "Cannot have more than one point in cshock"
            successFlag = -1
            RETURN
        END IF

        density = initialDens
        ! C-shock initialization

        IF (ALLOCATED(tn)) deallocate(tn, ti, tgc, tgr, tg)
        allocate(tn(points), ti(points), tgc(points), tgr(points), tg(points))
        
        currentTimeOld = 0.0
        driftVel       = 0.0
        zn0            = 0.0
        vn0            = 0.0
        mun            = 2*mh
        grainRadius5   = grainRadius/4.e-5
        dens6          = density(dstep)/1.e6

        ! maxTemp set by vs and pre-shock density
        ! Polynomial fits to values taken from Draine et al. 1983 have been made and coefficients placed here.
        ! Tested with 1e3 < dens < 1e6
        ! Fits only available for density of 1e4 and 1e6 so in between we average
        IF (initialDens .gt. 10**5.5) THEN
            maxTemp = (2.91731*vs*vs)-(23.78974*vs)+225.204167337
        ELSE IF (initialDens .gt. 10.0**4.5) THEN
            maxTemp = (3.38989*vs*vs)+(16.6519*vs)+96.569
            maxTemp = 0.5*maxTemp
        ELSE
            maxTemp = (0.47258*vs*vs)+(40.44161*vs)-128.635455216
        END IF    
    
        ! tsat proportional to 1/pre-shock density
        ! For more information on saturation time see:
        ! Sect. 4.2. in Jimenez-Serra et al. 2008
        tsat = (-15.38729*vs*vs*vs)+(2069.56962*vs*vs)-(90272.826991*vs)+1686858.54278
        tsat = tsat/initialDens

        ! First, we definre the initial parameters that define the C-shock structure
        ! 1. Length of the dissipation region and dissipation time
        !    [see Dopita & Sutherland 2003]
        dlength         = 12.0*pc*vs/initialDens
        dissipationTime = (dlength*1.0d-5/vs)/SECONDS_PER_YEAR

        ! 2. Parameters describing the decoupling between the ion and the neutral fluids. 
        coshinv1 = log((1/0.01)+sqrt((1/0.01)**2-1))
        coshinv2 = log((1/0.15)+sqrt((1/0.15)**2-1))

        ! 3. z2 (= zn) & z1 (= zi) 
        !    z2 is obtained by assuming that at z = dlength, 
        !    the velocity of the neutrals is 99% (vs-v0; see Draine et al. 1983);
        !    We assume that z2/z1=4.5 
        !    [for more see Sect. 4.1. in Jimenez-Serra et al. 2008]
        z2 = dlength/coshinv1
        z1 = z2/4.5

        ! 4. zmax is the distance at which Tn reaches its maximum. 
        !    This happens when the neutral fluid reaches velocities 
        !    that are almost 0.85% (vs-v0; see Draine et al. 1983); 
        !    z3 has to be 1/6 zmax
        !    [for more see Sect. 4.1. in Jimenez-Serra et al. 2008]
        zmax = dlength/coshinv2
        z3   = zmax/6

        ! 5. maxTemp is taken from Fig.9b in Draine et al. (1983) 
        !    and the at constant is derived as:
        a1 = 6.0
        at = (1/zmax)*((maxTemp-initialTemp)*(dexp(a1)-1.))**(1./6.)

        ! Second, we calculate v0 that depends on the alfven and the shock velocities.
        ! v0 is the final downstream velocity of ions and neutrals ()
        ! Magnetic field is in microGauss. 
        ! We assume strong magnetic field, i.e., 
        ! B0 = bm0*sqrt(initialDens) [Eq. (63); Draine, Roberge & Dalgarno 1983], 
        ! where the magnetic parameter bm0 = 1 microgauss
        ! For the general case, the Alfven velocity is calculated as 
        ! vA = B0/sqrt(4*pi*initialDens). 
        ! If we substitute the expression of B0, we obtain that:
        ! vA = bm0/sqrt(4*pi*mH)
        ! For historical reasons, we leave the calculation of the Alfven speed
        ! in neutrals that uses the magnetic parameter bm0
        bm0 = B0/sqrt(initialDens)
        bm0 = bm0*1D-06
        vA  = bm0/sqrt(4*pi*mh)
        vA  = vA/km

        !Calculation of v0, final velocity of ions/neutrals in the shock frame
        v0  = 2.
        v01 = 0
        DO WHILE (abs(v0-v01) .ge. 1e-6)
            v01 = v0
            g1  = -(vA**2*vs**2)/2
            g2  = v01**2-v01*vs-vA**2/2
            v0  = sqrt(g1/g2)
        END DO
        CALL sputteringSetup
    END SUBROUTINE initializePhysics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called every time loop in main.f90. Sets the timestep for the next output from   !
    ! UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    ! but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updateTargetTime
        IF (timeInYears .lt. 2.0*dissipationTime) THEN
            !get a nice sampling along the shock
            targetTime = (timeInYears+timestepFactor*dissipationTime)*SECONDS_PER_YEAR
        ELSE
            targetTime = (1.1*timeInYears)*SECONDS_PER_YEAR
        END IF
    END SUBROUTINE updateTargetTime

    ! Calculate shock properties for current time and set density, temperature and Av
    SUBROUTINE updatePhysics
        ! First calculate velocity of neutrals and position of shock front at currentTime
        call shst

        dzv = sqrt(2*vn*(vs-v0)-vn**2)
        dzv = 1./(dzv*km)
        dzv = z2*((vs-v0)/(vs-v0-vn))*dzv

        ! We introduce the gas temperature curve along the dissipation region of the
        ! C-shock. We also take into account that the gas and dust are decoupled. We
        ! use the equations for the collisional and radiative heating of grains of
        ! Draine, Roberge & Dalgarno (1983) and Hollenbach, Takahashi & Tielens (1991).
        ! For more on Tn and Ti see Sect. 2 in in Jimenez-Serra et al. (2008)
        tn(dstep) = initialTemp+((at*zn)**bt)/(dexp(zn/z3)-1)
        ti(dstep) = tn(dstep)+(mun*(driftVel*km)**2/(3*K_BOLTZ))
        ! Grain collisional heating; Eq. (58) in Draine et al. 1983
        tgc(dstep) = 15*(dens6/grainRadius5)**(0.1818)*(tn(dstep)/1000.0)**(0.2727)

        ! Grain radiative heating - currently not considered; 
        ! Based on Hollenbach, Takahashi & Tielens (1991)
        ! Frs        = 0.25*density(1)*mun*(vn*km)**3
        ! G0         = Frs
        ! trs0       = 12.2*G0**0.2
        ! tau100     = 2.7d2*G0/trs0**5
        ! tgr1       = 8.9d-11*nu0*G0*dexp(1.8*av(dstep))+2.7**5
        ! tgr2       = 3.4d-2*(0.42-log(3.5d-2*tau100*trs0))*tau100*trs0**6
        ! tgr(dstep) = (tgr1+tgr2)**0.2

        ! If we don't include the radiative heating that is characteristic
        ! of J-type shocks:
        tgr(dstep) = 0.0

        ! Total grain heating
        tg(dstep) = tgc(dstep)+tgr(dstep)

        ! Density change as shock evolves
        IF (timeInYears .gt. 0.0) THEN
            density = initialDens*vs/(vs-vn)
        END IF

        ! Temperature change as shock evolves
        IF (timeInYears .gt. 0.0) THEN
            tn(dstep)      = initialTemp+((at*zn)**bt)/(dexp(zn/z3)-1)
            gasTemp(dstep) = tn(dstep)
            ti(dstep)      = tn(dstep)+(mun*(driftVel*km)**2/(3*K_BOLTZ))
        ENDIF
        postShock = (timeInYears .gt. dissipationTime)

        IF ((gasTemp(dstep) .lt. minimumPostshockTemp) .AND. (postShock)) THEN
            gasTemp(dstep) = minimumPostshockTemp
        END IF
        dustTemp = gasTemp

    END SUBROUTINE updatePhysics

    ! For c-shock, sublimation is simply the sputtering subroutine
    SUBROUTINE sublimation(abund)
        REAL(dp)        :: abund(nspec+1, points)
        INTENT(INOUT)   :: abund
        REAL(dp)        :: timeDelta

        timeDelta = (currentTime-currentTimeOld)
        IF ((sum(abund(iceList,dstep)) .gt. 1d-25) .AND. (driftVel .gt. 0))&
        & CALL sputterIces(abund(:,dstep), driftVel, gasTemp(dstep), density(dstep), timeDelta)
        WHERE(abund.lt. 1.0d-50) abund = 0.0d-50
    END SUBROUTINE sublimation


    ! The subroutine below has been written by Izaskun Jimenez-Serra.
    ! It calculates the distance along the dissipation region (zn) 
    ! and the velocity of the gas as the shock evolves with time.
    SUBROUTINE shst
        REAL(dp) :: vn1,f1,f0,xcos,acosh
        INTEGER  :: loopCount
        ! We calculate the physical structure of the shock
        ! set vn1 arbitrarily high to ensure while loop is done at least once
        vn1       = 1d30
        vn        = vn0
        loopCount = 0
        DO WHILE ((abs(vn-vn1).ge.1.e-10) .and. (loopCount .lt. 100))
            vn1       = vn
            f1        = vs-vn1
            f0        = vs-vn0
            zn        = zn0+(currentTime-currentTimeOld)*km*(f1+f0)/2
            xcos      = zn/z2
            acosh     = 0.5*(dexp(xcos)+dexp(-xcos))
            vn        = (vs-v0)-((vs-v0)/acosh)
            loopCount = loopCount+1
        END  DO

        xcos  = zn/z1
        acosh = 0.5*(dexp(xcos)+dexp(-xcos))
        vi    = (vs-v0)-((vs-v0)/acosh)

        ! Store all variables as initial values for next iteration
        driftVel = vi-vn ! drift velocity
        zn0      = zn    ! distance into a dissipation region
        vn0      = vn    ! neutrals velocity
        v_out    = f0    ! final-shock velocity
    END SUBROUTINE shst

    SUBROUTINE ShockTypeCheck
        REAL(dp)   :: vAc        ! Alfven speed in charged fluid
        REAL(dp)   :: rhoc, rhon ! Density of charged and neutral fluid
        INTENT(IN) :: abund      ! Read the abundances from the 1st time step to evaluate the conditions


    END SUBROUTINE ShockTypeCheck
    
END MODULE cshock_mod
