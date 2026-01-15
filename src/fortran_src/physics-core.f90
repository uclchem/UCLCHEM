!This module contains all the must have physics variables. It is imported by both
!the other physics modules (so they can modify the physics) and the chemistry
!module (so it can use the physics for reaction rates).

MODULE physicscore
    USE constants
    USE DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    INTEGER :: dstep

    
    !Optional CR attentuation with column density and better H2 dissociation rates
    REAL(dp) :: h2CRPRate,zetaScale

    !variables either controlled by physics or that user may wish to change
    REAL(dp) :: timeInYears,targetTime,currentTimeold
    REAL(dp) ::  cloudSize
    REAL(dp), allocatable :: av(:),coldens(:),gasTemp(:),dustTemp(:),density(:),density_max(:)

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
    REAL(dp), PARAMETER :: ckLIon(10) = (/-3.331056497233d6,&
                                        1.207744586503d6,& 
                                        -1.913914106234d5,&
                                        1.731822350618d4,&
                                        -9.790557206178d2,&
                                        3.543830893824d1,&
                                        -8.034869454520d-1,& 
                                        1.04880859308d-2,&
                                        -6.188760100997d-5,&
                                        3.122820990797d-8/)
    
    REAL(dp), PARAMETER :: ckHion(10) = (/1.001098610761d7,&
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
    REAL(dp), PARAMETER :: ckLDiss(10)=(/1.582911005330d7,-6.465722684896d6, 1.172189025424d6, -1.237950798073d5, &
        8.393404654312d3, -3.788811358130d2, 1.138688455029d1, -2.197136304567d-1, &
        2.469841278950d-3, -1.232393620924d-5/)
    REAL(dp), PARAMETER :: ckHDiss(10)=(/1.217227462831d7,-4.989649250304d6, 9.079152156645d5, -9.624890825395d4, &
        6.551161486120d3, -2.968976216187d2, 8.959037875226d0, -1.735757324445d-1, &
        1.959267277734d-3, -9.816996707980d-6/)

CONTAINS
    !basic initialization of physics. All physics modules should call this and then
    !do their own initialization.
    SUBROUTINE coreInitializePhysics(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        timeInYears=currentTime/SECONDS_PER_YEAR

        ! Modules not restarted in python wraps so best to reset everything manually.
        IF (ALLOCATED(av)) DEALLOCATE(av,coldens,gasTemp,dustTemp,density,density_max)
        ALLOCATE(av(points),coldens(points),gasTemp(points),dustTemp(points),density(points),density_max(points))

        cloudSize = (rout-rin)*pc
        gasTemp=initialTemp
        dustTemp=gasTemp
        density=initialDens
        currentTimeOld=0.0
        IF (.not. ((ionModel .eq. 'L') .or. (ionModel .eq. 'H'))) THEN
            successFlag=-1
            write(*,*) "Error: ionModel must be either L or H"
            RETURN
        END IF
        IF ((improvedH2CRPDissociation) .and. (.not. cosmicRayAttenuation)) THEN
            successFlag=-1
            write(*,*) "Error: improvedH2CRPDissociation requires cosmicRayAttentuation to also be True"
            RETURN
        END IF

        !calculate initial column density as distance from core edge to current point * density
        DO dstep=1,points
            coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
        END DO
          !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av= baseAv + coldens/1.6d21
        zetaScale=zeta
    END SUBROUTINE coreInitializePhysics

    SUBROUTINE coreUpdatePhysics
        !calculate column density. Remember dstep counts from core centre to edge
        !and coldens should be amount of gas from edge to parcel.
        coldens(dstep)=cloudSize/real(points)*density(dstep)

        ! add previous column densities to current as we move into cloud to get total
        IF (dstep .lt. points) coldens(dstep)=coldens(dstep)+coldens(dstep-1)

        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv + coldens(dstep)/1.6d21
        if (.not. heatingFlag) then 
            dustTemp(dstep)=gasTemp(dstep)
        end if

        IF (cosmicRayAttenuation) CALL ionizationDependency
    END SUBROUTINE coreUpdatePhysics

    pure FUNCTION densdot(density)
    ! Returns the time derivative of the density.                                     
    ! Analytical function taken from Rawlings et al. 1992                             
    ! Called from chemistry.f90, density integrated with abundances so this gives ydot
    REAL(dp), INTENT(IN) :: density
    REAL(dp) :: densdot
    !Rawlings et al. 1992 freefall collapse. With freefallFactor for B-field etc
    IF ((density .lt. finalDens) .and. (freefall)) THEN
        densdot=freefallFactor*(density**4./initialDens)**0.33*&
        &(8.4d-30*initialDens*((density/initialDens)**0.33-1.))**0.5
    ELSE
        densdot=0.0
    ENDIF
    END FUNCTION densdot


    pure FUNCTION dByDnDensdot(density)
    !Defunct function which provides the necessary derivative d(dn/dt)/dn
    !in the case one uses a Jacobian.
    REAL(dp), INTENT(IN) :: density
    REAL(dp) :: dByDnDensdot
    !Rawlings et al. 1992 freefall collapse. With freefallFactor for B-field etc
    IF (density .lt. finalDens) THEN
        dByDndensdot=freefallFactor*8.4d-30*(density**3)*((9.0d0*((density/initialDens)**0.33))-8.0d0)
        dByDnDensdot=dByDnDensdot/(6.0d0*(((density**4.0)/initialDens)**0.66))
        dByDnDensdot=dByDnDensdot/dsqrt(initialDens*8.4d-30*(((density/initialDens)**0.33))-1.0d0)
    ELSE
        dByDnDensdot=0.0
    ENDIF
    END FUNCTION dByDnDensdot

    SUBROUTINE ionizationDependency
        REAL(dp) :: dissSum,dRate,zSum,ionRate
        INTEGER :: k
        !Attenuate CR by column density
        zeta = 1.0
        zSum = 0
        DO k=0,9,1
            IF (ionModel .eq. 'L') THEN
                ionRate=ckLIon(k+1)*log10(coldens(dstep))**k
            ELSEIF (ionModel .eq. 'H') THEN
                ionRate=ckHIon(k+1)*log10(coldens(dstep))**k
            ELSE 
                write(*,*) "WARNING: ionModel switch must be 0 or 1"
            ENDIF
            zSum=zSum+ionRate
        END DO

        ! update/overwrite zeta with attenuated value
        zeta = ((10**zSum)/1.3d-17)* zetaScale

        !rate calculation for H2 dissociation
        IF (improvedH2CRPDissociation) THEN
            dissSum = 0
            DO k=0,9,1
                IF (ionModel .eq. 'L') THEN
                    dRate=ckLDiss(k+1)*log10(coldens(dstep))**k
                ELSEIF (ionModel .eq. 'H') THEN
                    dRate=ckHDiss(k+1)*log10(coldens(dstep))**k
                ELSE 
                    write(*,*) "WARNING: ionModel switch must be L or H"
                ENDIF
                dissSum=dissSum+dRate
            END DO
            h2CRPRate=(10**dissSum)*zetaScale
        END IF
    END SUBROUTINE ionizationDependency

    ! Estimate the column density
    SUBROUTINE findcoldens_core2edge(coldens,rin,rho0,density_scale_radius,density_power_index,r)
      REAL(dp),intent(in) :: rin,r,rho0,density_scale_radius,density_power_index
      REAL(dp),intent(out) :: coldens
      INTEGER :: i,np
      REAL(dp) :: dr,drho,size,r1,r2

      np = 10000
      size = r-rin ![size] in pc
      dr = size/np ![dr] in pc
      coldens = 0.0d0
      IF (size .le. 0.0d0) return

      DO i=1,np
         r1 = rin + (i-1)*dr ![r1] in pc
         r2 = rin + i*dr ![r2] in pc
         drho = 0.5d0*(ngas_r(r2,rho0,density_scale_radius,density_power_index)+ngas_r(r1,rho0,density_scale_radius,density_power_index))
         coldens = coldens + drho*dr*pc
      END DO

    END SUBROUTINE findcoldens_core2edge

    SUBROUTINE findcoldens_edge2core(coldens,rho0,density_scale_radius,density_power_index,r)
        REAL(dp),intent(in) :: rho0,density_scale_radius,density_power_index,r
        REAL(dp),intent(out):: coldens
        if (r.gt.density_scale_radius) then
            coldens = rho0*density_scale_radius*pc/(density_power_index-1.d0) * (r/density_scale_radius)**(1.d0-density_power_index)
        else
            coldens = rho0*density_scale_radius*pc*(density_power_index/(density_power_index-1.d0)-r/density_scale_radius)
        end if
    END SUBROUTINE findcoldens_edge2core

    ! The profile of the gas volumn density
    ! REAL(dp) FUNCTION rhofit(r,rho0,r0,a)
    REAL(dp) FUNCTION ngas_r(r,rho0,density_scale_radius,density_power_index)
      REAL(dp) :: r,rho0,density_scale_radius,density_power_index
      ! [r] in pc, [density_scale_radius] in pc
      ngas_r = rho0/(1.d0 + (r/density_scale_radius)**density_power_index)

    END FUNCTION ngas_r

    REAL(dp) FUNCTION initialDens_r(r,p)
        REAL(dp) :: logn0, logr0,n0_init,r0_init
        REAL(dp) :: r,t,p
        t = 0.0d0
        logn0=61.8d0*(1.175d6-t)**(-0.01) - 49.4d0
        logr0=-28.5d0*(1.175e6-t)**(-0.01) + 28.93d0
        n0_init=10**(logn0)
        r0_init=10**(logr0) * aunit
        initialDens_r=1.0+(r/r0_init)**p
        initialDens_r = n0_init/initialDens_r
    END FUNCTION initialDens_r

END MODULE physicscore
