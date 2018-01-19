MODULE EXPLOSIONS
IMPLICIT NONE
EXTERNAL dvode

    !Explosion variables
    double precision :: initialDensity=1.0d20,density,initTemp=1000.0,temperature
    double precision :: expansionFactor,trappingFactor=1.0,adiabaticIndex=1
    double precision,parameter :: SOUND_SPEED=1.0d4,INITIAL_RADIUS=1.0d-5,FINAL_TIME=1.0d-8
    double precision, parameter :: HYDROGEN_THRESHOLD=0.05,JCR_FACTOR=1.0/5.41061d8,CORE_ATOMS=1.0d8
    integer,parameter :: nEquation=35, nReactions=18
    integer :: i

    double precision :: rate(nReactions)

    !integrator stuff
    double precision :: expCurrentTime,expTargetTime
    double precision :: h2oNormalization

    !do an explosion?
    logical, parameter :: EXPLOSION_RUN=.true.
CONTAINS
  SUBROUTINE explodeCycle(inputAbunds)
    double precision, INTENT(INOUT) :: inputAbunds(:)
    integer :: grainIndices(nEquation),gasIndices(nEquation)
    double precision :: Y(nEquation)
    integer :: ITOL,ITASK,ISTATE,IOPT,MESFLG,lrw,liw
    integer :: MXSTEP,MF
    integer,allocatable :: IWORK(:)
    double precision :: reltol,rpar,ipar
    double precision, allocatable :: RWORK(:),abstol(:)

    !Take the bits of main abundance array that is involved in explosion
    grainIndices=(/3,6,23,29,31,36,38,42,77,86,94,96,104,106,108,119,127 &
     &       ,129,132,169,171,177,189,191,193,198,200,205,207,227,229 &
     &       ,237,239,241,243/)

    !they will be put back into gas so keep indices of gas equivalents
    gasIndices=(/2,5,22,28,30,35,37,41,76,85,93,95,103,105,107,118,126,128 &
     &       ,131,168,170,176,188,190,192,197,199,204,206,226,228,236 &
     &       ,238,240,242/)

    Y=inputAbunds(grainIndices)
    h2oNormalization=Y(8)
    Y=Y/h2oNormalization
    inputAbunds(grainIndices)=0
    !Dvode settings
    ISTATE=1;MF=22;ITOL=2;ITASK=1;IOPT=1;MESFLG=1
    reltol=1e-6;MXSTEP=10000

    write(*,*) ISTATE

    LIW=30+nEquation
    LRW=22+(9*nEquation)+(2*nEquation*nEquation)
    ALLOCATE(IWORK(LIW),RWORK(LRW),abstol(nEquation))
    IWORK(6)=MXSTEP

    abstol=1.0d-16*Y
    WHERE(abstol<1d-30) abstol=1.0d-30
    rate=1.0d-28 !from Rawlings 2013/Coutens 2017

    write(*,*) "EXPLOSION!"
    expansionFactor=trappingFactor*SOUND_SPEED/INITIAL_RADIUS
    expCurrentTime=1d-20
    i=0

    DO WHILE (expCurrentTime .LT. FINAL_TIME)
        WHERE(Y<1.0d-30) Y=1.0d-30
        IF (expCurrentTime.lt.0.00001*FINAL_TIME) THEN
          expTargetTime=expCurrentTime*10.0
          write(99,*) expCurrentTime,density,temperature
          !CH3+OH+H2O->CH3OH abundances written in order
          write(89,*) Y(3),Y(7),Y(11),Y(16)
          !expCurrentTime=expTargetTime
        ELSE
          expTargetTime=expCurrentTime+0.00001*FINAL_TIME
          if (i.eq.10000)then
            write(99,*) expCurrentTime,density,temperature
            !CH3+OH+H2O->CH3OH abundances written in order
            write(89,*) Y(3),Y(7),Y(11),Y(16)
            !expCurrentTime=expTargetTime
            i=0
          END IF

        END IF
        CALL expansion
        CALL DVODE(F,nEquation,Y,expCurrentTime,expTargetTime,ITOL,reltol,abstol,ITASK,ISTATE,IOPT,&
                 &      RWORK,LRW,IWORK,LIW,JAC,MF,RPAR,IPAR)
        SELECT CASE(ISTATE)
            CASE(-1)
                !More steps required for this problem
                MXSTEP=MXSTEP*2    
                write(79,*)'Call to LSODE returned -1 meaning that MXSTEP exceeded'
                write(79,*)'but the integration was successful'
                write(79,*)'Doubling MXSTEP from:',MXSTEP,' to:',MXSTEP*2
                ISTATE=3
                IOPT=1
                IWORK(6)=MXSTEP
            CASE(-2)
                !Tolerances are too small for machine but succesful to current currentTime
                abstol=abstol*10.0
                ISTATE=3
            CASE(-3)
                write(79,*) "DVODE found invalid inputs"
                write(79,*) "abstol"
                write(79,*) abstol
                STOP
            CASE(-4)
                write(*,*) Y
                !Successful as far as currentTime but many errors.
                !Make targetTime smaller and just go again
                expTargetTime=(expCurrentTime+expTargetTime)/2.0
                ISTATE=1
            CASE DEFAULT
                IOPT=0
                ISTATE=3
        END SELECT
        i=i+1

    END DO

    Y=Y*h2oNormalization
    inputAbunds(gasIndices)=inputAbunds(gasIndices)+Y
    ISTATE=1
    write(*,*) "EXPLODED"
  END SUBROUTINE explodeCycle


  !update physical parameters of explosion,see Cecchi-Pestellini 2010.
  SUBROUTINE expansion
      !Gas expands at some fraction of sound speed

      density=(1.0+expansionFactor*expCurrentTime)**3
      density=initialDensity/density

      temperature=(1.0+expansionFactor*expCurrentTime)**adiabaticIndex
      temperature=initTemp/temperature
  END SUBROUTINE expansion

  !Function checks whether ices should now explode
  LOGICAL FUNCTION explosionLimit(hydrogenFractionalAbundace,totalMantleAbundance,dustAbundance)
      double precision hydrogenFractionalAbundace,totalMantleAbundance,dustAbundance
      double precision atomsPerGrain
      !Rawlings
      atomsPerGrain=(CORE_ATOMS*JCR_FACTOR)+totalMantleAbundance

      !me
      !atomsPerGrain=(CORE_ATOMS*dustAbundance)+totalMantleAbundance

      explosionLimit= (hydrogenFractionalAbundace/atomsPerGrain .GT. HYDROGEN_THRESHOLD)
      write(*,*) hydrogenFractionalAbundace/atomsPerGrain
      write(82,*) hydrogenFractionalAbundace,totalMantleAbundance,atomsPerGrain
  END FUNCTION explosionLimit


  SUBROUTINE F(nEquation, T, Y, YDOT)
      !DLSODE calls this subroutine to ask it what the RHS of the equations dy/dt=... are    

      INTEGER :: nEquation
      DOUBLE PRECISION :: T,Y(nEquation),YDOT(nEquation)
      DOUBLE PRECISION :: D,loss,prod

      D=density

      YDOT(1) = 0.0
      YDOT(2) = 0.0
      LOSS = -RATE(1)*Y(8)*Y(7)*D*D-2*RATE(6)*Y(8)*D*D-RATE(7)*Y(8) &
     &       *Y(5)*D*D-RATE(8)*Y(8)*Y(15)*D*D-RATE(9)*Y(13)*Y(8)*D*D &
     &       -RATE(10)*Y(8)*Y(10)*D*D
      YDOT(3) = Y(3)*LOSS
      YDOT(4) = 0.0
      LOSS = -RATE(2)*Y(8)*Y(7)*D*D-RATE(7)*Y(8)*Y(3)*D*D-RATE(11)*Y(8) &
     &       *Y(15)*D*D-RATE(12)*Y(13)*Y(8)*D*D-RATE(13)*Y(8)*Y(10)*D*D
      YDOT(5) = Y(5)*LOSS
      YDOT(6) = 0.0
      LOSS = -RATE(1)*Y(8)*Y(3)*D*D-RATE(2)*Y(8)*Y(5)*D*D-RATE(3)*Y(8) &
     &       *Y(15)*D*D-RATE(4)*Y(13)*Y(8)*D*D-RATE(5)*Y(8)*Y(10)*D*D
      YDOT(7) = Y(7)*LOSS
      LOSS = -RATE(1)*Y(3)*Y(7)*D*D-RATE(2)*Y(5)*Y(7)*D*D-RATE(3)*Y(15) &
     &       *Y(7)*D*D-RATE(4)*Y(13)*Y(7)*D*D-RATE(5)*Y(10)*Y(7)*D*D &
     &       -RATE(6)*Y(3)*Y(3)*D*D-RATE(7)*Y(3)*Y(5)*D*D-RATE(8)*Y(3) &
     &       *Y(15)*D*D-RATE(9)*Y(13)*Y(3)*D*D-RATE(10)*Y(3)*Y(10)*D*D &
     &       -RATE(11)*Y(15)*Y(5)*D*D-RATE(12)*Y(13)*Y(5)*D*D-RATE(13) &
     &       *Y(10)*Y(5)*D*D-RATE(14)*Y(15)*Y(15)*D*D-RATE(15)*Y(13) &
     &       *Y(15)*D*D-RATE(16)*Y(15)*Y(10)*D*D-RATE(17)*Y(13)*Y(13)*D &
     &       *D-RATE(18)*Y(13)*Y(10)*D*D
      PROD = +RATE(1)*Y(8)*Y(3)*Y(7)*D*D+RATE(2)*Y(8)*Y(5)*Y(7)*D*D &
     &       +RATE(3)*Y(8)*Y(15)*Y(7)*D*D+RATE(4)*Y(13)*Y(8)*Y(7)*D*D &
     &       +RATE(5)*Y(8)*Y(10)*Y(7)*D*D+RATE(6)*Y(8)*Y(3)*Y(3)*D*D &
     &       +RATE(7)*Y(8)*Y(3)*Y(5)*D*D+RATE(8)*Y(8)*Y(3)*Y(15)*D*D &
     &       +RATE(9)*Y(13)*Y(8)*Y(3)*D*D+RATE(10)*Y(8)*Y(3)*Y(10)*D*D &
     &       +RATE(11)*Y(8)*Y(15)*Y(5)*D*D+RATE(12)*Y(13)*Y(8)*Y(5)*D*D &
     &       +RATE(13)*Y(8)*Y(10)*Y(5)*D*D+RATE(14)*Y(8)*Y(15)*Y(15)*D &
     &       *D+RATE(15)*Y(13)*Y(8)*Y(15)*D*D+RATE(16)*Y(8)*Y(15)*Y(10) &
     &       *D*D+RATE(17)*Y(13)*Y(13)*Y(8)*D*D+RATE(18)*Y(13)*Y(8) &
     &       *Y(10)*D*D
      YDOT(8) = PROD+Y(8)*LOSS
      YDOT(9) = 0.0
      LOSS = -RATE(5)*Y(8)*Y(7)*D*D-RATE(10)*Y(8)*Y(3)*D*D-RATE(13) &
     &       *Y(8)*Y(5)*D*D-RATE(16)*Y(8)*Y(15)*D*D-RATE(18)*Y(13)*Y(8) &
     &       *D*D
      YDOT(10) = Y(10)*LOSS
      PROD = +RATE(6)*Y(8)*Y(3)*Y(3)*D*D
      YDOT(11) = PROD
      YDOT(12) = 0.0
      LOSS = -RATE(4)*Y(8)*Y(7)*D*D-RATE(9)*Y(8)*Y(3)*D*D-RATE(12)*Y(8) &
     &       *Y(5)*D*D-RATE(15)*Y(8)*Y(15)*D*D-2*RATE(17)*Y(8)*D*D &
     &       -RATE(18)*Y(8)*Y(10)*D*D
      YDOT(13) = Y(13)*LOSS
      PROD = +RATE(7)*Y(8)*Y(3)*Y(5)*D*D
      YDOT(14) = PROD
      LOSS = -RATE(3)*Y(8)*Y(7)*D*D-RATE(8)*Y(8)*Y(3)*D*D-RATE(11)*Y(8) &
     &       *Y(5)*D*D-2*RATE(14)*Y(8)*D*D-RATE(15)*Y(13)*Y(8)*D*D &
     &       -RATE(16)*Y(8)*Y(10)*D*D
      YDOT(15) = Y(15)*LOSS
      PROD = +RATE(1)*Y(8)*Y(3)*Y(7)*D*D
      YDOT(16) = PROD
      YDOT(17) = 0.0
      PROD = +RATE(2)*Y(8)*Y(5)*Y(7)*D*D
      YDOT(18) = PROD
      YDOT(19) = 0.0
      PROD = +RATE(10)*Y(8)*Y(3)*Y(10)*D*D
      YDOT(20) = PROD
      YDOT(21) = 0.0
      PROD = +RATE(13)*Y(8)*Y(10)*Y(5)*D*D
      YDOT(22) = PROD
      PROD = +RATE(9)*Y(13)*Y(8)*Y(3)*D*D
      YDOT(23) = PROD
      PROD = +RATE(8)*Y(8)*Y(3)*Y(15)*D*D
      YDOT(24) = PROD
      PROD = +RATE(5)*Y(8)*Y(10)*Y(7)*D*D
      YDOT(25) = PROD
      PROD = +RATE(12)*Y(13)*Y(8)*Y(5)*D*D
      YDOT(26) = PROD
      PROD = +RATE(11)*Y(8)*Y(15)*Y(5)*D*D
      YDOT(27) = PROD
      PROD = +RATE(3)*Y(8)*Y(15)*Y(7)*D*D
      YDOT(28) = PROD
      PROD = +RATE(4)*Y(13)*Y(8)*Y(7)*D*D
      YDOT(29) = PROD
      PROD = +RATE(18)*Y(13)*Y(8)*Y(10)*D*D
      YDOT(30) = PROD
      PROD = +RATE(16)*Y(8)*Y(15)*Y(10)*D*D
      YDOT(31) = PROD
      YDOT(32) = 0.0
      PROD = +RATE(14)*Y(8)*Y(15)*Y(15)*D*D
      YDOT(33) = PROD
      PROD = +RATE(15)*Y(13)*Y(8)*Y(15)*D*D
      YDOT(34) = PROD
      PROD = +RATE(17)*Y(13)*Y(13)*Y(8)*D*D
      YDOT(35) = PROD
  END SUBROUTINE F

  !This is a dummy for DLSODE, it has to call it but we do not use it.
  subroutine JAC (nEquation,T,Y,ML,MU,PD,NROWPD)
       INTEGER  nEquation, ML, MU, NROWPD
       DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
  END SUBROUTINE JAC    
END MODULE EXPLOSIONS