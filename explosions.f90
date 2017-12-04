MODULE EXPLOSIONS
IMPLICIT NONE
    !Explosion variables
    double precision :: initialDensity=1.0d20,density,initialTemp=500.0,temperature
    double precision :: expansionFactor,trappingFactor=1.0,adiabaticIndex=1
    double precision,parameter :: SOUND_SPEED=1.0d4,INITIAL_RADIUS=1.0d-5,HYDROGEN_THRESHOLD=0.05
    double precision, parameter :: FINAL_TIME=1.0d-9

    !integrator stuff
    double precision :: currentTime,targetTime
    
CONTAINS
    SUBROUTINE explodeCycle()
        expansionFactor=trappingFactor*SOUND_SPEED/INITIAL_RADIUS
        currentTime=0.0
        

        DO WHILE (currentTime .LT. FINAL_TIME)
            targetTime=currentTime+0.001*FINAL_TIME
            CALL expansion
            write(99,*) currentTime,density,temperature
            currentTime=targetTime
        END DO
    END SUBROUTINE explodeCycle


    !update physical parameters of explosion,see Cecchi-Pestellini 2010.
    SUBROUTINE expansion
        !Gas expands at some fraction of sound speed

        density=(1.0+expansionFactor*currentTime)**3
        density=initialDensity/density

        temperature=(1.0+expansionFactor*currentTime)**adiabaticIndex
        temperature=initialTemp/temperature
    END SUBROUTINE expansion

    !Function checks whether ices should now explode
    LOGICAL FUNCTION explosionLimit(hydrogenFractionalAbundace)
        double precision hydrogenFractionalAbundace
        explosionLimit= (hydrogenFractionalAbundace .GT. HYDROGEN_THRESHOLD)
    END FUNCTION explosionLimit

END MODULE EXPLOSIONS