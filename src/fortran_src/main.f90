! 2015 UCLCHEM by Serena Viti update by Jon Holdship
! Rewritten for FORTRAN 95 and modulised
PROGRAM uclchem
!Everything to do with physics should be stored in a physics module based on physics-template.f90
!UCLCHEM uses density and temperature in chemical calculations so these MUST be provided, everything else is
!user dependent
USE uclchemwrap, only: cloud,hot_core,cshock
USE io, only: inputId
USE constants, only: dp
IMPLICIT NONE
    CHARACTER (LEN=100):: modelType
    CHARACTER (LEN=100):: paramFile
    CHARACTER(:), ALLOCATABLE :: paramDict
    REAL(dp) :: abundances(500),dissipationResult
    INTEGER :: success,fileLength
    !Any subset of the parameters can be passed in a file on program start
    !see example.inp
    CALL GET_COMMAND_ARGUMENT(1, modelType)  
    CALL GET_COMMAND_ARGUMENT(2, paramFile)

    OPEN(unit=inputId, file=paramFile, action="read", &
    form="unformatted", access="stream")
    INQUIRE(unit=inputId, size=fileLength)

    ALLOCATE(character(fileLength) :: paramDict)
    READ(inputId) paramDict
    CLOSE(inputId)

    SELECT CASE(modelType)
    CASE("CLOUD")
        CALL cloud(paramDict,"",abundances,success)
    CASE("HOTCORE")
        CALL hot_core(3,3.0d2,paramDict,"",abundances,success)
    CASE("CSHOCK")
        CALL cshock(20.0d0,paramDict,"",abundances,dissipationResult,success)
    CASE default
        write(*,*) 'Model type not recognised'
        WRITE(*,*) 'Supported models are: CLOUD, CSHOCK, HOTCORE'
        STOP
    END SELECT
END PROGRAM uclchem