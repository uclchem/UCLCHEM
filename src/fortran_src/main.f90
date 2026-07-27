! 2022 UCLCHEM v3.0
! The canonical main file where the core code is written is actually found in wrap.f90
! main.f90 just provides a simple fortran interface to the core code so that a binary can
! be built and used directly from the command line.
program uclchem
use f2py_constants, only: nspec
use constants, only: dp
use io, only: inputId
use uclchemwrap, only: cloud, hot_core, cshock, postprocess
implicit none
    character (LEN=100) :: modelType
    character (LEN=100) :: paramFile
    character (LEN=32) :: model_arg1, model_arg2, model_arg3
    character(:), allocatable :: paramDict
    real(dp) :: abundances(nspec),dissipationResult
    integer :: success,fileLength,model_index
    real(8) :: max_temp, vshock, timestep_factor, minimum_temp
    !Any subset of the parameters can be passed in a file on program start
    !see example.inp
    call GET_COMMAND_ARGUMENT(1, modelType)
    call GET_COMMAND_ARGUMENT(2, paramFile)

    open(unit=inputId, file=paramFile, action="read", &
    form="unformatted", access="stream")
    inquire(unit=inputId, size=fileLength)

    allocate(character(fileLength) :: paramDict)
    read(inputId) paramDict
    close(inputId)
    select case(modelType)
    case("CLOUD")
        call cloud(paramDict,"",.false.,.false.,abundances,success)
    case("HOTCORE")
        !call hot_core (temp_indx,maxTemp,....)
        ! Read the strings from the CLI, then convert them to int/float
        call GET_COMMAND_ARGUMENT(3, model_arg1)
        call GET_COMMAND_ARGUMENT(4, model_arg2)
        read(model_arg1, *) model_index
        read(model_arg2, *) max_temp
        call hot_core(model_index,max_temp,paramDict,"",.false.,.false.,&
        &abundances,success)
    case("CSHOCK")
        !call cshock(vs,timestep_factor,minimum_temp,....)
        !sensible defaults in order: 20.0, 0.01d0, 10.0
        ! Read the strings from the CLI, then convert them to int/float
        call GET_COMMAND_ARGUMENT(3, model_arg1)
        call GET_COMMAND_ARGUMENT(4, model_arg2)
        call GET_COMMAND_ARGUMENT(5, model_arg3)
        read(model_arg1, *) vshock
        read(model_arg2, *) timestep_factor
        read(model_arg3, *) minimum_temp
       call cshock(vshock,timestep_factor,minimum_temp,paramDict,"",.false.,.false.,&
       &abundances,dissipationResult,success)
    case("POSTPROCESS")
       call postprocess(paramDict,"",.false.,.false.,abundances,success)
    case default
        write(*,*) "Model type not recognized"
        write(*,*) "Supported models are: CLOUD, CSHOCK, HOTCORE and POSTPROCESS"
        stop
    end select
end program uclchem
