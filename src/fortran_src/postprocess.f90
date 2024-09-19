! Module for post-processing data from hydrodynamical simulations
! Requires modifications to some core functions of the code
MODULE postprocess_mod
    USE constants
    USE DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp
    USE physicscore
    USE network
    USE f2py_constants
    IMPLICIT NONE

    ! character(len=100) :: trajecfile
    logical :: lusecoldens
    logical :: usepostprocess = .true.
    ! integer,parameter :: tfid=66
    integer :: ntime,tstep
    double precision, allocatable, dimension(:) :: ltime, ldens, lra, lzeta, lradfield, lgtemp, ldtemp
    double precision, allocatable, dimension(:) :: lnh, lnh2, lnco, lnc
    
CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called at start of UCLCHEM run
    ! Uses values in defaultparamters.f90 and any inputs to set initial values        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics(successFlag, timegrid, densgrid, radgrid, zetagrid, gtempgrid,&
        &dtempgrid, usecoldens, timepoints, nhgrid, nh2grid, ncogrid, ncgrid)
      INTEGER, INTENT(OUT) :: successFlag
      INTEGER, INTENT(IN) :: timepoints
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: timegrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: densgrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: radgrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: zetagrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: gtempgrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: dtempgrid
      LOGICAL, INTENT(IN) :: usecoldens
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: nhgrid
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: nh2grid
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: ncogrid
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: ncgrid
      ! write(*,*) 'Initialising postprocessing module'
      

      successFlag=0
      ! Create a local variable that can be recycled in updatePhysics
      lusecoldens = usecoldens
      ! Check if NHgrid is present, if so, we need to use the custom column densities
      if ( lusecoldens) then 
        cloudSize=0. ! Shielding column densities supplied separately
      end if
      endatfinaldensity = .false. ! Needs to end at final time to work properly
      freefall = .false. ! Won't work with freefall on
      
      if (.not. allocated(ltime)) then 
        allocate (ltime(timepoints), ldens(timepoints), lradfield(timepoints), lzeta(timepoints),&
       &lgtemp(timepoints), ldtemp(timepoints))
      end if
      ! Store the parameter grids into a local module contexts
      ltime(:) = timegrid
      ldens(:) = densgrid
      lradfield(:) = radgrid
      lzeta(:) = zetagrid
      lgtemp(:) = gtempgrid
      ldtemp(:) = dtempgrid
      ! If we have custom column densities, allocate and store them
      if (lusecoldens) then
      ! Only allocate the column densities if we need them:
        if (.not. allocated(lnh)) then
          allocate (lnh(timepoints), lnh2(timepoints), lnco(timepoints), lnc(timepoints))
        end if
        lnh(:) = nhgrid
        lnh2(:) = nh2grid
        lnco(:) = ncogrid
        lnc(:) = ncgrid
      end if 

      ! Initialise values to t=1 and overwrite them.
      tstep = 1
      targettime = ltime(tstep)
      density(dstep) = ldens(tstep)
      gastemp(dstep) = lgtemp(tstep)
      dusttemp(dstep) = ldtemp(tstep)
      radfield = lradfield(tstep)
      zeta = lzeta(tstep)
      if (lusecoldens) then
        coldens(dstep) = lnh(tstep)
        av(dstep) = 5.348e-22 * coldens(dstep)
      end if 
      ! Set final time to end of tracer histories
      finaltime = timegrid(timepoints)/seconds_per_year
    END SUBROUTINE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE updateTargetTime
      ! Set target time from postprocessing data
      targettime = ltime(tstep) + 1.*seconds_per_year

      ! write(*,"('Integrating chemistry to timestep ',I3,' ',ES10.3,' years')") tstep,targettime/seconds_per_year

    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics

      ! Update physical properties to values at tstep (== targettime)
      targettime = ltime(tstep)
      density(dstep) = ldens(tstep)
      gastemp(dstep) = lgtemp(tstep)
      dusttemp(dstep) = ldtemp(tstep) 
      radfield = lradfield(tstep)
      zeta = lzeta(tstep)
      
      if (lusecoldens) then
        coldens(dstep) = lnh(tstep)
        av(dstep) = 5.348e-22 * coldens(dstep)
      end if 
      ! If this is the last point, update tstep to move onto next time dump
      tstep = tstep + 1
      
    END SUBROUTINE updatePhysics

    SUBROUTINE sublimation(abund)
        ! This subroutine must be in every physics module so we dummy it here.
        !f2py integer, intent(aux) :: points
        DOUBLE PRECISION :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
    END SUBROUTINE sublimation    
END MODULE postprocess_mod


