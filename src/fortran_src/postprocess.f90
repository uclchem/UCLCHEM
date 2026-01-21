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
    logical :: luseav
    logical :: usepostprocess = .true.
    ! integer,parameter :: tfid=66
    integer :: tstep
    integer :: max_tstep  ! Last valid (non-zero) timestep
    integer :: postprocess_error  ! Non-zero if a fatal error occurred during updatePhysics
    double precision, allocatable, dimension(:) :: ltime, ldens, lra, lzeta, lradfield, lgtemp, ldtemp
    double precision, allocatable, dimension(:) :: lnh, lnh2, lnco, lnc, lav

    
CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called at start of UCLCHEM run
    ! Uses values in defaultparamters.f90 and any inputs to set initial values        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics(successFlag, timegrid, densgrid, radgrid, zetagrid, gtempgrid,&
        &dtempgrid, useav, avgrid, usecoldens, timepoints, nhgrid, nh2grid, ncogrid, ncgrid)
      INTEGER, INTENT(OUT) :: successFlag
      INTEGER, INTENT(IN) :: timepoints
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: timegrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: densgrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: radgrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: zetagrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: gtempgrid
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints) :: dtempgrid
      LOGICAL, INTENT(IN) :: useav
      DOUBLE PRECISION, INTENT(IN), DIMENSION(timePoints), OPTIONAL :: avgrid
      LOGICAL, INTENT(IN) :: usecoldens
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: nhgrid
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: nh2grid
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: ncogrid
      DOUBLE PRECISION, INTENT(IN), OPTIONAL, DIMENSION(timePoints) :: ncgrid
      

      successFlag=0
      postprocess_error = 0
      ! Create a local variable that can be recycled in updatePhysics
      lusecoldens = usecoldens
      luseav = useav
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
      
      ! Find last non-zero timestep (arrays may be zero-padded)
      max_tstep = timepoints
      do while (max_tstep > 1 .and. ltime(max_tstep) == 0.0d0)
        max_tstep = max_tstep - 1
      end do
      
      if (useav) then
        if (.not. allocated(lav)) then
          allocate (lav(timepoints))
        end if
        lav(:) = avgrid
      end if
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
      targettime = ltime(tstep) * SECONDS_PER_YEAR

      ! Use the first profile values as the starting physics for ALL depth points
      ! so that the initial output (time=0) reflects the provided tracer grids.
      density(:) = ldens(tstep)
      gastemp(:) = lgtemp(tstep)
      dusttemp(:) = ldtemp(tstep)
      radfield = lradfield(tstep)
      zeta = lzeta(tstep)
      if (luseav) then
        av(:) = lav(tstep)
      else if (lusecoldens) then
        coldens(:) = lnh(tstep)
        av(:) = 5.348e-22 * coldens(:)
      end if 

      ! Also ensure the current depth step has consistent values (dstep may not be 1)
      density(dstep) = ldens(tstep)
      gastemp(dstep) = lgtemp(tstep)
      dusttemp(dstep) = ldtemp(tstep)
      radfield = lradfield(tstep)
      zeta = lzeta(tstep)
      if (luseav) then
        av(dstep) = lav(tstep)
      else if (lusecoldens) then
        coldens(dstep) = lnh(tstep)
        av(dstep) = 5.348e-22 * coldens(dstep)
      end if
      ! Set final time to end of tracer histories (use last non-zero timestep)
      finaltime = timegrid(max_tstep)
    END SUBROUTINE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE updateTargetTime
      ! Set target time from postprocessing data (use exact dump time)
      ! Note: ltime is in years, but targettime must be in seconds
      targettime = ltime(tstep) * SECONDS_PER_YEAR
      
      ! If we've reached the final timestep, force end of simulation by setting finalTime to current time
      IF (tstep .ge. max_tstep) THEN
        finaltime = ltime(max_tstep)
      END IF

    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics

      ! Update physical properties to values at tstep (== targettime)
      targettime = ltime(tstep)
      density(dstep) = ldens(tstep)

      

      ! All checks passed - assign values to multi point arrays
      gastemp(dstep) = lgtemp(tstep)
      dusttemp(dstep) = ldtemp(tstep)
      density(dstep) = ldens(tstep)
      ! And single values:
      radfield = lradfield(tstep)
      zeta = lzeta(tstep)

      if (luseav) then
        av(dstep) = lav(tstep)
      else if (lusecoldens) then
        coldens(dstep) = lnh(tstep)
        av(dstep) = 5.348e-22 * coldens(dstep)
      end if 

      IF ((density(dstep) .ne. density(dstep)) .OR. (density(dstep) .le. 0.0d0) .OR. &
          (gastemp(dstep) .ne. gastemp(dstep)) .OR. (gastemp(dstep) < 1.0d0) .OR. &
          (dusttemp(dstep) .ne. dusttemp(dstep)) .OR. (dusttemp(dstep) < 1.0d0)) THEN
        write(*,*) 'POSTPROCESS_updatePhysics: FATAL invalid physics at tstep=', tstep, &
                  ' ldens=', ldens(tstep), ' lgtemp=', lgtemp(tstep), ' ldtemp=', ldtemp(tstep)
        postprocess_error = PHYSICS_UPDATE_ERROR
        RETURN
      END IF

      ! Increment tstep only if we haven't reached the final valid timestep
      IF (tstep < max_tstep) THEN
        tstep = tstep + 1
      END IF
    
      
    END SUBROUTINE updatePhysics

    SUBROUTINE sublimation(abund)
        ! This subroutine must be in every physics module so we dummy it here.
        !f2py integer, intent(aux) :: points
        DOUBLE PRECISION :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
    END SUBROUTINE sublimation    
END MODULE postprocess_mod


