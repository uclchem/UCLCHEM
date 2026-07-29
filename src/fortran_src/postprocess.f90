! Module for post-processing data from hydrodynamical simulations
! Requires modifications to some core functions of the code
module postprocess_mod
    use constants, only: SECONDS_PER_YEAR, PHYSICS_UPDATE_ERROR
    use DEFAULTPARAMETERS
    use f2py_constants, only: nSpec
    !f2py INTEGER, parameter :: dp
    use physicscore
    implicit none

    ! character(len=100) :: trajecfile
    logical :: lusecoldens
    logical :: luseav
    logical :: usepostprocess = .true.
    ! integer,parameter :: tfid=66
    integer :: tstep
    integer :: max_tstep  ! Last valid (non-zero) timestep
    integer :: postprocess_error  ! Non-zero if a fatal error occurred during updatePhysics
    real(dp), allocatable, dimension(:) :: ltime, ldens, lra, lzeta, lradfield, lgtemp, ldtemp
    real(dp), allocatable, dimension(:) :: lnh, lnh2, lnco, lnc, lav

    public

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called at start of UCLCHEM run
    ! Uses values in defaultparamters.f90 and any inputs to set initial values        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initializePhysics(successFlag, timegrid, densgrid, radgrid, zetagrid, gtempgrid,&
        &dtempgrid, useav, avgrid, usecoldens, timepoints, nhgrid, nh2grid, ncogrid, ncgrid)
      integer, intent(out) :: successFlag
      integer, intent(in) :: timepoints
      real(dp), intent(in), dimension(timePoints) :: timegrid
      real(dp), intent(in), dimension(timePoints) :: densgrid
      real(dp), intent(in), dimension(timePoints) :: radgrid
      real(dp), intent(in), dimension(timePoints) :: zetagrid
      real(dp), intent(in), dimension(timePoints) :: gtempgrid
      real(dp), intent(in), dimension(timePoints) :: dtempgrid
      logical, intent(in) :: useav
      real(dp), intent(in), dimension(timePoints), optional :: avgrid
      logical, intent(in) :: usecoldens
      real(dp), intent(in), optional, dimension(timePoints) :: nhgrid
      real(dp), intent(in), optional, dimension(timePoints) :: nh2grid
      real(dp), intent(in), optional, dimension(timePoints) :: ncogrid
      real(dp), intent(in), optional, dimension(timePoints) :: ncgrid


      successFlag=0
      postprocess_error = 0
      ! Create a local variable that can be recycled in updatePhysics
      lusecoldens = usecoldens
      luseav = useav
      ! Check if NHgrid is present, if so, we need to use the custom column densities
      if (lusecoldens) then
        cloudSize=0.0  ! Shielding column densities supplied separately
      end if
      freefall = .false.  ! Won't work with freefall on

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
      do while (max_tstep > 1 .and. ltime(max_tstep) == 0.0_dp)
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

      ! Initialize values to t=1 and overwrite them.
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
        av(:) = 5.348e-22_dp * coldens(:)
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
        av(dstep) = 5.348e-22_dp * coldens(dstep)
      end if
      ! Set final time to end of tracer histories (use last non-zero timestep)
      finaltime = timegrid(max_tstep)
    end subroutine initializePhysics

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine updateTargetTime
      ! Set target time from postprocessing data (use exact dump time)
      ! Note: ltime is in years, but targettime must be in seconds
      targettime = ltime(tstep) * SECONDS_PER_YEAR

      ! If we've reached the final timestep, force end of simulation by setting finalTime to current time
      if (tstep >= max_tstep) then
        finaltime = ltime(max_tstep)
      end if

    end subroutine updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine updatePhysics

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
        av(dstep) = 5.348e-22_dp * coldens(dstep)
      end if

      if ((density(dstep) /= density(dstep)) .OR. (density(dstep) <= 0.0_dp) .OR. &
          (gastemp(dstep) /= gastemp(dstep)) .OR. (gastemp(dstep) < 1.0_dp) .OR. &
          (dusttemp(dstep) /= dusttemp(dstep)) .OR. (dusttemp(dstep) < 1.0_dp)) then
        write(*,*) "POSTPROCESS_updatePhysics: FATAL invalid physics at tstep=", tstep, &
                  " ldens=", ldens(tstep), " lgtemp=", lgtemp(tstep), " ldtemp=", ldtemp(tstep)
        postprocess_error = PHYSICS_UPDATE_ERROR
        return
      end if

      ! Increment tstep only if we haven't reached the final valid timestep
      if (tstep < max_tstep) then
        tstep = tstep + 1
      end if


    end subroutine updatePhysics

    subroutine sublimation(abund, lpoints)
        ! This subroutine must be in every physics module so we dummy it here.
        integer, intent(in) :: lpoints
        real(dp), intent(inout) :: abund(nspec+1,lpoints)
    end subroutine sublimation
end module postprocess_mod


