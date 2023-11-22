! Module for post-processing data from hydrodynamical simulations
! Requires modifications to some core functions of the code
MODULE postprocess_mod
    USE physicscore
    USE network
    USE constants
    IMPLICIT NONE

    character(len=100) :: trajecfile
    logical :: lgpost=.false.
    integer,parameter :: tfid=66
    integer :: ntime,tstep
    double precision,allocatable :: timegrid(:),densgrid(:,:),tempgrid(:,:),nhgrid(:,:),nh2grid(:,:),ncogrid(:,:)
    
CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Called at start of UCLCHEM run
    ! Uses values in defaultparamters.f90 and any inputs to set initial values        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics(successFlag)
      INTEGER, INTENT(OUT) :: successFlag
      integer ::i,j,k
      double precision :: junk
      
      successFlag=1

      lgpost = .true.

        cloudSize=0. ! Shielding column densities supplied separately
        endatfinaldensity = .false. ! Needs to end at final time to work properly
        freefall = .false. ! Won't work with freefall on

        allocate(timegrid(ntime),densgrid(points,ntime),tempgrid(points,ntime),nhgrid(points,ntime),nh2grid(points,ntime),ncogrid(points,ntime))

        open(unit=tfid,file=trim(trajecfile),status='old')

        ! Read tracer particle histories from file
        do i=1,points
           do j=1,ntime
              read(tfid,*) timegrid(j),(junk,k=1,3),densgrid(i,j),(junk,k=1,3),tempgrid(i,j),(junk,k=1,2),nhgrid(i,j),nh2grid(i,j),ncogrid(i,j)
           end do
        end do

        close(unit=tfid)

        ! Initialise values to t=1
        tstep = 1

        density(:) = densgrid(:,tstep)
        gastemp(:) = tempgrid(:,tstep)
        dusttemp(:) = tempgrid(:,tstep) ! No separate dust temperature as yet...
        coldens(:) = nhgrid(:,tstep)
        av(:) = 5.348e-22 * coldens(:)

        ! Set final time to end of tracer histories
        finaltime = timegrid(ntime)/seconds_per_year

    END SUBROUTINE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE updateTargetTime

      ! Set target time from postprocessing data
      targettime = timegrid(tstep) + 1.*seconds_per_year

      write(*,"('Integrating chemistry to timestep ',I3,' ',ES10.3,' years')") tstep,targettime/seconds_per_year

    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics

      ! Update physical properties to values at tstep (== targettime)
      density(dstep) = densgrid(dstep,tstep)
      gastemp(dstep) = tempgrid(dstep,tstep)
      dusttemp(dstep) = tempgrid(dstep,tstep) ! No separate dust temperature as yet...
      coldens(dstep) = nhgrid(dstep,tstep)
      av(dstep) = 5.348e-22 * coldens(dstep)

      ! If this is the last point, update tstep to move onto next time dump
      if (dstep .eq. points) tstep = tstep + 1
      
    END SUBROUTINE updatePhysics

    SUBROUTINE sublimation(abund)
    ! This subroutine must be in every physics module so we dummy it here.
        DOUBLE PRECISION :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
    END SUBROUTINE sublimation    
END MODULE postprocess_mod


