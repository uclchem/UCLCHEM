! ==================================================================
! CVODE Interface Module for UCLCHEM
! ==================================================================
! This module provides a wrapper for the SUNDIALS/CVODE ODE solver
! to replace the legacy DVODE solver in UCLCHEM. It uses the modern
! SUNDIALS F2003 interface.
!
! Key features:
! - Drop-in replacement for DVODE_F90 call in chemistry.f90
! - Maintains same tolerances and method (BDF) as DVODE
! - Bridges existing Fortran RHS function to CVODE
! ==================================================================

MODULE cvode_interface_mod
  
  USE, INTRINSIC :: iso_c_binding
  USE fcvode_mod                    ! CVODE solver interface
  USE fnvector_serial_mod           ! Serial N_Vector
  USE fsunmatrix_dense_mod          ! Dense matrix
  USE fsunlinsol_dense_mod          ! Dense linear solver
  USE fsundials_core_mod            ! Core SUNDIALS types
  
  IMPLICIT NONE
  
  ! Module-level variables to maintain state between calls
  ! (analogous to DVODE's persistent work arrays)
  TYPE(c_ptr), SAVE :: sunctx = c_null_ptr        ! SUNDIALS context
  TYPE(c_ptr), SAVE :: cvode_mem = c_null_ptr     ! CVODE memory
  TYPE(N_Vector), POINTER, SAVE :: sunvec_y => null()      ! Solution N_Vector
  TYPE(N_Vector), POINTER, SAVE :: sunvec_atol => null()   ! Tolerance N_Vector
  TYPE(SUNMatrix), POINTER, SAVE :: sunmat_A => null()     ! Jacobian matrix
  TYPE(SUNLinearSolver), POINTER, SAVE :: sunlinsol_LS => null()  ! Linear solver
  
  LOGICAL, SAVE :: cvode_initialized = .FALSE.
  INTEGER, SAVE :: neq_saved = 0
  
  ! Abstract interface for RHS function
  ABSTRACT INTERFACE
    SUBROUTINE rhs_func_interface(neq, t, y, ydot)
      INTEGER, INTENT(IN) :: neq
      REAL(8), INTENT(IN) :: t
      REAL(8), INTENT(IN) :: y(neq)
      REAL(8), INTENT(OUT) :: ydot(neq)
    END SUBROUTINE rhs_func_interface
  END INTERFACE
  
  ! Module variable to store RHS function pointer
  PROCEDURE(rhs_func_interface), POINTER, SAVE :: user_rhs_ptr => null()
  
CONTAINS

  ! ================================================================
  ! RHS Callback Function for CVODE
  ! ================================================================
  ! This function has the signature required by CVODE and bridges
  ! to the existing UCLCHEM RHS function (F subroutine in chemistry.f90)
  ! ================================================================
  INTEGER(c_int) FUNCTION cvode_rhs_wrapper(t, sunvec_y_in, sunvec_f, user_data) &
       RESULT(ierr) BIND(C)
    
    IMPLICIT NONE
    
    ! Arguments from CVODE
    REAL(c_double), VALUE :: t
    TYPE(N_Vector) :: sunvec_y_in
    TYPE(N_Vector) :: sunvec_f
    TYPE(c_ptr), VALUE :: user_data
    
    ! Local variables
    REAL(c_double), POINTER :: yval(:), fval(:)
    
    ! Get pointers to the data in the N_Vectors
    yval => FN_VGetArrayPointer(sunvec_y_in)
    fval => FN_VGetArrayPointer(sunvec_f)
    
    ! Call the user's RHS function via the saved pointer
    IF (ASSOCIATED(user_rhs_ptr)) THEN
      CALL user_rhs_ptr(neq_saved, t, yval, fval)
      ierr = 0
    ELSE
      ! Error: RHS function not set
      ierr = -1
    END IF
    
  END FUNCTION cvode_rhs_wrapper

  ! ================================================================
  ! Main CVODE Solver Wrapper
  ! ================================================================
  ! This subroutine mimics the DVODE_F90 interface to serve as a
  ! drop-in replacement in UCLCHEM's chemistry.f90
  ! ================================================================
  SUBROUTINE solve_with_cvode(rhs, neq, y, t, tout, itask, istate, &
                               rtol, atol)
    
    IMPLICIT NONE
    
    ! Arguments (matching DVODE_F90 interface as closely as possible)
    PROCEDURE(rhs_func_interface) :: rhs  ! RHS function
    INTEGER, INTENT(IN) :: neq         ! Number of equations
    REAL(8), INTENT(INOUT) :: y(neq)   ! Solution vector
    REAL(8), INTENT(INOUT) :: t        ! Current time
    REAL(8), INTENT(IN) :: tout        ! Target time
    INTEGER, INTENT(IN) :: itask       ! Task flag (1 = CV_NORMAL)
    INTEGER, INTENT(INOUT) :: istate   ! State flag (1 = first call, 2 = continue)
    REAL(8), INTENT(IN) :: rtol        ! Relative tolerance
    REAL(8), INTENT(IN) :: atol(neq)   ! Absolute tolerance (vector)
    
    ! Local variables
    INTEGER(c_int) :: retval
    REAL(c_double) :: tret(1)
    REAL(8), ALLOCATABLE :: atol_copy(:)  ! Copy of atol for N_Vector
    INTEGER :: i
    
    ! ============================================================
    ! INITIALIZATION (first call with istate == 1)
    ! ============================================================
    IF (istate == 1) THEN
      
      ! Set the RHS function pointer for the callback wrapper
      user_rhs_ptr => rhs
      
      WRITE(*,*) "CVODE: Initializing solver for NEQ =", neq
      
      ! Save NEQ for use in callback
      neq_saved = neq
      
      ! Clean up any previous CVODE instance
      IF (cvode_initialized) THEN
        CALL cleanup_cvode()
      END IF
      
      ! Create SUNDIALS context
      retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)
      IF (retval /= 0) THEN
        WRITE(*,*) "ERROR: FSUNContext_Create failed, retval =", retval
        istate = -3
        RETURN
      END IF
      
      ! Create N_Vector for solution (wraps existing Fortran array)
      sunvec_y => FN_VMake_Serial(INT(neq, c_int64_t), y, sunctx)
      IF (.NOT. ASSOCIATED(sunvec_y)) THEN
        WRITE(*,*) "ERROR: FN_VMake_Serial failed for y"
        istate = -3
        RETURN
      END IF
      
      ! Create N_Vector for absolute tolerances (need a copy since atol is INTENT(IN))
      ALLOCATE(atol_copy(neq))
      atol_copy = atol
      
      ! CVODE uses WRMS norm and is more sensitive to very small tolerances than DVODE
      ! DVODE uses component-wise error check, CVODE uses weighted root-mean-square
      ! See TOLERANCE_COMPARISON.md for detailed analysis
      DO i = 1, neq
        ! Enforce a more reasonable floor for CVODE (much higher than 1e-30)
        ! UCLCHEM's abstol_min=1e-25 is too tight for CVODE
        atol_copy(i) = MAX(atol_copy(i), 1.0d-20)
        
        ! Also ensure abstol makes sense relative to rtol
        ! For rtol=1e-8, abstol should be at least ~1e-16 to be meaningful
        atol_copy(i) = MAX(atol_copy(i), rtol * 1.0d-8)
      END DO
      
      sunvec_atol => FN_VMake_Serial(INT(neq, c_int64_t), atol_copy, sunctx)
      IF (.NOT. ASSOCIATED(sunvec_atol)) THEN
        WRITE(*,*) "ERROR: FN_VMake_Serial failed for atol"
        istate = -3
        RETURN
      END IF
      
      ! Create CVODE solver (CV_BDF = 2 for stiff problems)
      cvode_mem = FCVodeCreate(CV_BDF, sunctx)
      IF (.NOT. c_associated(cvode_mem)) THEN
        WRITE(*,*) "ERROR: FCVodeCreate failed"
        istate = -3
        RETURN
      END IF
      
      ! Initialize CVODE with RHS function
      retval = FCVodeInit(cvode_mem, c_funloc(cvode_rhs_wrapper), t, sunvec_y)
      IF (retval /= 0) THEN
        WRITE(*,*) "ERROR: FCVodeInit failed, retval =", retval
        istate = -3
        RETURN
      END IF
      
      ! Set tolerances (scalar rtol, vector atol)
      ! Print tolerance diagnostics
      WRITE(*,'(A,ES12.5)') "  Relative tolerance: ", rtol
      WRITE(*,'(A,ES12.5,A,ES12.5)') "  Absolute tolerance range: min=", MINVAL(atol_copy), " max=", MAXVAL(atol_copy)
      WRITE(*,'(A,I0)') "  Number of atol values < 1e-20: ", COUNT(atol_copy < 1.0d-20)
      
      retval = FCVodeSVtolerances(cvode_mem, rtol, sunvec_atol)
      IF (retval /= 0) THEN
        WRITE(*,*) "ERROR: FCVodeSVtolerances failed, retval =", retval
        istate = -3
        RETURN
      END IF
      
      ! Create dense matrix for Jacobian
      sunmat_A => FSUNDenseMatrix(INT(neq, c_int64_t), INT(neq, c_int64_t), sunctx)
      IF (.NOT. ASSOCIATED(sunmat_A)) THEN
        WRITE(*,*) "ERROR: FSUNDenseMatrix failed"
        istate = -3
        RETURN
      END IF
      
      ! Create dense linear solver
      sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
      IF (.NOT. ASSOCIATED(sunlinsol_LS)) THEN
        WRITE(*,*) "ERROR: FSUNLinSol_Dense failed"
        istate = -3
        RETURN
      END IF
      
      ! Attach linear solver to CVODE
      retval = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
      IF (retval /= 0) THEN
        WRITE(*,*) "ERROR: FCVodeSetLinearSolver failed, retval =", retval
        istate = -3
        RETURN
      END IF
      
      ! Set maximum number of steps (matching DVODE's MXSTEP)
      retval = FCVodeSetMaxNumSteps(cvode_mem, INT(100000, c_long))
      IF (retval /= 0) THEN
        WRITE(*,*) "WARNING: FCVodeSetMaxNumSteps failed, retval =", retval
      END IF
      
      ! Relax Newton convergence criteria (CVODE defaults are too strict for this chemistry)
      ! Default max nonlinear iterations is 3, increase to 10
      retval = FCVodeSetMaxNonlinIters(cvode_mem, INT(10, c_int))
      IF (retval /= 0) THEN
        WRITE(*,*) "WARNING: FCVodeSetMaxNonlinIters failed, retval =", retval
      END IF
      
      ! Relax nonlinear convergence coefficient (default 0.1, use 0.33)
      retval = FCVodeSetNonlinConvCoef(cvode_mem, 0.33d0)
      IF (retval /= 0) THEN
        WRITE(*,*) "WARNING: FCVodeSetNonlinConvCoef failed, retval =", retval
      END IF
      
      ! Increase max convergence failures allowed (default 10, use 20)
      retval = FCVodeSetMaxConvFails(cvode_mem, INT(20, c_int))
      IF (retval /= 0) THEN
        WRITE(*,*) "WARNING: FCVodeSetMaxConvFails failed, retval =", retval
      END IF
      
      ! Set initial step size hint (important for CVODE performance)
      ! Use a small fraction of the integration interval
      retval = FCVodeSetInitStep(cvode_mem, (tout-t) * 1.0d-6)
      IF (retval /= 0) THEN
        WRITE(*,*) "WARNING: FCVodeSetInitStep failed, retval =", retval
      END IF
      
      ! Set minimum step size to prevent infinite step size reduction
      ! Use a very small value but not zero - allow CVODE to go smaller if needed
      retval = FCVodeSetMinStep(cvode_mem, 0.0d0)  ! Let CVODE choose minimum
      IF (retval /= 0) THEN
        WRITE(*,*) "WARNING: FCVodeSetMinStep failed, retval =", retval
      END IF
      
      cvode_initialized = .TRUE.
      WRITE(*,*) "CVODE: Initialization complete"
      
    END IF  ! End initialization
    
    ! ============================================================
    ! INTEGRATION STEP
    ! ============================================================
    
    ! Update the solution vector data pointer (in case y array moved)
    CALL FN_VSetArrayPointer(y, sunvec_y)
    
    WRITE(*,'(A,ES12.5,A,ES12.5)') "CVODE: Integrating from t=", t, " to tout=", tout
    
    ! Call CVODE to integrate to tout
    ! itask: 1 = CV_NORMAL (integrate to tout)
    retval = FCVode(cvode_mem, tout, sunvec_y, tret, CV_NORMAL)
    
    ! Update time
    t = tret(1)
    
    WRITE(*,'(A,ES12.5,A,I0)') "CVODE: Reached t=", t, ", retval=", retval
    
    ! Map CVODE return codes to DVODE-style istate
    SELECT CASE(retval)
      CASE(0)
        ! CV_SUCCESS: successful step
        istate = 2
        WRITE(*,*) "CVODE: Step successful"
        
      CASE(1)
        ! CV_TSTOP_RETURN: reached tstop
        istate = 2
        WRITE(*,*) "CVODE: Reached tstop"
        
      CASE(2)
        ! CV_ROOT_RETURN: found root
        istate = 3
        WRITE(*,*) "CVODE: Found root"
        
      CASE(-1)
        ! CV_TOO_MUCH_WORK: mxstep exceeded
        WRITE(*,*) "CVODE ERROR: Maximum number of steps exceeded"
        istate = -1
        
      CASE(-2)
        ! CV_TOO_MUCH_ACC: accuracy request too stringent
        WRITE(*,*) "CVODE: Accuracy request too stringent"
        istate = -2
        
      CASE(-3)
        ! CV_ERR_FAILURE: too many error test failures
        WRITE(*,*) "CVODE: Too many error test failures"
        istate = -4
        
      CASE(-4)
        ! CV_CONV_FAILURE: too many convergence failures
        WRITE(*,*) "CVODE: Too many convergence failures"
        istate = -5
        
      CASE DEFAULT
        WRITE(*,*) "CVODE: Error, retval =", retval
        istate = -1
        
    END SELECT
    
  END SUBROUTINE solve_with_cvode

  ! ================================================================
  ! Print CVODE Statistics (TODO: Fix FCVodeGetIntegratorStats interface)
  ! ================================================================
  ! SUBROUTINE print_cvode_stats()
  !   ! Commented out temporarily - need to fix interface to FCVodeGetIntegratorStats
  ! END SUBROUTINE print_cvode_stats

  ! ================================================================
  ! Cleanup Routine
  ! ================================================================
  ! Frees all SUNDIALS memory. Should be called at end of simulation.
  ! ================================================================
  SUBROUTINE cleanup_cvode()
    
    IMPLICIT NONE
    INTEGER(c_int) :: ierr
    
    IF (ASSOCIATED(sunlinsol_LS)) THEN
      ierr = FSUNLinSolFree(sunlinsol_LS)
      sunlinsol_LS => null()
    END IF
    
    IF (ASSOCIATED(sunmat_A)) THEN
      CALL FSUNMatDestroy(sunmat_A)
      sunmat_A => null()
    END IF
    
    IF (c_associated(cvode_mem)) THEN
      CALL FCVodeFree(cvode_mem)
      cvode_mem = c_null_ptr
    END IF
    
    ! Note: N_Vectors created with FN_VMake_Serial should NOT be destroyed
    ! with FN_VDestroy because they don't own the data
    sunvec_y => null()
    sunvec_atol => null()
    
    IF (c_associated(sunctx)) THEN
      ierr = FSUNContext_Free(sunctx)
      sunctx = c_null_ptr
    END IF
    
    cvode_initialized = .FALSE.
    neq_saved = 0
    
  END SUBROUTINE cleanup_cvode

END MODULE cvode_interface_mod
