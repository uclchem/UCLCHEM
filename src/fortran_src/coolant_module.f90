module COOLANT_MODULE
   use constants, only: dp, PI, K_BOLTZ, C, MIN_ABUND, HP, T_CMB, MH, &
       COOLANT_SOLVER_ERROR, COOLANT_POP_TOL_ERROR, COOLANT_CONFIG_ERROR, &
       COOLANT_DATA_ERROR, COOLANT_FILE_ERROR, COOLANT_FREQ_TOL_ERROR

   use defaultparameters, only: freq_rel_tol, pop_rel_tol, negative_abundance_tol
   use F2PY_CONSTANTS, only: NCOOLANTS, coolantFiles, coolantNames, coolantDataDir, &
                              coolantConversionFactors, coolantConversionMode, coolantParentNames, &
                              N_TOTAL_LEVELS
   use network, only: nH, nH2, nHe, nHx, nelec
   implicit none

   public

   integer, parameter :: populationsId = 53

   integer, parameter :: CACHE_SIZE = 10
   ! Tolerances are provided via the DEFAULTPARAMETERS module (freq_rel_tol, pop_rel_tol).

   !  Sparse per-partner collisional rate storage (replaces dense C_COEFF(7,NLEVEL,NLEVEL,1000))
   type partner_data_t
      integer  :: partner_id           ! original 1-7 LAMDA partner identifier
      integer  :: ntemp                ! actual number of temperatures from file
      integer  :: ncoll                ! number of stored transitions (forward + reverse)
      real(dp), allocatable :: temperature(:)   ! (ntemp)
      integer,  allocatable :: i_idx(:)         ! (ncoll) upper level index
      integer,  allocatable :: j_idx(:)         ! (ncoll) lower level index
      real(dp), allocatable :: c_coeff(:,:)     ! (ncoll, ntemp) rate coefficients
   end type partner_data_t

!  Specify the properties that define each coolant species
   type COOLANT_TYPE

      character(LEN=20) :: NAME  ! Coolant species name
      character(LEN=256) :: FILENAME  ! Name of the coolant data file

      integer :: INDEX  ! Index number of the coolant species
      integer :: NLEVEL  ! Number of levels in the system
      integer :: NPARTNER  ! Number of collision partners present in file
      integer :: partner_map(NCOOLANTS)  ! partner_map(id) -> index in partners(:); 0 if absent

      real(dp) :: MOLECULAR_MASS, previousCooling  ! Molecular mass of the coolant species

      real(dp), allocatable :: ENERGY(:),WEIGHT(:)  ! Energy (K) and statistical weight of each level
      real(dp), allocatable :: A_COEFF(:,:),B_COEFF(:,:)  ! Einstein A and B coefficients for each transition between levels
      real(dp), allocatable :: FREQUENCY(:,:)  ! Frequency (Hz) of each transition between levels
      type(partner_data_t), allocatable :: partners(:)  ! (npartner) sparse per-partner collisional data

      logical :: CONVERGED  ! Flag indicating whether the level populations of all particles have converged

      !Below were stolen from population type since I think we just want one coolant.
      real(dp) :: DENSITY  ! Total number density of the species (cm^-3)
      real(dp) :: LINEWIDTH  ! Doppler line width of the species (cm s^-1)


      real(dp), allocatable :: POPULATION(:)  ! Population density (cm^-3) of each level
      real(dp), allocatable :: PREVIOUS_POPULATION(:)  ! Population density calculated at the previous iteration step
      real(dp), allocatable :: EMISSIVITY(:,:)  ! Local emissivity (erg cm^-3 s^-1) of each transition
      real(dp), allocatable :: OPACITY(:,:)  ! Optical depth of each transition along each HEALPix ray to the PDR surface (or simulation boundary)
      real(dp), allocatable :: LAMBDA(:,:)  ! Lambda operator value for each transition

      ! Persistent workspace arrays to avoid repeated allocations (added for performance)
      real(dp), allocatable :: R(:,:)  ! Transition rate matrix (s^-1)
      real(dp), allocatable :: A(:,:)  ! Coefficient matrix for SE equations
      real(dp), allocatable :: B(:)  ! Right-hand-side vector for SE equations
      real(dp), allocatable :: RADIATION_FIELD(:,:)  ! Mean integrated radiation field
      real(dp), allocatable :: COLLISIONAL_RATE(:,:)  ! Collisional rate coefficients
      integer, allocatable :: IPIV(:), INDEX_ROW(:), INDEX_COL(:)  ! Gauss-Jordan workspace

      ! Collisional rate cache (multi-parameter, 1% relative tolerance)
      ! Cache has 10 entries per coolant, matched by all physical parameters
      real(dp) :: CACHE_TOLERANCE       ! Relative tolerance for cache matching, default 0.01 (1%)
      real(dp) :: CACHED_TEMPERATURE(CACHE_SIZE)  ! Cached temperature (K)
      real(dp) :: CACHED_N_H2(CACHE_SIZE)         ! Cached n(H2) = density * abundance(H2)
      real(dp) :: CACHED_N_ELEC(CACHE_SIZE)       ! Cached n(e-) = density * abundance(e-)
      real(dp) :: CACHED_N_H(CACHE_SIZE)          ! Cached n(H) = density * abundance(H)
      real(dp) :: CACHED_N_HE(CACHE_SIZE)         ! Cached n(He) = density * abundance(He)
      real(dp) :: CACHED_N_HPLUS(CACHE_SIZE)      ! Cached n(H+) = density * abundance(H+)
      real(dp), allocatable :: CACHED_COLLISIONAL_RATE(:,:,:)  ! (10, NLEVEL, NLEVEL)
      integer :: CACHE_INDEX  ! Current cache position (round-robin)

   end type COOLANT_TYPE



   type(COOLANT_TYPE), allocatable :: coolants(:)
   ! NCOOLANTS, coolantFiles, and coolantNames are imported from F2PY_CONSTANTS
   ! They are generated at build time by MakeRates based on user configuration
   integer :: coolantIndices(NCOOLANTS)

   ! Coolant population restart mode control
   ! 0 = WARM (default): Rescale existing populations when density changes, initialize to LTE on first call
   ! 1 = FORCE_LTE: Always reset to LTE before SE iteration
   ! 2 = FORCE_GROUND: Always reset to ground state before SE iteration
   !f2py INTEGER :: coolant_restart_mode
   integer :: coolant_restart_mode = 0

   ! Track whether populations have been initialized (per coolant)
   logical :: coolant_populations_initialized = .false.
   real(dp) :: CLOUD_DENSITY,CLOUD_COLUMN,CLOUD_SIZE

   ! Temperature at which level populations were last (re)solved (mode 0 only)
   real(dp) :: last_levpop_temperature = -1.0_dp
   ! Relative temperature change threshold to trigger a re-solve in mode 0 (default 1%)
   real(dp) :: coolant_temp_recompute_threshold = 1e-2_dp
   ! Flag: when .TRUE., the next call to MANAGE_COOLANT_POPULATIONS forces a re-solve
   logical :: coolant_levpop_force_recompute = .false.

   ! Module-level error state for runtime errors (checked in time loop)
   integer :: coolant_error_flag = 0
   character(LEN=256) :: coolant_error_message = ""


contains
   !=======================================================================
   !  STOLEN FROM UCL-PDR
   !  Read in the coolant species, their energy level structure, transition
   !  properties and collisional rate coefficients. The specified files are
   !  assumed to contain entries in the LAMDA/RADEX format, allowing files
   !  to be downloaded directly from the online database.
   !
   !-----------------------------------------------------------------------
   subroutine READ_COOLANTS(successFlag)
      integer, intent(inout) :: successFlag
      integer :: I,J,K,L,M,N,INDEX,IER
      integer :: NLEVEL,NLINE,NTEMP,NPARTNER,NCOLL,PARTNER_ID,NCNT
      real(dp), allocatable :: RATE_BUF(:)
      logical,  allocatable :: SEEN(:,:)
      integer :: actual_total_levels
      integer, parameter :: coolantID = 81

!      WRITE(*,*) 'DEBUG: READ_COOLANTS starting. NCOOLANTS=', NCOOLANTS
!      WRITE(*,*) 'DEBUG: coolantDataDir=', TRIM(coolantDataDir)

      if (ALLOCATED(coolants)) deallocate(coolants)
!      WRITE(*,*) 'DEBUG: Allocating coolants array of size ', NCOOLANTS
      allocate(coolants(NCOOLANTS))
      do N=1,NCOOLANTS  ! Loop over coolants
         coolants(N)%FILENAME=coolantFiles(N)
   !     Open the input file
!         WRITE(*,*) 'DEBUG: Attempting to open coolant file ', TRIM(coolants(N)%FILENAME)
         ! WRITE(*,*) "Trying to open: ", TRIM(dataDir)//coolants(N)%FILENAME
         ! OPEN(UNIT=1,FILE='Datafiles/Collisional-Rates/'//coolants(N)%FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')
         ! TODO: Remove magic number 1
         open(unit=coolantID, FILE=TRIM(coolantDataDir)//coolants(N)%FILENAME, IOSTAT=IER, ACTION="READ", STATUS="OLD")
         read(coolantID,*,IOSTAT=IER)  ! Skip the first comment line

   !     Produce an error message if the file does not exist (or cannot be opened for whatever reason)
         if(IER/=0) then
            write(*,*) "ERROR! Cannot open coolant data file ",TRIM(coolants(N)%FILENAME), &
                & " for input. I tried to open:", TRIM(coolantDataDir)//coolants(N)%FILENAME
            close(coolantID)
            successFlag = COOLANT_FILE_ERROR
            return
         end if

         read(coolantID,*,IOSTAT=IER) coolants(N)%NAME  ! Read the name of the coolant
!         WRITE(*,*) 'DEBUG: Read coolant name = ', TRIM(coolants(N)%NAME)
         read(coolantID,*,IOSTAT=IER)
!         WRITE(*,*) 'DEBUG: Skipped comment line (molecular weight header)'
         read(coolantID,*,IOSTAT=IER) coolants(N)%MOLECULAR_MASS  ! Read the molecular mass
!         WRITE(*,*) 'DEBUG: Read coolant molecular mass = ', coolants(N)%MOLECULAR_MASS
         read(coolantID,*,IOSTAT=IER)
!         WRITE(*,*) 'DEBUG: Skipped comment line (number of energy levels header)'
         coolants(N)%INDEX=0  ! Initialize the coolant species index (assigned later)

   !     Read the number of levels and allocate the energy, statistical weight,
   !     Einstein A & B coefficient and transition frequency arrays accordingly
!         WRITE(*,*) 'DEBUG: About to read NLEVEL'
         read(coolantID,*,IOSTAT=IER) NLEVEL
!         WRITE(*,*) 'DEBUG: Read NLEVEL =', NLEVEL
         read(coolantID,*,IOSTAT=IER)
         if(NLEVEL<2) then
            write(*,*) "ERROR! Incorrect number of energy levels in coolant data file ",&
                        &TRIM(coolants(N)%FILENAME)," (NLEVEL=",NLEVEL,")"
            close(coolantID)
            successFlag = COOLANT_DATA_ERROR
            return
         end if
         coolants(N)%NLEVEL=NLEVEL
         allocate(coolants(N)%ENERGY(1:NLEVEL))
         allocate(coolants(N)%WEIGHT(1:NLEVEL))
         allocate(coolants(N)%A_COEFF(1:NLEVEL,1:NLEVEL))
         allocate(coolants(N)%B_COEFF(1:NLEVEL,1:NLEVEL))
         allocate(coolants(N)%FREQUENCY(1:NLEVEL,1:NLEVEL))

   !     Initialize the level energies, statistical weights,
   !     Einstein coefficients and transition frequencies
         coolants(N)%ENERGY=0.0_dp
         coolants(N)%WEIGHT=0.0_dp
         coolants(N)%A_COEFF=0.0_dp
         coolants(N)%B_COEFF=0.0_dp
         coolants(N)%FREQUENCY=0.0_dp

   !     Read the energy (cm^-1) and statistical weight of each level
         do L=1,NLEVEL  ! Loop over levels
            read(coolantID,*,IOSTAT=IER) I,coolants(N)%ENERGY(I),coolants(N)%WEIGHT(I)
            coolants(N)%ENERGY(I)=coolants(N)%ENERGY(I)*C*HP  ! Convert from cm^-1 to erg
         end do  ! End of loop over levels
         read(coolantID,*,IOSTAT=IER)

   !     Read the Einstein A coefficient (s^-1) and frequency (GHz) of each radiative transition
         read(coolantID,*,IOSTAT=IER) NLINE
         read(coolantID,*,IOSTAT=IER)
         do L=1,NLINE  ! Loop over radiative transitions
            read(coolantID,*,IOSTAT=IER) INDEX,I,J,coolants(N)%A_COEFF(I,J),coolants(N)%FREQUENCY(I,J)
            coolants(N)%FREQUENCY(I,J)=coolants(N)%FREQUENCY(I,J)*1.0e9_dp  ! Convert from GHz to Hz
   !        Calculate the Einstein B coefficient using B_ij = A_ij/(2.h.nu^3/c^2)
            coolants(N)%B_COEFF(I,J)=coolants(N)%A_COEFF(I,J)/(2*HP*(coolants(N)%FREQUENCY(I,J)**3)/(C**2))
   !        Calculate the Einstein B coefficient for the reverse transition from detailed balance
            coolants(N)%B_COEFF(J,I)=coolants(N)%B_COEFF(I,J)*(coolants(N)%WEIGHT(I)/coolants(N)%WEIGHT(J))
         end do  ! End of loop over radiative transitions
         read(coolantID,*,IOSTAT=IER)

   !     Calculate the transition frequencies between all levels (even if forbidden)
         do I=1,NLEVEL
            do J=1,NLEVEL
   !           Check that the calculated and measured frequencies differ by <1.0%
   !           Produce an error message if the difference between them is greater
               if(coolants(N)%FREQUENCY(I,J)/=0.0_dp) then
                  if(ABS(coolants(N)%FREQUENCY(I,J)-ABS(coolants(N)%ENERGY(I)-coolants(N)%ENERGY(J))/HP) &
                      & /coolants(N)%FREQUENCY(I,J)>freq_rel_tol) then
                     write(*,*) "ERROR! Calculated frequency differs from measured frequency beyond configured tolerance."
                     write(*,"('Coolant: ',A)") TRIM(coolants(N)%NAME)
                     write(*,"('Tolerance (fraction)=',F10.6)") freq_rel_tol
                     write(*,"(1PD12.5,'Hz vs',1PD12.5,'Hz')") ABS(coolants(N)%ENERGY(I)-coolants(N)%ENERGY(J))/HP, &
                                                             & coolants(N)%FREQUENCY(I,J)
                     close(coolantID)
                     successFlag = COOLANT_FREQ_TOL_ERROR
                     return
                  end if
               else
                  coolants(N)%FREQUENCY(I,J)=ABS(coolants(N)%ENERGY(I)-coolants(N)%ENERGY(J))/HP
               end if
            end do
         end do

   !     Read the collisional rate coefficients (cm^3 s^-1) for each collision partner
   !     Sparse storage: only actual transitions and temperature points are allocated
         read(coolantID,*,IOSTAT=IER) NPARTNER
         coolants(N)%NPARTNER = NPARTNER
         allocate(coolants(N)%partners(NPARTNER))
         coolants(N)%partner_map = 0

         do L=1,NPARTNER  ! Loop over collision partners
            read(coolantID,*,IOSTAT=IER)           ! skip comment
            read(coolantID,*,IOSTAT=IER) PARTNER_ID
            if(PARTNER_ID<1 .OR. PARTNER_ID>7) then
               write(*,*) "ERROR! Unrecognized collision partner ID in coolant data file ",&
                           &TRIM(coolants(N)%FILENAME)," (ID=",PARTNER_ID,")"
               close(coolantID)
               successFlag = COOLANT_DATA_ERROR
               return
            end if
            read(coolantID,*,IOSTAT=IER)           ! skip comment
            read(coolantID,*,IOSTAT=IER) NCOLL
            read(coolantID,*,IOSTAT=IER)           ! skip comment
            read(coolantID,*,IOSTAT=IER) NTEMP
            read(coolantID,*,IOSTAT=IER)           ! skip temperature comment

            coolants(N)%partners(L)%partner_id = PARTNER_ID
            coolants(N)%partners(L)%ntemp      = NTEMP
            coolants(N)%partner_map(PARTNER_ID) = L

            allocate(coolants(N)%partners(L)%temperature(NTEMP))
            ! Allocate for up to 2*NCOLL transitions (forward + detailed-balance reverse)
            allocate(coolants(N)%partners(L)%i_idx(2*NCOLL))
            allocate(coolants(N)%partners(L)%j_idx(2*NCOLL))
            allocate(coolants(N)%partners(L)%c_coeff(2*NCOLL, NTEMP))
            coolants(N)%partners(L)%c_coeff = 0.0_dp

            read(coolantID,*,IOSTAT=IER) (coolants(N)%partners(L)%temperature(K),K=1,NTEMP)
            read(coolantID,*,IOSTAT=IER)           ! skip collisions comment

            allocate(RATE_BUF(NTEMP))
            allocate(SEEN(NLEVEL, NLEVEL))
            SEEN = .false.
            NCNT = 0

            do M=1,NCOLL  ! Loop over collisional transitions
               read(coolantID,*,IOSTAT=IER) INDEX,I,J,(RATE_BUF(K),K=1,NTEMP)
   !           Store forward transition (I->J) if not already present
               if (.NOT. SEEN(I,J)) then
                  NCNT = NCNT + 1
                  coolants(N)%partners(L)%i_idx(NCNT) = I
                  coolants(N)%partners(L)%j_idx(NCNT) = J
                  coolants(N)%partners(L)%c_coeff(NCNT,1:NTEMP) = RATE_BUF(1:NTEMP)
                  SEEN(I,J) = .true.
               end if
   !           Compute reverse (J->I) via detailed balance if not already present
               if (.NOT. SEEN(J,I)) then
                  SEEN(J,I) = .true.
                  NCNT = NCNT + 1
                  coolants(N)%partners(L)%i_idx(NCNT) = J
                  coolants(N)%partners(L)%j_idx(NCNT) = I
                  do K=1,NTEMP
                     if (RATE_BUF(K)/=0.0_dp) then
                        coolants(N)%partners(L)%c_coeff(NCNT,K) = RATE_BUF(K) &
                           * (coolants(N)%WEIGHT(I) / coolants(N)%WEIGHT(J)) &
                           * EXP(-(coolants(N)%ENERGY(I) - coolants(N)%ENERGY(J)) &
                                  / (K_BOLTZ * coolants(N)%partners(L)%temperature(K)))
                     end if
                  end do
               end if
            end do  ! End of loop over collisional transitions

            coolants(N)%partners(L)%ncoll = NCNT
            deallocate(RATE_BUF)
            deallocate(SEEN)

         end do  ! End of loop over collision partners
         coolants(N)%previousCooling=0.0_dp
         close(coolantID)

      end do  ! End of loop over coolants


       do N=1,NCOOLANTS
         allocate(coolants(N)%POPULATION(1:coolants(N)%NLEVEL))
         allocate(coolants(N)%PREVIOUS_POPULATION(1:coolants(N)%NLEVEL))
         allocate(coolants(N)%EMISSIVITY(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         !ALLOCATE(coolants(N)%OPACITY(0:NRAYS-1,1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         allocate(coolants(N)%OPACITY(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         allocate(coolants(N)%LAMBDA(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         coolants(N)%POPULATION=0.0
         coolants(N)%PREVIOUS_POPULATION=0.0
         coolants(N)%EMISSIVITY=0.0
         coolants(N)%OPACITY=0.0
         coolants(N)%LAMBDA=0.0

         ! Allocate persistent workspace arrays to avoid repeated allocations
         allocate(coolants(N)%R(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         allocate(coolants(N)%A(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         allocate(coolants(N)%B(1:coolants(N)%NLEVEL))
         allocate(coolants(N)%RADIATION_FIELD(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         allocate(coolants(N)%COLLISIONAL_RATE(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         allocate(coolants(N)%IPIV(1:coolants(N)%NLEVEL))
         allocate(coolants(N)%INDEX_ROW(1:coolants(N)%NLEVEL))
         allocate(coolants(N)%INDEX_COL(1:coolants(N)%NLEVEL))
         coolants(N)%R=0.0_dp
         coolants(N)%A=0.0_dp
         coolants(N)%B=0.0_dp
         coolants(N)%RADIATION_FIELD=0.0_dp
         coolants(N)%COLLISIONAL_RATE=0.0_dp
         coolants(N)%IPIV=0
         coolants(N)%INDEX_ROW=0
         coolants(N)%INDEX_COL=0

         ! Allocate and initialize collisional rate cache (10 entries per coolant)
         allocate(coolants(N)%CACHED_COLLISIONAL_RATE(10,1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         coolants(N)%CACHED_COLLISIONAL_RATE=0.0_dp
         coolants(N)%CACHED_TEMPERATURE=-1.0_dp  ! -1 indicates empty cache entry
         coolants(N)%CACHED_N_H2=-1.0_dp
         coolants(N)%CACHED_N_ELEC=-1.0_dp
         coolants(N)%CACHED_N_H=-1.0_dp
         coolants(N)%CACHED_N_HE=-1.0_dp
         coolants(N)%CACHED_N_HPLUS=-1.0_dp
         coolants(N)%CACHE_TOLERANCE=0.01_dp  ! Default tolerance: 1% relative precision
         coolants(N)%CACHE_INDEX=1  ! Start at first cache position
      end do

      ! Validate that the actual total levels matches the compile-time constant.
      ! A mismatch means the coolant data files changed since MakeRates was run.
      actual_total_levels = 0
      do N = 1, NCOOLANTS
         actual_total_levels = actual_total_levels + coolants(N)%NLEVEL
      end do
      if (actual_total_levels /= N_TOTAL_LEVELS) then
         write(*,*) "ERROR: Runtime total energy levels (", actual_total_levels, &
                    ") does not match compile-time N_TOTAL_LEVELS (", N_TOTAL_LEVELS, ")."
         write(*,*) "The coolant data files have changed since MakeRates was run."
         write(*,*) "Fix: Re-run MakeRates and reinstall UCLCHEM."
         successFlag = COOLANT_DATA_ERROR
         return
      end if

   end subroutine READ_COOLANTS

   !=======================================================================
   !
   !  Calculate the Doppler line width of each coolant species for all
   !  particles. Contributions from both thermal and turbulent motions
   !  are included.
   !
   !-----------------------------------------------------------------------
   subroutine UPDATE_COOLANT_LINEWIDTHS(gasTemperature,turbVel)
      real(dp), intent(in) :: gasTemperature, turbVel
      real(dp) :: thermVel
      integer :: N

      if (.NOT. ALLOCATED(coolants)) then
         write(*,*) "ERROR: coolants not allocated!"
      end if

      !calculate constant factors of thermal velocity
      !keep it squared since we'll be adding in quadrature with turbVel
      thermVel=2.0_dp*K_BOLTZ*gasTemperature/MH

      ! Loop over coolants
      do N=1,NCOOLANTS
         COOLANTS(N)%LINEWIDTH = SQRT((thermVel/coolants(N)%MOLECULAR_MASS)+(turbVel*turbVel))  ! v_thermal = (2kT/m)^1/2
      end do  ! End of loop over coolants
   end subroutine UPDATE_COOLANT_LINEWIDTHS


   subroutine UPDATE_COOLANT_ABUNDANCES(gasDensity,gasTemperature,abundances)
      real(dp), intent(in) :: gasDensity,gasTemperature
      real(dp), intent(in) :: abundances(:)
      real(dp) :: fraction
      integer :: N
      do N=1,NCOOLANTS
         if (coolantNames(N)=="p-H2") then
            fraction=1.0_dp/(1.0_dp+GET_ORTHO_PARA_RATIO(gasTemperature))
            ! write(*,*) gasTemperature,"p",fraction
            coolants(N)%density=abundances(coolantIndices(N))*gasDensity*fraction
         else if (coolantNames(N) == "o-H2") then
             fraction=1.0_dp/(1.0_dp+1.0_dp/GET_ORTHO_PARA_RATIO(gasTemperature))
             coolants(N)%density=abundances(coolantIndices(N))*gasDensity*fraction
             ! write(*,*) gasTemperature,"o",fraction
         else if (coolantNames(N) == "o-H2O") then
            coolants(N)%density=abundances(coolantIndices(N))*gasDensity*0.75_dp
         else if (coolantNames(N) == "p-H2O") then
            coolants(N)%density=abundances(coolantIndices(N))*gasDensity*0.25_dp

         else
             coolants(N)%density=abundances(coolantIndices(N))*gasDensity
             !write(*,*) coolantNames(N),abundances(coolantIndices(N))
         end if

         ! Clamp tiny negative densities (solver noise) to 1e-30 floor
         if (coolants(N)%density < MIN_ABUND) coolants(N)%density = MIN_ABUND

         ! Sanity check: only error on clearly unphysical values (memory corruption, etc.).
         ! Negative abundances during solver steps are normal solver noise; physical
         ! divergence is caught between solver steps in integrateODESystem.
         if (abundances(coolantIndices(N)) > 1.0e+10_dp) then
            write(*,"(A,I3,A,A,A,A,A,1PE12.4)") &
               "ERROR: Coolant #", N, " ('", TRIM(coolantNames(N)), &
               "') has unphysical abundance for parent species '", TRIM(coolantParentNames(N)), &
               "': ", abundances(coolantIndices(N))
            write(*,"(A,I4,A,I4)") &
               "  coolantIndices(N)=", coolantIndices(N), " (max allowed=", SIZE(abundances), ")"
            write(*,*) "  This likely indicates a configuration error (wrong parent_species) or memory corruption."
            coolant_error_flag = COOLANT_CONFIG_ERROR
            write(coolant_error_message, "(A,A,A)") "Unphysical abundance for coolant ", &
               TRIM(coolantNames(N)), ". Configuration error or memory corruption."
            return
         end if

         ! Info: Warn about very small densities (but don't stop - they'll be skipped by threshold check)
         if (coolants(N)%DENSITY < 1.0e-40_dp .AND. abundances(coolantIndices(N)) > 0.0_dp) then
            ! This is normal - species just hasn't formed yet or has negligible abundance
            ! The coolant will be skipped in MANAGE_COOLANT_POPULATIONS
         end if
      end do
   end subroutine UPDATE_COOLANT_ABUNDANCES

   !=======================================================================
   !
   !  Calculate the level populations at LTE for the given species.
   !
   !-----------------------------------------------------------------------
   subroutine CALCULATE_LTE_POPULATIONS(NLEVEL,ENERGY,WEIGHT,DENSITY, &
                                      & TEMPERATURE,POPULATION)


      integer, intent(in)  :: NLEVEL
      real(dp),     intent(in)  :: ENERGY(:),WEIGHT(:)
      real(dp),     intent(in)  :: DENSITY,TEMPERATURE
      real(dp),     intent(out) :: POPULATION(:)

      integer :: ILEVEL
      real(dp) :: PARTITION_FUNCTION

   !  Initialize the level populations
      POPULATION=0.0_dp

      PARTITION_FUNCTION=0.0_dp
      do ILEVEL=1,NLEVEL
         POPULATION(ILEVEL)=WEIGHT(ILEVEL)*EXP(-ENERGY(ILEVEL)/(K_BOLTZ*TEMPERATURE))
         if (isnan(population(ilevel))) write(*,*) "LTE",temperature
         PARTITION_FUNCTION=PARTITION_FUNCTION+POPULATION(ILEVEL)
         if (isnan(PARTITION_FUNCTION)) write(*,*) density,PARTITION_FUNCTION
      end do
      POPULATION=POPULATION*DENSITY/PARTITION_FUNCTION
      !WRITE(*,*) "LTE"
      !WRITE(*,*) population
      do ILEVEL=1,NLEVEL
         if (isnan(population(ilevel))) write(*,*) ilevel,PARTITION_FUNCTION,density,temperature
      end do
   !  Check that the sum of the level populations matches the total density within the configured tolerance
      if(ABS(DENSITY-SUM(POPULATION))/DENSITY>pop_rel_tol) then
         write(*,"('ERROR: Sum of LTE level populations differs from the total density by',F6.2,'%')") &
            & 1.0e2_dp*ABS(SUM(POPULATION)-DENSITY)/DENSITY
         coolant_error_flag = COOLANT_POP_TOL_ERROR
         write(coolant_error_message, "(A,F6.2,A)") &
            "LTE population sum deviates from density by ", &
            1.0e2_dp*ABS(SUM(POPULATION)-DENSITY)/DENSITY, "%"
         return
      end if

   end subroutine CALCULATE_LTE_POPULATIONS

   !=======================================================================
   !
   !  Calculate the line opacity of each coolant transition along each
   !  HEALPix ray to the PDR surface (or simulation boundary) for all
   !  particles. Calculations can be performed using either the Euler
   !  or Trapezium integration scheme (specified by a compiler flag).
   !
   !  UCLPDR major difference is we assume a single ray and instead of looking
   !  at particles along the ray, we assume homogenous medium at a distance
   !  of STEP_SIZE from edge
   !-----------------------------------------------------------------------
   subroutine CALCULATE_LINE_OPACITIES()
      !INTEGER    :: NRAYS=1 !hard coding 1 ray
      integer :: N,ILEVEL,JLEVEL
      real(dp) :: STEP_SIZE,FACTOR1,FACTOR2,FACTOR3

      coolant: do N=1,NCOOLANTS  ! Loop over coolants
         if (coolants(N)%CONVERGED) cycle coolant
         coolants(N)%OPACITY = 0.0_dp
         levels_i: do ILEVEL=1,coolants(N)%NLEVEL  ! Loop over levels (i)
            levels_j: do JLEVEL=1,coolants(N)%NLEVEL  ! Loop over levels (j)
               if(coolants(N)%A_COEFF(ILEVEL,JLEVEL)==0) cycle levels_j
               !Factor 1 combines constants: = A_ij*c^3/8*pi*nu_ij^3
               FACTOR1 = (coolants(N)%A_COEFF(ILEVEL,JLEVEL)*C**3)/(8*PI*coolants(N)%FREQUENCY(ILEVEL,JLEVEL)**3)


               !want to divide populations in factor3 by density of current and multiply by density of cloud
               !averaged over the column to surface then multiply by distance to cloud surface
               !that is (size*average_density/density) or column_density/density
               STEP_SIZE = CLOUD_COLUMN/CLOUD_DENSITY
               !STEP_SIZE = CLOUD_SIZE!
               FACTOR2 = 1.0_dp/coolants(N)%LINEWIDTH

               !Difference between average weight of ith level and average weight of jth level
               !divided by weight of jth level.
               ! = (n_j.g_i - n_i.g_j)/g_j = n_i.(n_j.g_i/n_i.g_j - 1)
               FACTOR3 = (coolants(N)%POPULATION(JLEVEL)*coolants(N)%WEIGHT(ILEVEL)  &
                     & -  coolants(N)%POPULATION(ILEVEL)*coolants(N)%WEIGHT(JLEVEL)) &
                     &  /coolants(N)%WEIGHT(JLEVEL)

               ! IF (TRIM(coolants(N)%NAME) == "C+") THEN
               !    write(*,*) "            ---->>>> [1]opacities = ", coolants(N)%OPACITY(ILEVEL,JLEVEL)
               !    write(*,*) "            ---->>>> Factor1 = ", FACTOR1
               !    write(*,*) "            ---->>>> Factor2 = ", FACTOR2
               !    write(*,*) "            ---->>>> Factor3 = ", FACTOR3
               !    write(*,*) "            ---->>>> Step_size = ", CLOUD_COLUMN, CLOUD_DENSITY, STEP_SIZE
               ! END IF

               coolants(N)%OPACITY(ILEVEL,JLEVEL) = coolants(N)%OPACITY(ILEVEL,JLEVEL) &
                                                       & + FACTOR1*FACTOR2*FACTOR3*STEP_SIZE
               ! IF (TRIM(coolants(N)%NAME) == "C+") THEN
                  ! write(*,*) "            ---->>>> [2]opacities = ", coolants(N)%OPACITY(ILEVEL,JLEVEL)
               ! END IF
               ! dtau_ij = A_ij*c^3/8*pi*nu_ij^3 * 1/(delta)v_D * n_i.(n_j.g_i/n_i.g_j - 1) * dr
            end do levels_j
         end do levels_i
      end do coolant
   end subroutine CALCULATE_LINE_OPACITIES

!=======================================================================
!
!  Calculate the lambda operator for each allowed coolant transition,
!  needed to solve the level populations using the Accelerated Lambda
!  Iteration (ALI) method.
!
!  UCLPDR lambda operator seems to only take into account change in opacity
!  between a point and surface. If we have homogenous medium then you get
!  zero. Do we need a different way to calculate this?
!-----------------------------------------------------------------------
subroutine CALCULATE_LAMBDA_OPERATOR()
   integer :: N,i,j
   real(dp) :: dTau_1, ALI_ij

   coolant: do N=1,NCOOLANTS  ! Loop over coolants
      if(coolants(N)%CONVERGED) cycle coolant

      levels_i: do i=1,coolants(N)%nLevel
         levels_j: do j=1,coolants(N)%nLevel
            if (coolants(N)%A_COEFF(i,j) == 0.0_dp) cycle levels_j

            !We want difference in opacity between current "particle" and next
            !But I'm dealing with 1 particle so let's just set full opacity as difference between here and surface
            dTau_1=coolants(N)%OPACITY(I,J)

            if (dTau_1 /= 0.0_dp) then
               ALI_ij=2.0_dp*(1.0_dp-exp(-dTau_1))/dTau_1
               ALI_ij=1.0_dp/(1.0_dp+(ALI_ij+2.0_dp)/(dTau_1)) - 1.0_dp
            else
               ALI_ij=0.0
            end if
            coolants(N)%lambda(i,j)=ALI_ij+1.0_dp
         end do levels_j
      end do levels_i
   end do coolant
end subroutine CALCULATE_LAMBDA_OPERATOR



!=======================================================================
!
!  Calculate the level populations of a given coolant by constructing
!  the matrix of transition rates and solving to find the populations
!  assuming statistical equilibrium. The resulting set of statistical
!  equilibrium equations take the form:
!
!     n_i.sum_j R_ij = sum_j n_j.R_ji
!  or:
!     n_i.sum_j R_ij - sum_j n_j.R_ji = 0
!
!  where n_i is the population density (cm^-3) of level i and R_ij is
!  the transition rate (s^-1) from level i to level j. By rearranging
!  these equations, they can then be put into the matrix form: A.n=0,
!  where A is a coefficient matrix and n is a vector containing the N
!  population densities of all the levels. The right-hand side of the
!  matrix equation is then a vector of length N, composed of zeroes.
!
!  The elements of the coefficient matrix A are specified as follows:
!
!     A_ij = -R_ji    (j==i)
!     A_ii = sum_j R_ij (j!=i)
!
!  Since this set of equilibrium equations is not independent, one of
!  the equations has to be replaced by the conservation equation:
!
!     sum_j n_j = n_tot
!
!  where n_tot is the density of the coolant species in all levels.
!
!  Therefore, the last row of the coefficient matrix is replaced with
!  this summation over all levels, and the last right-hand-side value
!  is set to the total density of the coolant species (cm^-3).
!
!  This system of linear equations is then solved using Gauss-Jordan
!  elimination to determine the population densities of all N levels.
!
!-----------------------------------------------------------------------
subroutine CALCULATE_LEVEL_POPULATIONS(COOLANT,gasTemperature,gasDensity,abundances,dustTemp)
   real(dp), intent(in) :: gasDensity,gasTemperature,dustTemp
   real(dp), intent(in) :: abundances(:)
   type(COOLANT_TYPE),  intent(inout)    :: COOLANT
   integer :: I,J,N,ierr
   real(dp)     :: SUM

!  Reset the persistent workspace arrays (no allocation needed)
   COOLANT%A=0.0_dp
   COOLANT%B=0.0_dp
   COOLANT%R=0.0_dp
   ! write(*,*) "temp",gasTemperature
   ! write(*,*) "pop",coolant%population

!  Debug: print dustTemp before entering CONSTRUCT_TRANSITION_MATRIX
!   WRITE(*,'(A,F12.4,A,F12.4,A,I3)') "DEBUG dustTemp=", dustTemp, &
!      " gasTemp=", gasTemperature, " NLEVEL=", COOLANT%NLEVEL
   if (ISNAN(dustTemp)) write(*,*) "WARNING: dustTemp is NaN!"
   if (dustTemp <= 0.0_dp) write(*,*) "WARNING: dustTemp <= 0!"

!  Construct the matrix of transition rates R_ij (s^-1)
   call CONSTRUCT_TRANSITION_MATRIX(COOLANT,COOLANT%R,gasTemperature,gasDensity&
      &,abundances,dustTemp)


!  Fill the coefficient matrix A and the right-hand-side vector b
   levels_i: do I=1,COOLANT%NLEVEL
      SUM=0.0_dp
      levels_j: do J=1,COOLANT%NLEVEL
         if(J==I) cycle levels_j
         if(ISNAN(COOLANT%R(J,I))) then
            write(*,*) "ERROR: NaN in transition matrix for coolant ", TRIM(coolant%name)
            do N=1,COOLANT%NLEVEL
               write(*,*) COOLANT%R(N,:)
            end do
            coolant_error_flag = COOLANT_SOLVER_ERROR
            write(coolant_error_message, "(A,A)") "NaN in transition matrix for coolant ", &
               TRIM(coolant%name)
            return
         end if
         ! IF (COOLANT%R(J,I) .eq. 0.0) write(*,*) coolant%name,I,J
         COOLANT%A(I,J) = -COOLANT%R(J,I)
         SUM=SUM+COOLANT%R(I,J)
      end do levels_j
      COOLANT%A(I,I)=SUM
      if (COOLANT%A(I,I) == 0.0_dp) write(*,*) coolant%name,I
   end do levels_i
   COOLANT%B=0.0_dp

!  Replace the last equilibrium equation in the transition matrix with
!  the conservation equation (i.e. the sum of the population densities
!  over all levels), and replace the last entry in the right-hand-side
!  vector with the total density of the coolant species.
   COOLANT%A(COOLANT%NLEVEL,:)=1  ! Sum over all levels, sum_j n_j
   COOLANT%B(COOLANT%NLEVEL)=coolant%density

!  Call the Gauss-Jordan solver (the solution is returned in vector b)
   call GAUSS_JORDAN(COOLANT%NLEVEL,COOLANT%A,COOLANT%B,COOLANT%IPIV,COOLANT%INDEX_ROW,COOLANT%INDEX_COL,ierr)
   if (ierr /= 0) then
      coolant_error_flag = COOLANT_SOLVER_ERROR
      write(coolant_error_message, "(A,A)") "Singular matrix in Gauss-Jordan for coolant ", &
         TRIM(COOLANT%NAME)
      return
   end if

!  Replace negative or NaN values caused by numerical noise around zero
   do I=1,COOLANT%NLEVEL
      if(.NOT.COOLANT%B(I)>=0) COOLANT%B(I)=0.0_dp
      if (ISNAN(COOLANT%B(I))) then
         write(*,*) "NaN in population: ", I, COOLANT%B(I),COOLANT%NLEVEL
         write(*,*) COOLANT%A
      end if
   end do

!  Store the previously-calculated population densities for comparison
   coolant%PREVIOUS_POPULATION=coolant%POPULATION

!  Store the new population densities
   coolant%POPULATION=COOLANT%B

!  Note: Arrays are persistent and reused, no deallocation needed

end subroutine CALCULATE_LEVEL_POPULATIONS
!=======================================================================

!=======================================================================
!
!  Standard linear equation solver using Gauss-Jordon elimination taken
!  directly from Numerical Recipes (Chapter B2).
!
!  A is an NxN input coefficient matrix, B is an input vector of size N
!  containing the right-hand-side values. On output, A is replaced by its
!  matrix inverse and B is replaced by the corresponding set of solution
!  values.
!
!-----------------------------------------------------------------------
subroutine GAUSS_JORDAN(N,A,B,IPIV,INDEX_ROW,INDEX_COL,IERR)
   use numerics, only: SWAP

   integer, intent(in)    :: N
   real(dp),     intent(inout) :: A(:,:)
   real(dp),     intent(inout) :: B(:)
   integer, intent(inout) :: IPIV(:),INDEX_ROW(:),INDEX_COL(:)
   integer, intent(out)   :: IERR

   integer :: I,J,K,L,IROW,ICOL
   real(dp) :: MAX,DUMMY,PIVINV

!  Reset the persistent workspace arrays (no allocation needed)
   IERR=0
   ICOL=0
   IROW=0
   IPIV=0

   do I=1,N  ! Main loop over columns to be reduced
      MAX=0.0_dp
      do J=1,N
         if(IPIV(J)/=1) then
            do K=1,N
               if(IPIV(K)==0) then
                  if(ISNAN(A(J,K))) then
                     write(*,*)
                     write(*,*) "ERROR! NaN found in coefficient matrix A of Gauss-Jordan routine"
                     write(*,"(' A(',I3,',',I3,') =',F4.1)") J,K,A(J,K)
                     write(*,*) "A ="
                     do L=1,N
                        write(*,"(100(0PD9.1))") A(L,:)
                     end do
                     write(*,*)
                     IERR = -1
                     return
                  end if
                  if(ABS(A(J,K))>=MAX) then
                     MAX=ABS(A(J,K))
                     IROW=J
                     ICOL=K
                  end if
               else if(IPIV(K)>1) then
                  write(*,*)
                  write(*,*) "ERROR! Singular matrix found in Gauss-Jordan routine (#1)"
                  write(*,*)
                  IERR = -1
                  return
               end if
            end do
         end if
      end do
      IPIV(ICOL)=IPIV(ICOL)+1
      if(IROW/=ICOL) then
         call SWAP(A(IROW,:),A(ICOL,:))
         call SWAP(B(IROW),B(ICOL))
      end if
      INDEX_ROW(I)=IROW
      INDEX_COL(I)=ICOL
      if(A(ICOL,ICOL)==0.0_dp) then
         write(*,*)
         write(*,*) "ERROR! Singular matrix found in Gauss-Jordan routine (#2)"
         write(*,*)
         IERR = -1
         return
      end if
      PIVINV=1.0_dp/A(ICOL,ICOL)
      A(ICOL,ICOL)=1.0_dp
      A(ICOL,:)=A(ICOL,:)*PIVINV
      B(ICOL)=B(ICOL)*PIVINV
      do L=1,N
         if(L/=ICOL) then
            DUMMY=A(L,ICOL)
            A(L,ICOL)=0.0_dp
            A(L,:)=A(L,:)-A(ICOL,:)*DUMMY
            B(L)=B(L)-B(ICOL)*DUMMY
         end if
      end do
   end do  ! End of main loop over columns

!  Unscramble the solution by interchanging pairs of columns
!  in the reverse order to which the permutation was built up
   do L=N,1,-1
      call SWAP(A(:,INDEX_ROW(L)),A(:,INDEX_COL(L)))
   end do

!  Note: Arrays are persistent and reused, no deallocation needed

end subroutine GAUSS_JORDAN

!=======================================================================
!=======================================================================
!
!  Construct the matrix of transition rates for a given coolant species
!  using the escape probability formalism to determine the local field.
!
!-----------------------------------------------------------------------
subroutine CONSTRUCT_TRANSITION_MATRIX(COOLANT,TRANSITION_MATRIX,gasTemperature,gasDensity,abundances,dustTemp)
   type(COOLANT_TYPE), intent(inout)    :: COOLANT
   real(dp), intent(in)    :: gasTemperature,gasDensity,dustTemp
   real(dp), intent(in) :: abundances(:)
   real(dp), intent(out)   :: TRANSITION_MATRIX(:,:)

   integer :: I,J,K,cache_hit
   real(dp) :: dustEmissivity
   real(dp) :: S_ij,B_ij,B_ij_CMB,B_ij_DUST,BETA_ij
   real(dp) :: S_ij_PREVIOUS,LAMBDA_ij,rhoGrain,dustDensity
   real(dp) :: n_H2,n_elec,n_H,n_He,n_Hplus
   logical :: found_in_cache

   ! Debug: Check for NaN inputs
   if (ISNAN(dustTemp)) then
      write(*,*) "ERROR: dustTemp is NaN in CONSTRUCT_TRANSITION_MATRIX for ", TRIM(COOLANT%NAME)
      return
   end if
   if (ISNAN(gasTemperature)) then
      write(*,*) "ERROR: gasTemperature is NaN in CONSTRUCT_TRANSITION_MATRIX for ", TRIM(COOLANT%NAME)
      return
   end if

!  Reset the persistent workspace arrays (no allocation needed)
   COOLANT%RADIATION_FIELD=0.0_dp
   ! TRANSITION_MATRIX=0.0_dp
!  Initialize the local emissivities
   COOLANT%EMISSIVITY=0.0_dp

   levels_i: do I=1,COOLANT%NLEVEL
      levels_j: do J=1,COOLANT%NLEVEL
         if(J>=I) exit levels_j

!        Calculate the background radiation field including contributions
!        from CMB blackbody emission and dust modified blackbody emission
         if(COOLANT%FREQUENCY(I,J)<1.0e15_dp) then
            B_ij_CMB=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2)/ &
                & (EXP(HP*COOLANT%FREQUENCY(I,J)/(K_BOLTZ*T_CMB))-1.0_dp)
!#ifdef DUST
           rhoGrain=2.0_dp  ! Grain mass density (g cm^-3)
           dustDensity=1.0e-12_dp*gasDensity
           dustEmissivity=(rhoGrain*dustDensity)*(0.01_dp*(1.3_dp*COOLANT%FREQUENCY(I,J)/3.0e11_dp))

           ! Protect against invalid dustTemp
           if (dustTemp <= 0.0_dp .OR. ISNAN(dustTemp)) then
              B_ij_DUST = 0.0_dp
           else
              B_ij_DUST=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2)/ &
                  & (EXP(HP*COOLANT%FREQUENCY(I,J)/(K_BOLTZ*dustTemp))-1.0_dp)*dustEmissivity
              if (ISNAN(B_ij_DUST)) B_ij_DUST = 0.0_dp
           end if
!#else
 !           B_ij_DUST=0.T
!#endif
            B_ij=B_ij_CMB+B_ij_DUST
         else
            B_ij=0.0_dp
         end if

!        If the population n_i is zero then the source function is zero and the
!        mean integrated radiation field is just the background radiation field
         if(COOLANT%POPULATION(I)< 1.0e-20_dp) then
            S_ij=0.0_dp
            BETA_ij=1.0_dp

!        If the difference between n_i.g_j and n_j.g_i is vanishingly small, then
!        calculate the source function directly and set the escape probability to 1
         else if(ABS(COOLANT%POPULATION(I)*COOLANT%WEIGHT(J)-COOLANT%POPULATION(J)*COOLANT%WEIGHT(I))==0) then
            S_ij=HP*COOLANT%FREQUENCY(I,J)*COOLANT%POPULATION(I)*COOLANT%A_COEFF(I,J)/(4.0_dp*PI)
            BETA_ij=1.0_dp
           ! write(*,*) "vanishingly",HP,COOLANT%FREQUENCY(I,J),COOLANT%POPULATION(I),COOLANT%A_COEFF(I,J),(4.0*PI)
         else
!           Calculate the source function
            S_ij=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2) &
              & /((COOLANT%POPULATION(J)*COOLANT%WEIGHT(I)) &
              &  /(COOLANT%POPULATION(I)*COOLANT%WEIGHT(J))-1.0_dp)

!           Calculate the escape probability
            BETA_ij=ESCAPE_PROBABILITY(COOLANT%OPACITY(I,J))
            !write(*,*) "regular",HP,COOLANT%FREQUENCY(I,J),COOLANT%POPULATION(J),COOLANT%WEIGHT(I),COOLANT%POPULATION(I),COOLANT%WEIGHT(J)
         end if

!        Calculate the local emissivity (erg cm^-3 s^-1) for the transition line
         if(COOLANT%POPULATION(I)>0) then
            COOLANT%EMISSIVITY(I,J)=COOLANT%POPULATION(I)*COOLANT%A_COEFF(I,J) &
                                             & *HP*COOLANT%FREQUENCY(I,J)*BETA_ij*(S_ij-B_ij)/S_ij
         end if
!        Calculate the mean integrated radiation field <J_ij>
         COOLANT%RADIATION_FIELD(I,J)=(1.0_dp-BETA_ij)*S_ij+BETA_ij*B_ij

         ! Protect against NaN in radiation field
         if (ISNAN(COOLANT%RADIATION_FIELD(I,J))) then
            COOLANT%RADIATION_FIELD(I,J) = 0.0_dp
         end if

!        Lambda operator keeps breaking so lets not use ALI for now
! !        Calculate the source function in the same manner for
! !        the populations determined on the previous iteration
!          IF(COOLANT%PREVIOUS_POPULATION(I).lt.1.e-20_dp) THEN
!             S_ij_PREVIOUS=0.0_dp
!          ELSE IF(ABS(COOLANT%PREVIOUS_POPULATION(I)*COOLANT%WEIGHT(J)-COOLANT%PREVIOUS_POPULATION(J)*COOLANT%WEIGHT(I)).EQ.0) THEN
!             S_ij_PREVIOUS=HP*COOLANT%FREQUENCY(I,J)*COOLANT%PREVIOUS_POPULATION(I)*COOLANT%A_COEFF(I,J)/(4.0*PI)
!          ELSE
!             S_ij_PREVIOUS=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2) &
!               & /((COOLANT%PREVIOUS_POPULATION(J)*COOLANT%WEIGHT(I)) &
!               &  /(COOLANT%PREVIOUS_POPULATION(I)*COOLANT%WEIGHT(J))-1.0_dp)
!          END IF

! !        Use the Accelerated Lambda Iteration method to speed up convergence
! !        by amplifying the incremental difference of the new source function
!          LAMBDA_ij=COOLANT%LAMBDA(I,J)
!          COOLANT%RADIATION_FIELD(I,J)=(1.0_dp-BETA_ij)*(LAMBDA_ij*S_ij+(1.0_dp-LAMBDA_ij)*S_ij_PREVIOUS)+BETA_ij*B_ij

         COOLANT%RADIATION_FIELD(J,I)=COOLANT%RADIATION_FIELD(I,J)
      end do levels_j  ! End of loop over levels (j)
   end do levels_i  ! End of loop over levels (i)

!  Calculate or retrieve cached collisional rates
!  Calculate number densities of collision partners for cache key
   n_H2 = gasDensity * abundances(nH2)
   n_elec = gasDensity * abundances(nelec)
   n_H = gasDensity * abundances(nH)
   n_He = gasDensity * abundances(nHe)
   n_Hplus = gasDensity * abundances(nHx)

   found_in_cache = .false.
   cache_hit = -1

!  Search cache for matching parameters (within tolerance)
   cache_index: do K=1,10
      if (COOLANT%CACHED_TEMPERATURE(K) < 0.0_dp) cycle cache_index  ! Skip empty slots
      if (IS_WITHIN_TOLERANCE(COOLANT%CACHED_TEMPERATURE(K), gasTemperature, COOLANT%CACHE_TOLERANCE) .AND. &
          IS_WITHIN_TOLERANCE(COOLANT%CACHED_N_H2(K), n_H2, COOLANT%CACHE_TOLERANCE) .AND. &
          IS_WITHIN_TOLERANCE(COOLANT%CACHED_N_ELEC(K), n_elec, COOLANT%CACHE_TOLERANCE) .AND. &
          IS_WITHIN_TOLERANCE(COOLANT%CACHED_N_H(K), n_H, COOLANT%CACHE_TOLERANCE) .AND. &
          IS_WITHIN_TOLERANCE(COOLANT%CACHED_N_HE(K), n_He, COOLANT%CACHE_TOLERANCE) .AND. &
          IS_WITHIN_TOLERANCE(COOLANT%CACHED_N_HPLUS(K), n_Hplus, COOLANT%CACHE_TOLERANCE)) then
         found_in_cache = .true.
         cache_hit = K
         exit cache_index
      end if
   end do cache_index

   if (found_in_cache) then
      ! Use cached collisional rates
      COOLANT%COLLISIONAL_RATE = COOLANT%CACHED_COLLISIONAL_RATE(cache_hit,:,:)
   else
      ! Calculate new collisional rates
      COOLANT%COLLISIONAL_RATE=0.0_dp
      call CALCULATE_COLLISIONAL_RATES(COOLANT,gasDensity,gasTemperature,abundances,COOLANT%COLLISIONAL_RATE)

      ! Store in cache (round-robin replacement)
      COOLANT%CACHED_TEMPERATURE(COOLANT%CACHE_INDEX) = gasTemperature
      COOLANT%CACHED_N_H2(COOLANT%CACHE_INDEX) = n_H2
      COOLANT%CACHED_N_ELEC(COOLANT%CACHE_INDEX) = n_elec
      COOLANT%CACHED_N_H(COOLANT%CACHE_INDEX) = n_H
      COOLANT%CACHED_N_HE(COOLANT%CACHE_INDEX) = n_He
      COOLANT%CACHED_N_HPLUS(COOLANT%CACHE_INDEX) = n_Hplus
      COOLANT%CACHED_COLLISIONAL_RATE(COOLANT%CACHE_INDEX,:,:) = COOLANT%COLLISIONAL_RATE

      ! Advance cache index (circular)
      COOLANT%CACHE_INDEX = COOLANT%CACHE_INDEX + 1
      if (COOLANT%CACHE_INDEX > 10) COOLANT%CACHE_INDEX = 1
   end if
   ! IF (TRIM(coolant%name) == "C+") THEN
   !    write(*,*) "      **** opaticy: ", coolant%opacity
   !    write(*,*) "      **** coefficients", BETA_ij, S_ij, B_ij
   !    write(*,*) "      **** coolant: ", TRIM(coolant%name), "' radiation_field=", COOLANT%RADIATION_FIELD
   ! END IF
!  Construct the transition matrix: R_ij = A_ij + B_ij.<J> + C_ij
   do I=1,COOLANT%NLEVEL
      do J=1,COOLANT%NLEVEL
         if (isnan(coolant%A_COEFF(I,J))) write(*,*) "detected NaN A_COEFF!"
         if (isnan(coolant%b_COEFF(I,J))) write(*,*) "detected NaN b_COEFF!"
         if (isnan(COOLANT%RADIATION_FIELD(I,J))) write(*,*) "detected NaN RADFIELD!"
         if (isnan(COOLANT%COLLISIONAL_RATE(I,J))) write(*,*) "Detection NaN collisional rates!"
         TRANSITION_MATRIX(I,J)=COOLANT%A_COEFF(I,J) &
                            & + COOLANT%B_COEFF(I,J)*COOLANT%RADIATION_FIELD(I,J) &
                            & + COOLANT%COLLISIONAL_RATE(I,J)
      end do
   end do

!  Note: Arrays are persistent and reused, no deallocation needed

end subroutine CONSTRUCT_TRANSITION_MATRIX
!=======================================================================


!=======================================================================
!
!  Calculate the total collisional rates (s^-1) for a given coolant at
!  the specified temperature by summing the individual rates from each
!  of the available collision partners by linear interpolation between
!  the available rate coefficients at specific temperature values.
!
!-----------------------------------------------------------------------
subroutine CALCULATE_COLLISIONAL_RATES(COOLANT,DENSITY,TEMPERATURE, &
                                     & ABUNDANCE,COLLISIONAL_RATE)
  ! USE FUNCTIONS_MODULE


   type(COOLANT_TYPE), intent(in)  :: COOLANT
   real(dp),      intent(in)  :: DENSITY,TEMPERATURE
   real(dp),      intent(in)  :: ABUNDANCE(:)
   real(dp),      intent(out) :: COLLISIONAL_RATE(:,:)

   integer :: I,J,K,L,M,KLO,KHI,PARTNER_ID,NTEMP_P
   real(dp) :: PARA_FRACTION,ORTHO_FRACTION
   real(dp) :: STEP,C_COEFF,ABUND_FACTOR

!  Initialize the collisional rates
   COLLISIONAL_RATE=0.0_dp

!  Calculate the H2 ortho/para ratio at equilibrium for the specified
!  temperature and the resulting fractions of H2 in para & ortho form
   PARA_FRACTION=1.0_dp/(1.0_dp+GET_ORTHO_PARA_RATIO(TEMPERATURE))
   ORTHO_FRACTION=1.0_dp-PARA_FRACTION

   do L=1,COOLANT%NPARTNER  ! Loop over collision partners (sparse)
      PARTNER_ID = COOLANT%partners(L)%partner_id
      NTEMP_P    = COOLANT%partners(L)%ntemp

!     Determine the abundance factor for this collision partner
      if(PARTNER_ID==1) then
         ABUND_FACTOR = DENSITY * ABUNDANCE(nH2)
      else if(PARTNER_ID==2) then
         ABUND_FACTOR = DENSITY * ABUNDANCE(nH2) * PARA_FRACTION
      else if(PARTNER_ID==3) then
         ABUND_FACTOR = DENSITY * ABUNDANCE(nH2) * ORTHO_FRACTION
      else if(PARTNER_ID==4) then
         ABUND_FACTOR = DENSITY * ABUNDANCE(nelec)
      else if(PARTNER_ID==5) then
         ABUND_FACTOR = DENSITY * ABUNDANCE(nH)
      else if(PARTNER_ID==6) then
         ABUND_FACTOR = DENSITY * ABUNDANCE(nHe)
      else  ! PARTNER_ID.EQ.7: protons
         ABUND_FACTOR = DENSITY * ABUNDANCE(nHx)
      end if

!     Determine the two nearest temperature values in the partner's temperature grid
      KLO=0
      KHI=0
      temperatures: do K=1,NTEMP_P
         if(COOLANT%partners(L)%temperature(K)>TEMPERATURE) then
            KLO=K-1
            KHI=K
            exit temperatures
         end if
      end do temperatures

!     Clamp to valid range
      if(KHI==0) then
         KLO=NTEMP_P
         KHI=NTEMP_P
      else if(KHI==1) then
         KLO=1
         KHI=1
      end if

!     Linear interpolation step fraction
      if(KLO==KHI) then
         STEP=0.0_dp
      else
         STEP=(TEMPERATURE-COOLANT%partners(L)%temperature(KLO)) &
           & /(COOLANT%partners(L)%temperature(KHI)-COOLANT%partners(L)%temperature(KLO))
      end if

!     Accumulate rates over stored (I,J) pairs only (sparse iteration)
      do M=1,COOLANT%partners(L)%ncoll
         I = COOLANT%partners(L)%i_idx(M)
         J = COOLANT%partners(L)%j_idx(M)
         C_COEFF = COOLANT%partners(L)%c_coeff(M,KLO) &
                 + (COOLANT%partners(L)%c_coeff(M,KHI) - COOLANT%partners(L)%c_coeff(M,KLO)) * STEP
         COLLISIONAL_RATE(I,J) = COLLISIONAL_RATE(I,J) + C_COEFF * ABUND_FACTOR
      end do

   end do  ! End of loop over collision partners

end subroutine CALCULATE_COLLISIONAL_RATES

!=======================================================================
!
!  Check for convergence in the population densities of all levels of
!  each coolant of each particle. Set the relevant convergence flags.
!  Calculate the percentage of particles that have converged for each
!  coolant.
!
!-----------------------------------------------------------------------
function CHECK_CONVERGENCE() result(all_coolants_converged)
   logical :: all_coolants_converged

   integer :: I,N
   real(dp) :: RELATIVE_CHANGE
   real(dp), parameter :: POPULATION_LIMIT=1.0e-14_dp
   real(dp), parameter :: POPULATION_CONVERGENCE_CRITERION=1.0e-2_dp
   logical :: convergence(NCOOLANTS)

   coolant: do N=1,NCOOLANTS  ! Loop over coolants
      coolants(N)%CONVERGED=.true.
      levels: do I=1,coolants(N)%NLEVEL  ! Loop over levels

!        Skip this level if its population density is below the cut-off
         if(coolants(N)%POPULATION(I)<POPULATION_LIMIT*coolants(N)%DENSITY) cycle levels

!        Skip this level if its population density has not changed
         if(coolants(N)%POPULATION(I)==coolants(N)%PREVIOUS_POPULATION(I)) cycle levels

!        Calculate the relative change in population density between this iteration and the previous
         RELATIVE_CHANGE=ABS(coolants(N)%POPULATION(I)-coolants(N)%PREVIOUS_POPULATION(I)) &
                       & *2/(coolants(N)%POPULATION(I)+coolants(N)%PREVIOUS_POPULATION(I))

!        If the relative change is greater than the criterion for convergence, set the flag to false
         if(RELATIVE_CHANGE>POPULATION_CONVERGENCE_CRITERION) then
            coolants(N)%CONVERGED=.false.
            exit levels
         end if
      end do levels
      convergence(N)=coolants(N)%converged
   end do coolant
   all_coolants_converged=ALL(convergence)
end function CHECK_CONVERGENCE


!=======================================================================
!
!  Warm restart: rescale existing coolant populations to new densities.
!  This preserves the population distribution while adjusting to new
!  total densities. The old populations (which sum to old density) are
!  rescaled proportionally so they sum to the new density.
!
!  This is useful when physical conditions change slightly - the relative
!  population distribution should be similar, just scaled to the new total.
!
!-----------------------------------------------------------------------
subroutine WARMSTART_COOLANT_POPULATIONS()
   integer :: N
   real(dp) :: old_total, scale_factor

   do N=1,NCOOLANTS
      ! Calculate the current total population
      old_total = SUM(coolants(N)%POPULATION)

      ! Avoid division by zero - if no populations exist, skip rescaling
      if (old_total > 0.0_dp) then
         ! Calculate scale factor to match new density
         scale_factor = coolants(N)%DENSITY / old_total

         ! Rescale all level populations
         coolants(N)%POPULATION = coolants(N)%POPULATION * scale_factor
         coolants(N)%PREVIOUS_POPULATION = coolants(N)%PREVIOUS_POPULATION * scale_factor

         ! Mark as not converged to trigger refinement
         coolants(N)%CONVERGED = .false.
      end if
   end do

end subroutine WARMSTART_COOLANT_POPULATIONS


!=======================================================================
!
!  Manage coolant populations with automatic initialization and warm restart.
!  This should be called BEFORE the SE iteration loop to ensure populations
!  are properly initialized or restarted based on the coolant_restart_mode.
!
!  Modes:
!    0 = WARM (default): Initialize to LTE on first call, then rescale when density changes
!    1 = FORCE_LTE: Always reset to LTE before each SE iteration (original behavior)
!    2 = FORCE_GROUND: Always reset to ground state before each SE iteration
!
!-----------------------------------------------------------------------
subroutine MANAGE_COOLANT_POPULATIONS(gasTemperature)
   real(dp), intent(in) :: gasTemperature
   integer :: N
   real(dp) :: old_total, scale_factor

   ! Mode 1: FORCE_LTE - always reset to LTE (original behavior)
   if (coolant_restart_mode == 1) then
      coolant: do N=1,NCOOLANTS
         ! Skip coolants with negligible density
         if (coolants(N)%DENSITY < 1.0e-40_dp) then
            coolants(N)%POPULATION = 0.0_dp
            coolants(N)%PREVIOUS_POPULATION = 0.0_dp
            coolants(N)%EMISSIVITY = 0.0_dp
            coolants(N)%CONVERGED = .true.
            cycle coolant
         end if

         call CALCULATE_LTE_POPULATIONS(coolants(N)%NLEVEL, &
                                       coolants(N)%ENERGY, &
                                       coolants(N)%WEIGHT, &
                                       coolants(N)%DENSITY, &
                                       gasTemperature, &
                                       coolants(N)%POPULATION)
         coolants(N)%PREVIOUS_POPULATION = coolants(N)%POPULATION
         coolants(N)%CONVERGED = .false.
      end do coolant
      return
   end if

   ! Mode 2: FORCE_GROUND - always reset to ground state
   if (coolant_restart_mode == 2) then
      coolant_2: do N=1,NCOOLANTS
         ! Skip coolants with negligible density
         if (coolants(N)%DENSITY < 1.0e-40_dp) then
            coolants(N)%POPULATION = 0.0_dp
            coolants(N)%PREVIOUS_POPULATION = 0.0_dp
            coolants(N)%EMISSIVITY = 0.0_dp
            coolants(N)%CONVERGED = .true.
            cycle coolant_2
         end if

         coolants(N)%POPULATION = 0.0_dp
         coolants(N)%POPULATION(1) = coolants(N)%DENSITY
         coolants(N)%PREVIOUS_POPULATION = coolants(N)%POPULATION
         coolants(N)%CONVERGED = .false.
      end do coolant_2
      return
   end if

   ! Mode 0: WARM (default) - smart initialization and warm restart
   if (.NOT. coolant_populations_initialized) then
      ! First call: Initialize populations to LTE
      coolant_3: do N=1,NCOOLANTS
         ! Skip coolants with negligible density (< 1e-40 cm^-3)
         ! This avoids wasting computation on absent species
         if (coolants(N)%DENSITY < 1.0e-40_dp) then
            coolants(N)%POPULATION = 0.0_dp
            coolants(N)%PREVIOUS_POPULATION = 0.0_dp
            coolants(N)%EMISSIVITY = 0.0_dp
            coolants(N)%CONVERGED = .true.  ! Mark as converged to skip SE solver
            cycle coolant_3
         end if

         call CALCULATE_LTE_POPULATIONS(coolants(N)%NLEVEL, &
                                       coolants(N)%ENERGY, &
                                       coolants(N)%WEIGHT, &
                                       coolants(N)%DENSITY, &
                                       gasTemperature, &
                                       coolants(N)%POPULATION)
         coolants(N)%PREVIOUS_POPULATION = coolants(N)%POPULATION
         coolants(N)%CONVERGED = .false.
      end do coolant_3
      coolant_populations_initialized = .true.
      last_levpop_temperature = gasTemperature
   else
      ! Subsequent calls: Warm restart - rescale if density changed
      coolant_4: do N=1,NCOOLANTS
         ! Skip coolants with negligible density (< 1e-40 cm^-3)
         if (coolants(N)%DENSITY < 1.0e-40_dp) then
            coolants(N)%POPULATION = 0.0_dp
            coolants(N)%PREVIOUS_POPULATION = 0.0_dp
            coolants(N)%EMISSIVITY = 0.0_dp
            coolants(N)%CONVERGED = .true.
            cycle coolant_4
         end if

         old_total = SUM(coolants(N)%POPULATION)

         ! Check if density changed significantly (> 0.1% relative change)
         if (old_total > 0.0_dp) then
            if (ABS(coolants(N)%DENSITY - old_total) / old_total > 1.0e-3_dp) then
               ! Density changed: rescale populations proportionally
               scale_factor = coolants(N)%DENSITY / old_total
               coolants(N)%POPULATION = coolants(N)%POPULATION * scale_factor
               coolants(N)%PREVIOUS_POPULATION = coolants(N)%PREVIOUS_POPULATION * scale_factor
               coolants(N)%CONVERGED = .false.
            end if
         else
            ! No existing populations: initialize to ground state
            coolants(N)%POPULATION = 0.0_dp
            coolants(N)%POPULATION(1) = coolants(N)%DENSITY
            coolants(N)%PREVIOUS_POPULATION = coolants(N)%POPULATION
            coolants(N)%CONVERGED = .false.
         end if
      end do coolant_4

      ! Also reset convergence if temperature changed significantly or a recompute was forced
      if (coolant_levpop_force_recompute .OR. &
          ABS(gasTemperature - last_levpop_temperature) / MAX(last_levpop_temperature, 1.0_dp) &
              .GT. coolant_temp_recompute_threshold) then
         do N = 1, NCOOLANTS
            if (coolants(N)%DENSITY >= 1.0e-40_dp) coolants(N)%CONVERGED = .false.
         end do
         last_levpop_temperature = gasTemperature
         coolant_levpop_force_recompute = .false.
      end if
   end if

end subroutine MANAGE_COOLANT_POPULATIONS


!=======================================================================
!
!  Get the current coolant restart mode
!
!-----------------------------------------------------------------------
function GET_COOLANT_RESTART_MODE() result(current_coolant_restart_mode)
    integer :: current_coolant_restart_mode
   current_coolant_restart_mode = coolant_restart_mode
end function GET_COOLANT_RESTART_MODE


!=======================================================================
!
!  Set the coolant restart mode
!
!  mode values:
!    0 = WARM (default): Initialize to LTE on first call, rescale on density change
!    1 = FORCE_LTE: Always reset to LTE before SE iteration
!    2 = FORCE_GROUND: Always reset to ground state before SE iteration
!
!-----------------------------------------------------------------------
subroutine SET_COOLANT_RESTART_MODE(mode)
   integer, intent(in) :: mode
   if (mode < 0 .OR. mode > 2) then
      write(*,*) "WARNING: Invalid coolant_restart_mode ", mode
      write(*,*) "Valid modes: 0 (WARM), 1 (FORCE_LTE), 2 (FORCE_GROUND)"
      write(*,*) "Keeping current mode: ", coolant_restart_mode
      return
   end if
   coolant_restart_mode = mode
   ! Reset initialization flag when mode changes to force reinitialization
   coolant_populations_initialized = .false.
end subroutine SET_COOLANT_RESTART_MODE


!=======================================================================
!
!  Calculate the ortho-to-para ratio of H2 at thermal equilibrium
!  for the specified temperature, making use of energy level data
!  if available, or an approximation if not.
!
!-----------------------------------------------------------------------
function GET_ORTHO_PARA_RATIO(TEMPERATURE) result(ORTHO_PARA_RATIO)
   real(dp), intent(in) :: TEMPERATURE
   real(dp) :: ORTHO_PARA_RATIO


   integer :: I,J,N,ORTHO_INDEX,PARA_INDEX
   real(dp) :: I_ORTHO,I_PARA,ORTHO_FRACTION,PARA_FRACTION
!  Check if coolant data is available for the ortho and para forms
   ORTHO_INDEX=0
   PARA_INDEX=0
   do N=1,NCOOLANTS
      if(COOLANTS(N)%NAME=="o-H2") ORTHO_INDEX=N
      if(COOLANTS(N)%NAME=="p-H2") PARA_INDEX=N
   end do

   !  Calculate the exact ortho/para ratio if molecular data is available
   if(ORTHO_INDEX/=0 .AND. PARA_INDEX/=0) then
      I_ORTHO=1.0_dp
      I_PARA=0.0_dp  ! Total nuclear spins of the two forms

      ! Calculate the ortho/para ratio of H2 using the expression
      ! from Poelman & Spaans (2005, A&A, 440, 559, equation 11)
      ORTHO_FRACTION=0.0_dp
      do I=1,COOLANTS(ORTHO_INDEX)%NLEVEL
         ORTHO_FRACTION=ORTHO_FRACTION+COOLANTS(ORTHO_INDEX)%WEIGHT(I) &
                     & *EXP(-COOLANTS(ORTHO_INDEX)%ENERGY(I)/(K_BOLTZ*TEMPERATURE))
      end do

      PARA_FRACTION=0.0_dp
      do I=1,COOLANTS(PARA_INDEX)%NLEVEL
         PARA_FRACTION=PARA_FRACTION+COOLANTS(PARA_INDEX)%WEIGHT(I) &
                    & *EXP(-COOLANTS(PARA_INDEX)%ENERGY(I)/(K_BOLTZ*TEMPERATURE))
      end do

      ORTHO_PARA_RATIO=(2*I_ORTHO+1)/(2*I_PARA+1)*(ORTHO_FRACTION/PARA_FRACTION)

   else
     ! Approximate the ortho/para ratio of H2 using the expression
     ! given by Flower & Watt (1985, MNRAS, 213, 991, equation 2)
     ORTHO_PARA_RATIO=9.0_dp*EXP(-170.5_dp/TEMPERATURE)

     ! Limit the ortho/para ratio to its statistical limit
     if(ORTHO_PARA_RATIO>3.0_dp) ORTHO_PARA_RATIO=3.0_dp

   end if

end function GET_ORTHO_PARA_RATIO

function ESCAPE_PROBABILITY(TAU) result(BETA)
   real(dp),  intent(in) :: TAU
   integer :: K
   real(dp)  :: BETA

!  Initialize the escape probability values along each ray
   BETA=0.0_dp

!     Limit the escape probability to unity for masing transitions (tau <= 0)
   if(TAU<=0) then
      BETA=1.0_dp

!     Prevent floating point overflow caused by very low opacity (tau < 1E-8)
   else if(ABS(TAU)<1.0e-8_dp) then
      BETA=1.0_dp

!     For all other cases use the standard escape probability formalism
   else
      BETA=(1.0_dp-EXP(-TAU))/TAU
   end if

!  The total escape probability must be divided by the number of rays to
!  account for the fraction of the total solid angle covered by each ray
!  (assuming that each ray covers the same fraction of the total 4*pi sr).
!  In the case of only 1 ray (i.e., semi-infinite slab geometry) the ray
!  subtends a solid angle of 2*pi sr, since the photons escape through the
!  hemisphere in the outward direction, so its escape probability should
!  be divided by two.
   !BETA=0.5*BETA
end function ESCAPE_PROBABILITY

!=======================================================================
!
!  Check if two values match within relative tolerance for cache lookup.
!  Uses symmetric relative difference: 2|a-b|/(|a|+|b|) < tolerance
!
!-----------------------------------------------------------------------
function IS_WITHIN_TOLERANCE(cached_val, current_val, tol) result(within_tolerance)
   real(dp), intent(in) :: cached_val, current_val, tol
   logical :: within_tolerance
   real(dp) :: rel_diff

!  Handle near-zero values (both must be negligible)
   if (ABS(cached_val) < MIN_ABUND .AND. ABS(current_val) < MIN_ABUND) then
      within_tolerance = .true.
      return
   end if

!  One value near-zero but not the other - no match
   if (ABS(cached_val) < MIN_ABUND .OR. ABS(current_val) < MIN_ABUND) then
      within_tolerance = .false.
      return
   end if

!  Calculate symmetric relative difference
   rel_diff = 2.0_dp * ABS(cached_val - current_val) / (ABS(cached_val) + ABS(current_val))
   within_tolerance = (rel_diff < tol)

end function IS_WITHIN_TOLERANCE

subroutine writePopulations(fileName,modelNumber)
   character(*), intent(in) :: fileName, modelNumber
   character(LEN=11), allocatable :: populationLabels(:)
   integer :: i,n,p

   allocate(populationLabels(1:SUM((/(coolants(N)%nLevel,N=1,NCOOLANTS)/))))
   p=1
   do n=1,NCOOLANTS
      do I=1,coolants(n)%nLevel
         write(populationLabels(p),"(I3)") I-1
         populationLabels(p)="n("//TRIM(ADJUSTL(coolants(N)%NAME))//","//&
                                 &TRIM(ADJUSTL(populationLabels(p)))//")"
          p=p+1
      end do
   end do

   open(populationsId, file=fileName, status="unknown", action="write")
   write(populationsId,"(A5,999(':',A11))") "MODEL",populationLabels
   write(populationsId,"(A2,999(':',E13.5))") modelNumber,(coolants(N)%POPULATION,N=1,NCOOLANTS)
   close(populationsId)
end subroutine writePopulations


subroutine writeOpacities(fileName, modelNumber)
   character(*), intent(in) :: fileName, modelNumber

   character(LEN=11), allocatable :: LINE_LABELS(:)

   integer  :: I,J,N,P

!  Create the transition label for each emission line
   allocate(LINE_LABELS(1:COUNT((/(coolants(N)%A_COEFF,N=1,NCOOLANTS)/)>0)))
   P=1
   coolant: do N=1,NCOOLANTS  ! Loop over coolants
      levels_i: do I=1,coolants(N)%NLEVEL  ! Loop over levels (i)
         levels_j: do J=1,coolants(N)%NLEVEL  ! Loop over levels (j)
            if(coolants(N)%A_COEFF(I,J)==0) cycle levels_j

!           Assume all coolants with < 6 levels are atoms with fine-structure lines
            if(coolants(N)%NLEVEL<6) then

!              Label fine-structure lines with their wavelength in micron (or in Angstrom, if appropriate)
               if(coolants(N)%FREQUENCY(I,J)>=3.0e14_dp) then
                  write(LINE_LABELS(P),"(I4)") NINT(C/coolants(N)%FREQUENCY(I,J)*1.0e8_dp)  ! Wavelength in Angstrom (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//"A"
               else if(coolants(N)%FREQUENCY(I,J)>=3.0e13_dp) then
                  write(LINE_LABELS(P),"(F5.2)") C/coolants(N)%FREQUENCY(I,J)*1.0e4_dp  ! Wavelength in micron (2 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//"um"
               else if(coolants(N)%FREQUENCY(I,J)>=1.0e13_dp) then
                  write(LINE_LABELS(P),"(F5.1)") C/coolants(N)%FREQUENCY(I,J)*1.0e4_dp  ! Wavelength in micron (1 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//"um"
               else
                  write(LINE_LABELS(P),"(I4)") NINT(C/coolants(N)%FREQUENCY(I,J)*1.0e4_dp)  ! Wavelength in micron (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//"um"
               end if

!              Handle labeling of both neutral and ionized atomic species
               LINE_LABELS(P)="["//coolants(N)%NAME(1:VERIFY(coolants(N)%NAME,"+ ",.true.))  &
                               & //REPEAT("I",COUNT_SUBSTRING(coolants(N)%NAME,"+"))//"I] " &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))

!           Assume all coolants with more levels are molecules with pure rotational or ro-vibrational lines
            else

!              If the line frequency corresponds to a wavelength in the micron range
!              or shorter, then label the line with its wavelength in micron instead
               if(coolants(N)%FREQUENCY(I,J)>=3.0e14_dp) then
                  write(LINE_LABELS(P),"(I4)") NINT(C/coolants(N)%FREQUENCY(I,J)*1.0e8_dp)  ! Wavelength in Angstrom (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(coolants(N)%NAME))//" "//TRIM(ADJUSTL(LINE_LABELS(P)))//"A"
               else if(coolants(N)%FREQUENCY(I,J)>=3.0e13_dp) then
                  write(LINE_LABELS(P),"(F5.2)") C/coolants(N)%FREQUENCY(I,J)*1.0e4_dp  ! Wavelength in micron (2 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(coolants(N)%NAME))//" "//TRIM(ADJUSTL(LINE_LABELS(P)))//"um"
               else if(coolants(N)%FREQUENCY(I,J)>=1.0e13_dp) then
                  write(LINE_LABELS(P),"(F5.1)") C/coolants(N)%FREQUENCY(I,J)*1.0e4_dp  ! Wavelength in micron (1 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(coolants(N)%NAME))//" "//TRIM(ADJUSTL(LINE_LABELS(P)))//"um"
               else if(coolants(N)%FREQUENCY(I,J)>=6.0e12_dp) then
                  write(LINE_LABELS(P),"(I4)") NINT(C/coolants(N)%FREQUENCY(I,J)*1.0e4_dp)  ! Wavelength in micron (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(coolants(N)%NAME))//" "//TRIM(ADJUSTL(LINE_LABELS(P)))//"um"
               else
                  if(J>100) then
                     write(LINE_LABELS(P),"(I4,'-',I3)") I-1,J-1
                  else if(J>10) then
                     write(LINE_LABELS(P),"(I4,'-',I2)") I-1,J-1
                  else
                     write(LINE_LABELS(P),"(I4,'-',I1)") I-1,J-1
                  end if
                  LINE_LABELS(P)=TRIM(ADJUSTL(coolants(N)%NAME))//"(" &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))//")"
                  LINE_LABELS(P)=REPEAT(" ",(LEN(LINE_LABELS(P))-LEN_TRIM(LINE_LABELS(P)))/2) &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))
               end if

            end if
            P=P+1
         end do levels_j  ! End of loop over levels (j)
      end do levels_i  ! End of loop over levels (i)
   end do coolant  ! End of loop over coolants

!  Open and write to the output file
   open(populationsId, file=fileName, status="replace", action="write")
   write(populationsId, "('Particle',999(2X,A11))") LINE_LABELS
   do N=1,NCOOLANTS
      do I=1,coolants(N)%NLEVEL
         levels_j_2: do J=1,coolants(N)%NLEVEL
            if(coolants(N)%A_COEFF(I,J)==0) cycle levels_j_2
            write(populationsId,"(999ES13.5)",advance="NO") coolants(N)%OPACITY(I,J)
         end do levels_j_2
      end do
   end do
   close(populationsId)

   deallocate(LINE_LABELS)

end subroutine writeOpacities

!=======================================================================
!
!  Return the number of non-overlapping occurrences
!  of a given substring within the specified string.
!
!-----------------------------------------------------------------------
function COUNT_SUBSTRING(STRING,SUBSTRING) result(COUNT)


   character(LEN=*), intent(in) :: STRING,SUBSTRING
   integer :: I,COUNT,POSITION

   COUNT=0
   if(LEN(SUBSTRING)==0) return

   I=1
   do
      POSITION=INDEX(STRING(I:),SUBSTRING)
      if(POSITION==0) return
      COUNT=COUNT+1
      I=I+POSITION+LEN(SUBSTRING)
   end do

end function COUNT_SUBSTRING
end module COOLANT_MODULE

