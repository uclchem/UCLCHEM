MODULE CONSTANTS
   INTEGER, parameter :: dp = KIND(1.0D+0)
   REAL(dp), parameter :: C  = 2.99792458D+10 !Speed of light in cgs
   REAL(dp), PARAMETER :: K_BOLTZ = 1.38065040D-16 ! Boltzmann constant cgs
   REAL(dp), PARAMETER :: HP = 6.62606896D-27 !Planck constant in cgs
   REAL(dp), PARAMETER :: HP_SI = 6.62607015D-34 !Planck constant in SI
   REAL(dp), PARAMETER :: REDUCED_PLANCK=1.054571628d-27
   REAL(dp), PARAMETER :: MH = 1.67262164D-24 !H nucleus mass in cgs
   REAL(dp), PARAMETER :: AMU=1.66053892d-24 !atomic mass unit in cgs
   REAL(dp), PARAMETER :: PI = 3.141592654
   REAL(dp), PARAMETER :: K_BOLTZ_SI=1.38d-23 !Boltzmann constant SI
   REAL(dp), PARAMETER :: PC=3.086d18 !parsec in cgs
   REAL(dp), PARAMETER :: au=2.063d5 !1 AU in cgs
   REAL(dp), PARAMETER :: KM=1.d5 !kilometre in cgs
   REAL(dp), PARAMETER :: SECONDS_PER_YEAR=3.16d7
   REAL(dp), PARAMETER :: T_CMB=2.73
   REAL(dp), PARAMETER :: EV = 1.60217646D-12 ! electron volt in erg
   REAL(dp), PARAMETER :: GRAV_G = 6.674d-8 !gravitational constant in cgs
   REAL(dp), PARAMETER :: SB_CONST=5.6704d-5 !Stefan Boltzmann constant in cgs
   REAL(dp), PARAMETER :: HABING_TO_DRAINE = 1 / 1.7 !conversion factor from Habing to Draine field

   !Error codes for python wrap
   INTEGER, PARAMETER :: PARAMETER_READ_ERROR=-1
   INTEGER, PARAMETER :: PHYSICS_INIT_ERROR=-2
   INTEGER, PARAMETER :: CHEM_INIT_ERROR=-3
   INTEGER, PARAMETER :: INT_UNRECOVERABLE_ERROR=-4
   INTEGER, PARAMETER :: INT_TOO_MANY_FAILS_ERROR=-5
   INTEGER, PARAMETER :: NOT_ENOUGH_TIMEPOINTS_ERROR=-6
   INTEGER, PARAMETER :: PHYSICS_UPDATE_ERROR=-7
   
CONTAINS
SUBROUTINE DUMMY_THREE(dummy_three_output)
   integer, intent(out) :: dummy_three_output
   dummy_three_output = 1
END SUBROUTINE DUMMY_THREE

END MODULE CONSTANTS