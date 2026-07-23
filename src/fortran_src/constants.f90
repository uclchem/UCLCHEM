module CONSTANTS
   integer, parameter :: dp = KIND(1.0D+0)
   real(dp), parameter :: C  = 2.99792458D+10  !Speed of light in cgs
   real(dp), parameter :: K_BOLTZ = 1.38065040D-16  ! Boltzmann constant cgs
   real(dp), parameter :: HP = 6.62606896D-27  !Planck constant in cgs
   real(dp), parameter :: HP_SI = 6.62607015D-34  !Planck constant in SI
   real(dp), parameter :: REDUCED_PLANCK=1.054571628d-27
   real(dp), parameter :: MH = 1.67262164D-24  !H nucleus mass in cgs
   real(dp), parameter :: AMU=1.66053892d-24  !atomic mass unit in cgs
   real(dp), parameter :: PI = 3.141592654
   real(dp), parameter :: K_BOLTZ_SI=1.38d-23  !Boltzmann constant SI
   real(dp), parameter :: PC=3.086d18  !parsec in cgs
   real(dp), parameter :: au=2.063d5  !1 AU in cgs
   real(dp), parameter :: KM=1.0d5  !kilometer in cgs
   real(dp), parameter :: SECONDS_PER_YEAR=3.16d7
   real(dp), parameter :: T_CMB=2.73
   real(dp), parameter :: EV = 1.60217646D-12  ! electron volt in erg
   real(dp), parameter :: GRAV_G = 6.674d-8  !gravitational constant in cgs
   real(dp), parameter :: SB_CONST=5.6704d-5  !Stefan Boltzmann constant in cgs
   real(dp), parameter :: HABING_TO_DRAINE = 1 / 1.7  !conversion factor from Habing to Draine field
   real(dp), parameter :: Lsun = 3.828d+33  ! Sun luminosity in cgs
   real(dp), parameter :: aunit = 1.495978d13  ! AU in cm
   real(dp), parameter :: N_AVOGADRO=6.022140857d23  !Avogadro constant
   real(dp), parameter :: KCAL_TO_JOULE=4.184d3  !Constant to convert kcal to J
   real(dp), parameter :: uISRF = 8.64d-13  !Energy density of the interstellar radiation field in cgs
   real(dp), parameter :: uISRF_UV = 5.29e-14  !Energy density of the UV interstellar radiation field in cgs

   !Sentinel value used in network.f90 for absent reaction types
   integer, parameter :: REAC_NOT_PRESENT = 99999

   !Error codes for python wrap
   integer, parameter :: PARAMETER_READ_ERROR=-1
   integer, parameter :: PHYSICS_INIT_ERROR=-2
   integer, parameter :: CHEM_INIT_ERROR=-3
   integer, parameter :: INT_UNRECOVERABLE_ERROR=-4
   integer, parameter :: INT_TOO_MANY_FAILS_ERROR=-5
   integer, parameter :: NOT_ENOUGH_TIMEPOINTS_ERROR=-6
   integer, parameter :: PHYSICS_UPDATE_ERROR=-7
   integer, parameter :: SOLVER_STATS_OVERFLOW_ERROR=-8
   integer, parameter :: COOLANT_FILE_ERROR=-9
   integer, parameter :: COOLANT_DATA_ERROR=-10
   integer, parameter :: COOLANT_FREQ_TOL_ERROR=-11
   integer, parameter :: COOLANT_POP_TOL_ERROR=-12
   integer, parameter :: COOLANT_SOLVER_ERROR=-13
   integer, parameter :: COOLANT_CONFIG_ERROR=-14
   integer, parameter :: NEGATIVE_ABUNDANCE_ERROR=-15
   integer, parameter :: CONSERVATION_ERROR=-16
   integer, parameter :: ZERO_INNER_RADIUS_ERROR=-17

contains
subroutine DUMMY_THREE(dummy_three_output)
   integer, intent(out) :: dummy_three_output
   dummy_three_output = 1
end subroutine DUMMY_THREE

end module CONSTANTS
