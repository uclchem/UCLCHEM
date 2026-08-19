module CONSTANTS
    implicit none

    public

    integer, parameter :: dp = KIND(1.0D+0)
    real(dp), parameter :: C  = 2.99792458e+10_dp  !Speed of light in cgs
    real(dp), parameter :: K_BOLTZ = 1.38065040e-16_dp  ! Boltzmann constant cgs
    real(dp), parameter :: HP = 6.62606896e-27_dp  !Planck constant in cgs
    real(dp), parameter :: HP_SI = 6.62607015e-34_dp  !Planck constant in SI
    real(dp), parameter :: REDUCED_PLANCK=1.054571628e-27_dp
    real(dp), parameter :: MH = 1.67262164e-24_dp  !H nucleus mass in cgs
    real(dp), parameter :: AMU=1.66053892e-24_dp  !atomic mass unit in cgs
    real(dp), parameter :: PI = 3.141592654_dp
    real(dp), parameter :: K_BOLTZ_SI=1.38e-23_dp  !Boltzmann constant SI
    real(dp), parameter :: PC=3.086e18_dp  !parsec in cgs
    real(dp), parameter :: au=2.063e5_dp  !1 AU in cgs
    real(dp), parameter :: KM=1.0e5_dp  !kilometer in cgs
    real(dp), parameter :: SECONDS_PER_YEAR=3.16e7_dp
    real(dp), parameter :: T_CMB=2.73_dp  ! Cosmic microwave background temperature in K
    real(dp), parameter :: ZETA_0=1.3e-17_dp  ! Standard cosmic ray ionization rate in s-1
    real(dp), parameter :: EV = 1.60217646e-12_dp  ! electron volt in erg
    real(dp), parameter :: GRAV_G = 6.674e-8_dp  !gravitational constant in cgs
    real(dp), parameter :: SB_CONST=5.6704e-5_dp  !Stefan Boltzmann constant in cgs
    real(dp), parameter :: HABING_TO_DRAINE = 1.0_dp / 1.7_dp  !conversion factor from Habing to Draine field
    real(dp), parameter :: Lsun = 3.828e+33_dp  ! Sun luminosity in cgs
    real(dp), parameter :: aunit = 1.495978e13_dp  ! AU in cm
    real(dp), parameter :: N_AVOGADRO=6.022140857e23_dp  !Avogadro constant
    real(dp), parameter :: KCAL_TO_JOULE=4.184e3_dp  !Constant to convert kcal to J
    real(dp), parameter :: uISRF = 8.64e-13_dp  !Energy density of the interstellar radiation field in cgs
    real(dp), parameter :: uISRF_UV = 5.29e-14_dp  !Energy density of the UV interstellar radiation field in cgs
    real(dp), parameter :: MIN_ABUND = 1.0e-30_dp  !Minimum abundance allowed

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
