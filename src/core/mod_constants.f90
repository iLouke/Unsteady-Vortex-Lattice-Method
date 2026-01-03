!> =============================================================================
!> \module   mod_constants
!> \brief    Physical and Mathematical Constants for ROS
!> \details  Provides fundamental constants used throughout the simulation:
!>           - Mathematical constants
!>           - Air properties at sea level
!>           - Gas constants
!>           - Gravitational acceleration
!>
!> \sources  Mark Drella, "Flight Vehicle Aerodynamics", 2014
!>
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2025
!> =============================================================================
module mod_constants
    use mod_precision
    implicit none
    private

    ! Mathematical Constants
    public :: PI, TWO_PI, DEG2RAD, RAD2DEG, ZERO
    ! Air Properties at Sea Level
    public :: RHO_SL, P_SL, T_SL, SPEED_OF_SOUND_SL, MU_SL, T_s, R_GAS, GAMMA
    ! Gravitational Acceleration
    public :: G_ACCEL

    real(wp), parameter :: PI       = 4.0_wp * atan(1.0_wp)
    real(wp), parameter :: TWO_PI   = 2.0_wp * PI
    real(wp), parameter :: DEG2RAD  = PI / 180.0_wp
    real(wp), parameter :: RAD2DEG  = 180.0_wp / PI
    real(wp), parameter :: ZERO  = 1.0e-12_wp

    ! Air Properties at Sea Level - Mark Drella, Equation (1.1)
    real(wp), parameter :: RHO_SL               = 1.225_wp     ! Density (kg/m3)
    real(wp), parameter :: P_SL                 = 101325.0_wp  ! Pressure (Pa)
    real(wp), parameter :: T_SL                 = 288.15_wp    ! Temperature (K)
    real(wp), parameter :: SPEED_OF_SOUND_SL    = 340.29_wp    ! Speed of Sound (m/s)
    real(wp), parameter :: MU_SL                = 1.7894e-5_wp ! Dynamic Viscosity (kg/m-s)

    ! Specific Gas Constant for Air (J/kgK) - Mark Drella, Equation (1.7)
    real(wp), parameter :: R_GAS = 287.04_wp
    ! Sutherland's Constant for Air (K) - Mark Drella, Equation (1.4)
    real(wp), parameter :: T_s = 110.4_wp

    ! Ratio of Specific Heats for Air
    real(wp), parameter :: GAMMA = 1.4_wp

    ! Gravitational Acceleration (m/s^2)
    real(wp), parameter :: G_ACCEL = 9.80665_wp

end module mod_constants