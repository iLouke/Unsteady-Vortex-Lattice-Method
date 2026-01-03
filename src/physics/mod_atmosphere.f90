!> =============================================================================
!> \module   mod_atmosphere
!> \brief    ISA Standard Atmosphere Model (Object-Oriented)
!> \details  Calculates atmospheric properties up to 47km (Troposphere/Stratosphere).
!>           Stores state in a derived type for easy access by other modules.
!> 
!> \sources  Mark Drella, "Flight Vehicle Aerodynamics", 2014
!>
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2025
!> =============================================================================
module mod_atmosphere
    use mod_precision
    use mod_constants
    implicit none
    private

    ! Public Interface
    public :: t_atmosphere

    !> \type t_atmosphere
    !> \brief Container for atmospheric state and methods
    type :: t_atmosphere
        ! --- State Variables ---
        real(wp) :: altitude       = 0.0_wp  ! [m]
        real(wp) :: temperature    = 0.0_wp  ! [K]
        real(wp) :: pressure       = 0.0_wp  ! [Pa]
        real(wp) :: density        = 0.0_wp  ! [kg/m^3]
        real(wp) :: viscosity      = 0.0_wp  ! [Pa.s]
        real(wp) :: speed_of_sound = 0.0_wp  ! [m/s]

    contains
        ! --- Methods ---
        !> Updates the state based on a new altitude
        procedure :: update  => atm_update
        !> Prints current state to console or file
        procedure :: display => atm_display
    end type t_atmosphere

contains

    !> =========================================================================
    !> \brief   Update atmospheric properties for a given altitude.
    !> \details Populates the components of the derived type.
    !> \param   altitude  Input altitude in meters.
    !> =========================================================================
    subroutine atm_update(self, altitude)
        class(t_atmosphere), intent(inout) :: self
        real(wp),            intent(in)    :: altitude

        ! 1. Store Altitude
        self%altitude = altitude

        ! 2. Compute Properties using private pure functions
        self%temperature    = calc_temperature(altitude)
        self%pressure       = calc_pressure(altitude)
        self%density        = calc_density(self%pressure, self%temperature)
        self%viscosity      = calc_viscosity(self%temperature)
        
        ! 3. Compute Speed of Sound (Adiabatic assumption)
        self%speed_of_sound = sqrt(GAMMA * self%pressure / self%density)

    end subroutine atm_update

    !> =========================================================================
    !> \brief   Display the current state of the atmosphere object.
    !> \param   unit  (Optional) File unit number to write to. Default is stdout.
    !> =========================================================================
    subroutine atm_display(self, unit)
        class(t_atmosphere), intent(in) :: self
        integer, intent(in), optional   :: unit
        integer :: u

        ! Default to standard output (6) if unit not provided
        u = 6
        if (present(unit)) u = unit

        write(u, '(A)') " "
        write(u, '(A)') "--- ISA Atmospheric Properties ---"
        write(u, '(A, T30, F12.3, 1X, A)') "  Altitude:",       self%altitude,       "m"
        write(u, '(A, T30, F12.3, 1X, A)') "  Temperature:",    self%temperature,    "K"
        write(u, '(A, T30, ES12.4, 1X, A)') "  Pressure:",       self%pressure,       "Pa"
        write(u, '(A, T30, F12.4, 1X, A)') "  Density:",        self%density,        "kg/m^3"
        write(u, '(A, T30, ES12.4, 1X, A)') "  Dyn. Viscosity:", self%viscosity,      "Pa.s"
        write(u, '(A, T30, F12.3, 1X, A)') "  Speed of Sound:", self%speed_of_sound, "m/s"
        write(u, '(A)') "----------------------------------"
        write(u, '(A)') " "
    end subroutine atm_display


    ! ==========================================================================
    ! PRIVATE PURE HELPER FUNCTIONS
    ! ==========================================================================

    !> \brief Equations (1.3) - Temperature Model
    pure function calc_temperature(alt) result(T)
        real(wp), intent(in) :: alt
        real(wp) :: T, z, term1, term2
        
        z = alt * 1.0e-3_wp ! Convert to km

        if (z < 47.0_wp) then
            term1 = exp(35.75_wp - 3.25_wp * z)
            term2 = exp(-3.0_wp + 0.0003_wp * z**3)
            T = 216.65_wp + 2.0_wp * log(1.0_wp + term1 + term2)
        else
            T = 216.65_wp ! Stratosphere constant
        end if
    end function calc_temperature

    !> \brief Equation (1.2) - Pressure Model
    pure function calc_pressure(alt) result(P)
        real(wp), intent(in) :: alt
        real(wp) :: P, z, frac, expo
        
        z = alt * 1.0e-3_wp ! Convert to km
        
        frac = (0.0015_wp * z**2) / (1.0_wp - 0.018_wp*z + 0.0011_wp*z**2)
        expo = exp(-0.118_wp * z - frac)
        P    = P_SL * expo
    end function calc_pressure

    !> \brief Equation (1.7) - Ideal Gas Law
    pure function calc_density(P, T) result(rho)
        real(wp), intent(in) :: P, T
        real(wp) :: rho
        rho = P / (R_GAS * T)
    end function calc_density

    !> \brief Equation (1.4) - Sutherland's Law
    pure function calc_viscosity(T) result(mu)
        real(wp), intent(in) :: T
        real(wp) :: mu
        mu = MU_SL * ((T / T_s)**1.5_wp) * ((T_SL + T_s) / (T + T_s))
    end function calc_viscosity

end module mod_atmosphere