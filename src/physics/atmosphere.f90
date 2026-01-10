module atmosphere
    use constants
    USE logger_io
    implicit none
    private

    ! Public Interface
    public :: Air

    !> \type Air
    type :: Air
        real(wp) :: altitude       = 0.0_wp  ! [m]
        real(wp) :: temperature    = 0.0_wp  ! [K]
        real(wp) :: pressure       = 0.0_wp  ! [Pa]
        real(wp) :: density        = 0.0_wp  ! [kg/m^3]
        real(wp) :: viscosity      = 0.0_wp  ! [Pa.s]
        real(wp) :: speed_of_sound = 0.0_wp  ! [m/s]
    contains
        procedure :: update  => atm_update
        procedure :: display => atm_display
    end type Air

contains

    !> Update atmospheric properties
    subroutine atm_update(self, altitude)
        class(Air), intent(inout) :: self
        real(wp),            intent(in)    :: altitude

        self%altitude = altitude
        self%temperature    = calc_temperature(altitude)
        self%pressure       = calc_pressure(altitude)
        self%density        = calc_density(self%pressure, self%temperature)
        self%viscosity      = calc_viscosity(self%temperature)
        self%speed_of_sound = sqrt(GAMMA * self%pressure / self%density)
    end subroutine atm_update

    !> Display state (Formatted to match Config_IO)
    subroutine atm_display(self)
        class(Air), intent(in) :: self
        CHARACTER(LEN=256) :: Line

        ! template WRITE(Line, '(A30, F15.4, 1X, A)') "   Ref Area (S):",  Cfg%RefArea,   "m^2"
        CALL Write_Log(" [ATMOSPHERE (ISA)]")
        !Altitude
        WRITE(Line, '(A30, F15.4, 1X, A)') "   Altitude:",       self%altitude,       "m"
        CALL Write_Log(Line)
        ! Temperature
        WRITE(Line, '(A30, F15.4, 1X, A)') "   Temperature:",    self%temperature,    "K"
        CALL Write_Log(Line)
        ! Pressure
        WRITE(Line, '(A30, ES15.4, 1X, A)') "   Pressure:",       self%pressure,       "Pa"
        CALL Write_Log(Line)
        ! Density
        WRITE(Line, '(A30, F15.4, 1X, A)') "   Density:",        self%density,        "kg/m^3"
        CALL Write_Log(Line)
        ! Dynamic Viscosity
        WRITE(Line, '(A30, ES15.4, 1X, A)') "   Dyn. Viscosity:", self%viscosity,      "Pa.s"
        CALL Write_Log(Line)
        ! Speed of Sound
        WRITE(Line, '(A30, F15.4, 1X, A)') "   Speed of Sound:", self%speed_of_sound, "m/s"
        CALL Write_Log(Line)
        CALL Write_Log(" ")
    end subroutine atm_display

    ! -- Private Helper Functions --
    pure function calc_temperature(alt) result(T)
        real(wp), intent(in) :: alt
        real(wp) :: T, z, term1, term2
        z = alt * 1.0e-3_wp 
        if (z < 47.0_wp) then
            term1 = exp(35.75_wp - 3.25_wp * z)
            term2 = exp(-3.0_wp + 0.0003_wp * z**3)
            T = 216.65_wp + 2.0_wp * log(1.0_wp + term1 + term2)
        else
            T = 216.65_wp 
        end if
    end function calc_temperature

    pure function calc_pressure(alt) result(P)
        real(wp), intent(in) :: alt
        real(wp) :: P, z, frac, expo
        z = alt * 1.0e-3_wp
        frac = (0.0015_wp * z**2) / (1.0_wp - 0.018_wp*z + 0.0011_wp*z**2)
        expo = exp(-0.118_wp * z - frac)
        P    = P_SL * expo
    end function calc_pressure

    pure function calc_density(P, T) result(rho)
        real(wp), intent(in) :: P, T
        real(wp) :: rho
        rho = P / (R_GAS * T)
    end function calc_density

    pure function calc_viscosity(T) result(mu)
        real(wp), intent(in) :: T
        real(wp) :: mu
        mu = MU_SL * ((T / T_s)**1.5_wp) * ((T_SL + T_s) / (T + T_s))
    end function calc_viscosity

end module atmosphere