module mod_atmosphere
    ! --------------------------------------------------
    ! Module for Atmospheric Properties and Calculations
    ! Sources :
    !   Mark Drella, "Flight Vehicle Aerodynamics", 2014
    ! --------------------------------------------------
    use mod_precision
    use mod_constants
    implicit none
    private

    public :: compute_isa_properties

contains

    ! ------------------------------------------------------------------
    ! Standard Atmosphere Model (Troposphere < 47km)
    ! ------------------------------------------------------------------
    subroutine compute_isa_properties(altitude, temp, pressure, rho, mu, speed_of_sound, verbose, log_unit)
        ! Inputs:
        ! altitude  - Altitude in meters
        ! Outputs:
        ! rho       - Air density (kg/m3)
        ! temp      - Air temperature (K)
        ! pressure  - Air pressure (Pa)
        real(wp), intent(in)            :: altitude
        logical , intent(in), optional  :: verbose
        integer , intent(in), optional  :: log_unit
        real(wp), intent(out)           :: temp, pressure, rho, mu, speed_of_sound
        
        ! 1. Temperature Model
        temp = get_temperature(altitude)
        ! 2. Pressure Model
        pressure = get_pressure(altitude)
        ! 3. Density (Ideal Gas Law)
        rho = get_density(pressure, temp)
        ! 4. Viscosity (Sutherland's Law)
        mu = get_viscosity(temp)
        ! 5. Speed of Sound
        ! Assuming gamma = 1.4 for air
        speed_of_sound = get_speed_of_sound(p = pressure, rho = rho)

        ! Optional Verbose Output
        if (present(verbose)) then
            if (verbose) then
                call write_isa_data(altitude, temp, pressure, rho, mu, speed_of_sound)
            end if
        end if
        ! Optional Log Output
        if (present(log_unit)) then
            call log_isa_data(altitude, temp, pressure, rho, mu, speed_of_sound, log_unit)
        end if
    end subroutine compute_isa_properties

    ! --------------------------------------------------
    ! Compute Temperature at given altitude
    ! --------------------------------------------------
    real(wp) function get_temperature(altitude)
        ! Equation (1.3)
        ! Altitude in meters
        real(wp), intent(in) :: altitude
        
        real(wp) :: z           ! Altitude in km
        real(wp) :: logarithm                    ! Helper variable
        real(wp) :: exponential_1, exponential_2 ! Helper variables
        
        z = altitude / 1000.0_wp  ! Convert to km for calculation
        if (z < 47.0_wp) then


            exponential_1   = exp(35.75_wp - 3.25_wp*z)
            exponential_2   = exp( -3.0_wp + 0.0003_wp*z**3)
            logarithm       = log(1.0_wp + exponential_1 + exponential_2)

            get_temperature = 216.65_wp + 2.0_wp * logarithm
        else
            ! Above 47 km, constant temperature (Stratosphere)
            get_temperature = 216.65_wp
        end if
    end function get_temperature

    ! --------------------------------------------------
    ! Compute Pressure at given altitude
    ! --------------------------------------------------
    real(wp) function get_pressure(altitude)
        ! Equation (1.2)
        ! Altitude in meters
        real(wp), intent(in) :: altitude
        
        real(wp) :: z, exponential, fraction

        z = altitude / 1000.0_wp  ! Convert to km for calculation

        ! Split the expression for numerical stability
        fraction =  (0.0015_wp * z**2)/(1.0_wp - 0.018_wp*z + 0.0011_wp*z**2)
        exponential =  exp(-0.118_wp * z - fraction)
        ! Final Pressure Calculation
        get_pressure = P_SL * exponential
    end function get_pressure

    ! --------------------------------------------------
    ! Compute Density using Ideal Gas Law
    ! --------------------------------------------------
    real(wp) function get_density(pressure, temperature)
        ! Equation (1.7)
        ! rho = P / (R * T)
        real(wp), intent(in) :: pressure, temperature

        get_density = pressure / (R_GAS * temperature)
    end function get_density

    ! --------------------------------------------------
    ! Compute Dynamic Viscosity using Sutherland's Law
    ! --------------------------------------------------
    real(wp) function get_viscosity(temperature)
        ! Equation (1.4)
        real(wp), intent(in) :: temperature
        real(wp) :: T                       ! Helper variable
        real(wp) :: fraction_1, fraction_2  ! Helper variables


        ! Change the variable name for clarity and readability
        T = temperature

        fraction_1 = (T / T_s)**1.5_wp
        fraction_2 = (T_SL + T_s)/(T + T_s)
        ! Final Viscosity Calculation
        get_viscosity = MU_SL * fraction_1 * fraction_2
    end function get_viscosity

    ! --------------------------------------------------
    ! Compute Speed of Sound
    ! --------------------------------------------------
    real(wp) function get_speed_of_sound(p, rho, T, h)
        ! Equation (1.70)
        ! GAMMA - Ratio of Specific Heats, is defined in mod_constants
        ! R_GAS - Specific Gas Constant  , is defined in mod_constants
        real(wp), intent(in), optional :: p, rho, T, h

        if (present(p) .and. present(rho)) then
            get_speed_of_sound = sqrt(GAMMA * p / rho)
            return
        end if
        if (present(T)) then
            get_speed_of_sound = sqrt(GAMMA * R_GAS * T)
            return
        end if
        if (present(h)) then
            get_speed_of_sound = sqrt((GAMMA - 1.0_wp) * h)
        end if

        get_speed_of_sound = sqrt(GAMMA * R_GAS * T)
    end function get_speed_of_sound


    ! --------------------------------------------------
    ! HELPER SUBROUTINES
    ! Logging and Writing ISA Data
    ! --------------------------------------------------

    subroutine log_isa_data(altitude, temp, pressure, rho, mu, speed_of_sound, log_unit)
        real(wp), intent(in) :: altitude, temp, pressure, rho, mu, speed_of_sound
        integer , intent(in) :: log_unit

        write(log_unit,'(a)') ' '
        write(log_unit,'(a)') '--- ISA Atmospheric Properties ---'
        write(log_unit,'(a)') ' '
        write(log_unit,"(a,f10.3)") 'Altitude          (m)    : ', altitude
        write(log_unit,"(a,f10.3)") 'Temperature       (K)    : ', temp
        write(log_unit,"(a,f10.3)") 'Pressure          (Pa)   : ', pressure
        write(log_unit,"(a,f10.3)") 'Density           (kg/m3): ', rho
        write(log_unit,"(a,f10.8)") 'Dynamic Viscosity (Pa.s) : ', mu
        write(log_unit,"(a,f10.3)") 'Speed of Sound    (m/s)  : ', speed_of_sound
    end subroutine log_isa_data

    subroutine write_isa_data(altitude, temp, pressure, rho, mu, speed_of_sound)
        real(wp), intent(in) :: altitude, temp, pressure, rho, mu, speed_of_sound

        write(*,'(a)') ' '
        write(*,'(a)') '--- ISA Atmospheric Properties ---'
        write(*,'(a)') ' '
        write(*,"(a,f10.3)") 'Altitude (m): ', altitude
        write(*,"(a,f10.3)") 'Temperature (K): ', temp
        write(*,"(a,f10.3)") 'Pressure (Pa): ', pressure
        write(*,"(a,f10.3)") 'Density (kg/m3): ', rho
        write(*,"(a,f10.8)") 'Dynamic Viscosity (Pa.s): ', mu
        write(*,"(a,f10.3)") 'Speed of Sound (m/s): ', speed_of_sound
    end subroutine write_isa_data
    
end module mod_atmosphere