! module mod_atmosphere
!     use mod_precision
!     use mod_constants ! Assuming you have PI, R_GAS, GAMMA here
!     implicit none
!     private

!     public :: compute_isa_properties, compute_reynolds, compute_dynamic_viscosity

!     ! Standard Sea Level Constants (ISA)
!     real(wp), parameter :: T0_SL = 288.15_wp      ! Temp (K)
!     real(wp), parameter :: P0_SL = 101325.0_wp    ! Pressure (Pa)
!     real(wp), parameter :: RHO0_SL = 1.225_wp     ! Density (kg/m3)
!     real(wp), parameter :: G_ACCEL = 9.80665_wp   ! Gravity (m/s2)
!     real(wp), parameter :: R_AIR = 287.05_wp      ! Gas Constant (J/kgK)
!     real(wp), parameter :: LAPSE_RATE = -0.0065_wp ! Temp lapse rate (K/m)

! contains

!     ! ------------------------------------------------------------------
!     ! Standard Atmosphere Model (Troposphere < 11km)
!     ! ------------------------------------------------------------------
!     subroutine compute_isa_properties(altitude, rho, temp, pressure, speed_of_sound)
!         real(wp), intent(in) :: altitude           ! Meters
!         real(wp), intent(out) :: rho, temp, pressure, speed_of_sound
        
!         real(wp) :: theta ! Temperature ratio

!         ! 1. Temperature Model (T = T0 + L * h)
!         temp = T0_SL + (LAPSE_RATE * altitude)

!         ! 2. Pressure Model
!         theta = temp / T0_SL
!         pressure = P0_SL * (theta ** (-G_ACCEL / (LAPSE_RATE * R_AIR)))

!         ! 3. Density (Ideal Gas Law: rho = P / RT)
!         rho = pressure / (R_AIR * temp)

!         ! 4. Speed of Sound (a = sqrt(gamma * R * T))
!         ! Assuming gamma = 1.4 for air
!         speed_of_sound = sqrt(1.4_wp * R_AIR * temp)

!     end subroutine compute_isa_properties

!     ! ------------------------------------------------------------------
!     ! Sutherland's Law for Dynamic Viscosity
!     ! ------------------------------------------------------------------
!     function compute_dynamic_viscosity(temp) result(mu)
!         real(wp), intent(in) :: temp ! Kelvin
!         real(wp) :: mu
        
!         real(wp), parameter :: MU_REF = 1.716e-5_wp
!         real(wp), parameter :: T_REF = 273.15_wp
!         real(wp), parameter :: S_CONST = 110.4_wp

!         mu = MU_REF * ((temp / T_REF)**1.5_wp) * ( (T_REF + S_CONST) / (temp + S_CONST) )
!     end function compute_dynamic_viscosity

!     ! ------------------------------------------------------------------
!     ! Reynolds Number Calculator
!     ! ------------------------------------------------------------------
!     function compute_reynolds(rho, V, L, mu) result(Re)
!         real(wp), intent(in) :: rho, V, L, mu
!         real(wp) :: Re
        
!         if (mu < 1.0e-15_wp) then
!             Re = 0.0_wp ! Avoid division by zero
!         else
!             Re = (rho * V * L) / mu
!         end if
!     end function compute_reynolds

! end module mod_atmosphere