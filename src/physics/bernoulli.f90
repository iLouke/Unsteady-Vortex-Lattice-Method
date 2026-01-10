!> =============================================================================
!> \module   bernoulli
!> \brief    Pressure and Aerodynamic Coefficient Utilities
!> \details  Implements the Unsteady Bernoulli equation for both incompressible 
!>           and compressible flows (isentropic relations).
!> 
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2026
!> =============================================================================
module bernoulli
    use constants
    implicit none

    private
    public :: compute_bernoulli_pressure, compute_pressure_coefficient

contains

    !> =========================================================================
    !> \brief   Computes local static pressure using Bernoulli's Law.
    !> \details Defaults to Incompressible (Drella Eq 1.107). 
    !>          If compressible is true, uses Isentropic relations (Drella Eq 1.112).
    !>          Marked as elemental to support array-based calculations.
    !> 
    !> \param   P_inf         Freestream static pressure [Pa]
    !> \param   Rho_inf       Freestream density [kg/m^3]
    !> \param   V_inf         Freestream velocity magnitude [m/s]
    !> \param   V_local       Local velocity magnitude [m/s]
    !> \param   dphi_dt       (Optional) Unsteady term d(phi)/dt. Default = 0.0
    !> \param   compressible  (Optional) Use isentropic compressibility. Default = F
    !> \param   gamma_val     (Optional) Heat capacity ratio. Defaults to GAMMA
    !> \return  P_local       Calculated local static pressure [Pa]
    !> =========================================================================
    elemental function compute_bernoulli_pressure(P_inf, Rho_inf, V_inf, &
                                                  V_local, dphi_dt,     &
                                                  compressible, gamma_val) result(P_local)
        real(wp), intent(in) :: P_inf, Rho_inf, V_inf, V_local
        real(wp), intent(in), optional :: dphi_dt, gamma_val
        logical,  intent(in), optional :: compressible
        real(wp) :: P_local

        ! Internal local variables
        real(wp) :: t_unsteady, g_ratio
        logical  :: is_comp
        real(wp) :: a_inf_sq, M_inf_sq, bracket, base, expo

        ! 1. Resolve Optional Arguments
        t_unsteady = 0.0_wp; if (present(dphi_dt)) t_unsteady = dphi_dt
        is_comp = .false.  ; if (present(compressible)) is_comp = compressible
        g_ratio = GAMMA    ; if (present(gamma_val)) g_ratio = gamma_val

        ! 2. Calculation Dispatch
        if (.not. is_comp) then
            ! --- INCOMPRESSIBLE (Drella Eq 1.107) ---
            ! p = p_inf + 0.5*rho*(V_inf^2 - V^2) - rho*(dphi/dt)
            P_local = P_inf + (0.5_wp * Rho_inf * (V_inf**2 - V_local**2)) - &
                      (Rho_inf * t_unsteady)
        else
            ! --- COMPRESSIBLE / ISENTROPIC (Drella Eq 1.112) ---
            ! Calculate local speed of sound and Mach number
            a_inf_sq = (g_ratio * P_inf) / Rho_inf
            M_inf_sq = (V_inf**2) / max(a_inf_sq, ZERO)
            
            ! Inner bracket term: [ 1 - (V/V_inf)^2 - (2/V_inf^2)*dphi_dt ]
            bracket = 1.0_wp
            if (V_inf > ZERO) then
                bracket = 1.0_wp - (V_local**2 / V_inf**2) - (2.0_wp * t_unsteady / V_inf**2)
            end if

            ! Final relation: P = P_inf * [ 1 + (g-1)/2 * M^2 * bracket ] ^ (g/(g-1))
            base = 1.0_wp + ((g_ratio - 1.0_wp) / 2.0_wp) * M_inf_sq * bracket
            expo = g_ratio / max(g_ratio - 1.0_wp, ZERO)
            P_local = P_inf * (base ** expo)
        end if
    end function compute_bernoulli_pressure

    !> =========================================================================
    !> \brief   Computes the non-dimensional Pressure Coefficient (Cp).
    !> \details Supports two calculation paths:
    !>          1. Precise: Uses P_local, P_inf, Rho_inf (Valid for all flows).
    !>          2. Velocity Approx: Uses 1 - (V/V_inf)^2 (Incompressible only).
    !> =========================================================================
    elemental function compute_pressure_coefficient(V_inf, V_local, P_local, P_inf, Rho_inf) result(Cp)
        real(wp), intent(in) :: V_inf
        real(wp), intent(in), optional :: V_local, P_local, P_inf, Rho_inf
        real(wp) :: Cp, q_inf

        ! Option 1: Precise Pressure Definition (Highest Priority)
        if (present(P_local) .and. present(P_inf) .and. present(Rho_inf)) then
            q_inf = 0.5_wp * Rho_inf * (V_inf**2)
            Cp = merge((P_local - P_inf) / q_inf, 0.0_wp, q_inf > ZERO)

        ! Option 2: Incompressible Velocity-Squared Approximation
        elseif (present(V_local)) then
            Cp = merge(1.0_wp - (V_local**2 / V_inf**2), 0.0_wp, V_inf > ZERO)

        else
            ! Fallback if insufficient arguments are provided
            Cp = 0.0_wp
        end if
    end function compute_pressure_coefficient

end module bernoulli