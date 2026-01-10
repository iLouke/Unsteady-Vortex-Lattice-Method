!> =============================================================================
!> \module   conservation_law
!> \brief    Fluid Conservation Laws and Flux Utilities
!> \details  Calculates mass, momentum, energy, and enthalpy fluxes across 
!>           control volume surfaces. Implements the Newtonian stress tensor.
!>
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2026
!> =============================================================================
module conservation_law
    use constants
    implicit none

    private
    public :: compute_mass_flux, calculate_momentum_flux
    public :: calculate_internal_energy_flux, compute_enthalpy_flux
    public :: compute_viscous_stress_vector, compute_stress_tensor

contains

    !> =========================================================================
    !> \brief   Computes Mass Flux across a surface.
    !> \details mass_flux = rho * (V . n)
    !> =========================================================================
    pure function compute_mass_flux(density, velocity, normal) result(mass_flux)
        real(wp), intent(in) :: density, velocity(3), normal(3)
        real(wp)             :: mass_flux
        mass_flux = density * dot_product(velocity, normal)
    end function compute_mass_flux

    !> =========================================================================
    !> \brief   Computes Momentum Flux vector across a surface.
    !> \details momentum_flux = [rho * (V . n)] * V
    !> =========================================================================
    pure function calculate_momentum_flux(density, velocity, normal) result(momentum_flux)
        real(wp), intent(in) :: density, velocity(3), normal(3)
        real(wp)             :: momentum_flux(3)
        momentum_flux = (density * dot_product(velocity, normal)) * velocity
    end function calculate_momentum_flux

    !> =========================================================================
    !> \brief   Computes Internal Energy Flux across a surface.
    !> \details energy_flux = rho * (V . n) * e_total
    !> =========================================================================
    pure function calculate_internal_energy_flux(density, velocity, normal, &
                                                 internal_e, total_e) result(energy_flux)
        real(wp), intent(in)           :: density, velocity(3), normal(3), internal_e
        real(wp), intent(in), optional :: total_e
        real(wp)                       :: energy_flux, local_total_e

        if (present(total_e)) then
            local_total_e = total_e
        else
            ! e_total = e_internal + 0.5 * V^2
            local_total_e = internal_e + 0.5_wp * dot_product(velocity, velocity)
        end if

        energy_flux = density * dot_product(velocity, normal) * local_total_e
    end function calculate_internal_energy_flux

    !> =========================================================================
    !> \brief   Computes Enthalpy Flux across a surface.
    !> \details enthalpy_flux = rho * (V . n) * h_total
    !> =========================================================================
    pure function compute_enthalpy_flux(density, velocity, normal, &
                                        internal_h, total_h) result(enthalpy_flux)
        real(wp), intent(in)           :: density, velocity(3), normal(3), internal_h
        real(wp), intent(in), optional :: total_h
        real(wp)                       :: enthalpy_flux, local_total_h

        if (present(total_h)) then
            local_total_h = total_h
        else
            ! h_total = h_static + 0.5 * V^2
            local_total_h = internal_h + 0.5_wp * dot_product(velocity, velocity)
        end if

        enthalpy_flux = density * dot_product(velocity, normal) * local_total_h
    end function compute_enthalpy_flux

    !> =========================================================================
    !> \brief   Computes the Viscous Stress Vector.
    !> \details vector = [Tau] . n
    !> =========================================================================
    pure function compute_viscous_stress_vector(tau_tensor, normal) result(stress_vector)
        real(wp), intent(in) :: tau_tensor(3,3), normal(3)
        real(wp)             :: stress_vector(3)
        stress_vector = matmul(tau_tensor, normal)
    end function compute_viscous_stress_vector

    !> =========================================================================
    !> \brief   Computes the Newtonian Viscous Stress Tensor.
    !> \details Implements Tau_ij = mu * ( (du_i/dx_j + du_j/dx_i) - 2/3 * delta_ij * div(V) )
    !> =========================================================================
    pure function compute_stress_tensor(mu, vel_grad) result(tau)
        real(wp), intent(in) :: mu, vel_grad(3,3)
        real(wp)             :: tau(3,3)
        real(wp)             :: div_v, bulk_term
        real(wp), parameter  :: TWO_THIRDS = 2.0_wp / 3.0_wp

        ! Divergence of velocity field (Trace of velocity gradient)
        div_v = vel_grad(1,1) + vel_grad(2,2) + vel_grad(3,3)
        bulk_term = TWO_THIRDS * div_v

        ! --- Normal Stresses (Diagonal) ---
        tau(1,1) = mu * (2.0_wp * vel_grad(1,1) - bulk_term)
        tau(2,2) = mu * (2.0_wp * vel_grad(2,2) - bulk_term)
        tau(3,3) = mu * (2.0_wp * vel_grad(3,3) - bulk_term)

        ! --- Shear Stresses (Off-Diagonal, Symmetric) ---
        tau(1,2) = mu * (vel_grad(1,2) + vel_grad(2,1))
        tau(1,3) = mu * (vel_grad(1,3) + vel_grad(3,1))
        tau(2,3) = mu * (vel_grad(2,3) + vel_grad(3,2))

        ! Symmetric properties
        tau(2,1) = tau(1,2)
        tau(3,1) = tau(1,3)
        tau(3,2) = tau(2,3)

    end function compute_stress_tensor

end module conservation_law