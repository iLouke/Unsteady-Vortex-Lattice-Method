program test_conservation_law
    use mod_precision
    use mod_conservation_law
    implicit none

    real(wp) :: rho, mu
    real(wp) :: vel(3), normal(3)
    real(wp) :: mass_flux_res
    real(wp) :: gradV(3,3), stress(3,3)

    print *, "=========================================="
    print *, "       TEST: CONSERVATION LAWS            "
    print *, "=========================================="

    ! -----------------------------------------------------------
    ! TEST 1: Mass Flux
    ! -----------------------------------------------------------
    rho = 1.225_wp
    vel = [10.0_wp, 0.0_wp, 0.0_wp]
    normal = [1.0_wp, 0.0_wp, 0.0_wp]
    
    mass_flux_res = compute_mass_flux(rho, vel, normal)
    call assert_approx(mass_flux_res, 12.25_wp, "Mass Flux (Aligned)")

    ! -----------------------------------------------------------
    ! TEST 2: Viscous Stress
    ! -----------------------------------------------------------
    mu = 0.01_wp
    gradV = 0.0_wp
    gradV(1,2) = 1.0_wp ! du/dy = 1.0
    
    stress = compute_stress_tensor(mu, gradV)
    
    call assert_approx(stress(1,2), 0.01_wp, "Stress Tau_xy")
    call assert_approx(stress(2,1), 0.01_wp, "Stress Tau_yx (Symmetry)")

    print *, "------------------------------------------"
    print *, " Conservation Laws Tests Finished "
    print *, "------------------------------------------"

contains

    subroutine assert_approx(actual, expected, test_name)
        real(wp), intent(in) :: actual, expected
        character(len=*), intent(in) :: test_name
        real(wp), parameter :: tol = 1.0e-5_wp

        if (abs(actual - expected) < tol) then
            print *, "[PASS] ", test_name
        else
            print *, "[FAIL] ", test_name
            stop 1
        end if
    end subroutine assert_approx

end program test_conservation_law