program test_bernoulli
    use mod_precision
    use mod_constants
    use mod_bernoulli
    implicit none

    real(wp) :: P_inf, V_inf, rho
    real(wp) :: V_loc, P_loc, P_expected
    real(wp) :: Cp_res, q_inf
    
    ! Array testing (Verification of Elemental attribute)
    real(wp) :: V_arr(3), Cp_arr(3)
    integer  :: fail_count = 0

    print *, "=========================================="
    print *, "          TEST: MOD_BERNOULLI             "
    print *, "=========================================="

    ! Setup Standard Sea Level Conditions
    P_inf = P_SL
    V_inf = 100.0_wp
    rho   = RHO_SL

    ! -----------------------------------------------------------
    ! TEST 1: Incompressible Stagnation Pressure
    ! -----------------------------------------------------------
    print *, "[TEST 1] Incompressible Stagnation Point..."
    V_loc = 0.0_wp
    P_expected = P_inf + 0.5_wp * rho * (V_inf**2)
    
    P_loc = compute_bernoulli_pressure(P_inf, rho, V_inf, V_loc, compressible=.false.)
    call assert_approx(P_loc, P_expected, "Stagnation Pressure", fail_count)

    ! -----------------------------------------------------------
    ! TEST 2: Compressible Pressure (Mach ~ 0.3)
    ! -----------------------------------------------------------
    print *, "[TEST 2] Compressible Isentropic Pressure..."
    ! At Mach 0.29, compressibility effects should be small but non-zero
    P_loc = compute_bernoulli_pressure(P_inf, rho, V_inf, V_loc, compressible=.true.)
    
    ! Check that compressible stagnation pressure is higher than incompressible
    if (P_loc > P_expected) then
        print *, "   [PASS] Compressibility correction detected (+)"
    else
        print *, "   [FAIL] Compressibility correction missing"
        fail_count = fail_count + 1
    end if

    ! -----------------------------------------------------------
    ! TEST 3: Vectorized Pressure Coefficient
    ! -----------------------------------------------------------
    print *, "[TEST 3] Vectorized Cp (Elemental check)..."
    V_arr  = [0.0_wp, V_inf, V_inf * 0.5_wp]
    Cp_arr = compute_pressure_coefficient(V_inf, V_local=V_arr)
    
    ! Expected: [1.0, 0.0, 0.75] since Cp = 1 - (V/Vinf)^2
    call assert_approx(Cp_arr(1), 1.0_wp,  "Cp Stagnation (Array)", fail_count)
    call assert_approx(Cp_arr(2), 0.0_wp,  "Cp Freestream (Array)", fail_count)
    call assert_approx(Cp_arr(3), 0.75_wp, "Cp Mid-velocity (Array)", fail_count)

    print *, "------------------------------------------"
    if (fail_count == 0) then
        print *, " ALL BERNOULLI TESTS PASSED "
    else
        print *, " TESTS FAILED: ", fail_count
        call exit(1)
    end if

contains

    subroutine assert_approx(actual, expected, label, failures)
        real(wp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        if (abs(actual - expected) < 1.0e-5_wp) then
            print *, "   [PASS] ", label
        else
            print *, "   [FAIL] ", label
            print *, "          Expected:", expected, " Got:", actual
            failures = failures + 1
        end if
    end subroutine assert_approx

end program test_bernoulli