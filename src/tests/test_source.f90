program test_source
    use mod_precision
    use mod_constants
    use mod_source
    implicit none

    real(wp) :: p(2), p1(2), p2(2), sigma
    real(wp) :: vel(2)
    integer  :: fail_count = 0

    print *, "=========================================="
    print *, "          TEST: MOD_SOURCE                "
    print *, "=========================================="

    ! --- Test 1: Point Source Radial Flow ---
    print *, "[TEST 1] Point Source Induction..."
    sigma = TWO_PI
    vel = point_source_2D(p=[1.0_wp, 0.0_wp], p0=[0.0_wp, 0.0_wp], sigma=sigma)
    
    ! V = (sigma/2pi*r) * r_hat. For r=1, V=1.0
    call assert_approx(vel(1), 1.0_wp, "Point Source u-velocity", fail_count)

    ! --- Test 2: 2D Panel Jump Condition ---
    print *, "[TEST 2] Source Panel Jump Condition..."
    p1 = [-1.0_wp, 0.0_wp]
    p2 = [ 1.0_wp, 0.0_wp]
    sigma = 1.0_wp
    
    ! Evaluate just above the center of the panel
    vel = compute_source_2D_panel_velocity(p=[0.0_wp, 0.0001_wp], p1=p1, p2=p2, sigma=sigma)
    
    ! Normal velocity w should be exactly sigma/2 = 0.5
    call assert_approx(vel(2), 0.5_wp, "Source Panel Normal Jump (w=0.5)", fail_count)

    print *, "------------------------------------------"
    if (fail_count == 0) print *, " ALL SOURCE TESTS PASSED "
    if (fail_count /= 0) call exit(1)

contains
    subroutine assert_approx(actual, expected, label, failures)
        real(wp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        if (abs(actual - expected) < 1.0e-4_wp) then
            print *, "   [PASS] ", label
        else
            print *, "   [FAIL] ", label
            print *, "          Expected:", expected, " Got:", actual
            failures = failures + 1
        end if
    end subroutine assert_approx
end program test_source