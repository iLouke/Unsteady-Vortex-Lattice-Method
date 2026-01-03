program test_doublet
    use mod_precision
    use mod_constants
    use mod_doublet
    implicit none

    real(wp) :: p(2), p0(2), mu_val
    real(wp) :: vel(2)
    integer  :: fail_count = 0

    print *, "=========================================="
    print *, "          TEST: MOD_DOUBLET               "
    print *, "=========================================="

    ! --- Test 1: Point Doublet ---
    print *, "[TEST 1] Point Doublet Induction..."
    p0 = [0.0_wp, 0.0_wp]
    p  = [1.0_wp, 0.0_wp] 
    mu_val = TWO_PI          
    
    vel = point_doublet_2D(p, p0, mu_val)
    ! Magnitude should be 1.0 at r=1 for mu=2pi
    call assert_approx(norm2(vel), 1.0_wp, "Doublet Point Magnitude", fail_count)

    ! --- Test 2: Doublet Panel Symmetry ---
    print *, "[TEST 2] 2D Doublet Panel Symmetry..."
    ! Evaluation point: exactly above the center of a panel from -1 to 1
    p = [0.0_wp, 0.5_wp]
    vel = compute_doublet_panel2D_velocity(p, p1=[-1.0_wp, 0.0_wp], &
                                              p2=[ 1.0_wp, 0.0_wp], mu=1.0_wp)
    
    ! At the symmetry plane (x=0), u must be exactly 0.0
    call assert_approx(vel(1), 0.0_wp, "Panel Symmetry (u=0)", fail_count)

    ! --- Test 3: Off-Center induction ---
    ! If u is 0 at center, let's check a point where it shouldn't be 0
    p = [0.5_wp, 0.5_wp]
    vel = compute_doublet_panel2D_velocity(p, p1=[-1.0_wp, 0.0_wp], &
                                              p2=[ 1.0_wp, 0.0_wp], mu=1.0_wp)
    if (abs(vel(1)) > ZERO) then
        print *, "   [PASS] Off-center induction detected"
    else
        print *, "   [FAIL] No induction detected off-center"
        fail_count = fail_count + 1
    end if

    print *, "------------------------------------------"
    if (fail_count /= 0) call exit(1)

contains

    subroutine assert_approx(actual, expected, label, failures)
        real(wp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        if (abs(actual - expected) < 1.0e-5_wp) then
            print *, "   [PASS] ", label
        else
            print *, "   [FAIL] ", label, " Got: ", actual, " Exp: ", expected
            failures = failures + 1
        end if
    end subroutine assert_approx

end program test_doublet