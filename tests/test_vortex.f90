program test_vortex
    use constants
    use vortex
    implicit none

    real(wp) :: p_f(3), p1(3), p2(3), circ
    real(wp) :: vel3d(3)
    integer  :: fail_count = 0

    print *, "=========================================="
    print *, "          TEST: MOD_VORTEX                "
    print *, "=========================================="

    print *, "[TEST 1] 3D Vortex Line Magnitude..."
    p1 = [0.0_wp, -1.0_wp, 0.0_wp]
    p2 = [0.0_wp,  1.0_wp, 0.0_wp]
    p_f = [1.0_wp,  0.0_wp, 0.0_wp]
    circ = 4.0_wp * PI
    
    vel3d = compute_vortex_line_velocity(p_f, p1, p2, circ)
    
    ! The magnitude at this point should be exactly sqrt(2)
    call assert_approx(norm2(vel3d), sqrt(2.0_wp), "Vortex Line Magnitude", fail_count)

    print *, "[TEST 2] Horseshoe Vortex (Far Field)..."
    ! Evaluation point far away should have near-zero velocity
    vel3d = compute_horseshoe_velocity([100.0_wp, 0.0_wp, 0.0_wp], &
                                       [0.0_wp, -1.0_wp, 0.0_wp], &
                                       [0.0_wp, 1.0_wp, 0.0_wp], 1.0_wp)
    
    if (norm2(vel3d) < 1.0_wp) then
         print *, "   [PASS] Horseshoe Far-Field decay"
    else
         print *, "   [FAIL] Horseshoe unexpected magnitude"
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

end program test_vortex