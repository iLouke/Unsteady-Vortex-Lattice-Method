program test_atmosphere
    use mod_precision
    use mod_constants
    use mod_atmosphere
    implicit none

    type(t_atmosphere) :: atm
    integer :: fail_count = 0
    integer :: unit_log, io_stat

    print *, "======================================="
    print *, "      TEST: MOD_ATMOSPHERE (OOP)       "
    print *, "======================================="

    ! -----------------------------------------------------------
    ! TEST 1: Sea Level Constants
    ! -----------------------------------------------------------
    print *, "[TEST 1] Verifying Sea Level (0m)..."
    call atm%update(0.0_wp)

    call assert_approx(atm%temperature,    T_SL,              "Temperature SL",    fail_count)
    call assert_approx(atm%pressure,       P_SL,              "Pressure SL",       fail_count)
    call assert_approx(atm%density,        RHO_SL,            "Density SL",        fail_count)
    call assert_approx(atm%speed_of_sound, SPEED_OF_SOUND_SL, "Speed of Sound SL", fail_count)

    ! -----------------------------------------------------------
    ! TEST 2: High Altitude (10km / Cruise)
    ! -----------------------------------------------------------
    print *, "[TEST 2] Verifying High Altitude (10,000m)..."
    call atm%update(10000.0_wp)

    ! Reference values for 10km ISA: T ~ 223.25K, P ~ 26436 Pa
    if (atm%temperature < T_SL .and. atm%pressure < P_SL) then
        print *, "   [PASS] Correct lapse rate trends detected."
    else
        print *, "   [FAIL] Physics violation at altitude."
        fail_count = fail_count + 1
    end if

    ! -----------------------------------------------------------
    ! TEST 3: Logging to File
    ! -----------------------------------------------------------
    print *, "[TEST 3] Verifying Output Method..."
    open(newunit=unit_log, file='atm_test_output.txt', status='replace', iostat=io_stat)
    if (io_stat == 0) then
        call atm%display(unit_log)
        close(unit_log)
        print *, "   [PASS] Display method wrote to unit successfully."
    else
        print *, "   [FAIL] File I/O Error."
        fail_count = fail_count + 1
    end if

    print *, "---------------------------------------"
    if (fail_count == 0) then
        print *, " ALL ATMOSPHERE TESTS PASSED "
    else
        print *, " TESTS FAILED: ", fail_count
        call exit(1)
    end if

contains

    subroutine assert_approx(actual, expected, label, failures)
        real(wp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        if (abs(actual - expected) < 1.0e-2_wp) then
            print *, "   [PASS] ", label
        else
            print *, "   [FAIL] ", label
            print *, "          Expected:", expected, " Got:", actual
            failures = failures + 1
        end if
    end subroutine assert_approx

end program test_atmosphere