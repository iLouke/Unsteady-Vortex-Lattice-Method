program test_atmosphere
    use mod_precision
    use mod_constants
    use mod_atmosphere
    implicit none

    ! Variables
    real(wp) :: alt, temp, p, rho, mu, sos
    real(wp), parameter :: TOL = 1.0e-2_wp ! Tolerance for comparisons
    integer :: unit_log
    integer :: fail_count = 0  ! Track number of failures

    print *, "======================================="
    print *, "    RUNNING ATMOSPHERE MODULE TESTS    "
    print *, "======================================="

    ! --- Test Case 1: Sea Level (0 m) ---
    print *, " "
    print *, "[TEST 1] Verifying Sea Level Conditions..."
    
    alt = 0.0_wp
    call compute_isa_properties(alt, temp, p, rho, mu, sos)

    ! Use the helper subroutine to check values
    ! Note: Ensure T_SL, P_SL, etc. match the names in your mod_constants
    call assert_close(temp, T_SL             , "Temperature"   , fail_count)
    call assert_close(p   , P_SL             , "Pressure"      , fail_count)
    call assert_close(rho , RHO_SL           , "Density"       , fail_count)
    call assert_close(mu  , MU_SL            , "Viscosity"     , fail_count)
    call assert_close(sos , SPEED_OF_SOUND_SL, "Speed of Sound", fail_count)


    ! --- Test Case 2: Logging ---
    print *, " "
    print *, "[TEST 2] Verifying File I/O..."
    
    open(newunit=unit_log, file='test_log.txt', status='replace', iostat=fail_count)
    
    if (fail_count == 0) then
        call compute_isa_properties(1000.0_wp, temp, p, rho, mu, sos, log_unit=unit_log)
        close(unit_log)
        print *, "   [PASS] Log file created successfully."
    else
        print *, "   [FAIL] Could not open log file."
        fail_count = fail_count + 1
    end if

    ! --- FINAL SUMMARY ---
    print *, "---------------------------------------"
    if (fail_count == 0) then
        print *, " ALL TESTS PASSED."
    else
        print *, " TESTS FAILED: ", fail_count
        ! Exit with error code 1 so CTest marks this as "Failed"
        call exit(1) 
    end if

! ------------------------------------------------------------------
! Internal Subroutines (The "Contains" block)
! ------------------------------------------------------------------
contains

    ! Helper to assert that two reals are close enough
    subroutine assert_close(actual, expected, label, failures)
        real(wp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        integer, intent(inout) :: failures
        
        real(wp) :: diff

        diff = abs(actual - expected)

        if (diff < TOL) then
            ! Optional: Use generic formatting (g0) to make it cleaner
            print *, "   [PASS] ", label
        else
            print *, "   [FAIL] ", label
            print *, "          Expected: ", expected
            print *, "          Got:      ", actual
            print *, "          Diff:     ", diff
            failures = failures + 1
        end if
    end subroutine assert_close

end program test_atmosphere