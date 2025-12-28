! Defines wp (real64) constants
module mod_precision
    implicit none
    private

    public :: wp

    integer, parameter :: wp = selected_real_kind(15, 307) ! At least 15 decimal digits, exponent up to 10^307

end module mod_precision