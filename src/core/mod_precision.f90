!> =============================================================================
!> \module   mod_precision
!> \brief    Defines the working precision for real numbers.
!> \details  Contains the definition of the working precision used throughout the code.
!>           
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2025
!> =============================================================================
module mod_precision
    implicit none
    private

    public :: wp

    integer, parameter :: wp = selected_real_kind(15, 307) ! At least 15 decimal digits, exponent up to 10^307

end module mod_precision