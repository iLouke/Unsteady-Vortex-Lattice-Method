module mod_math
    use mod_precision
    use mod_constants
    implicit none

    private
    
    public :: cross_product, vector_norm, unit_vector, rotate_2d, is_in_bounds, linspace
    public :: interpolation, extrapolation
    
    ! Overload evaluate_1d to handle both scalar and array inputs
    public :: evaluate_1d
    interface evaluate_1d
        module procedure evaluate_1d_scalar
        module procedure evaluate_1d_array
    end interface

contains

    !> =========================================================================
    !> \brief   Scalar version of 1D Lookup
    !> =========================================================================
    pure function evaluate_1d_scalar(x, x_data, y_data, method) result(y)
        real(wp), intent(in) :: x
        real(wp), intent(in) :: x_data(:), y_data(:)
        character(len=*), intent(in), optional :: method
        real(wp) :: y

        if (is_in_bounds(x, x_data)) then
            y = interpolation(x, x_data, y_data, method)
        else
            y = extrapolation(x, x_data, y_data, method)
        end if
    end function evaluate_1d_scalar

    !> =========================================================================
    !> \brief   Array version of 1D Lookup
    !> \details This allows you to pass an entire array of X values.
    !> =========================================================================
    pure function evaluate_1d_array(x_arr, x_data, y_data, method) result(y_arr)
        real(wp), intent(in) :: x_arr(:)
        real(wp), intent(in) :: x_data(:), y_data(:)
        character(len=*), intent(in), optional :: method
        real(wp) :: y_arr(size(x_arr))
        integer :: i

        do i = 1, size(x_arr)
            y_arr(i) = evaluate_1d_scalar(x_arr(i), x_data, y_data, method)
        end do
    end function evaluate_1d_array

    !> =========================================================================
    !> \brief   1D Interpolation (Pure, not Elemental)
    !> =========================================================================
    pure function interpolation(x, x_data, y_data, method) result(y)
        real(wp), intent(in) :: x
        real(wp), intent(in) :: x_data(:), y_data(:)
        character(len=*), intent(in), optional :: method
        real(wp) :: y, term, slope
        integer :: n, i, j
        character(len=16) :: method_local

        method_local = 'linear' 
        if (present(method)) method_local = trim(adjustl(method))
        n = size(x_data)
        y = 0.0_wp

        select case (trim(method_local))
        case ('linear')
            if (abs(x - x_data(n)) < ZERO) then
                y = y_data(n)
                return
            end if

            do i = 1, n-1
                if (x >= x_data(i) .and. x <= x_data(i+1)) then
                    slope = (y_data(i+1) - y_data(i)) / (x_data(i+1) - x_data(i))
                    y = y_data(i) + slope * (x - x_data(i))
                    return
                end if
            end do
            
            y = merge(y_data(1), y_data(n), x < x_data(1))

        case ('lagrange')
            do i = 1, n
                term = y_data(i)
                do j = 1, n
                    if (i /= j) term = term * (x - x_data(j)) / (x_data(i) - x_data(j))
                end do
                y = y + term
            end do
        end select
    end function interpolation

    !> =========================================================================
    !> \brief   1D Extrapolation (Pure, not Elemental)
    !> =========================================================================
    pure function extrapolation(x, x_data, y_data, method) result(y)
        real(wp), intent(in) :: x
        real(wp), intent(in) :: x_data(:), y_data(:)
        character(len=*), intent(in), optional :: method
        real(wp) :: y, slope, term
        integer :: n, i, j
        character(len=16) :: method_local

        method_local = 'linear' 
        if (present(method)) method_local = trim(adjustl(method))
        n = size(x_data)

        select case (trim(method_local))
        case ('linear')
            if (x < x_data(1)) then
                slope = (y_data(2) - y_data(1)) / (x_data(2) - x_data(1))
                y = y_data(1) + slope * (x - x_data(1))
            else
                slope = (y_data(n) - y_data(n-1)) / (x_data(n) - x_data(n-1))
                y = y_data(n) + slope * (x - x_data(n))
            end if
        case ('lagrange')
            y = 0.0_wp
            do i = 1, n
                term = y_data(i)
                do j = 1, n
                    if (i /= j) term = term * (x - x_data(j)) / (x_data(i) - x_data(j))
                end do
                y = y + term
            end do
        case default
            y = 0.0_wp
        end select
    end function extrapolation

    pure function is_in_bounds(x, x_data) result(inside)
        real(wp), intent(in) :: x
        real(wp), intent(in) :: x_data(:)
        logical :: inside
        inside = (x >= x_data(1) .and. x <= x_data(size(x_data)))
    end function is_in_bounds

    ! --- Vector Operations ---
    pure function cross_product(v1, v2) result(v_cross)
        real(wp), intent(in) :: v1(3), v2(3)
        real(wp) :: v_cross(3)
        v_cross = [v1(2)*v2(3)-v1(3)*v2(2), v1(3)*v2(1)-v1(1)*v2(3), v1(1)*v2(2)-v1(2)*v2(1)]
    end function cross_product

    pure function vector_norm(v) result(mag)
        real(wp), intent(in) :: v(3)
        real(wp) :: mag
        mag = sqrt(dot_product(v, v))
    end function vector_norm

    pure function unit_vector(v) result(u_vec)
        real(wp), intent(in) :: v(3)
        real(wp) :: u_vec(3)
        real(wp) :: mag
        mag = vector_norm(v)
        u_vec = merge(v / mag, [0.0_wp, 0.0_wp, 0.0_wp], mag > ZERO)
    end function unit_vector

    pure subroutine rotate_2d(point_in, angle, point_out, angle_type)
        real(wp), intent(in) :: point_in(2), angle
        real(wp), intent(out) :: point_out(2)
        character(len=*), intent(in), optional :: angle_type
        real(wp) :: ang, c, s
        ang = merge(angle * DEG2RAD, angle, present(angle_type) .and. trim(adjustl(angle_type)) == 'degrees')
        c = cos(ang); s = sin(ang)
        point_out = [point_in(1)*c - point_in(2)*s, point_in(1)*s + point_in(2)*c]
    end subroutine rotate_2d

    function linspace(start_val, end_val, n) result(arr)
        real(wp), intent(in) :: start_val, end_val
        integer, intent(in) :: n
        real(wp), allocatable :: arr(:)
        integer :: i
        allocate(arr(n))
        if (n == 1) then
            arr(1) = start_val
        else
            do i = 1, n
                arr(i) = start_val + (i-1) * (end_val - start_val) / (n - 1)
            end do
        end if
    end function linspace

end module mod_math