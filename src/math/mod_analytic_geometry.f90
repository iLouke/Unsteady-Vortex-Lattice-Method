! Dot, Cross, Norm, Coordinate Transforms
module analytic_geometry
    use mod_precision
    use mod_constants
    implicit none

    private
    
    public :: cross_product, rotate_2d

    contains

    ! --------------------------------------------------
    ! Cross Product of two 3D vectors
    ! --------------------------------------------------
    function cross_product(v1, v2) result(v_cross)
        real(wp), intent(in) :: v1(3), v2(3)
        real(wp) :: v_cross(3)

        v_cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v_cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v_cross(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end function cross_product

    !--------------------------------------------------
    ! Rotation In Two Dimensions
    !--------------------------------------------------
    subroutine rotate_2d(point_in, angle, point_out, angle_type)
        ! CCW Rotation of a 2D point by a given angle
        real(wp), intent(in)  :: point_in(2)
        real(wp), intent(in)  :: angle
        real(wp), intent(out) :: point_out(2)
        character(len=*), intent(in), optional :: angle_type
        
        ! .. Local Variables ..
        real(wp) :: angle_rad
        real(wp) :: cos_angle, sin_angle

        angle_rad = angle

        if (present(angle_type)) then
            if (angle_type == 'degrees') then
                angle_rad = angle * DEG2RAD
            else if (angle_type == 'radians') then
                angle_rad = angle
            end if
        end if

        cos_angle = cos(angle_rad)
        sin_angle = sin(angle_rad)

        point_out(1) = point_in(1) * cos_angle - point_in(2) * sin_angle
        point_out(2) = point_in(1) * sin_angle + point_in(2) * cos_angle
    end subroutine rotate_2d

    


end module analytic_geometry
