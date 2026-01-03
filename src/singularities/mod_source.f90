!> =============================================================================
!> \module   mod_source
!> \brief    2D/3D Source Singularities (sigma)
!> \details  Implements point sources, 2D constant-strength panels, and 3D 
!>           quadrilateral source panels (Hess & Smith method).
!> 
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2026
!> =============================================================================
module mod_source
    use mod_precision
    use mod_constants
    use mod_math
    implicit none
    private

    ! Public Interface
    public :: point_source_2D
    public :: compute_source_2D_panel_velocity
    public :: compute_source_quad_velocity

contains

    !> =========================================================================
    !> \brief   Induced velocity from a 2D Point Source.
    !> \param   p      Field Point (x,y)
    !> \param   p0     Source Location (x0,y0)
    !> \param   sigma  Source Strength
    !> \return  velocity Vector (u,v)
    !> =========================================================================
    pure function point_source_2D(p, p0, sigma) result(velocity)
        real(wp), intent(in) :: p(2), p0(2)
        real(wp), intent(in) :: sigma
        real(wp) :: velocity(2)
        
        real(wp) :: r2, r_vec(2)

        r_vec = p - p0
        r2 = dot_product(r_vec, r_vec)

        if (r2 < ZERO) then
            velocity = 0.0_wp
        else
            velocity = (sigma / (TWO_PI * r2)) * r_vec
        end if
    end function point_source_2D

    !> =========================================================================
    !> \brief   Constant Strength 2D Source Panel.
    !> \details Based on Katz & Plotkin, Eqs 10.17 - 10.21.
    !> \param   p      Field point (x, z)
    !> \param   p1, p2 Panel nodes
    !> \param   sigma  Strength
    !> =========================================================================
    pure function compute_source_2D_panel_velocity(p, p1, p2, sigma) result(vel_gl)
        real(wp), intent(in) :: p(2), p1(2), p2(2)
        real(wp), intent(in) :: sigma
        real(wp) :: vel_gl(2)

        real(wp) :: p_vec(2), tang(2), norm(2)
        real(wp) :: L, d(2), loc(2)
        real(wp) :: u_loc, w_loc, r1_sq, r2_sq, th1, th2

        ! 1. Geometry
        p_vec = p2 - p1
        L     = norm2(p_vec)
        
        if (L < ZERO) then
            vel_gl = 0.0_wp
            return
        end if

        tang = p_vec / L
        norm = [-tang(2), tang(1)]

        ! 2. Transformation
        d = p - p1
        loc = [dot_product(d, tang), dot_product(d, norm)]
        r1_sq = loc(1)**2 + loc(2)**2
        r2_sq = (loc(1) - L)**2 + loc(2)**2

        ! 3. Local Velocity (u=Axial, w=Normal)
        if (r1_sq < ZERO .or. r2_sq < ZERO) then
            u_loc = 0.0_wp
        else
            u_loc = (sigma / (4.0_wp * PI)) * log(r1_sq / r2_sq)
        end if

        ! Surface singularity check
        if (abs(loc(2)) < ZERO .and. (loc(1) > 0.0_wp .and. loc(1) < L)) then
            w_loc = merge(0.5_wp*sigma, -0.5_wp*sigma, loc(2) >= 0.0_wp)
        else
            th1 = atan2(loc(2), loc(1))
            th2 = atan2(loc(2), loc(1) - L)
            w_loc = (sigma / TWO_PI) * (th2 - th1)
        end if

        ! 4. Local -> Global
        vel_gl = u_loc * tang + w_loc * norm
    end function compute_source_2D_panel_velocity

    !> =========================================================================
    !> \brief   3D Constant Strength Source Quadrilateral (Hess & Smith).
    !> \details Includes Far-Field point source approximation for efficiency.
    !> \param   p       Field point (3)
    !> \param   corners Quad nodes (4,3)
    !> \param   sigma   Strength
    !> =========================================================================
    function compute_source_quad_velocity(p, corners, sigma) result(vel_gl)
        real(wp), intent(in) :: p(3), corners(4,3), sigma
        real(wp) :: vel_gl(3)

        real(wp) :: center(3), normal(3), tx(3), ty(3), r_vec(3)
        real(wp) :: loc_p(3), loc_c(4,2), d_s(4), r_s(4)
        real(wp) :: area, dist_sq, u_l, v_l, w_l, term
        integer  :: k, kn
        real(wp), parameter :: FAR_FIELD_LIMIT = 5.0_wp

        ! 1. Geometry & Far-Field Check
        center = sum(corners, dim=1) * 0.25_wp
        r_vec  = p - center
        dist_sq = dot_product(r_vec, r_vec)
        
        ! Calculate Area and Normal
        r_vec = cross_product(corners(3,:) - corners(1,:), corners(4,:) - corners(2,:))
        area  = 0.5_wp * norm2(r_vec)
        if (area < ZERO) return

        if (dist_sq > (FAR_FIELD_LIMIT * sqrt(area))**2) then
            ! Point Source Approximation
            vel_gl = (sigma * area / (4.0_wp * PI * dist_sq * sqrt(dist_sq))) * (p - center)
            return
        end if

        ! 2. Coordinate System
        normal = r_vec / norm2(r_vec)
        tx     = unit_vector(corners(2,:) - corners(1,:))
        ty     = cross_product(normal, tx)

        ! 3. Local Transformation
        call project_to_panel(p, corners(1,:), tx, ty, normal, loc_p)
        do k = 1, 4
            call project_to_panel(corners(k,:), corners(1,:), tx, ty, normal, r_vec)
            loc_c(k,:) = r_vec(1:2)
            r_s(k) = sqrt((loc_p(1)-loc_c(k,1))**2 + (loc_p(2)-loc_c(k,2))**2 + loc_p(3)**2)
        end do

        ! 4. Hess & Smith Components
        u_l = 0.0_wp; v_l = 0.0_wp; w_l = 0.0_wp
        do k = 1, 4
            kn = merge(1, k + 1, k == 4)
            d_s(k) = sqrt((loc_c(kn,1)-loc_c(k,1))**2 + (loc_c(kn,2)-loc_c(k,2))**2)
            
            ! Induced velocity kernels
            term = log(abs((r_s(k) + r_s(kn) + d_s(k)) / (r_s(k) + r_s(kn) - d_s(k) + ZERO)))
            u_l = u_l + ((loc_c(kn,2) - loc_c(k,2)) / d_s(k)) * term
            v_l = v_l + ((loc_c(k,1) - loc_c(kn,1)) / d_s(k)) * term
            
            ! Solid angle component for w (Simplified robust implementation)
            if (abs(loc_p(3)) > ZERO) then
                w_l = w_l + atan2(loc_p(3)*d_s(k)*(r_s(k)+r_s(kn)), &
                                  r_s(k)*r_s(kn)*d_s(k)**2) ! Placeholder for full H&S atan
            end if
        end do

        vel_gl = (u_l * tx + v_l * ty + w_l * normal) * (sigma / (4.0_wp * PI))
    end function compute_source_quad_velocity

    subroutine project_to_panel(p_glob, orig, tx, ty, tz, p_loc)
        real(wp), intent(in) :: p_glob(3), orig(3), tx(3), ty(3), tz(3)
        real(wp), intent(out) :: p_loc(3)
        real(wp) :: tmp(3)
        tmp = p_glob - orig
        p_loc = [dot_product(tmp, tx), dot_product(tmp, ty), dot_product(tmp, tz)]
    end subroutine project_to_panel

end module mod_source