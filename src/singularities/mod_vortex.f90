!> =============================================================================
!> \module   mod_vortex
!> \brief    Vortex Singularity Models (2D/3D)
!> \details  Implements point vortices, constant-strength vortex panels, 
!>           3D vortex lines, rings, and horseshoe vortices.
!>
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2026
!> =============================================================================
module mod_vortex
    use mod_precision
    use mod_constants
    implicit none
    private

    ! Public Interface
    public :: point_vortex_2D
    public :: compute_vortex_panel2D_velocity
    public :: compute_vortex_line_velocity
    public :: compute_vortex_ring_velocity
    public :: compute_horseshoe_velocity

contains

    !> =========================================================================
    !> \brief   Induced velocity from a 2D Point Vortex.
    !> \param   p      Field point (x, y)
    !> \param   p0     Vortex location (x0, y0)
    !> \param   gamma  Circulation strength (CCW positive)
    !> \return  vel    Induced velocity vector (u, v)
    !> =========================================================================
    pure function point_vortex_2D(p, p0, gamma) result(vel)
        real(wp), intent(in) :: p(2), p0(2)
        real(wp), intent(in) :: gamma
        real(wp) :: vel(2)

        real(wp) :: r2, r_vec(2)

        r_vec = p - p0
        r2 = dot_product(r_vec, r_vec)

        if (r2 < ZERO) then
            vel = 0.0_wp
        else
            ! CCW Convention: u = -Gamma/2pi * dy/r^2, v = Gamma/2pi * dx/r^2
            vel(1) = -(gamma / (TWO_PI * r2)) * r_vec(2)
            vel(2) =  (gamma / (TWO_PI * r2)) * r_vec(1)
        end if
    end function point_vortex_2D

    !> =========================================================================
    !> \brief   2D Constant Strength Vortex Panel.
    !> \details Based on Katz & Plotkin, Eqs 10.39 - 10.40.
    !> \param   p      Field point (x, z)
    !> \param   p1, p2 Panel nodes
    !> \param   gamma  Circulation per unit length
    !> =========================================================================
    pure function compute_vortex_panel2D_velocity(p, p1, p2, gamma) result(vel_gl)
        real(wp), intent(in) :: p(2), p1(2), p2(2)
        real(wp), intent(in) :: gamma
        real(wp) :: vel_gl(2)

        real(wp) :: p_vec(2), tang(2), norm_v(2)
        real(wp) :: L, d(2), loc(2)
        real(wp) :: u_l, w_l, r1_sq, r2_sq, th1, th2

        p_vec = p2 - p1
        L     = norm2(p_vec)
        
        if (L < ZERO) then
            vel_gl = 0.0_wp
            return
        end if

        tang   = p_vec / L
        norm_v = [-tang(2), tang(1)]

        d   = p - p1
        loc = [dot_product(d, tang), dot_product(d, norm_v)]
        r1_sq = loc(1)**2 + loc(2)**2
        r2_sq = (loc(1) - L)**2 + loc(2)**2

        ! Tangential Velocity (Eq 10.39)
        if (abs(loc(2)) < ZERO .and. (loc(1) > 0.0_wp .and. loc(1) < L)) then
            u_l = merge(0.5_wp*gamma, -0.5_wp*gamma, loc(2) >= 0.0_wp)
        else
            th1 = atan2(loc(2), loc(1))
            th2 = atan2(loc(2), loc(1) - L)
            u_l = (gamma / TWO_PI) * (th2 - th1)
        end if

        ! Normal Velocity (Eq 10.40)
        w_l = merge((gamma / (4.0_wp * PI)) * log(r2_sq / r1_sq), 0.0_wp, &
                    r1_sq > ZERO .and. r2_sq > ZERO)

        vel_gl = u_l * tang + w_l * norm_v
    end function compute_vortex_panel2D_velocity

    !> =========================================================================
    !> \brief   Induced velocity from a 3D Vortex Line Segment.
    !> \details Based on Katz & Plotkin, Eq 10.115 (Biot-Savart Law).
    !> \param   p      Field point (x, y, z)
    !> \param   p1, p2 Segment start/end nodes
    !> \param   gamma  Circulation strength
    !> =========================================================================
    pure function compute_vortex_line_velocity(p, p1, p2, gamma) result(vel)
        real(wp), intent(in) :: p(3), p1(3), p2(3)
        real(wp), intent(in) :: gamma
        real(wp) :: vel(3)

        real(wp) :: r1(3), r2(3), r0(3), cp(3)
        real(wp) :: r1_m, r2_m, cp_sq, k_fact

        r1 = p - p1
        r2 = p - p2
        r0 = p2 - p1
        
        ! Cross product of diagonals for direction
        cp = [r1(2)*r2(3)-r1(3)*r2(2), r1(3)*r2(1)-r1(1)*r2(3), r1(1)*r2(2)-r1(2)*r2(1)]
        cp_sq = dot_product(cp, cp)
        r1_m  = norm2(r1)
        r2_m  = norm2(r2)

        if (r1_m < ZERO .or. r2_m < ZERO .or. cp_sq < ZERO) then
            vel = 0.0_wp
        else
            ! Biot-Savart Kernel
            k_fact = (gamma / (4.0_wp * PI * cp_sq)) * &
                     (dot_product(r0, r1)/r1_m - dot_product(r0, r2)/r2_m)
            vel = k_fact * cp
        end if
    end function compute_vortex_line_velocity

    !> =========================================================================
    !> \brief   Induced velocity from a 3D Constant Strength Vortex Ring.
    !> \details Superposition of 4 vortex line segments.
    !> \param   p        Field point
    !> \param   corners  4 corners of the quad (CCW)
    !> \param   gamma    Circulation
    !> =========================================================================
    pure function compute_vortex_ring_velocity(p, corners, gamma) result(vel)
        real(wp), intent(in) :: p(3), corners(4,3), gamma
        real(wp) :: vel(3)
        integer  :: k, kn

        vel = 0.0_wp
        do k = 1, 4
            kn = merge(1, k + 1, k == 4)
            vel = vel + compute_vortex_line_velocity(p, corners(k,:), corners(kn,:), gamma)
        end do
    end function compute_vortex_ring_velocity

    !> =========================================================================
    !> \brief   Induced velocity from a 3D Horseshoe Vortex.
    !> \details Consists of a bound segment and two semi-infinite trailing legs.
    !> \param   p        Field point
    !> \param   pa, pb   Nodes of the bound vortex (Left to Right)
    !> \param   gamma    Circulation
    !> \param   wake_dir (Optional) Downstream direction
    !> =========================================================================
    pure function compute_horseshoe_velocity(p, pa, pb, gamma, wake_dir) result(vel)
        real(wp), intent(in)           :: p(3), pa(3), pb(3), gamma
        real(wp), intent(in), optional :: wake_dir(3)
        real(wp) :: vel(3)

        real(wp) :: w_vec(3), p_inf_a(3), p_inf_b(3), L_wake

        ! 1. Define wake direction and length (numerical infinity)
        w_vec = [1.0_wp, 0.0_wp, 0.0_wp]
        if (present(wake_dir)) w_vec = wake_dir / norm2(wake_dir)
        
        L_wake = 1000.0_wp * max(norm2(pb - pa), 1.0_wp)

        ! 2. Define semi-infinite endpoints
        p_inf_a = pa + (L_wake * w_vec)
        p_inf_b = pb + (L_wake * w_vec)

        ! 3. Sum segments: Leg A (Inf -> A) + Bound (A -> B) + Leg B (B -> Inf)
        vel = compute_vortex_line_velocity(p, p_inf_a, pa, gamma) + &
              compute_vortex_line_velocity(p, pa, pb, gamma)       + &
              compute_vortex_line_velocity(p, pb, p_inf_b, gamma)

    end function compute_horseshoe_velocity

end module mod_vortex