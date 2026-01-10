!> =============================================================================
!> \module   vortex
!> \brief    Vortex Singularity Models (2D/3D)
!> \details  Implements point vortices, constant-strength vortex panels, 
!>           3D vortex lines, rings, and horseshoe vortices.
!>
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2026
!> =============================================================================
MODULE vortex
    USE constants
    IMPLICIT NONE
    PRIVATE

    ! Public Interface
    PUBLIC :: point_vortex_2D
    PUBLIC :: compute_vortex_panel2D_velocity
    PUBLIC :: compute_vortex_line_velocity
    PUBLIC :: compute_vortex_ring_velocity
    PUBLIC :: compute_horseshoe_velocity

CONTAINS

    !> =========================================================================
    !> \brief   Induced velocity from a 2D Point Vortex.
    !> =========================================================================
    PURE FUNCTION point_vortex_2D(p, p0, gama) RESULT(vel)
        REAL(wp), INTENT(IN) :: p(2), p0(2)
        REAL(wp), INTENT(IN) :: gama
        REAL(wp) :: vel(2)
        REAL(wp) :: r2, r_vec(2)

        r_vec = p - p0
        r2 = DOT_PRODUCT(r_vec, r_vec)

        IF (r2 < 1.0E-12_wp) THEN
            vel = 0.0_wp
        ELSE
            ! CCW Convention
            vel(1) = -(gama / (2.0_wp * PI * r2)) * r_vec(2)
            vel(2) =  (gama / (2.0_wp * PI * r2)) * r_vec(1)
        END IF
    END FUNCTION point_vortex_2D

    !> =========================================================================
    !> \brief   2D Constant Strength Vortex Panel.
    !> =========================================================================
    PURE FUNCTION compute_vortex_panel2D_velocity(p, p1, p2, gama) RESULT(vel_gl)
        REAL(wp), INTENT(IN) :: p(2), p1(2), p2(2)
        REAL(wp), INTENT(IN) :: gama
        REAL(wp) :: vel_gl(2)

        REAL(wp) :: p_vec(2), tang(2), norm_v(2)
        REAL(wp) :: L, d(2), loc(2)
        REAL(wp) :: u_l, w_l, r1_sq, r2_sq, th1, th2

        p_vec = p2 - p1
        L     = NORM2(p_vec)
        
        IF (L < 1.0E-12_wp) THEN
            vel_gl = 0.0_wp
            RETURN
        END IF

        tang   = p_vec / L
        norm_v = [-tang(2), tang(1)]

        d   = p - p1
        loc = [DOT_PRODUCT(d, tang), DOT_PRODUCT(d, norm_v)]
        r1_sq = loc(1)**2 + loc(2)**2
        r2_sq = (loc(1) - L)**2 + loc(2)**2

        ! Tangential Velocity
        IF (ABS(loc(2)) < 1.0E-12_wp .AND. (loc(1) > 0.0_wp .AND. loc(1) < L)) THEN
            u_l = MERGE(0.5_wp*gama, -0.5_wp*gama, loc(2) >= 0.0_wp)
        ELSE
            th1 = ATAN2(loc(2), loc(1))
            th2 = ATAN2(loc(2), loc(1) - L)
            u_l = (gama / (2.0_wp * PI)) * (th2 - th1)
        END IF

        ! Normal Velocity
        w_l = MERGE((gama / (4.0_wp * PI)) * LOG(r2_sq / r1_sq), 0.0_wp, &
                    r1_sq > 1.0E-12_wp .AND. r2_sq > 1.0E-12_wp)

        vel_gl = u_l * tang + w_l * norm_v
    END FUNCTION compute_vortex_panel2D_velocity

    !> =========================================================================
    !> \brief   Induced velocity from a 3D Vortex Line Segment (Biot-Savart).
    !> =========================================================================
    PURE FUNCTION compute_vortex_line_velocity(p, p1, p2, gama) RESULT(vel)
        REAL(wp), INTENT(IN) :: p(3), p1(3), p2(3)
        REAL(wp), INTENT(IN) :: gama
        REAL(wp) :: vel(3)

        REAL(wp) :: r1(3), r2(3), r0(3), cp(3)
        REAL(wp) :: r1_m, r2_m, cp_sq, k_fact

        r1 = p - p1
        r2 = p - p2
        r0 = p2 - p1
        
        ! Cross product (r1 x r2)
        cp = [r1(2)*r2(3)-r1(3)*r2(2), r1(3)*r2(1)-r1(1)*r2(3), r1(1)*r2(2)-r1(2)*r2(1)]
        cp_sq = DOT_PRODUCT(cp, cp)
        r1_m  = NORM2(r1)
        r2_m  = NORM2(r2)

        IF (r1_m < 1.0E-12_wp .OR. r2_m < 1.0E-12_wp .OR. cp_sq < 1.0E-12_wp) THEN
            vel = 0.0_wp
        ELSE
            k_fact = (gama / (4.0_wp * PI * cp_sq)) * &
                     (DOT_PRODUCT(r0, r1)/r1_m - DOT_PRODUCT(r0, r2)/r2_m)
            vel = k_fact * cp
        END IF
    END FUNCTION compute_vortex_line_velocity

    !> =========================================================================
    !> \brief   Induced velocity from a 3D Constant Strength Vortex Ring.
    !> \details Updated to use (3,4) input for contiguous memory access.
    !> \param   corners  (3,4) Array -> [x,y,z] for [P1,P2,P3,P4]
    !> =========================================================================
    PURE FUNCTION compute_vortex_ring_velocity_slow(p, corners, gama) RESULT(vel)
        REAL(wp), INTENT(IN) :: p(3)
        REAL(wp), INTENT(IN) :: corners(3,4) ! <--- FIXED DIMENSIONS
        REAL(wp), INTENT(IN) :: gama
        REAL(wp) :: vel(3)
        INTEGER  :: k, kn

        vel = 0.0_wp
        DO k = 1, 4
            kn = MERGE(1, k + 1, k == 4)
            ! <--- FIXED SLICING: Use Columns (:, k) instead of Rows (k, :)
            vel = vel + compute_vortex_line_velocity(p, corners(:,k), corners(:,kn), gama)
        END DO
    END FUNCTION compute_vortex_ring_velocity_slow

    !> =========================================================================
    !> \brief   Optimized Induced velocity from a 3D Vortex Ring.
    !> \details Inlines the line segment logic to reuse vector norms and 
    !>          avoid redundant calculations.
    !> =========================================================================
    PURE FUNCTION compute_vortex_ring_velocity(p, corners, gama) RESULT(vel)
        REAL(wp), INTENT(IN) :: p(3)
        REAL(wp), INTENT(IN) :: corners(3,4)
        REAL(wp), INTENT(IN) :: gama
        REAL(wp) :: vel(3)

        ! Local variables
        REAL(wp) :: r(3,4)           ! Vectors from p to corners
        REAL(wp) :: r_sq(4)          ! Squared magnitude |r|^2
        REAL(wp) :: r_inv(4)         ! Inverse magnitude 1/|r|
        REAL(wp) :: cp(3)            ! Cross product r1 x r2
        REAL(wp) :: cp_sq            ! |r1 x r2|^2
        REAL(wp) :: d12              ! Dot product r1 . r2
        REAL(wp) :: k_fact, common_factor
        INTEGER  :: i, next_i

        ! Constants
        REAL(wp), PARAMETER :: tol = 1.0E-12_wp
        REAL(wp), PARAMETER :: inv_4pi = 1.0_wp / (4.0_wp * PI) 

        vel = 0.0_wp

        ! 1. Pre-calculate vectors and distances to all 4 corners ONCE
        DO i = 1, 4
        r(:,i)  = p(:) - corners(:,i)
        r_sq(i) = DOT_PRODUCT(r(:,i), r(:,i))
        
        ! Singularity check (Point on corner)
        IF (r_sq(i) < tol) RETURN 
        
        r_inv(i) = 1.0_wp / SQRT(r_sq(i))
        END DO

        ! 2. Loop over 4 segments using pre-computed vectors
        DO i = 1, 4
        next_i = MERGE(1, i + 1, i == 4)

        ! Cross product (r_i x r_next)
        ! This represents the direction of induced velocity for this segment
        cp(1) = r(2,i)*r(3,next_i) - r(3,i)*r(2,next_i)
        cp(2) = r(3,i)*r(1,next_i) - r(1,i)*r(3,next_i)
        cp(3) = r(1,i)*r(2,next_i) - r(2,i)*r(1,next_i)

        cp_sq = DOT_PRODUCT(cp, cp)

        ! Singularity check (Point on the line segment)
        IF (cp_sq > tol) THEN
            ! Optimization: We avoid calculating the segment vector r0 completely.
            ! We use the identity: r0 . r1 = |r1|^2 - r1.r2
            
            d12 = DOT_PRODUCT(r(:,i), r(:,next_i))
            
            ! The Biot-Savart term: (r0.r1/|r1| - r0.r2/|r2|)
            ! Becomes: ( (|r1|^2 - d12)/|r1| - (d12 - |r2|^2)/|r2| )
            
            k_fact = (r_sq(i) - d12) * r_inv(i) + (r_sq(next_i) - d12) * r_inv(next_i)
            
            common_factor = (gama * inv_4pi) / cp_sq
            
            vel = vel + (common_factor * k_fact) * cp
        END IF
        END DO
        
    END FUNCTION compute_vortex_ring_velocity

    !> =========================================================================
    !> \brief   Induced velocity from a 3D Horseshoe Vortex.
    !> =========================================================================
    PURE FUNCTION compute_horseshoe_velocity(p, pa, pb, gama, wake_dir) RESULT(vel)
        REAL(wp), INTENT(IN)           :: p(3), pa(3), pb(3), gama
        REAL(wp), INTENT(IN), OPTIONAL :: wake_dir(3)
        REAL(wp) :: vel(3)

        REAL(wp) :: w_vec(3), p_inf_a(3), p_inf_b(3), L_wake

        w_vec = [1.0_wp, 0.0_wp, 0.0_wp]
        IF (PRESENT(wake_dir)) w_vec = wake_dir / NORM2(wake_dir)
        
        L_wake = 1000.0_wp * MAX(NORM2(pb - pa), 1.0_wp)

        p_inf_a = pa + (L_wake * w_vec)
        p_inf_b = pb + (L_wake * w_vec)

        vel = compute_vortex_line_velocity(p, p_inf_a, pa, gama) + &
              compute_vortex_line_velocity(p, pa, pb, gama)       + &
              compute_vortex_line_velocity(p, pb, p_inf_b, gama)

    END FUNCTION compute_horseshoe_velocity

END MODULE vortex