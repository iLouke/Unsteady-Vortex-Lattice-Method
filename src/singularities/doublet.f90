! !> =============================================================================
! !> \module   doublet
! !> \brief    2D/3D Doublet Singularities (mu)
! !> \details  Provides velocity induction for point doublets and constant-strength
! !>           doublet panels.
! !> 
! !> \author   Georgios Loukas / PhD Candidate NTUA
! !> \date     2026
! !> =============================================================================
! module doublet
!     use constants
!     implicit none
!     private

!     ! Public Interface
!     public :: point_doublet_2D
!     public :: compute_doublet_panel2D_velocity

! contains

!     !> =========================================================================
!     !> \brief   Calculate the velocity induced by a 2D Point Doublet.
!     !> \details The doublet is oriented along the X-axis by default.
!     !> 
!     !> \param   p        Field Point (x,y) where velocity is evaluated
!     !> \param   p0       Location of the doublet (x0,y0)
!     !> \param   mu       Doublet strength (Circulation * distance)
!     !> \param   beta     (Optional) Orientation angle in radians (CCW from X-axis)
!     !> \return  velocity Induced (u,v) velocity vector
!     !> =========================================================================
!     pure function point_doublet_2D(p, p0, mu, beta) result(velocity)
!         real(wp), intent(in)           :: p(2), p0(2)
!         real(wp), intent(in)           :: mu
!         real(wp), intent(in), optional :: beta
!         real(wp) :: velocity(2)
        
!         real(wp) :: d(2), d_loc(2)
!         real(wp) :: r2, inv_r4, u_loc, v_loc
!         real(wp) :: c, s

!         ! 1. Calculate relative distance vector
!         d = p - p0
!         r2 = dot_product(d, d)
        
!         ! 2. Check for singularity
!         if (r2 < ZERO) then
!             velocity = 0.0_wp
!             return
!         end if

!         ! 3. Coordinate Transformation to Local Frame
!         if (present(beta)) then
!             c = cos(beta)
!             s = sin(beta)
!             ! Rotate field point relative to doublet orientation
!             d_loc(1) =  d(1) * c + d(2) * s
!             d_loc(2) = -d(1) * s + d(2) * c
!         else
!             d_loc = d
!         end if

!         ! 4. Velocity in Local Frame (Standard X-axis Doublet)
!         ! Based on gradient of Potential: Phi = -mu/(2*PI) * (x/r^2)
!         inv_r4 = 1.0_wp / (r2 * r2)
!         u_loc = (mu / TWO_PI) * (d_loc(2)**2 - d_loc(1)**2) * inv_r4
!         v_loc = (mu / TWO_PI) * (-2.0_wp * d_loc(1) * d_loc(2)) * inv_r4

!         ! 5. Transform Velocity back to Global Frame
!         if (present(beta)) then
!             velocity(1) = u_loc * c - v_loc * s
!             velocity(2) = u_loc * s + v_loc * c
!         else
!             velocity = [u_loc, v_loc]
!         end if

!     end function point_doublet_2D

!     !> =========================================================================
!     !> \brief   Calculate the velocity induced by a Constant Strength Doublet Panel.
!     !> \details Based on Katz & Plotkin, Eqs 10.29 - 10.30.
!     !>          Equivalence: A 2D doublet panel is equivalent to two point vortices 
!     !>          at the edges.
!     !> 
!     !> \param   p        Field point (x, z)
!     !> \param   p1       Panel start node (x1, z1)
!     !> \param   p2       Panel end node (x2, z2)
!     !> \param   mu       Doublet strength (constant across panel)
!     !> \return  vel_global Induced (u, w) velocity vector
!     !> =========================================================================
!     pure function compute_doublet_panel2D_velocity(p, p1, p2, mu) result(vel_global)
!         real(wp), intent(in) :: p(2), p1(2), p2(2)
!         real(wp), intent(in) :: mu
!         real(wp) :: vel_global(2)

!         real(wp) :: p_vec(2), tang(2), norm(2)
!         real(wp) :: L, d(2), loc(2)
!         real(wp) :: u_loc, w_loc, r1_sq, r2_sq

!         ! 1. Geometric Properties
!         p_vec = p2 - p1
!         L     = norm2(p_vec)
        
!         if (L < ZERO) then
!             vel_global = 0.0_wp
!             return
!         end if

!         tang = p_vec / L
!         norm = [-tang(2), tang(1)]

!         ! 2. Coordinate Transformation (Global -> Local)
!         d = p - p1
!         loc(1) = dot_product(d, tang)
!         loc(2) = dot_product(d, norm)

!         ! 3. Distances from edges
!         r1_sq = loc(1)**2 + loc(2)**2
!         r2_sq = (loc(1) - L)**2 + loc(2)**2

!         ! 4. Local Velocity Calculation
!         if (r1_sq < ZERO .or. r2_sq < ZERO) then
!             vel_global = 0.0_wp
!             return
!         end if

!         ! Eq 10.29 & 10.30
!         u_loc = -(mu / TWO_PI) * ( (loc(2) / r1_sq) - (loc(2) / r2_sq) )
!         w_loc =  (mu / TWO_PI) * ( (loc(1) / r1_sq) - ((loc(1) - L) / r2_sq) )

!         ! Handling self-induced surface singularity (On the panel)
!         if (abs(loc(2)) < ZERO .and. (loc(1) > 0.0_wp .and. loc(1) < L)) then
!              u_loc = 0.0_wp
!         end if

!         ! 5. Local -> Global Transformation
!         vel_global = u_loc * tang + w_loc * norm

!     end function compute_doublet_panel2D_velocity

! end module doublet