MODULE vector_utils
  USE constants, only: wp, ZERO, DEG2RAD
  IMPLICIT NONE

  PRIVATE
    PUBLIC :: Create_Rotation_Matrix, Cross_Product, Normalize

CONTAINS

  !> Returns 3x3 Rotation Matrix from Euler Angles (Degrees)
  PURE FUNCTION Create_Rotation_Matrix(Alpha, Beta, Gamma) RESULT(R)
    REAL(wp), INTENT(IN) :: Alpha, Beta, Gamma ! Degrees
    REAL(wp) :: R(3,3)
    REAL(wp) :: a, b, g ! Radians
    REAL(wp) :: ca, sa, cb, sb, cg, sg

    ! Convert to Radians
    a = Alpha * DEG2RAD
    b = Beta  * DEG2RAD
    g = Gamma * DEG2RAD

    ca = COS(a); sa = SIN(a)
    cb = COS(b); sb = SIN(b)
    cg = COS(g); sg = SIN(g)

    ! Legacy Logic
    R(1,1) = ca*cb
    R(1,2) = ca*sb*sg - sa*cg
    R(1,3) = ca*sb*cg + sa*sg

    R(2,1) = sa*cb
    R(2,2) = sa*sb*sg + ca*cg
    R(2,3) = sa*sb*cg - ca*sg

    R(3,1) = -sb
    R(3,2) = cb*sg
    R(3,3) = cb*cg
  END FUNCTION Create_Rotation_Matrix

  !> Standard Cross Product
  PURE FUNCTION Cross_Product(u, v) RESULT(w)
    REAL(wp), INTENT(IN) :: u(3), v(3)
    REAL(wp) :: w(3)
    w(1) = u(2)*v(3) - u(3)*v(2)
    w(2) = u(3)*v(1) - u(1)*v(3)
    w(3) = u(1)*v(2) - u(2)*v(1)
  END FUNCTION Cross_Product

  !> Safe Normalize
  PURE FUNCTION Normalize(v) RESULT(n)
    REAL(wp), INTENT(IN) :: v(3)
    REAL(wp) :: n(3)
    REAL(wp) :: mag
    mag = NORM2(v)
    IF (mag > ZERO) THEN
       n = v / mag
    ELSE
       n = 0.0_wp
    END IF
  END FUNCTION Normalize

END MODULE vector_utils