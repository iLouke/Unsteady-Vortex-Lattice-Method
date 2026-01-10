MODULE linalg
  USE constants
  USE logger_io
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Factorize_Matrix, Solve_Linear_System_LU

  ! Define Interfaces for LAPACK routines
  INTERFACE
     ! DGETRF: LU Factorization
     SUBROUTINE DGETRF(M, N, A, LDA, IPIV, INFO)
       USE constants
       INTEGER, INTENT(IN) :: M, N, LDA
       REAL(wp), INTENT(INOUT) :: A(LDA, *)
       INTEGER, INTENT(OUT) :: IPIV(*), INFO
     END SUBROUTINE DGETRF

     ! DGETRS: Solve using LU Factors
     SUBROUTINE DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       USE constants
       CHARACTER(LEN=1), INTENT(IN) :: TRANS
       INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
       REAL(wp), INTENT(IN) :: A(LDA, *)
       INTEGER, INTENT(IN) :: IPIV(*)
       REAL(wp), INTENT(INOUT) :: B(LDB, *)
       INTEGER, INTENT(OUT) :: INFO
     END SUBROUTINE DGETRS
  END INTERFACE

CONTAINS

  !> ===========================================================================
  !> \brief Performs LU Decomposition on Matrix A.
  !> \details A is overwritten with the LU factors.
  !> ===========================================================================
  SUBROUTINE Factorize_Matrix(N, A, IPIV)
    INTEGER, INTENT(IN)    :: N
    REAL(wp), INTENT(INOUT):: A(N,N)  ! On exit, contains L and U factors
    INTEGER, INTENT(OUT)   :: IPIV(N) ! Pivot indices
    
    INTEGER :: INFO

    CALL DGETRF(N, N, A, N, IPIV, INFO)

    IF (INFO < 0) THEN
       CALL Write_Log("ERROR: DGETRF Argument illegal.")
       STOP
    ELSE IF (INFO > 0) THEN
       CALL Write_Log("ERROR: Matrix is Singular (Zero pivot).")
       STOP
    END IF
  END SUBROUTINE Factorize_Matrix

  !> ===========================================================================
  !> \brief Solves Ax=b using pre-computed LU factors from Factorize_Matrix.
  !> ===========================================================================
  SUBROUTINE Solve_Linear_System_LU(N, LU, IPIV, b, x)
    INTEGER, INTENT(IN)    :: N
    REAL(wp), INTENT(IN)   :: LU(N,N) ! The factorized matrix
    INTEGER, INTENT(IN)    :: IPIV(N) ! The pivot indices
    REAL(wp), INTENT(IN)   :: b(N)    ! The RHS vector
    REAL(wp), INTENT(OUT)  :: x(N)    ! The Solution
    
    INTEGER :: INFO
    REAL(wp) :: B_Mat(N, 1)

    ! Copy b to B_Mat (LAPACK expects a matrix for RHS)
    B_Mat(:, 1) = b(:)

    ! 'N' means No transpose of A
    CALL DGETRS('N', N, 1, LU, N, IPIV, B_Mat, N, INFO)

    IF (INFO /= 0) THEN
       CALL Write_Log("ERROR: DGETRS failed.")
       STOP
    END IF

    x(:) = B_Mat(:, 1)

  END SUBROUTINE Solve_Linear_System_LU

END MODULE linalg