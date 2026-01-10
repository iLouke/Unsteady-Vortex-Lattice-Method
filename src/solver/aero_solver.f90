MODULE aero_solver
  USE constants
  USE mesh_types
  USE wake_types
  USE wake_solver
  USE config_io
  USE logger_io
  USE linalg        ! Now using the split routines
  USE vector_utils
  USE vtk_io 
  USE vortex, ONLY: compute_vortex_ring_velocity
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Run_Unsteady_Simulation

CONTAINS

   SUBROUTINE INITIALISE_AIC(Body, AIC, Config)
      TYPE(BodyMesh), INTENT(IN) :: Body
      TYPE(SimSettings), INTENT(IN) :: Config
      REAL(wp), INTENT(OUT) :: AIC(Body%NPanels, Body%NPanels)

      INTEGER :: i, j, k, panels, node_idx
      INTEGER :: SYM, GRND
      REAL(wp) :: Corners(3,4), Normal(3), ControlPoint(3), VelInduced(3)
      REAL(wp) :: VortexStrength, MultY, MultZ

      AIC = 0.0_wp
      panels = Body%NPanels

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP PRIVATE(i, j, k, sym, grnd, ControlPoint, Normal, VelInduced, Corners, VortexStrength, MultY, MultZ, node_idx) &
      !$OMP SCHEDULE(DYNAMIC)
      do i = 1, panels
         ControlPoint = Body%Panels(i)%Center
         Normal       = Body%Panels(i)%LocalZ

         do j = 1, panels
            VelInduced = 0.0_wp
            
            do sym = 0, Config%Symmetry
               do grnd = 0, Config%GroundEffect
                  
                  VortexStrength = 1.0_wp * ((-1.0_wp)**sym) * ((-1.0_wp)**grnd)
                  MultY = (-1.0_wp)**sym
                  MultZ = (-1.0_wp)**grnd

                  ! --- Geometry Transformation (Inlined) ---
                  DO k = 1, 4
                     node_idx = Body%Panels(j)%Vertices(k)
                     Corners(1,k) = Body%Nodes(node_idx)%r(1)
                     Corners(2,k) = Body%Nodes(node_idx)%r(2) * MultY
                     Corners(3,k) = Body%Nodes(node_idx)%r(3) * MultZ
                  END DO
                  ! -----------------------------------------

                  IF (.NOT. CHECK_FAR_FIELD(ControlPoint, Corners(:,1), Corners(:,2), Corners(:,3), Corners(:,4))) THEN
                      VelInduced = VelInduced + compute_vortex_ring_velocity(ControlPoint, Corners, VortexStrength)
                  END IF
               end do
            end do
            
            AIC(i,j) = AIC(i,j) + DOT_PRODUCT(VelInduced, Normal)
         end do
      end do
      !$OMP END PARALLEL DO
   END SUBROUTINE INITIALISE_AIC

   ! Helper FUNCTION to calculate far field
   PURE FUNCTION CHECK_FAR_FIELD(p, p1, p2, p3, p4) RESULT(is_far)
      REAL(wp), INTENT(IN) :: p(3), p1(3), p2(3), p3(3), p4(3)
      REAL(wp) :: diagonal, distance, center(3)
      LOGICAL :: is_far

      if (NORM2(p3 - p4) < 1.0E-12_wp) then
         diagonal = MAX(NORM2(p3 - p1), NORM2(p3 - p2))
         center   = (1.0_wp/3.0_wp) * (p1 + p2 + p3)
      else
         diagonal = MAX(NORM2(p3 - p1), NORM2(p4 - p2))
         center   = 0.25_wp * (p1 + p2 + p3 + p4)
      end if
      distance = NORM2(p - center)

      is_far = distance > 5.0_wp * diagonal
   END FUNCTION CHECK_FAR_FIELD

  SUBROUTINE Run_Unsteady_Simulation(Body, Cfg)
    TYPE(BodyMesh), INTENT(INOUT) :: Body
    TYPE(SimSettings), INTENT(IN) :: Cfg
    TYPE(WakeMesh) :: Wake
    
    INTEGER :: Step, i
    REAL(wp), ALLOCATABLE :: AIC(:,:), RHS(:), BodyGamma(:)
    INTEGER, ALLOCATABLE  :: IPIV(:) ! <--- Added for LU Pivots
    
    REAL(wp) :: V_wake_ind(3), V_inf_vec(3)
    CHARACTER(LEN=64) :: LogBuf
    
    CALL Write_Log(" ")
    CALL Write_Log("--- Starting Unsteady VLM Simulation ---")

    ALLOCATE(AIC(Body%NPanels, Body%NPanels))
    ALLOCATE(RHS(Body%NPanels))
    ALLOCATE(BodyGamma(Body%NPanels))
    ALLOCATE(IPIV(Body%NPanels)) ! <--- Allocate Pivots

    ! ---------------------------------------------------------
    ! 1. PRE-PROCESS: Build and Factorize AIC ONCE
    ! ---------------------------------------------------------
    CALL Write_Log("Building AIC Matrix...")
    CALL INITIALISE_AIC(Body, AIC, Cfg)
    
    CALL Write_Log("Factorizing AIC Matrix (LU Decomposition)...")
    ! AIC is now overwritten with LU factors. We keep this "inverted" form.
    CALL Factorize_Matrix(Body%NPanels, AIC, IPIV) 
    
    ! ---------------------------------------------------------
    ! 2. Initialize Wake
    ! ---------------------------------------------------------
    CALL Initialize_Unsteady_Wake(Body, Wake, Cfg)
    V_inf_vec = [Cfg%Vinf, 0.0_wp, 0.0_wp]

    ! ---------------------------------------------------------
    ! 3. Time Loop
    ! ---------------------------------------------------------
    DO Step = 1, Cfg%MaxIter
       WRITE(LogBuf, '(A, I4)') " Time Step: ", Step
       CALL Write_Log(TRIM(LogBuf))
       
       ! A. Convect Wake
       CALL Convect_Wake(Wake, Cfg)
       
       ! B. Shed New Row
       CALL Shed_Wake_Row(Body, Wake)
       
       ! C. Update RHS (Must be done every step)
       !    RHS = -(V_inf + V_wake) . n
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, V_wake_ind)
       DO i = 1, Body%NPanels
          V_wake_ind = Calc_Wake_Induced_Velocity(Wake, Body%Panels(i)%Center)
          RHS(i) = -DOT_PRODUCT(V_wake_ind + V_inf_vec, Body%Panels(i)%LocalZ)
       END DO
       !$OMP END PARALLEL DO
       
       ! D. Solve using PRE-COMPUTED Factors
       !    AIC contains LU factors, IPIV contains pivots
       CALL Solve_Linear_System_LU(Body%NPanels, AIC, IPIV, RHS, BodyGamma)
       
       ! E. Update Circulation
       DO i = 1, Body%NPanels
          Body%Panels(i)%Circulation = BodyGamma(i)
       END DO
       
       ! F. Update Wake Gamma
       CALL Update_Wake_Circulation(Wake, Body)
       
       ! G. Exports
       CALL Export_Unsteady_Wake_VTK(Wake, Step)
       
    END DO

    CALL Export_Results_To_VTK(Body, Cfg)
    
    DEALLOCATE(AIC, RHS, BodyGamma, IPIV)
    CALL Write_Log("--- Simulation Complete ---")
    
  END SUBROUTINE Run_Unsteady_Simulation

END MODULE aero_solver