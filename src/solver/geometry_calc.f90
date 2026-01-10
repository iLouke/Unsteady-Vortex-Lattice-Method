MODULE geometry_calc
  USE constants
  USE mesh_types
  USE config_io
  USE vector_utils
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Transform_Mesh_Coordinates, Compute_Panel_Properties

CONTAINS

  !-------------------------------------------------------------------------
  ! Rotates and Translates the Mesh Nodes
  !-------------------------------------------------------------------------
  SUBROUTINE Transform_Mesh_Coordinates(Body, Cfg)
    TYPE(BodyMesh), INTENT(INOUT) :: Body
    TYPE(SimSettings), INTENT(IN) :: Cfg
    
    REAL(wp) :: RotMat(3,3)
    REAL(wp) :: P_Old(3), P_New(3)
    INTEGER  :: i

    ! 1. Build Matrix
    RotMat = Create_Rotation_Matrix(Cfg%Alpha, Cfg%Beta, Cfg%Gamma)

    PRINT *, " "
    PRINT *, "--- Geometry Transformation ---"
    PRINT *, "   Applying Rotation (A/B/G):", Cfg%Alpha, Cfg%Beta, Cfg%Gamma
    PRINT *, "   Center of Rotation (CG):  ", Cfg%CG

    ! 2. Transform every node
    DO i = 1, Body%NNodes
       P_Old = Body%Nodes(i)%r
       
       ! Translate to Origin (Relative to CG)
       P_Old = P_Old - Cfg%CG
       
       ! Rotate
       P_New = MATMUL(RotMat, P_Old)
       
       ! Translate Back
       Body%Nodes(i)%r = P_New + Cfg%CG
    END DO

  END SUBROUTINE Transform_Mesh_Coordinates


  !-------------------------------------------------------------------------
  ! Calculates Normal, Tangent, Center for every panel
  !-------------------------------------------------------------------------
  SUBROUTINE Compute_Panel_Properties(Body)
    TYPE(BodyMesh), INTENT(INOUT) :: Body
    
    INTEGER :: i, n1, n2, n3, n4
    REAL(wp) :: P1(3), P2(3), P3(3), P4(3)
    REAL(wp) :: Diag1(3), Diag2(3), CrossP(3), VecToCenter(3)

    PRINT *, "--- Computing Panel Properties ---"

    DO i = 1, Body%NPanels
       ! Node Indices
       n1 = Body%Panels(i)%Vertices(1)
       n2 = Body%Panels(i)%Vertices(2)
       n3 = Body%Panels(i)%Vertices(3)
       n4 = Body%Panels(i)%Vertices(4)
       
       ! Coordinates
       P1 = Body%Nodes(n1)%r
       P2 = Body%Nodes(n2)%r
       P3 = Body%Nodes(n3)%r
       
       ! 1. CENTER Calculation
       IF (Body%Panels(i)%NVertices == 3) THEN
          Body%Panels(i)%Center = (P1 + P2 + P3) / 3.0_wp
       ELSE
          P4 = Body%Nodes(n4)%r
          Body%Panels(i)%Center = (P1 + P2 + P3 + P4) / 4.0_wp
       END IF

       ! 2. NORMAL (Local Z)
       IF (Body%Panels(i)%NVertices == 4) THEN
          Diag1 = P3 - P1
          Diag2 = P4 - P2
          CrossP = Cross_Product(Diag1, Diag2)
       ELSE
          CrossP = Cross_Product(P2 - P1, P3 - P1)
       END IF
       
       Body%Panels(i)%LocalZ = Normalize(CrossP)
       Body%Panels(i)%Area   = 0.5_wp * NORM2(CrossP)

       ! 3. LONGITUDINAL / TANGENT (Local X)
       ! Vector from Center to Node 1 (Legacy Logic)
       VecToCenter = P1 - Body%Panels(i)%Center
       Body%Panels(i)%LocalX = Normalize(VecToCenter)

       ! 4. TRANSVERSE (Local Y)
       ! Normal x Longitudinal
       Body%Panels(i)%LocalY = Cross_Product(Body%Panels(i)%LocalZ, &
                                             Body%Panels(i)%LocalX)

    END DO
    
    PRINT *, "   Processed", Body%NPanels, "panels."

  END SUBROUTINE Compute_Panel_Properties

END MODULE geometry_calc