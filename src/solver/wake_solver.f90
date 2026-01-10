MODULE wake_solver
  USE constants
  USE mesh_types
  USE wake_types
  USE config_io
  USE logger_io
  USE vortex, ONLY: compute_vortex_ring_velocity
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Initialize_Unsteady_Wake, Shed_Wake_Row, Convect_Wake, &
            Calc_Wake_Induced_Velocity, Update_Wake_Circulation

CONTAINS

  !> Identify TE nodes and setup the Wake Mesh structure
  SUBROUTINE Initialize_Unsteady_Wake(Body, Wake, Cfg)
    TYPE(BodyMesh), INTENT(INOUT) :: Body
    TYPE(WakeMesh), INTENT(OUT) :: Wake
    TYPE(SimSettings), INTENT(IN) :: Cfg
    
    INTEGER :: i, k, n_te, n1, n2
    REAL(wp) :: max_x
    CHARACTER(LEN=64) :: Buf
    
    CALL Write_Log(" ")
    CALL Write_Log("   --- Detecting Trailing Edge ---")

    ! =========================================================
    ! 1. DETECT TRAILING EDGES (Restored Logic)
    ! =========================================================
    n_te = 0
    
    ! Find Max X for geometric fallback
    max_x = -1.0E9_wp
    DO i = 1, Body%NNodes
       IF (Body%Nodes(i)%r(1) > max_x) max_x = Body%Nodes(i)%r(1)
    END DO
    
    DO i = 1, Body%NPanels
       Body%Panels(i)%IsTrailingEdge = .FALSE.
       
       ! Standard VLM Winding: Nodes 2 and 3 are the rear edge
       n1 = Body%Panels(i)%Vertices(2)
       n2 = Body%Panels(i)%Vertices(3)
       
       ! Check A: Markers (from legacy .DAT file)
       IF (Body%Nodes(n1)%Marker /= 0 .AND. Body%Nodes(n2)%Marker /= 0) THEN
          Body%Panels(i)%IsTrailingEdge = .TRUE.
          Body%Panels(i)%WakeNodes = [n2, n1] ! Winding order for wake shedding
          n_te = n_te + 1
          
       ! Check B: Geometric (At the back of the wing)
       ELSE IF (ABS(Body%Nodes(n1)%r(1) - max_x) < 0.01_wp .AND. &
                ABS(Body%Nodes(n2)%r(1) - max_x) < 0.01_wp) THEN
          Body%Panels(i)%IsTrailingEdge = .TRUE.
          Body%Panels(i)%WakeNodes = [n2, n1]
          n_te = n_te + 1
       END IF
    END DO

    IF (n_te == 0) THEN 
       CALL Write_Log("ERROR: No Trailing Edge detected.")
       STOP
    END IF

    WRITE(Buf, '(I4)') n_te
    CALL Write_Log("   Trailing Edge Panels Found: "//TRIM(ADJUSTL(Buf)))

    ! =========================================================
    ! 2. ALLOCATE WAKE
    ! =========================================================
    CALL Init_Wake_Mesh(Wake, n_te, Cfg%MaxIter)
    ALLOCATE(Wake%TE_Node_Map(Wake%N_Span_Nodes)) 
    
    ! =========================================================
    ! 3. MAP NODES
    ! =========================================================
    ! We must map the TE nodes to the Wake Row 1 indices.
    ! Simple Strip Assumption: Panels are ordered sequentially along span.
    ! Wake Node k   = TE Node 2 of Panel k
    ! Wake Node k+1 = TE Node 1 of Panel k
    
    k = 1 
    DO i = 1, Body%NPanels
       IF (Body%Panels(i)%IsTrailingEdge) THEN
          n1 = Body%Panels(i)%WakeNodes(1) 
          n2 = Body%Panels(i)%WakeNodes(2) 
          
          ! This simplistic mapping assumes panels are ordered 1..N across the span.
          ! For a complex mesh, you need a connectivity search here.
          ! For now, we trust the strip order.
          Wake%TE_Node_Map(k)   = n1
          Wake%TE_Node_Map(k+1) = n2
          k = k + 1
       END IF
    END DO
    
    ! Initialize "Row 1" (Attached to Body)
    DO k = 1, Wake%N_Span_Nodes
       Wake%Nodes(1, k)%r = Body%Nodes(Wake%TE_Node_Map(k))%r
    END DO

  END SUBROUTINE Initialize_Unsteady_Wake

  !> Create a new row of wake panels
  SUBROUTINE Shed_Wake_Row(Body, Wake)
    TYPE(BodyMesh), INTENT(IN) :: Body
    TYPE(WakeMesh), INTENT(INOUT) :: Wake
    INTEGER :: j, k
    
    Wake%CurrentStep = Wake%CurrentStep + 1
    
    ! Shift Node Rows (2 -> 3, 1 -> 2)
    DO j = Wake%CurrentStep + 1, 2, -1
       DO k = 1, Wake%N_Span_Nodes
          Wake%Nodes(j, k)%r = Wake%Nodes(j-1, k)%r
       END DO
    END DO
    
    ! Reset Row 1 to Body Position
    DO k = 1, Wake%N_Span_Nodes
       Wake%Nodes(1, k)%r = Body%Nodes(Wake%TE_Node_Map(k))%r
    END DO
    
    ! Shift Panel Circulations
    DO j = Wake%CurrentStep, 2, -1
       DO k = 1, Wake%N_Span_Panels
          Wake%Panels(j, k)%Circulation = Wake%Panels(j-1, k)%Circulation
       END DO
    END DO
  END SUBROUTINE Shed_Wake_Row

  !> Move all wake nodes by V_inf * dt
  SUBROUTINE Convect_Wake(Wake, Cfg)
    TYPE(WakeMesh), INTENT(INOUT) :: Wake
    TYPE(SimSettings), INTENT(IN) :: Cfg
    INTEGER :: i, j
    REAL(wp) :: V_conv(3)
    
    V_conv = [Cfg%Vinf, 0.0_wp, 0.0_wp]
    
    ! Start from Row 2 (Row 1 is pinned to the body)
    DO i = 2, Wake%CurrentStep + 1
       DO j = 1, Wake%N_Span_Nodes
          Wake%Nodes(i,j)%r = Wake%Nodes(i,j)%r + V_conv * Cfg%dt
       END DO
    END DO
  END SUBROUTINE Convect_Wake

  !> Calculate induced velocity from ALL wake rings
  FUNCTION Calc_Wake_Induced_Velocity(Wake, Target) RESULT(V_ind)
    TYPE(WakeMesh), INTENT(IN) :: Wake
    REAL(wp), INTENT(IN) :: Target(3)
    REAL(wp) :: V_ind(3), Corners(3,4)
    INTEGER :: i, j
    
    V_ind = 0.0_wp
    
    IF (Wake%CurrentStep == 0) RETURN
    
    DO i = 1, Wake%CurrentStep
       DO j = 1, Wake%N_Span_Panels
          ! Quad Points:
          ! P1(i+1, j)   P2(i+1, j+1)
          !    ^            |
          !    |            v
          ! P4(i,   j)   P3(i,   j+1)
          
          Corners(:, 1) = Wake%Nodes(i+1, j)%r
          Corners(:, 2) = Wake%Nodes(i+1, j+1)%r
          Corners(:, 3) = Wake%Nodes(i,   j+1)%r
          Corners(:, 4) = Wake%Nodes(i,   j)%r
          
          V_ind = V_ind + compute_vortex_ring_velocity( &
             Target, Corners, Wake%Panels(i,j)%Circulation)
       END DO
    END DO
  END FUNCTION Calc_Wake_Induced_Velocity

  !> Assign Gamma to the newly shed wake row (Row 1)
  SUBROUTINE Update_Wake_Circulation(Wake, Body)
    TYPE(WakeMesh), INTENT(INOUT) :: Wake
    TYPE(BodyMesh), INTENT(IN) :: Body
    INTEGER :: i, counter
    
    counter = 1
    DO i = 1, Body%NPanels
       IF (Body%Panels(i)%IsTrailingEdge) THEN
          ! Wake Gamma = Body Gamma (Kelvin's condition in simplest form)
          ! Ideally: Gamma_Wake = Gamma_Body_New - Gamma_Body_Old
          ! For steady start, usually Gamma_Wake = Gamma_Body
          Wake%Panels(1, counter)%Circulation = Body%Panels(i)%Circulation
          counter = counter + 1
       END IF
    END DO
  END SUBROUTINE Update_Wake_Circulation

END MODULE wake_solver