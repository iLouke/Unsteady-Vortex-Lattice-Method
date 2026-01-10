MODULE wake_types
  USE constants
  IMPLICIT NONE

  TYPE :: WakePoint
     REAL(wp) :: r(3)
  END TYPE WakePoint

  TYPE :: WakePanel
     INTEGER :: Vertices(4)
     REAL(wp) :: Circulation
  END TYPE WakePanel

  TYPE :: WakeMesh
     TYPE(WakePoint), ALLOCATABLE :: Nodes(:,:) 
     TYPE(WakePanel), ALLOCATABLE :: Panels(:,:)
     
     INTEGER :: N_Span_Nodes
     INTEGER :: N_Span_Panels
     INTEGER :: CurrentStep
     
     INTEGER, ALLOCATABLE :: TE_Node_Map(:)
  END TYPE WakeMesh

CONTAINS

  SUBROUTINE Init_Wake_Mesh(W, n_span_panels, max_steps)
    TYPE(WakeMesh), INTENT(OUT) :: W
    INTEGER, INTENT(IN) :: n_span_panels, max_steps
    INTEGER :: i, j
    
    W%N_Span_Panels = n_span_panels
    W%N_Span_Nodes  = n_span_panels + 1 
    W%CurrentStep   = 0
    
    ! Allocate Max Dimensions
    ALLOCATE(W%Nodes(max_steps + 1, W%N_Span_Nodes))
    ALLOCATE(W%Panels(max_steps, W%N_Span_Panels))
    
    ! --- FIX: Initialize using explicit loops ---
    DO i = 1, max_steps + 1
       DO j = 1, W%N_Span_Nodes
          W%Nodes(i,j)%r = 0.0_wp
       END DO
    END DO

    DO i = 1, max_steps
       DO j = 1, W%N_Span_Panels
          W%Panels(i,j)%Circulation = 0.0_wp
          W%Panels(i,j)%Vertices = 0
       END DO
    END DO
    
  END SUBROUTINE Init_Wake_Mesh

END MODULE wake_types