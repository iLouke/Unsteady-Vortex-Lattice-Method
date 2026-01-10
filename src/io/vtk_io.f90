MODULE vtk_io
  USE constants
  USE mesh_types
  USE wake_types  ! <--- Needed for WakeMesh
  USE config_io
  USE logger_io
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Export_Results_To_VTK, Export_Unsteady_Wake_VTK

CONTAINS

  !> ===========================================================================
  !> \brief Writes the Unsteady Wake Mesh to a numbered VTK file
  !> ===========================================================================
  SUBROUTINE Export_Unsteady_Wake_VTK(Wake, Step)
    TYPE(WakeMesh), INTENT(IN) :: Wake
    INTEGER, INTENT(IN) :: Step
    
    INTEGER :: u, i, j, n_pts, n_cells
    CHARACTER(LEN=64) :: Filename
    CHARACTER(LEN=10) :: StepStr
    INTEGER :: idx1, idx2, idx3, idx4
    INTEGER :: row_offset_curr, row_offset_next
    
    ! 1. Format Filename: wake_0001.vtk, wake_0002.vtk ...
    WRITE(StepStr, '(I4.4)') Step
    Filename = "wake_" // TRIM(StepStr) // ".vtk"
    
    OPEN(NEWUNIT=u, FILE=Filename, STATUS='REPLACE', ACTION='WRITE')
    
    WRITE(u, '(A)') "# vtk DataFile Version 3.0"
    WRITE(u, '(A)') "Unsteady VLM Wake"
    WRITE(u, '(A)') "ASCII"
    WRITE(u, '(A)') "DATASET POLYDATA"
    
    ! 2. Points
    ! Wake has (CurrentStep + 1) Rows and (N_Span_Nodes) Columns
    n_pts = (Wake%CurrentStep + 1) * Wake%N_Span_Nodes
    
    WRITE(u, '(A, I8, A)') "POINTS ", n_pts, " float"
    
    ! Write row by row
    DO i = 1, Wake%CurrentStep + 1
       DO j = 1, Wake%N_Span_Nodes
          WRITE(u, '(3F12.6)') REAL(Wake%Nodes(i, j)%r, 4)
       END DO
    END DO
    
    ! 3. Cells (Connectivity)
    n_cells = Wake%CurrentStep * Wake%N_Span_Panels
    
    WRITE(u, '(A, I8, I8)') "POLYGONS ", n_cells, n_cells * 5
    
    DO i = 1, Wake%CurrentStep
       DO j = 1, Wake%N_Span_Panels
          ! Calculate 0-based indices for VTK
          ! Row i, Col j. 
          ! P1=TopRight, P2=TopLeft, P3=BotLeft, P4=BotRight (Viewed from top?)
          ! Wake Nodes map: Index = (Row-1)*N_Span_Nodes + (Col-1)
          
          ! Let's follow the standard Ring winding: 
          ! Node(i+1, j) -> Node(i+1, j+1) -> Node(i, j+1) -> Node(i, j)
          
          
          
          row_offset_curr = (i - 1) * Wake%N_Span_Nodes
          row_offset_next = (i)     * Wake%N_Span_Nodes
          
          ! VTK Indices (0-based)
          idx1 = row_offset_next + (j - 1)      ! Node(i+1, j)
          idx2 = row_offset_next + (j)          ! Node(i+1, j+1)
          idx3 = row_offset_curr + (j)          ! Node(i,   j+1)
          idx4 = row_offset_curr + (j - 1)      ! Node(i,   j)
          
          WRITE(u, '(I2, 4(1X,I8))') 4, idx1, idx2, idx3, idx4
       END DO
    END DO
    
    ! 4. Data (Circulation)
    WRITE(u, '(A, I8)') "CELL_DATA ", n_cells
    WRITE(u, '(A)') "SCALARS Circulation float 1"
    WRITE(u, '(A)') "LOOKUP_TABLE default"
    
    DO i = 1, Wake%CurrentStep
       DO j = 1, Wake%N_Span_Panels
          WRITE(u, '(E12.5)') REAL(Wake%Panels(i, j)%Circulation, 4)
       END DO
    END DO
    
    CLOSE(u)
  END SUBROUTINE Export_Unsteady_Wake_VTK

  ! ... (Keep Export_Results_To_VTK and Write_Body_VTK as they were) ...
  SUBROUTINE Export_Results_To_VTK(Body, Cfg)
      TYPE(BodyMesh), INTENT(IN) :: Body
      TYPE(SimSettings), INTENT(IN) :: Cfg
      CALL Write_Body_VTK(Body)
      ! For steady call, we ignore wake or write rigid wake
  END SUBROUTINE Export_Results_To_VTK

  SUBROUTINE Write_Body_VTK(Body)
    TYPE(BodyMesh), INTENT(IN) :: Body
    INTEGER :: u, i, total_pts
    CHARACTER(LEN=64) :: Filename = "body.vtk"
    
    OPEN(NEWUNIT=u, FILE=Filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(u, '(A)') "# vtk DataFile Version 3.0"
    WRITE(u, '(A)') "VLM Body Surface"
    WRITE(u, '(A)') "ASCII"
    WRITE(u, '(A)') "DATASET POLYDATA"
    
    WRITE(u, '(A, I8, A)') "POINTS ", Body%NNodes, " float"
    DO i = 1, Body%NNodes
       WRITE(u, '(3F12.6)') REAL(Body%Nodes(i)%r, 4)
    END DO
    
    total_pts = 0
    DO i = 1, Body%NPanels
       total_pts = total_pts + 1 + Body%Panels(i)%NVertices
    END DO
    
    WRITE(u, '(A, I8, I8)') "POLYGONS ", Body%NPanels, total_pts
    DO i = 1, Body%NPanels
       IF (Body%Panels(i)%NVertices == 4) THEN
          WRITE(u, '(I2, 4(1X,I8))') 4, Body%Panels(i)%Vertices(1:4) - 1
       ELSE
          WRITE(u, '(I2, 3(1X,I8))') 3, Body%Panels(i)%Vertices(1:3) - 1
       END IF
    END DO
    
    WRITE(u, '(A, I8)') "CELL_DATA ", Body%NPanels
    WRITE(u, '(A)') "SCALARS Circulation float 1"
    WRITE(u, '(A)') "LOOKUP_TABLE default"
    DO i = 1, Body%NPanels
       WRITE(u, '(E12.5)') REAL(Body%Panels(i)%Circulation, 4)
    END DO
    
    WRITE(u, '(A)') "VECTORS Forces float"
    DO i = 1, Body%NPanels
       WRITE(u, '(3E12.5)') REAL(Body%Panels(i)%Force, 4)
    END DO
    CLOSE(u)
  END SUBROUTINE Write_Body_VTK

END MODULE vtk_io