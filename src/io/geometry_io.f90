MODULE geometry_io
  USE constants
  USE mesh_types
  USE config_io
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Auto_Detect_Mesh_Size, Load_Legacy_Geometry

CONTAINS

  !-------------------------------------------------------------------------
  ! Check file lengths and update Config counts
  !-------------------------------------------------------------------------
  SUBROUTINE Auto_Detect_Mesh_Size(Cfg)
    TYPE(SimSettings), INTENT(INOUT) :: Cfg
    INTEGER :: RealNodes, RealPanels
    
    PRINT *, " "
    PRINT *, "--- Verifying Geometry Files ---"
    
    RealNodes = Count_File_Lines(Cfg%DatFile)
    IF (RealNodes == 0) STOP "ERROR: Node file empty/missing."
    PRINT *, "  File '", TRIM(Cfg%DatFile), "': Found", RealNodes, "nodes."
    
    RealPanels = Count_File_Lines(Cfg%PanFile)
    IF (RealPanels == 0) STOP "ERROR: Panel file empty/missing."
    PRINT *, "  File '", TRIM(Cfg%PanFile), "': Found", RealPanels, "panels."

    Cfg%NGrid   = RealNodes
    Cfg%NPanels = RealPanels
  END SUBROUTINE Auto_Detect_Mesh_Size

  !-------------------------------------------------------------------------
  ! Parse Files
  !-------------------------------------------------------------------------
  SUBROUTINE Load_Legacy_Geometry(Cfg, Body)
    TYPE(SimSettings), INTENT(IN) :: Cfg
    TYPE(BodyMesh), INTENT(OUT)   :: Body
    INTEGER :: u, i, dummy, mark
    REAL(wp) :: x, y, z
    
    CALL Init_Mesh(Body, Cfg%NGrid, Cfg%NPanels)

    ! Nodes
    PRINT *, "Loading Nodes..."
    OPEN(NEWUNIT=u, FILE=TRIM(Cfg%DatFile), STATUS='OLD', ACTION='READ')
    DO i = 1, Cfg%NGrid
       READ(u, *) dummy, x, y, z, mark
       Body%Nodes(i)%r = [x, y, z + Cfg%HFlight]
       Body%Nodes(i)%Marker = mark
    END DO
    CLOSE(u)

    ! Panels
    PRINT *, "Loading Panels..."
    OPEN(NEWUNIT=u, FILE=TRIM(Cfg%PanFile), STATUS='OLD', ACTION='READ')
    DO i = 1, Cfg%NPanels
       READ(u, *) dummy, Body%Panels(i)%Vertices(1:4), Body%Panels(i)%PropertyID
       Body%Panels(i)%ID = i
       IF (Body%Panels(i)%Vertices(3) == Body%Panels(i)%Vertices(4)) THEN
          Body%Panels(i)%NVertices = 3
       ELSE
          Body%Panels(i)%NVertices = 4
       END IF
    END DO
    CLOSE(u)
  END SUBROUTINE Load_Legacy_Geometry

  !-------------------------------------------------------------------------
  ! Internal Helper
  !-------------------------------------------------------------------------
  FUNCTION Count_File_Lines(Filename) RESULT(Cnt)
    CHARACTER(LEN=*), INTENT(IN) :: Filename
    INTEGER :: Cnt, u, io
    LOGICAL :: ex
    Cnt = 0
    INQUIRE(FILE=Filename, EXIST=ex)
    IF (.NOT. ex) RETURN
    OPEN(NEWUNIT=u, FILE=Filename, STATUS='OLD', ACTION='READ', IOSTAT=io)
    IF (io /= 0) RETURN
    DO
       READ(u, *, IOSTAT=io)
       IF (io < 0) EXIT
       Cnt = Cnt + 1
    END DO
    CLOSE(u)
  END FUNCTION Count_File_Lines

END MODULE geometry_io