!> =============================================================================
!> \module   mesh_types
!> \brief    Defines Panel, Wing, Nodes
!> \details  
!>
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2025
!> =============================================================================
MODULE mesh_types
  USE constants
  IMPLICIT NONE

  TYPE :: Point3D
     REAL(wp) :: r(3)   ! Coordinates
     INTEGER  :: Marker ! Track Trailing Edges (1=TE, 0=Normal)
  END TYPE Point3D

  TYPE :: Panel
     INTEGER  :: ID
     INTEGER  :: Vertices(4)   
     INTEGER  :: NVertices     
     INTEGER  :: PropertyID    
     
     ! --- Geometry ---
     REAL(wp) :: Center(3)
     REAL(wp) :: Area
     
     ! --- Local Coordinate System (Basis Vectors) ---
     ! LocalX: Longitudinal (Streamwise)
     ! LocalY: Transverse (Spanwise)
     ! LocalZ: Normal
     REAL(wp) :: LocalX(3)
     REAL(wp) :: LocalY(3)
     REAL(wp) :: LocalZ(3)
     
     ! --- Physics ---
     REAL(wp) :: Circulation   
     REAL(wp) :: Force(3)
     REAL(wp) :: Moment(3)

     ! Topology
     LOGICAL :: IsTrailingEdge ! True if this panel sheds a wake
     INTEGER :: WakeNodes(2)   ! Indices of the two nodes forming the TE
  END TYPE Panel

  TYPE :: BodyMesh
     TYPE(Point3D), ALLOCATABLE   :: Nodes(:)
     TYPE(Panel)  , ALLOCATABLE   :: Panels(:)
     INTEGER  :: NPanels, NNodes
  END TYPE BodyMesh

CONTAINS

  ! Constructor to allocate the mesh dynamically
  SUBROUTINE Init_Mesh(Mesh, n_nodes, n_panels)
    TYPE(BodyMesh), INTENT(OUT) :: Mesh
    INTEGER, INTENT(IN) :: n_nodes, n_panels
    ALLOCATE(Mesh%Nodes(n_nodes))
    ALLOCATE(Mesh%Panels(n_panels))
    Mesh%NNodes = n_nodes
    Mesh%NPanels = n_panels
    Mesh%Nodes(:)%Marker = 0
  END SUBROUTINE Init_Mesh

END MODULE mesh_types