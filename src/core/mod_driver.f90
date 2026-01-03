!> =============================================================================
!> \module   mod_driver
!> \brief    High-Level API / Orchestrator for ROS
!> \details  Manages the lifecycle of a simulation:
!>           1. Configuration Load
!>           2. Mesh Generation
!>           3. Solver Dispatch (LLT, VLM, Panel)
!>           4. Data Export (VTK)
!>
!> \author   Georgios Loukas / PhD Candidate NTUA
!> \date     2025
!> =============================================================================
module mod_driver
!     use mod_precision
!     use mod_constants
!     use mod_mesh
!     use mod_vtk
!     use mod_atmosphere
    
!     ! --- Import Solvers ---
!     use mod_lifting_line
!     ! use mod_vlm_solver    ! <--- Future: Uncomment when ready
!     ! use mod_panel_solver  ! <--- Future: Uncomment when ready

!     ! --- Import Configs ---
!     use mod_cfg_LLT
!     ! use mod_cfg_VLM       ! <--- Future

!     implicit none
!     private

!     public :: ros_simulation

!     !> \type ros_simulation
!     !> \brief The main container for the simulation state
!     type :: ros_simulation
!         ! State Objects
!         type(t_mesh)     :: wing
!         type(vtk_writer) :: writer
        
!         ! Results Container (Generic)
!         real(wp) :: CL, CD, CY
!         real(wp), allocatable :: solution_field(:) ! Can hold Gamma or Cp
        
!     contains
!         procedure :: run_analysis
!         procedure :: save_results
!         procedure :: clean
!     end type ros_simulation

! contains

!     !> =========================================================================
!     !> \brief   Run a specific analysis method
!     !> \param   method    String: "LLT", "VLM", "PANEL"
!     !> \param   conf_file Path to .nml file
!     !> =========================================================================
!     subroutine run_analysis(self, method, conf_file)
!         class(ros_simulation), intent(inout) :: self
!         character(len=*), intent(in) :: method
!         character(len=*), intent(in) :: conf_file

!         ! Local Config Objects
!         type(config_llt) :: cfg_llt
!         ! type(config_vlm) :: cfg_vlm  ! <--- Future

!         print *, "================================================="
!         print *, "       ROS: Research on Open Solvers             "
!         print *, "       Method: ", trim(method)
!         print *, "================================================="

!         select case (trim(method))
        
!         ! ---------------------------------------------------------------------
!         ! CASE: LIFTING LINE THEORY
!         ! ---------------------------------------------------------------------
!         case ("LLT")
!             ! 1. Setup
!             print *, "--> [Driver] Initializing LLT Configuration..."
!             ! Check if file exists, otherwise use defaults inside config
!             call cfg_llt%initialize(input_file=conf_file)
!             call cfg_llt%display()

!             ! 2. Geometry
!             print *, "--> [Driver] Generating Lifting Line Mesh..."
!             call self%wing%generate(cfg_llt)

!             ! 3. Solve
!             print *, "--> [Driver] Running LLT Solver..."
!             call solve_llt_system(self%wing, cfg_llt%alpha, &
!                                   self%CL, self%CD, self%solution_field)

!         ! ---------------------------------------------------------------------
!         ! CASE: VORTEX LATTICE METHOD (Future Placeholder)
!         ! ---------------------------------------------------------------------
!         case ("VLM")
!             print *, "--> [Driver] VLM not yet implemented."
!             ! call cfg_vlm%init(conf_file)
!             ! call self%wing%generate_surface(cfg_vlm)
!             ! call solve_vlm(...)
        
!         ! ---------------------------------------------------------------------
!         ! CASE: PANEL METHOD (Future Placeholder)
!         ! ---------------------------------------------------------------------
!         case ("PANEL")
!              print *, "--> [Driver] Panel Method not yet implemented."

!         ! ---------------------------------------------------------------------
!         ! DEFAULT
!         ! ---------------------------------------------------------------------
!         case default
!             print *, "[ERROR] Unknown method: ", trim(method)
!             stop

!         end select

!         ! 4. Report to Console
!         print *, " "
!         print *, "--- Summary Results ---"
!         print '(A, F10.5)', "  CL : ", self%CL
!         print '(A, F10.5)', "  CD : ", self%CD
!         print *, "-----------------------"

!     end subroutine run_analysis

!     !> =========================================================================
!     !> \brief   Export results to VTK
!     !> =========================================================================
!     subroutine save_results(self, filename)
!         class(ros_simulation), intent(inout) :: self
!         character(len=*), intent(in) :: filename

!         print *, "--> [Driver] Exporting to ", trim(filename)
        
!         call self%writer%initialize()
        
!         ! 1. Save Grid
!         call self%writer%set_mesh(self%wing%pts, self%wing%conn)

!         ! 2. Save Scalar Data
!         if (allocated(self%solution_field)) then
!             call self%writer%add_property("Solution_Data", self%solution_field)
!         end if
        
!         ! 3. Save Geometry Data
!         if (allocated(self%wing%chords)) then
!             call self%writer%add_property("Chord", self%wing%chords)
!         end if

!         call self%writer%write(filename)

!     end subroutine save_results

!     !> =========================================================================
!     !> \brief   Destructor
!     !> =========================================================================
!     subroutine clean(self)
!         class(ros_simulation), intent(inout) :: self
!         call self%wing%free()
!         call self%writer%free()
!         if (allocated(self%solution_field)) deallocate(self%solution_field)
!     end subroutine clean

end module mod_driver