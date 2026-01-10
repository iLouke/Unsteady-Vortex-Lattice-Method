MODULE uvlm
  USE constants
  USE mesh_types
  USE atmosphere
  USE config_io
  USE geometry_io
  USE logger_io
  USE geometry_calc
  USE aero_solver
  USE vtk_io
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Initialize_System, Load_Geometry, Finalize_System, Run_Solver
  
  TYPE(SimSettings)  :: Config
  TYPE(BodyMesh)     :: Body

CONTAINS

  SUBROUTINE Initialize_System()
    CHARACTER(LEN=128) :: InputFilename
    CHARACTER(LEN=1)   :: UserResp
    LOGICAL            :: FileExists
    INTEGER            :: NumArgs

    ! 1. CLI Args
    NumArgs = COMMAND_ARGUMENT_COUNT()
    IF (NumArgs > 0) THEN
       CALL GET_COMMAND_ARGUMENT(1, InputFilename)
    ELSE
       PRINT *, "Please enter the name of the configuration file:"
       READ(*, '(A)') InputFilename
       IF (TRIM(InputFilename) == "") InputFilename = "config.nml"
    END IF

    ! 2. Start Logging
    CALL Init_Logging() 
    CALL Write_Log("System Initialized.")

    ! 3. Load Config
    INQUIRE(FILE=TRIM(InputFilename), EXIST=FileExists)
    IF (FileExists) THEN
       CALL Read_Settings_File(TRIM(InputFilename), Config)
    ELSE
       PRINT *, "Config not found. Create one? (y/n)"
       READ(*, '(A)') UserResp
       IF (UserResp == 'y' .OR. UserResp == 'Y') THEN
          CALL Run_Wizard_And_Save(TRIM(InputFilename), Config)
       ELSE
          STOP
       END IF
    END IF

    ! 4. Auto-Detect & Display
    CALL Auto_Detect_Mesh_Size(Config)
   !  CALL Display_Config(Config)

  END SUBROUTINE Initialize_System

  SUBROUTINE Load_Geometry()
    CHARACTER(LEN=128) :: Buf
    
    CALL Write_Log("Loading Geometry...")
    CALL Load_Legacy_Geometry(Config, Body)
    
    CALL Write_Log("Transforming Coordinates...")
    CALL Transform_Mesh_Coordinates(Body, Config)
    
    CALL Write_Log("Computing Panel Properties...")
    CALL Compute_Panel_Properties(Body)
    
    CALL Write_Log("Geometry Loaded.")
    
    CALL Write_Log(" ")
    CALL Write_Log("--- Geometric Parameters ---")
    WRITE(Buf, '(T30, F12.4)') Config%RefArea;   CALL Write_Log("  Wing Ref Area (S):"//TRIM(Buf))
    WRITE(Buf, '(T30, F12.4)') Config%WingSpan;  CALL Write_Log("  Wing Span (b):    "//TRIM(Buf))
    WRITE(Buf, '(T30, F12.4)') Config%AspectRatio; CALL Write_Log("  Aspect Ratio:     "//TRIM(Buf))
    CALL Write_Log("--------------------------")
  END SUBROUTINE Load_Geometry

  ! =========================================================
  !  UPDATED DRIVER ROUTINE
  ! =========================================================
  SUBROUTINE Run_Solver()
      CALL Write_Log(" ")
      CALL Write_Log("--- Executing Unsteady Solver ---")
      
      ! The new engine handles the Time Loop, Solving, and Loads internally
      CALL Run_Unsteady_Simulation(Body, Config)
      
      ! Note: Visualization export is now handled inside the unsteady loop 
      ! (or you can move it here if you only want the final state)
      
  END SUBROUTINE Run_Solver

  SUBROUTINE Finalize_System()
    CALL Finalize_IO()
    IF (ALLOCATED(Body%Nodes))  DEALLOCATE(Body%Nodes)
    IF (ALLOCATED(Body%Panels)) DEALLOCATE(Body%Panels)
    PRINT *, "Program Finished."
  END SUBROUTINE Finalize_System

END MODULE uvlm