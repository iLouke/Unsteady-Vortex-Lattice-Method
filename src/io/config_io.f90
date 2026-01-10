MODULE config_io
  USE constants, only: wp
  USE atmosphere
  USE logger_io
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: SimSettings, Read_Settings_File, Run_Wizard_And_Save, Display_Config

  TYPE :: SimSettings
    ! --- Embedded Sub-systems ---
    TYPE(Air) :: Air

    ! --- File Paths ---
    CHARACTER(LEN=64) :: DatFile    ! Nodes file
    CHARACTER(LEN=64) :: PanFile    ! Panels file

    ! --- Mesh Counts ---
    INTEGER :: NPanels, NGrid

    ! --- Geometry ---
    REAL(wp) :: RefArea, WingSpan, AspectRatio, MAC
    REAL(wp) :: CG(3)
    INTEGER  :: Symmetry
    INTEGER  :: GroundEffect

    ! --- Flow ---
    REAL(wp) :: Alpha, Beta, Gamma
    REAL(wp) :: Vinf, HFlight

    ! --- Solver ---
    REAL(wp) :: dt, Epsilon
    INTEGER  :: MaxIter
    REAL(wp) :: TotalTime
    REAL(wp) :: Time


    ! Additional Non printable settings can go here
    INTEGER :: Iteration = 0
  END TYPE SimSettings


CONTAINS

  !-------------------------------------------------------------------------
  ! Read Config File
  !-------------------------------------------------------------------------
  SUBROUTINE Read_Settings_File(Filename, Cfg)
    CHARACTER(LEN=*), INTENT(IN) :: Filename
    TYPE(SimSettings), INTENT(OUT) :: Cfg
    INTEGER :: u, io_stat
    CHARACTER(LEN=256) :: io_msg ! To store the specific error message
    
    ! Local vars for Namelist
    CHARACTER(LEN=64) :: DatFile, PanFile
    INTEGER :: Symmetry, MaxIter
    REAL(wp) :: Alpha, Beta, Gamma, Vinf, HFlight, dt, Epsilon, MAC, RefArea, WingSpan
    REAL(wp) :: CG(3)
    INTEGER :: GroundEffect
    
    NAMELIST /CONFIG/ DatFile, PanFile, Symmetry, RefArea, WingSpan, MAC, CG, &
                      Alpha, Beta, Gamma, Vinf, HFlight, dt, MaxIter, Epsilon

    ! Defaults
    DatFile="DELTA.DAT"; PanFile="DELTA.PAN"
    Symmetry=0; GroundEffect=0
    RefArea=1.0_wp; WingSpan=1.0_wp; MAC=1.0_wp; CG=0.0_wp
    Alpha=0.0_wp; Beta=0.0_wp; Gamma=0.0_wp; Vinf=100.0_wp
    HFlight=0.0_wp; dt=0.1_wp; MaxIter=100; Epsilon=1.0E-5_wp

    OPEN(NEWUNIT=u, FILE=Filename, STATUS='OLD', ACTION='READ', IOSTAT=io_stat)
    IF (io_stat /= 0) STOP "Error opening config file."
    
    ! --- CRITICAL FIX: Check IOSTAT and IOMSG ---
    READ(u, NML=CONFIG, IOSTAT=io_stat, IOMSG=io_msg)
    
    IF (io_stat /= 0) THEN
       CALL Write_Log(" ")
       CALL Write_Log("==========================================================")
       CALL Write_Log("           CRITICAL ERROR READING CONFIGURATION           ")
       CALL Write_Log("==========================================================")
       CALL Write_Log(" The file '" // TRIM(Filename) // "' could not be read.")
       CALL Write_Log(" ")
       CALL Write_Log(" Fortran Error Message:")
       CALL Write_Log(" >> " // TRIM(io_msg))
       CALL Write_Log(" ")
       CALL Write_Log(" Common Causes:")
       CALL Write_Log(" 1. Typo in variable name or value.")
       CALL Write_Log(" 2. Missing comma between values.")
       CALL Write_Log("==========================================================")
       STOP
    END IF
    
    CLOSE(u)

    ! Map to Object
    Cfg%DatFile=DatFile; Cfg%PanFile=PanFile

    Cfg%Symmetry=Symmetry

    Cfg%RefArea=RefArea; Cfg%WingSpan=WingSpan; Cfg%MAC=MAC; Cfg%CG=CG

    Cfg%Alpha=Alpha; Cfg%Beta=Beta; Cfg%Gamma=Gamma

    Cfg%Vinf=Vinf
    Cfg%HFlight=HFlight

    Cfg%dt=dt
    Cfg%MaxIter=MaxIter
    Cfg%TotalTime=dt*REAL(MaxIter,wp)
    Cfg%Time=0.0_wp

    Cfg%Epsilon=Epsilon

    ! Calculations
    IF (Cfg%RefArea > 1.0E-6_wp) THEN
        Cfg%AspectRatio = (Cfg%WingSpan**2) / Cfg%RefArea
    ELSE
        Cfg%AspectRatio = 0.0_wp
    END IF

    if (Cfg%HFlight < Cfg%WingSpan * 0.5_wp) then
      Cfg%GroundEffect = 1
    else
      Cfg%GroundEffect = 0
    end if

    CALL Cfg%Air%update(Cfg%HFlight)
    
  END SUBROUTINE Read_Settings_File

  !-------------------------------------------------------------------------
  ! Wizard Entry Point
  !-------------------------------------------------------------------------
  SUBROUTINE Run_Wizard_And_Save(Filename, Cfg)
    CHARACTER(LEN=*), INTENT(IN) :: Filename
    TYPE(SimSettings), INTENT(OUT) :: Cfg
    
    CALL Write_Log("-------------------------------------------------------")
    CALL Write_Log(" LAUNCHING SETUP WIZARD")
    CALL Write_Log("-------------------------------------------------------")
    
    CALL Interactive_Setup(Cfg)
    CALL Save_Config(Filename, Cfg)
    CALL Cfg%Air%update(Cfg%HFlight)
    
    PRINT *, "Settings saved to '", TRIM(Filename), "'."
  END SUBROUTINE Run_Wizard_And_Save

  !-------------------------------------------------------------------------
  ! Internal: Interactive Inputs
  !-------------------------------------------------------------------------
  SUBROUTINE Interactive_Setup(Cfg)
    TYPE(SimSettings), INTENT(OUT) :: Cfg
    
    Cfg%MaxIter = 100; Cfg%Epsilon = 1.0E-5_wp; Cfg%CG = 0.0_wp

    PRINT *, ">> Geometry Nodes Filename (.DAT) [Default: DELTA.DAT]:"
    READ(*, '(A)') Cfg%DatFile
    IF (TRIM(Cfg%DatFile) == "") Cfg%DatFile = "DELTA.DAT"

    PRINT *, ">> Geometry Panels Filename (.PAN) [Default: DELTA.PAN]:"
    READ(*, '(A)') Cfg%PanFile
    IF (TRIM(Cfg%PanFile) == "") Cfg%PanFile = "DELTA.PAN"
    
    PRINT *, ">> Is the Mesh Symmetric? (1=Yes, 0=No):"
    READ(*, *) Cfg%Symmetry

    PRINT *, ">> Wing Reference Area (S) [m^2]:"
    READ(*, *) Cfg%RefArea
    PRINT *, ">> Wing Span (b) [m]:"
    READ(*, *) Cfg%WingSpan

    IF (Cfg%RefArea > 0.0_wp) THEN
        Cfg%AspectRatio = (Cfg%WingSpan**2) / Cfg%RefArea
        PRINT *, "   -> Calculated Aspect Ratio:", Cfg%AspectRatio
    END IF

    PRINT *, ">> Mean Aerodynamic Chord (MAC) [m]:"
    READ(*, *) Cfg%MAC
    PRINT *, ">> Center of Gravity (CG) X Y Z [m]:"
    READ(*, *) Cfg%CG(1), Cfg%CG(2), Cfg%CG(3)
    !  yaw, pitch, and roll angles are α, β and γ
    PRINT *, ">> Yaw (Alpha), Pitch (Beta), Roll (Gamma) [deg]:"
    READ(*, *) Cfg%Alpha, Cfg%Beta, Cfg%Gamma
    
    PRINT *, ">> Vinf [m/s], HFlight [m]:"
    READ(*, *) Cfg%Vinf, Cfg%HFlight

    PRINT *, ">> dt [s]:"
    READ(*, *) Cfg%dt

    PRINT *, ">> MaxIter:"
    READ(*, *) Cfg%MaxIter
        PRINT *, "   -> Calculated Simulated Time [s]", Cfg%dt * REAL(Cfg%MaxIter, wp)

    PRINT *, ">> Epsilon:"
    READ(*, *) Cfg%Epsilon

  END SUBROUTINE Interactive_Setup

  !-------------------------------------------------------------------------
  ! Internal: Save to File
  !-------------------------------------------------------------------------
  SUBROUTINE Save_Config(Filename, Cfg)
    CHARACTER(LEN=*), INTENT(IN) :: Filename
    TYPE(SimSettings), INTENT(IN) :: Cfg
    INTEGER :: u
    
    OPEN(NEWUNIT=u, FILE=Filename, STATUS='UNKNOWN', ACTION='WRITE')
    WRITE(u, '(A)') "&CONFIG"
    WRITE(u, '(A)') " !--- Files ---"
    WRITE(u, '(A,A,A)')     " DatFile   = '", TRIM(Cfg%DatFile), "',"
    WRITE(u, '(A,A,A)')     " PanFile   = '", TRIM(Cfg%PanFile), "',"
    WRITE(u, '(A)') " !--- Geometry ---"
    WRITE(u, '(A,I0,A)')    " Symmetry = ", Cfg%Symmetry, ","
    WRITE(u, '(A,F10.3,A)') " RefArea   = ", Cfg%RefArea,   ", ! S"
    WRITE(u, '(A,F10.3,A)') " WingSpan  = ", Cfg%WingSpan,  ", ! b"
    WRITE(u, '(A,F10.3,A)') " MAC       = ", Cfg%MAC,       ", ! Mean Aerodynamic Chord [m]"
    WRITE(u, '(A,3(F10.3,1X),A)') " CG        = ", Cfg%CG, ","
    WRITE(u, '(A)') " !--- Flow ---"
    WRITE(u, '(A,F10.3,A)') " Alpha     = ", Cfg%Alpha,   ", ! Yaw [deg]"
    WRITE(u, '(A,F10.3,A)') " Beta      = ", Cfg%Beta,    ", ! Pitch [deg]"
    WRITE(u, '(A,F10.3,A)') " Gamma     = ", Cfg%Gamma,   ", ! Roll [deg]"
    WRITE(u, '(A,F10.3,A)') " Vinf      = ", Cfg%Vinf,    ", ! Velocity [m/s]"
    WRITE(u, '(A,F10.3,A)') " HFlight   = ", Cfg%HFlight, ", ! Height [m]"
    WRITE(u, '(A)') " !--- Solver ---"
    WRITE(u, '(A,F10.6,A)') " dt        = ", Cfg%dt,      ", ! Time Step [s]"
    WRITE(u, '(A,I0,A)')    " MaxIter   = ", Cfg%MaxIter, ", ! Maximum Iterations"
    WRITE(u, '(A,E12.4,A)') " Epsilon   = ", Cfg%Epsilon, ", ! Convergence Criterion"
    WRITE(u, '(A)') "/"
    CLOSE(u)
  END SUBROUTINE Save_Config

  !-------------------------------------------------------------------------
  ! Print the Full Configuration to Console
  !-------------------------------------------------------------------------
  SUBROUTINE Display_Config(Cfg)
    TYPE(SimSettings), INTENT(IN) :: Cfg
    CHARACTER(LEN=256) :: Line

    CALL Write_Log(" ")
    CALL Write_Log("==========================================================")
    CALL Write_Log("               SIMULATION CONFIGURATION                   ")
    CALL Write_Log("==========================================================")
    
    CALL Write_Log(" [FILES]")
    ! A30 ensures the label takes exactly 30 spaces. F15.4 takes 15 spaces.
    WRITE(Line, '(A30, A)') "   Nodes File:",    TRIM(Cfg%DatFile)
    CALL Write_Log(Line)
    WRITE(Line, '(A30, A)') "   Panels File:",   TRIM(Cfg%PanFile)
    CALL Write_Log(Line)
    CALL Write_Log(" ")
    
    CALL Write_Log(" [MESH STATS]")
    WRITE(Line, '(A30, I8)') "   Grid Points:",   Cfg%NGrid
    CALL Write_Log(Line)
    WRITE(Line, '(A30, I8)') "   Panels:",        Cfg%NPanels
    CALL Write_Log(Line)
    WRITE(Line, '(A30, A)')  "   Symmetry:",      MERGE("YES", "NO ", Cfg%Symmetry == 1)
    CALL Write_Log(Line)
    CALL Write_Log(" ")

    CALL Write_Log(" [GEOMETRY]")
    WRITE(Line, '(A30, F15.4, 1X, A)') "   Ref Area (S):",  Cfg%RefArea,   "m^2"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, F15.4, 1X, A)') "   Span (b):",      Cfg%WingSpan,  "m"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, F15.4)')        "   Aspect Ratio:",  Cfg%AspectRatio
    CALL Write_Log(Line)
    WRITE(Line, '(A30, F15.4, 1X, A)') "   MAC:",           Cfg%MAC,       "m"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, 3F10.4, 1X, A)')"   CG (x,y,z):",    Cfg%CG,        "m"
    CALL Write_Log(Line)
    CALL Write_Log(" ")

    CALL Write_Log(" [FLOW CONDITIONS]")
    WRITE(Line, '(A30, F15.3, 1X, A)') "   Velocity:",      Cfg%Vinf,    "m/s"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, F15.3, 1X, A)') "   Height:",        Cfg%HFlight, "m"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, F15.3, 1X, A)') "   Alpha:",         Cfg%Alpha,   "deg"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, F15.3, 1X, A)') "   Beta:",          Cfg%Beta,    "deg"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, F15.3, 1X, A)') "   Gamma:",         Cfg%Gamma,   "deg"
    CALL Write_Log(Line)
    
    CALL Cfg%Air%display()

    CALL Write_Log(" [SOLVER]")
    WRITE(Line, '(A30, F15.5, 1X, A)') "   Time Step:",     Cfg%dt,      "s"
    CALL Write_Log(Line)
    WRITE(Line, '(A30, I8)')           "   Max Iter:",      Cfg%MaxIter
    CALL Write_Log(Line)
    WRITE(Line, '(A30, ES15.4)')       "   Convergence:",   Cfg%Epsilon
    CALL Write_Log(Line)
    
    CALL Write_Log("==========================================================")
    CALL Write_Log(" ")
  END SUBROUTINE Display_Config

END MODULE config_io