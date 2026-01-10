MODULE logger_io
  USE constants
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Init_Logging, Write_Log, Write_Loads, Finalize_IO, Get_Log_Unit

  ! State Flag: Simple True/False. No more checking unit numbers.
  LOGICAL :: System_Logs_Open = .FALSE.

CONTAINS

  SUBROUTINE Init_Logging()
    ! Open fixed units defined in constants.f90
    OPEN(UNIT=LU_LOG,   FILE='log.txt',   STATUS='REPLACE', ACTION='WRITE')
    OPEN(UNIT=LU_LOADS, FILE='loads.txt', STATUS='REPLACE', ACTION='WRITE')
    
    ! Write Header to Loads
    WRITE(LU_LOADS, '(A)') "# Iter      Lift           Drag           Side           Roll           Pitch          Yaw"
    
    ! Set Flag
    System_Logs_Open = .TRUE.
  END SUBROUTINE Init_Logging

  ! Helper for Atmosphere module
  FUNCTION Get_Log_Unit() RESULT(u)
    INTEGER :: u
    u = LU_LOG
  END FUNCTION Get_Log_Unit

  SUBROUTINE Write_Log(Message)
    CHARACTER(LEN=*), INTENT(IN) :: Message
    
    ! 1. Console Output
    WRITE(*, '(A)') TRIM(Message)
    
    ! 2. File Output (Only if initialized)
    IF (System_Logs_Open) THEN
        WRITE(LU_LOG, '(A)') TRIM(Message)
        FLUSH(LU_LOG) ! Force write to disk
    END IF
  END SUBROUTINE Write_Log

  SUBROUTINE Write_Loads(Iter, Forces, Moments)
    INTEGER, INTENT(IN) :: Iter
    REAL(wp), INTENT(IN) :: Forces(3), Moments(3)
    
    ! Console Echo (Brief)
    WRITE(*, '(A,I3, A,F8.2, A,F8.2)') "Iter:", Iter, " L:", Forces(1), " D:", Forces(2)
    
    IF (System_Logs_Open) THEN
        ! Full Data to loads.txt
        WRITE(LU_LOADS, '(I4, 6F15.5)') Iter, Forces, Moments
        FLUSH(LU_LOADS)

        ! Echo Brief to log.txt as well
        WRITE(LU_LOG, '(A,I3, A,F8.2, A,F8.2)') "Iter:", Iter, " L:", Forces(1), " D:", Forces(2)
        FLUSH(LU_LOG)
    END IF
  END SUBROUTINE Write_Loads

  SUBROUTINE Finalize_IO()
    IF (System_Logs_Open) THEN
        CLOSE(LU_LOG)
        CLOSE(LU_LOADS)
        System_Logs_Open = .FALSE.
    END IF
  END SUBROUTINE Finalize_IO

END MODULE logger_io