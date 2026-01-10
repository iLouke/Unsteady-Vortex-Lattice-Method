PROGRAM ROS_Main
  USE uvlm
  IMPLICIT NONE

  CALL Initialize_System()
  CALL Load_Geometry()
  
  CALL Run_Solver()

  CALL Finalize_System()
END PROGRAM ROS_Main