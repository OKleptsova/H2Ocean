SUBROUTINE H2Ocean_Finalize
  USE Mod_All_Variable
  IMPLICIT NONE
  INTEGER :: i
  
  IF(enable_station_output) THEN
      CLOSE(station_fid)    
  ENDIF
  
  
  ! Close time file
  CLOSE(fid_time)
  
  
  IF(enable_zeta_max_output) THEN
      CALL H2Ocean_Zeta_Max_Output_Once
  ENDIF

  IF(enable_altim_output) THEN
      CALL H2Ocean_altim_Output_Once
  ENDIF
END SUBROUTINE H2Ocean_Finalize