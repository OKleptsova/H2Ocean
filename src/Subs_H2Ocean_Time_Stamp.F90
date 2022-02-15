SUBROUTINE H2Ocean_Time_Stamp
  USE Mod_All_Variable
  IMPLICIT NONE
  ! Find number of seconds since 1/1/1970
  !walltime_begin = TIME()
  walltime_begin = secnds(0.0)
  ! Regsiter the begin cpu time of the Programme. 
  CALL CPU_TIME (cputime_begin)
END SUBROUTINE H2Ocean_Time_Stamp


SUBROUTINE H2Ocean_Time_Report
  USE Mod_All_Variable
  IMPLICIT NONE
 ! Find number of seconds since 1/1/70
  walltime_end= secnds( walltime_begin )  
  ! Regsiter the begin cpu time of the Programme. 
  CALL CPU_TIME ( cputime_end ) 
  WRITE (*,*) 'CPU Time of computation was ', cputime_end - cputime_begin, ' seconds'
  WRITE (*,*) 'The Wall Clock Time of computation was ', walltime_end, ' seconds'  
 
END SUBROUTINE H2Ocean_Time_Report
