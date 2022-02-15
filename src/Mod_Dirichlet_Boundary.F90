MODULE Mod_Dirichlet_Boundary 
  USE Mod_Precision
  
  IMPLICIT NONE
  
  LOGICAL                        :: dirichlet_boundary = .FALSE.
  CHARACTER(len=400)             :: Dirichlet_file = ' '
  
  LOGICAL                        :: q_dirichlet_boundary = .FALSE.
  CHARACTER(len=400)             :: q_Dirichlet_file = ' '  
  
  LOGICAL                        :: u_dirichlet_boundary = .FALSE.
  CHARACTER(len=400)             :: u_Dirichlet_file = ' '

  REAL(kind=REAL_DP)             :: D_bo
  REAL(kind=REAL_DP)             :: q_D_bo, u_D_bo
  REAL(kind=REAL_DP),PARAMETER   :: D_bo_Inf = 999.d0
  
  
  !Store the Direchlet boundary
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:)       ::  D_boundary
    
  !Store the u Direchlet boundary
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:)       ::  u_D_boundary    

 !Store the q Direchlet boundary
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:)       ::  q_D_boundary  
  
    
  CONTAINS
  SUBROUTINE H2Ocean_Read_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nodenmb,n
    INTEGER(KIND=INT_KIND)  :: file_Dirichlet=100
    CHARACTER(len=400)      :: filename
    filename = TRIM(Dirichlet_file)  
    OPEN(unit=file_Dirichlet,file=filename,status="old",action='read') 
    READ(file_Dirichlet,*) nodenmb
    ALLOCATE(D_boundary(nodenmb,2))
    DO n=1,nodenmb
      READ(file_Dirichlet,*)D_boundary(n,:)
    ENDDO
    CLOSE(100)
  END SUBROUTINE H2Ocean_Read_Dirichlet_Boundary


  SUBROUTINE H2Ocean_Read_q_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nodenmb,n
    INTEGER(KIND=INT_KIND)  :: file_Dirichlet=100
    CHARACTER(len=400)      :: filename
    filename = TRIM(q_Dirichlet_file)  
    OPEN(unit=file_Dirichlet,file=filename,status="old",action='read') 
    READ(file_Dirichlet,*) nodenmb
    ALLOCATE(q_D_boundary(nodenmb,2))
    DO n=1,nodenmb
      READ(file_Dirichlet,*)q_D_boundary(n,:)
    ENDDO
    CLOSE(100)
  END SUBROUTINE H2Ocean_Read_q_Dirichlet_Boundary



  SUBROUTINE H2Ocean_Read_u_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nodenmb,n
    INTEGER(KIND=INT_KIND)  :: file_Dirichlet=100
    CHARACTER(len=400)      :: filename
    filename = TRIM(u_Dirichlet_file)  
    OPEN(unit=file_Dirichlet,file=filename,status="old",action='read') 
    READ(file_Dirichlet,*) nodenmb
    ALLOCATE(u_D_boundary(nodenmb,2))
    DO n=1,nodenmb
      READ(file_Dirichlet,*)u_D_boundary(n,:)
    ENDDO
    CLOSE(100)
  END SUBROUTINE H2Ocean_Read_u_Dirichlet_Boundary


  
  SUBROUTINE H2Ocean_Interp_D_bo(time_in) 
  IMPLICIT NONE
  REAL(kind=REAL_DP), INTENT(in)  :: time_in
  INTEGER                         :: HBD,n
  REAL(kind=REAL_DP)              :: t_low,t_high,D_low,D_high
  
  
  HBD =  SIZE(D_boundary,1)
  D_bo = D_bo_Inf + 1.0d0;
  IF(time_in>= MINVAL(D_boundary(:,1)) .and. time_in <= MAXVAL(D_boundary(:,1))) THEN
      DO n=1, HBD-1
          IF(time_in>= D_boundary(n,1) .and. time_in<= D_boundary(n+1,1)) THEN
              t_low  = D_boundary(n,1)
              t_high = D_boundary(n+1,1)
              D_low  = D_boundary(n,2)
              D_high = D_boundary(n+1,2)
              D_bo   = D_low + (D_high-D_low)/(t_high-t_low)*(time_in-t_low)
              EXIT
          ENDIF
      ENDDO
  ENDIF 

END SUBROUTINE H2Ocean_Interp_D_bo



SUBROUTINE H2Ocean_Interp_q_D_bo(time_in) 
  IMPLICIT NONE
  REAL(kind=REAL_DP), INTENT(in)  :: time_in
  INTEGER                         :: HBD,n
  REAL(kind=REAL_DP)              :: t_low,t_high,D_low,D_high
  
  
  HBD =  SIZE(q_D_boundary,1)
  q_D_bo = D_bo_Inf + 1.0d0;
  IF(time_in>= MINVAL(q_D_boundary(:,1)) .and. time_in <= MAXVAL(q_D_boundary(:,1))) THEN
      DO n=1, HBD-1
          IF(time_in>= q_D_boundary(n,1) .and. time_in<= q_D_boundary(n+1,1)) THEN
              t_low  = q_D_boundary(n,1)
              t_high = q_D_boundary(n+1,1)
              D_low  = q_D_boundary(n,2)
              D_high = q_D_boundary(n+1,2)
              q_D_bo = D_low + (D_high-D_low)/(t_high-t_low)*(time_in-t_low)
              EXIT
          ENDIF
      ENDDO
  ENDIF 

END SUBROUTINE H2Ocean_Interp_q_D_bo



SUBROUTINE H2Ocean_Interp_u_D_bo(time_in) 
  IMPLICIT NONE
  REAL(kind=REAL_DP), INTENT(in)  :: time_in
  INTEGER                         :: HBD,n
  REAL(kind=REAL_DP)              :: t_low,t_high,D_low,D_high
  
  
  HBD =  SIZE(u_D_boundary,1)
  u_D_bo = D_bo_Inf + 1.0d0;
  IF(time_in>= MINVAL(u_D_boundary(:,1)) .and. time_in <= MAXVAL(u_D_boundary(:,1))) THEN
      DO n=1, HBD-1
          IF(time_in>= u_D_boundary(n,1) .and. time_in<= u_D_boundary(n+1,1)) THEN
              t_low  = u_D_boundary(n,1)
              t_high = u_D_boundary(n+1,1)
              D_low  = u_D_boundary(n,2)
              D_high = u_D_boundary(n+1,2)
              u_D_bo = D_low + (D_high-D_low)/(t_high-t_low)*(time_in-t_low)
              EXIT
          ENDIF
      ENDDO
  ENDIF 

END SUBROUTINE H2Ocean_Interp_u_D_bo


  
END MODULE Mod_Dirichlet_Boundary