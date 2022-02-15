MODULE Mod_TriMesh_Interp
  USE Mod_Precision
  USE Mod_TriMesh
  IMPLICIT NONE
  
  CONTAINS
  
   SUBROUTINE Point_Nearest_Element(x0,y0,el_out)
    IMPLICIT NONE
    REAL(KIND=REAL_DP), INTENT(IN)       :: x0, y0
    INTEGER(KIND=INT_KIND), INTENT(OUT)  :: el_out
    REAL(KIND=REAL_DP),ALLOCATABLE, DIMENSION(:) ::  distance 
    INTEGER(KIND=INT_KIND) ::  el,min_pos(1) , j      
    el_out =0
    ALLOCATE(distance(TriMesh_elementnmb))
    DO el = 1,TriMesh_elementnmb
        distance(el) = (x0-element(el)%x)**2 + (y0-element(el)%y)**2
    ENDDO
    min_pos = MINLOC(distance)
    el_out = min_pos(1)   
    DEALLOCATE( distance)
  END SUBROUTINE Point_Nearest_Element


  SUBROUTINE TriMesh_Interp_Element(xin,yin,el,z_in,z_out)
    IMPLICIT NONE
    REAL(KIND=REAL_DP),      INTENT(IN)      :: xin,yin
    INTEGER(KIND=INT_KIND),  INTENT(IN)      :: el 
    REAL(KIND=REAL_DP),      INTENT(IN)      :: z_in(TriMesh_nodenmb)   
    REAL(KIND=REAL_DP),      INTENT(OUT)     :: z_out
    REAL(KIND=REAL_DP)                        :: xd(3),yd(3),zd(3)
    REAL(KIND=REAL_DP)                        :: A,dzdx, dzdy
 

    xd =  node(element(el)%nodeinds(:)) %x 
    yd =  node(element(el)%nodeinds(:)) %y 
    zd =  z_in(element(el)%nodeinds(:))  
    
    A = 0.5d0*(xd(1)*yd(2)+xd(2)*yd(3)+xd(3)*yd(1) &
              -xd(1)*yd(3)-xd(2)*yd(1)-xd(3)*yd(2))
    dzdx = 1./(2.0d0*A) *( (yd(2)-yd(3))*zd(1) +(yd(3)-yd(1))*zd(2) &
                          +(yd(1)-yd(2))*zd(3)  )
    dzdy = 1./(2.0d0*A) *( (xd(3)-xd(2))*zd(1) +(xd(1)-xd(3))*zd(2) &
                          +(xd(2)-xd(1))*zd(3)  )
                         
    z_out = zd(1) + dzdx*(xin-xd(1))  + dzdy*(yin-yd(1)) 
    
                    
  END SUBROUTINE TriMesh_Interp_Element 
  
  
END MODULE Mod_TriMesh_Interp


FUNCTION isPointInsideTriangle(x0,y0,x1,y1,x2,y2,x3,y3)
  USE Mod_Precision
  IMPLICIT NONE
  REAL(kind=REAL_DP), INTENT(IN):: x0,y0,x1,y1,x2,y2,x3,y3
  LOGICAL :: isPointInsideTriangle
  LOGICAL,EXTERNAL  :: isOnSameSide
  
  IF(       isOnSameSide(x0,y0,x1,y1,x2,y2,x3,y3) &
     .AND.  isOnSameSide(x0,y0,x2,y2,x3,y3,x1,y1) & 
     .AND.  isOnSameSide(x0,y0,x3,y3,x1,y1,x2,y2) ) THEN
      isPointInsideTriangle = .TRUE.
  ELSE
      isPointInsideTriangle = .FALSE.
  ENDIF    
END FUNCTION


FUNCTION isOnSameSide(x0,y0,x1,y1,x2,y2,x3,y3)
  USE Mod_Precision
  IMPLICIT NONE
  REAL(kind=REAL_DP), INTENT(IN):: x0,y0,x1,y1,x2,y2,x3,y3
  LOGICAL :: isOnSameSide
  REAL(kind=REAL_DP)  :: a, b, c
    a = y0 - y1; 
    b = x1 - x0; 
    c = x0 * y1 - x1 * y0;
    IF( (a * x2 + b * y2 + c) * (a * x3 + b * y3 + c) <= 0) THEN
       isOnSameSide = .TRUE.
    ELSE
       isOnSameSide = .FALSE.
    ENDIF
END FUNCTION isOnSameSide