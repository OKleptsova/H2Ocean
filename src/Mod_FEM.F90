MODULE Mod_FEM
!------------------------------------------------------------------------
!
!  Copyright (C) 2009-2011 Technische Universiteit Delft,
!  Haiyang Cui, Guus Stelling and Julie Pietrzak
!  
!  
!  Programmers:
!  Haiyang Cui
!  PhD student
!  Environmental Fluid Mechanics section
!  Faculty of Civil Engineering and Geosciences
!  Delft University of Technology
!  Delft, The Netherlands
!  Tel:   +31 15 278 5433
!  Email: H.Cui@tudelft.nl   cuiocean@gmail.com
!	
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA
!
!---Authors--------------------
!
!   1.00: Haiyang Cui
!------------------------------

!---Updates--------------------
!
!------------------------------

!---Purpose--------------------
!
!------------------------------
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_Parameter
  IMPLICIT NONE
  ! Finite Element 
  REAL(kind=REAL_DP)         :: Vol2D
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:) :: derivative_zeta_stdbf, derivative_vel_stdbf
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:) :: bafuzeta_x, bafuzeta_y
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:) :: bafuvel_x, bafuvel_y
  CONTAINS

  SUBROUTINE FEM_Standard_Element_Definition
!---------------------------------------------------------------
!  DEFINITION OF :
! - 2D STANDARD ELEMENTS
! - BASISFUNCTIONS ON 2D STANDARD ELEMENTS (stdbafu)
!  Zeta : stdbafu(1)=1-x-y   stdbafu(2)=x    stdbafu(3)=y
!  Velocity: stdbafu(1)=1-2y stdbafu(2)=2x+2y-1 stdbafu(3)=1-2x
! - SCALARPRODUCTS
! - DERIVATIVES OF STD.BASISFUNCTIONS!
!---------------------------------------------------------------
    IMPLICIT NONE
    Vol2D = 1./2.
 
    ALLOCATE(derivative_vel_stdbf(2, 3), derivative_zeta_stdbf(2,3))
 
    derivative_zeta_stdbf      =  0.
    derivative_zeta_stdbf(:,1) = -1.
    derivative_zeta_stdbf(1,2) =  1.
    derivative_zeta_stdbf(2,3) =  1.

    derivative_vel_stdbf       =  0.
    derivative_vel_stdbf(2,1)  = -2.
    derivative_vel_stdbf(:,2)  =  2.
    derivative_vel_stdbf(1,3)  = -2.
  END SUBROUTINE FEM_Standard_Element_Definition
  
  SUBROUTINE FEM_Basis_Functions 
    IMPLICIT NONE
    REAL(kind=REAL_DP)                     :: DET2D
    REAL(kind=REAL_DP), DIMENSION(2, 3)    :: derivative_loczeta, derivative_locvel
    REAL(kind=REAL_DP), DIMENSION(2,2)     :: jacobian2D, jacobian2D_inv
    INTEGER        :: elem, i

!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 
  
    ALLOCATE(bafuzeta_x(3,TriMesh_elementnmb), bafuzeta_y(3,TriMesh_elementnmb))
    ALLOCATE(bafuvel_x(3,TriMesh_elementnmb), bafuvel_y(3,TriMesh_elementnmb))
 
    !
    bafuzeta_x  = 0.0d0
    bafuzeta_y  = 0.0d0
    bafuvel_x   = 0.0d0
    bafuvel_y   = 0.0d0
    
!$OMP   PARALLEL PRIVATE(elem,i,jacobian2D,jacobian2D_inv,DET2D,derivative_loczeta,derivative_locvel,mythread)  &
!$OMP&           SHARED(TriMesh_elementnmb,bafuzeta_x,bafuzeta_y,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_elementnmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    
    DO elem=1,TriMesh_elementnmb
        CALL Local_Element_Def(elem,   &
        jacobian2D, jacobian2D_inv, DET2D,  derivative_loczeta, derivative_locvel)
        DO i=1,3
            bafuzeta_x(i,elem) = derivative_loczeta(1,i)
            bafuzeta_y(i,elem) = derivative_loczeta(2,i)
            bafuvel_x(i,elem)  = derivative_locvel(1,i)
            bafuvel_y(i,elem)  = derivative_locvel(2,i)
        ENDDO
    ENDDO ! end for elem
    
!$OMP END DO   
!$OMP END PARALLEL     
  END SUBROUTINE FEM_Basis_Functions
  


!----------------------------------------------------
! subroutine local_element_def
!
! 2D LOCAL ELEMENT DEFINITIONS
!
!  INPUT:
!   derivative_stabafu_x_2D(2, 3)
!
!  OUTPUT:
!   jacobian2D(2,2), jacobian2D_inv(2,2), DET
!   derivative_locbafu_x_2D(2, 3)
!
!  EXTERNAL SUBROUTINE: call matrix_inverse_2x2
!----------------------------------------------------
SUBROUTINE Local_Element_Def(elem, &
    jacobian2D, jacobian2D_inv, DET, derivative_loczeta, derivative_locvel) 
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)               :: elem 

      
    INTEGER          :: i, nd
    REAL(kind=REAL_DP), DIMENSION(2,2),  INTENT(OUT)  :: jacobian2D
    REAL(kind=REAL_DP), DIMENSION(2,2),  INTENT(OUT)  :: jacobian2D_inv
    REAL(kind=REAL_DP), intent(OUT)       :: DET
    REAL(kind=REAL_DP), DIMENSION(2, 3), INTENT(OUT)  :: derivative_loczeta
    REAL(kind=REAL_DP), DIMENSION(2, 3), INTENT(OUT)  :: derivative_locvel
    REAL(kind=REAL_DP), DIMENSION(2, 3)      :: local_cart
    REAL(kind=REAL_DP), DIMENSION(3, 2)      :: der_transp
    REAL(kind=REAL_DP)        :: meancos
    

!        print *, elem, element(elem)%nodeinds(1:3)
    DO i=1, 3
        nd=element(elem)%nodeinds(i)
        !
        ! cartesian coordinates
        IF( Coordinate_type == 1) THEN
            local_cart(1,i)= node(nd)%x
            local_cart(2,i)= node(nd)%y        
        ELSEIF( Coordinate_type == 2) THEN
            local_cart(1,i)=r_earth * node(nd)%x
            local_cart(2,i)=r_earth * node(nd)%y
        ENDIF
    ENDDO
    !
    ! TRANSFORMATION - MATRIX "jacobian"
    !
    DO i=1,2
        jacobian2D(1,i) = local_cart(1,i+1) - local_cart(1,1)
        IF(Coordinate_type ==2) THEN
            meancos = cos(0.5d0*(local_cart(2,i+1)/r_earth+local_cart(2,1)/r_earth))
            jacobian2D(1,i) = jacobian2D(1,i) * meancos
        ENDIF
        jacobian2D(2,i) = local_cart(2,i+1) - local_cart(2,1)
    ENDDO
    !
    !  INVERSE OF jacobian
    !
    CALL Matrix_Inverse_2x2(jacobian2D, jacobian2D_inv, DET)
    der_transp = MATMUL(TRANSPOSE(derivative_zeta_stdbf), jacobian2D_inv)
    derivative_loczeta = TRANSPOSE(der_transp)
    der_transp = MATMUL(TRANSPOSE(derivative_vel_stdbf), jacobian2D_inv)
    derivative_locvel = TRANSPOSE(der_transp)
END SUBROUTINE Local_Element_Def
!
!----------------------------------------------------------------------------
!
SUBROUTINE  Matrix_Inverse_2x2 (A, AINV, DET)
    !
    !  * * * * * * * * * * * * * * * * * * * * * * * * * *
    !  CALCULATE THE DETERMINATE AND INVERSE OF A(2,2)
    !  * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !  A = ORIGINAL MATRIX
    !  AINV = INVERSE OF MATRIX A
    !  DET  = DETERMINANT OF A
    !
    USE Mod_Precision
    IMPLICIT NONE
    !
    REAL(kind=REAL_DP), DIMENSION(2,2), INTENT(IN)  :: A
    REAL(kind=REAL_DP), DIMENSION(2,2), INTENT(OUT) :: AINV
    REAL(kind=REAL_DP), INTENT(OUT)     :: DET
    !
    INTEGER           :: i,j
    !
    DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    IF ( DET .eq. 0.0 )  THEN
        DO j=1,2
            WRITE(*,*) (A(i,j),i=1,2)
        END DO
        WRITE(*,*) 'SINGULAR 2X2 MATRIX'
        STOP
    ELSE
        AINV(1,1) =  A(2,2)/DET
        AINV(1,2) = -A(1,2)/DET
        AINV(2,1) = -A(2,1)/DET
        AINV(2,2) =  A(1,1)/DET
    ENDIF
END SUBROUTINE Matrix_Inverse_2x2


END MODULE Mod_FEM
