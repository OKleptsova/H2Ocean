!------------------------------------------------------------------------
! Main program for the unstructured finite volume model, H2Ocean.
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



MODULE Mod_Sparse_Matrix
 USE Mod_Precision
 IMPLICIT NONE
 SAVE
 
!Description:
!  For sparse matrix,  csr format
!  dim     -----  the number of unknows, e.g. the number of variables
!  nza     -----  the nubmer of entries in the sparse matrix specified in compressed-row-storage(CSR) by rowptr, colind, values
!  values  -----  array of matrix entry values in the same order as column indices colind
!  rowptr  -----  array of row-POINTERs, rowptr(1)=1; rowptr(i+1)-rowptr(i)=number of entries in row i
!  colind  -----  array of column indices(Fortran indexing, lowest indes 1, highest index dim) of matrix entries
!----------------------------------------
  TYPE sparse_matrix
      INTEGER  :: nza     
      INTEGER  :: dim
      REAL(kind=REAL_DP) ,   POINTER, DIMENSION(:)      :: values
      INTEGER(KIND=INT_KIND), POINTER, DIMENSION(:)     :: colind
      INTEGER(KIND=INT_KIND), POINTER, DIMENSION(:)     :: rowptr
 END TYPE sparse_matrix
 
 TYPE(sparse_matrix)                                    :: A_s
 REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)        :: RHS_s

 
 CONTAINS
 
 SUBROUTINE sparse_matrix_locate(nrow,ncol,A_in, offset)
   ! Find the locate of entry (nrow,ncol) in the sparse matrix
   IMPLICIT NONE
   INTEGER, INTENT(IN)                :: nrow, ncol     ! nod: the center node, ndin: the node needed to be located the position
   TYPE(sparse_matrix), INTENT(IN)    :: A_in
   INTEGER, INTENT(OUT)               :: offset             ! the postion in vector q_neighbour_nodes(nod)%addresses(i)
   INTEGER                            :: i,pos
   
   offset = 0 
   DO i = A_in%rowptr(nrow), A_in%rowptr(nrow+1)-1
      IF( ncol == A_in%colind(i)) THEN
          offset = i
          EXIT     ! as soon as the postion is found, jump out the loop, more efficient
      END IF     
   ENDDO
   
   IF(offset == 0) THEN
       WRITE(*,*) "There is no such entry (" , nrow, ",", ncol,") in this sparse matrix."
   ENDIF    
 END SUBROUTINE sparse_matrix_locate


! function adj_bandwidth ( node_num, adj_num, adj_row, adj )
FUNCTION adj_bandwidth ( A_in)
!*****************************************************************************80
!
!! ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about 
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) ADJ_BANDWIDTH, the bandwidth of the adjacency
!    matrix.
!
  IMPLICIT NONE
  TYPE(sparse_matrix), INTENT(IN)    :: A_in
 
  INTEGER(KIND=INT_KIND):: band_hi,band_lo
  INTEGER(KIND=INT_KIND):: col,i,j 
  INTEGER(KIND=INT_KIND):: adj_bandwidth
  
  band_lo = 0
  band_hi = 0

  DO i = 1, A_in%dim

    DO j = A_in%rowptr(i), A_in%rowptr(i+1)-1   
      col = A_in%colind(j)
      band_lo = MAX ( band_lo, i - col )
      band_hi = max ( band_hi, col - i )
    END DO

  END DO
  adj_bandwidth = band_lo + 1 + band_hi
  RETURN  
 END FUNCTION

 SUBROUTINE sparse_matrix_display_bandwidth(A_in)
   IMPLICIT NONE
   TYPE(sparse_matrix), INTENT(IN)    :: A_in
   INTEGER  :: bandwidth
   bandwidth = adj_bandwidth (A_in)
   WRITE(*,*) "The bandwith of this sparse matrix is: ", bandwidth
   
 END SUBROUTINE sparse_matrix_display_bandwidth

 
 
SUBROUTINE SparseSolver(A_in,x,rhs,maxits,iter_tol1,iter_tol2,method)
  IMPLICIT NONE
  TYPE(sparse_matrix), INTENT(IN)      :: A_in 
  REAL(KIND=REAL_DP) , INTENT(IN)  :: rhs(:)
  REAL(KIND=REAL_DP) , INTENT(INOUT) :: x(:)
  INTEGER, INTENT(IN)    :: maxits
  REAL(KIND=REAL_DP),  INTENT(IN)    :: iter_tol1,iter_tol2
  CHARACTER*(*),  INTENT(IN)    :: method

  INTEGER(KIND=INT_KIND)    :: n,nnz
 
  ! Other parameters
 
  INTEGER  lwk,nwk
  INTEGER,ALLOCATABLE,DIMENSION(:)  ::  jau,ju,iw
  INTEGER ipar(16),lfil,ierr 
  REAL*8,ALLOCATABLE,DIMENSION(:)   ::  au,wk
  REAL*8,ALLOCATABLE,DIMENSION(:)   ::  xran, fpar, al
  REAL*8  tol

  INTEGER :: i,j

  external cg,bcg,dbcg,bcgstab,tfqmr,gmres,fgmres,dqgmres
  external cgnr, fom, runrc, ilut,ilud 


  n=A_in%dim
  nnz = A_in%nza
  

  ALLOCATE(jau(nnz*2),ju(nnz*2),iw(n*3))
  ALLOCATE(au(nnz*2),wk(n*40))
  ALLOCATE(xran(n+1), fpar(16), al(n+1))

  lwk=n*40 
  ipar(2) = 2
  !ipar(3) = 1
  ipar(3) = 2
  ipar(4) = lwk
  ipar(5) = 10 
  ipar(6) = maxits
  fpar(1) = iter_tol1
  fpar(2) = iter_tol2

  lfil = 3
  tol = iter_tol1
  nwk = nnz*3
  
  CALL ilut(n,A_in%values,A_in%colind,A_in%rowptr,lfil,tol,au,jau,ju,nwk,wk,iw,ierr)



  ipar(2) = 2
  xran = x;
  SELECT CASE(TRIM(ADJUSTL(method))) 
  CASE('cg')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,cg)
  CASE('bcg')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,bcg)  
  CASE('dbcg')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,dbcg)
  CASE('bcgstab')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,bcgstab)
  CASE('tfqmr')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,tfqmr)
  CASE('gmres')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,gmres)
  CASE('fgmres')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,fgmres)
  CASE('dqgmres')
    CALL runrc(n,rhs,x,ipar,fpar,wk,xran,A_in%values,A_in%colind,A_in%rowptr,au,jau,ju,dqgmres)
  CASE DEFAULT
      WRITE(*,*) "At this moment, the method ", TRIM(ADJUSTL(method)), ' is not supported'
  END SELECT  
  
 
  DEALLOCATE(jau,ju,iw)
  DEALLOCATE(au,wk)
  DEALLOCATE(xran, fpar, al)
END SUBROUTINE SparseSolver

 
END MODULE Mod_Sparse_Matrix


 