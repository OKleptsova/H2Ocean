!------------------------------------------------------------------------ 
! Copyright (C) 2009-2011 Technische Universiteit Delft,
!  Haiyang Cui, Guus Stelling and Julie Pietrzak
!  --|-----------------------------------------------------------|--
!    |  H.Cui@tudelft.nl	                                     |
!    |  Haiyang Cui                                              |
!    |  Environmental Fluid Mechanics section                    |
!    |  Faculty of Civil Engineering and Geosciences             |
!    |  Delft University of Technology                           |
!    |  Room 2.92,Building 23,Stevinweg 1,2628CN                 |
!    |  Delft, The Netherlands                                   |  
!    |  Tel:   +31 15 278 5433                                   |
!  --|-----------------------------------------------------------|--

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
!------------------------------------------------------------------------
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Haiyang Cui                               |
!   --|-----------------------------------------------------------|--
 
!---Authors--------------------
!
!   1.00: Haiyang Cui
!------------------------------

!---Updates--------------------
!
!------------------------------
 
  
MODULE Mod_Mesh_DataTypes
 USE Mod_Precision
 IMPLICIT NONE
!---------------------------------------------------------------------
! node_type
! x   --- x coordinate
! y   --- y coordiante
! ind --- 0 interior node
!         1 solid boundary node, 
!         2 semi open semi Dirichlet node
!         3 open boundary node
!---------------------------------------------------------------------  
 TYPE node_type
    REAL(KIND=REAL_DP)            :: x
    REAL(KIND=REAL_DP)            :: y
    INTEGER(KIND=INT_KIND)        :: ind
 END TYPE node_type
 
 
!---------------------------------------------------------------------
! element_type
 
! nodeinds   --- store the three nodes' number
! edgeinds   --- store the three edges' number
! x          --- the x coordinate of the centroid
! y          --- the y coordinate of the centroid
! area       --- the surface area of this element
!--------------------------------------------------------------------- 
 TYPE element_type
    INTEGER(KIND=INT_KIND)            :: nodeinds(3)
    INTEGER(KIND=INT_KIND)            :: edgeinds(3)
    REAL(KIND=REAL_DP)                :: x
    REAL(KIND=REAL_DP)                :: y 
    REAL(KIND=REAL_DP)                :: area
 END TYPE element_type
 

!---------------------------------------------------------------------
! edge_type

! nodeinds   --- store the two nodes' ind
! edgeinds   --- store the two elements' ind
! dx, dy     --- the vector direction, point from n1 to n2, n1<n2
! ndx,ndy    --- the normal direction, point from left to right
! length     --- the length of the edge
!--------------------------------------------------------------------- 
 
 TYPE edge_type
    INTEGER(KIND=INT_KIND)               :: nodeinds(2)
    INTEGER(KIND=INT_KIND)               :: elementinds(2)
    REAL(KIND=REAL_DP)                   :: x
    REAL(KIND=REAL_DP)                   :: y
    REAL(KIND=REAL_DP)                   :: dx
    REAL(KIND=REAL_DP)                   :: dy
    REAL(KIND=REAL_DP)                   :: ndx
    REAL(KIND=REAL_DP)                   :: ndy  
    REAL(KIND=REAL_DP)                   :: length   
 END TYPE edge_type
END MODULE Mod_Mesh_DataTypes




!---Purpose--------------------
!   The TriMesh module. 
!   Main Subroutine is : TriMesh_Initialize
!   If you provide the node and element files, this subroutine will return
!   the derived data type Trimesh,  which inluces all the mesh information.
!------------------------------
  
MODULE MOD_TriMesh
  USE Mod_Precision
  USE Mod_Derived_DataType
  USE Mod_Mesh_DataTypes
  IMPLICIT NONE
  ! Type of mesh
  INTEGER(KIND=INT_KIND)                    :: Coordinate_type = 1
  ! Coordinate_type = 1, cartician coordinate, unit, [m]
  ! Coordinate_type = 2, radian coordinate,    unit, [degree]
       
  ! The two files, contain the nodes, and element information
  CHARACTER(len=400)                        :: node2Dfile=''
  CHARACTER(len=400)                        :: element2Dfile=''  
  
  INTEGER(KIND=INT_KIND)                    :: TriMesh_nodenmb
  INTEGER(KIND=INT_KIND)                    :: TriMesh_elementnmb
  INTEGER(KIND=INT_KIND)                    :: TriMesh_edgenmb
  
  TYPE(node_type),   ALLOCATABLE,DIMENSION(:)        :: node 
  TYPE(element_type),ALLOCATABLE,DIMENSION(:)        :: element 
  TYPE(edge_type),   ALLOCATABLE,DIMENSION(:)        :: edge 
  
  ! Link Information
  ! The neighbour elements of one node
  TYPE(address_type), ALLOCATABLE, DIMENSION(:)      :: node_ngb_elements    
  ! The neighbour nodes of one node
  TYPE(address_type), ALLOCATABLE, DIMENSION(:)      :: node_ngb_nodes
  
  ! The open boundary edge with index 2
  TYPE(address_type)                                 :: ob_edge_2
  ! The open boundary edge with index 3
  TYPE(address_type)                                 :: ob_edge_3  
  
  CONTAINS
  
  SUBROUTINE TriMesh_Initialize
    IMPLICIT NONE
  !Step1: TriMesh_ReadFile
  
  CALL TriMesh_ReadFile  

  ! Step2: Mesh scalling
  CALL TriMesh_Scaling
 
  ! Step3: Build the information of neighbour elements for a specified node
  CALL TriMesh_Build_Node_Ngb_Elements
 
  ! Step4: Build the information of neighbour nodes for a specified node
  CALL TriMesh_Build_Node_Ngb_Nodes
 
  ! Step5: Build edge arrays
  CALL TriMesh_Build_Edge_Arrays
 
  
  ! Step6: Compute element area
  CALL TriMesh_Compute_Element_Area
 
  ! Step7: Build open boundary edge, ind=2, ind=3
  
  CALL TriMesh_Build_Open_Boundary_Edge
  

  WRITE(*,'(a30,i8.6)') 'Total node      number is: ', TriMesh_nodenmb  
  WRITE(*,'(a30,i8.6)') 'Total element   number is: ', TriMesh_elementnmb  
  WRITE(*,'(a30,i8.6)') 'Total edge      number is: ', TriMesh_edgenmb  
  WRITE(*,'(a30,i8.6)') 'Total ob_edge_2 number is: ', ob_edge_2%nmb  
  WRITE(*,'(a30,i8.6)') 'Total ob_edge_3 number is: ', ob_edge_3%nmb   
 
END SUBROUTINE TriMesh_Initialize
  
  
SUBROUTINE TriMesh_ReadFile
    IMPLICIT NONE
    INTEGER                           :: fid = 100
    INTEGER                           :: IOErr
    LOGICAL                           :: alive
    
    
    INTEGER(KIND=INT_KIND)            :: n,nq,ind,n1,n2,n3
    REAL(KIND=REAL_DP)                :: x1,y1
        
    WRITE(*,*) 'Reading the mesh files....'

! Step 1 : Read the 2D triangle mesh files
! 1. Read the node2Dfile,     
    INQUIRE(file=node2Dfile, exist=alive)
    IF(alive) THEN
        OPEN (fid, file=TRIM(node2Dfile), status='old', IOSTAT=IOErr)
        IF( IOErr .NE. 0) THEN
            WRITE(*,*) 'Something wrong with the node2Dfile!'
            STOP
        ENDIF
        READ(fid,*) TriMesh_nodenmb   
        
          
        ALLOCATE (node(TriMesh_nodenmb))
        DO n = 1,TriMesh_nodenmb
            READ(fid,*) nq, x1, y1, ind
            node(n)%x       = x1
            node(n)%y       = y1
            node(n)%ind     = ind   
        ENDDO
        CLOSE(fid)
    ELSE
        WRITE (*, *) trim(node2Dfile), "doesn't exist."
    END IF
    
! 2. Read the elem2Dfile,     
    INQUIRE(file=element2Dfile, exist=alive)
    IF(alive) THEN
        OPEN (fid, file=TRIM(element2Dfile), status='old', IOSTAT=IOErr)
        IF( IOErr .NE. 0) THEN
            WRITE(*,*) 'Something wrong with the element2Dfile!'
            STOP
        ENDIF
        READ(fid,*) TriMesh_elementnmb     
        ALLOCATE (element(TriMesh_elementnmb))
        
        DO n = 1, TriMesh_elementnmb
            READ(fid,*) nq,n1,n2,n3
            element(n)%nodeinds(1) = n1
            element(n)%nodeinds(2) = n2
            element(n)%nodeinds(3) = n3
            element(n)%x = (node(n1)%x+node(n2)%x+node(n3)%x)/3.0d0
            element(n)%y = (node(n1)%y+node(n2)%y+node(n3)%y)/3.0d0
        ENDDO    
      
        CLOSE(fid)
    ELSE
        WRITE (*, *) trim(element2Dfile), "doesn't exist."
    END IF
    
  

  END SUBROUTINE TriMesh_ReadFile
  
  
  
  SUBROUTINE TriMesh_Scaling
    USE Mod_Parameter
    IMPLICIT NONE  
    INTEGER(KIND=INT_KIND) :: nd
    !  TRANSFORMATION OF 2D-NODE-ARRAYS, FROM DEGREE TO RADIAN     
    IF( Coordinate_type == 2) THEN
        DO nd=1,TriMesh_nodenmb
           node(nd)%x       = node(nd)%x *rad
           node(nd)%y       = node(nd)%y *rad
        ENDDO
       
        DO   nd=1,TriMesh_elementnmb 
           element(nd)%x    = element(nd)%x *rad
           element(nd)%y    = element(nd)%y *rad
        ENDDO
    ENDIF    
  END SUBROUTINE TriMesh_Scaling
  
  
  
  SUBROUTINE TriMesh_Build_Node_Ngb_Elements
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)                              :: j, a, b, c
    INTEGER(KIND=INT_KIND) , ALLOCATABLE, DIMENSION(:)  :: ind   
    !
    !--------------- 2D TriMesh:
    ! Builds Node_Ngb_Elements
    !       
    ALLOCATE(ind(TriMesh_nodenmb))
    ind = 0


    DO j = 1,TriMesh_elementnmb
        a = element(j)%nodeinds(1)
        b = element(j)%nodeinds(2)
        c = element(j)%nodeinds(3)
        ind(a) = ind(a) + 1
        ind(b) = ind(b) + 1
        ind(c) = ind(c) + 1
    ENDDO    
 
    IF (.NOT. ALLOCATED(node_ngb_elements)) THEN
        ALLOCATE (node_ngb_elements(TriMesh_nodenmb))
    ENDIF
    
    node_ngb_elements(1:TriMesh_nodenmb)%nmb = ind(1:TriMesh_nodenmb)
    
    DO j = 1,TriMesh_nodenmb
        ALLOCATE(node_ngb_elements(j)%addresses(ind(j)))
    ENDDO
    
    ind = 0   
     
    DO j=1,TriMesh_elementnmb
        a = element(j)%nodeinds(1)
        b = element(j)%nodeinds(2)
        c = element(j)%nodeinds(3)
        ind(a) = ind(a) + 1
        ind(b) = ind(b) + 1
        ind(c) = ind(c) + 1
        node_ngb_elements(a)%addresses(ind(a)) = j
        node_ngb_elements(b)%addresses(ind(b)) = j
        node_ngb_elements(c)%addresses(ind(c)) = j   
    ENDDO   
 
    DEALLOCATE(ind)
 
  END SUBROUTINE TriMesh_Build_Node_Ngb_Elements
  


  SUBROUTINE TriMesh_Build_Node_Ngb_Nodes
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)        :: j, a, b, c,el,k,m,ml(1)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: check
    INTEGER       :: count
    INTEGER(KIND=INT_KIND) , DIMENSION(100)   :: AUX=0
    !
    !--------------- 2D TriMesh:
    ! Builds  Node_Ngb_Nodes
    !       
    ALLOCATE( check(TriMesh_nodenmb))
    check = 0
    IF (.NOT. ALLOCATED(node_ngb_nodes)) THEN
        ALLOCATE (node_ngb_nodes(TriMesh_nodenmb))
    ENDIF    
   
    DO j = 1,TriMesh_nodenmb
        count = 0
        DO m = 1,node_ngb_elements(j)%nmb
            el = node_ngb_elements(j)%addresses(m)
            DO k=1, 3
                a = element(el)%nodeinds(k)
                IF (check(a) == 0) THEN
                    check(a) = 1
                    count = count + 1
                    aux(count) = a
                ENDIF
            ENDDO
        ENDDO
        node_ngb_nodes(j)%nmb = count
        
        ALLOCATE(node_ngb_nodes(j)%addresses(count))
        
        
        ! we need to sort array aux(1:count)
        DO m = count,1,-1
            ml = MAXLOC(aux(1:count))
            b = ml(1)
            node_ngb_nodes(j)%addresses(m) = aux(b)
            check(aux(b)) = 0
            aux(b) = -999
        ENDDO
    ENDDO ! end for nodenmb
    DEALLOCATE( check)
  END SUBROUTINE TriMesh_Build_Node_Ngb_Nodes
  
  SUBROUTINE TriMesh_Build_Edge_Arrays
    USE Mod_Parameter
    IMPLICIT NONE  
    
    INTEGER(KIND=INT_KIND)        :: n, el, nmb, cnt
    INTEGER(KIND=INT_KIND)        :: elnmb, nghbrnmb
    INTEGER(KIND=INT_KIND)        :: nghbrnodes(100)=0
    
    INTEGER(KIND=INT_KIND)        :: i, j, k, elnodes(3)
    REAL(kind=REAL_DP)                 :: vec1(3), vec2(3), nvec(3), r
    REAL(kind=REAL_DP)                 :: cosine
    INTEGER(KIND=INT_KIND)        :: edge2D,edg
    
 
    edge2D = 0
    cnt = 1
    vec1 = 0.
    vec2 = 0.
    nvec = 0.
  
    DO n = 1,TriMesh_nodenmb
        nghbrnmb = node_ngb_nodes(n)%nmb
        edge2D = edge2D + COUNT(node_ngb_nodes(n)%addresses>n)
    ENDDO
     
    TriMesh_edgenmb = edge2D
    IF (.NOT. ALLOCATED(edge)) THEN
            ALLOCATE (edge(TriMesh_edgenmb))
            edge(:)%nodeinds(1) = 0
            edge(:)%nodeinds(2) = 0
            edge(:)%elementinds(1) = 0
            edge(:)%elementinds(2) = 0
    ENDIF
               
     
    
    ! Build up the information about the edges
    DO n = 1,TriMesh_nodenmb
        elnmb = node_ngb_elements(n)%nmb
        nghbrnmb = node_ngb_nodes(n)%nmb
        
        nghbrnodes = 0
        nghbrnodes(1:nghbrnmb) = node_ngb_nodes(n)%addresses
        DO i = 1,nghbrnmb
            IF (nghbrnodes(i) <= n) CYCLE

            edge(cnt)%nodeinds(:) = (/n, nghbrnodes(i)/)   
               
            
            
            vec1(1) = node(nghbrnodes(i))%x - node(n)%x
            vec1(2) = node(nghbrnodes(i))%y - node(n)%y
            
            DO j = 1,elnmb
                el = node_ngb_elements(n)%addresses(j)
                elnodes = element(el)%nodeinds(:)
                
                !
                ! Finds whether edge cnt belongs to the element el or not
                !
                IF (COUNT(elnodes == n .or. elnodes == nghbrnodes(i)) == 2) THEN
                   
                    DO k = 1,3
                        IF (elnodes(k) /= n .and. elnodes(k) /= nghbrnodes(i)) THEN
                            vec2(1) = node(elnodes(k))%x - node(n)%x
                            vec2(2) = node(elnodes(k))%y - node(n)%y
                            EXIT
                        ENDIF
                    ENDDO
                    nvec(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
                    !
                    ! Finds whether element el is on the left or on the right of vec1
                    !
                    IF (nvec(3) > 0) THEN
                        edge(cnt)%elementinds(1) = el
                    ELSE
                        edge(cnt)%elementinds(2) = el
                    ENDIF
                    !
                    ! Finds which nodes of the triangle form the edge cnt
                    !
                    
                    IF (COUNT(elnodes(1:2) == n .or. elnodes(1:2) == nghbrnodes(i)) == 2) THEN
                        element(el)%edgeinds(1) = cnt
                    ELSEIF (COUNT(elnodes(2:3) == n .or. elnodes(2:3) == nghbrnodes(i)) == 2) THEN
                        element(el)%edgeinds(2) = cnt
                    ELSE
                        element(el)%edgeinds(3) = cnt
                    ENDIF            
                ENDIF
                
            ENDDO

            !
            ! Now construct the normal to the edge
            !
            
            IF(Coordinate_type == 2) THEN
                 cosine  = cos(0.5d0*(node(nghbrnodes(i))%y + node(n)%y))
                 vec1(1) = vec1(1) * cosine
                 vec1(1) = vec1(1) * r_earth
                 vec1(2) = vec1(2) * r_earth
            ENDIF
            
            edge(cnt)%dx =  vec1(1)
            edge(cnt)%dy =  vec1(2)

            edge(cnt)%ndx =  vec1(2)
            edge(cnt)%ndy = -vec1(1)
            
            edge(cnt)%length = SQRT(SUM(vec1*vec1))  
            r = edge(cnt)%length 
            edge(cnt)%ndx = edge(cnt)%ndx/r
            edge(cnt)%ndy = edge(cnt)%ndy/r
            
            
            IF (MINVAL(edge(cnt)%elementinds(:)) ==0 .and. nvec(3) < 0) THEN
                edge(cnt)%ndx  = -edge(cnt)%ndx
                edge(cnt)%ndy  = -edge(cnt)%ndy 
                edge(cnt)%elementinds(1) = edge(cnt)%elementinds(2)
                edge(cnt)%elementinds(2) = 0
                edge(cnt)%nodeinds(:) = (/edge(cnt)%nodeinds(2),edge(cnt)%nodeinds(1)/)
                
            ENDIF

            edge(cnt)%x = SUM(node(edge(cnt)%nodeinds(:))%x)/2.0d0
            edge(cnt)%y = SUM(node(edge(cnt)%nodeinds(:))%y)/2.0d0                   
            
            
            ! Come to the next edge
            cnt = cnt + 1
            
            
            
        ENDDO

    ENDDO

  END SUBROUTINE TriMesh_Build_Edge_Arrays   
  
  
  SUBROUTINE TriMesh_Compute_Element_Area
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)              :: el,eledges(3)
    REAL(KIND=REAL_DP)                  :: p,a,b,c
!Varibles for OpenMP
!$  integer :: chunk, chunksize  
!$  integer :: nthreads,OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 


!$OMP   PARALLEL PRIVATE(el,eledges,a,b,c,p,mythread)  &
!$OMP&         SHARED(TriMesh_elementnmb,element,edge,  &
!$OMP&         CHUNK,nthreads)
!$ NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_elementnmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    
    DO el=1,TriMesh_elementnmb
        eledges  = element(el)%edgeinds
        a = edge(eledges(1))%length
        b = edge(eledges(2))%length
        c = edge(eledges(3))%length
        p = 0.5d0*(a+b+c)
        element(el)%area = SQRT(p*(p-a)*(p-b)*(p-c))
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE TriMesh_Compute_Element_Area
  
  
  SUBROUTINE TriMesh_Build_Open_Boundary_Edge
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)         :: edg,nd1,nd2,edge_nmb2, edge_nmb3
    
    
    edge_nmb2 = 0
    edge_nmb3 = 0
    DO edg= 1, TriMesh_edgenmb
        nd1 = edge(edg)%nodeinds(1)
        nd2 = edge(edg)%nodeinds(2)
        
        IF( node(nd1)%ind == 2 .AND. node(nd2)%ind == 2) THEN
            edge_nmb2 = edge_nmb2 + 1
        ENDIF
        IF( node(nd1)%ind == 3 .AND. node(nd2)%ind == 3) THEN
            edge_nmb3 = edge_nmb3 + 1
        ENDIF        
    ENDDO
    
    
    ob_edge_2%nmb = edge_nmb2
    ob_edge_3%nmb = edge_nmb3

    ALLOCATE(ob_edge_2%addresses(ob_edge_2%nmb))
    ALLOCATE(ob_edge_3%addresses(ob_edge_3%nmb))
    
    edge_nmb2 = 0
    edge_nmb3 = 0
    DO edg= 1, TriMesh_edgenmb
        nd1 = edge(edg)%nodeinds(1)
        nd2 = edge(edg)%nodeinds(2)
        IF( node(nd1)%ind == 2 .AND. node(nd2)%ind == 2) THEN
            edge_nmb2 = edge_nmb2 + 1
            ob_edge_2%addresses(edge_nmb2) = edg
        ENDIF
        IF( node(nd1)%ind == 3 .AND. node(nd2)%ind == 3) THEN
            edge_nmb3 = edge_nmb3 + 1
            ob_edge_3%addresses(edge_nmb3) = edg
        ENDIF        
    ENDDO    
  END SUBROUTINE TriMesh_Build_Open_Boundary_Edge   
END MODULE MOD_TriMesh


 

