SUBROUTINE H2Ocean_NonHydro_Initialize
  USE Mod_All_Variable
  IMPLICIT NONE
  
  IF(stencil_method == 2) THEN
      CALL H2Ocean_Build_Node_Ngb_Nodes_Extended
  ENDIF
  
  CALL H2Ocean_Sparse_Matrix_Initialize
  CALL H2Ocean_Allocate_Nonhydro_Arrays
  CALL H2Ocean_Nonhydro_Build_grad_d
END SUBROUTINE H2Ocean_NonHydro_Initialize


SUBROUTINE H2Ocean_Build_Node_Ngb_Nodes_Extended
    USE Mod_TriMesh
    USE Mod_H2Ocean
    IMPLICIT NONE
    INTEGER                   :: j, k,l, m,n, a, b, c, nd_count, el, ml(1),edg, el_edg
    INTEGER, DIMENSION(100)   :: AUX=0
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ind, check
    
    ALLOCATE(ind(TriMesh_nodenmb), check(TriMesh_nodenmb))
    ! Builds node_extend_ngb_nodes
    !
    ALLOCATE(node_ngb_nodes_extended(TriMesh_nodenmb))
    check = 0
    DO j = 1,TriMesh_nodenmb
        nd_count = 0
        DO m = 1,node_ngb_elements(j)%nmb
            el = node_ngb_elements(j)%addresses(m)

            ! go through the three nodes of the element
            DO k=1, 3 
                a = element(el)%nodeinds(k)
                IF (check(a) == 0) THEN
                    check(a) = 1
                    nd_count = nd_count + 1
                    aux(nd_count) = a
                ENDIF
            ENDDO
            
            !go through the forth node, the outter one
            DO k=1,3
                edg = element(el)%edgeinds(k)
                IF( COUNT(edge(edg)%nodeinds(:) == j) == 0) THEN
                    DO l=1,2
                        el_edg = edge(edg)%elementinds(l)
                        IF(el_edg .ne. el .and. el_edg .ne. 0) THEN
                            DO n=1,3
                                a = element(el_edg)%nodeinds(n)
                                IF (check(a) == 0) THEN
                                    check(a) = 1
                                    nd_count = nd_count + 1
                                    aux(nd_count) = a
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDDO
                ENDIF
            ENDDO
        ENDDO
        node_ngb_nodes_extended(j)%nmb = nd_count
        
        ALLOCATE(node_ngb_nodes_extended(j)%addresses(nd_count))
        
        !
        ! we need to sort array aux(1:count)
        !
        DO m = nd_count,1,-1
            ml = MAXLOC(aux(1:nd_count))
            b = ml(1)
            node_ngb_nodes_extended(j)%addresses(m) = aux(b)
            check(aux(b)) = 0
            aux(b) = -999
        ENDDO
        
    ENDDO ! end for node2D
    DEALLOCATE(ind, check)

END SUBROUTINE H2Ocean_Build_Node_Ngb_Nodes_Extended


SUBROUTINE H2Ocean_Sparse_Matrix_Initialize
  USE Mod_TriMesh 
  USE Mod_Sparse_Matrix
  USE Mod_Precision
  USE Mod_H2Ocean
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND)   :: nd, i
  ALLOCATE(RHS_s(TriMesh_nodenmb))
  
  
  SELECT CASE(stencil_method)
  
  CASE(1)
    ! The dimension of the Matrix equals node2D
    A_s%dim = TriMesh_nodenmb
    ALLOCATE(A_s%rowptr(TriMesh_nodenmb + 1))
    A_s%rowptr(1) = 1
    DO nd = 1, TriMesh_nodenmb
        A_s%rowptr(nd+1) = A_s%rowptr(nd) + node_ngb_nodes(nd)%nmb
    END DO
  
    A_s%nza = A_s%rowptr(TriMesh_nodenmb+1) - 1
    WRITE(*,*) 'The number of nonzeros in A_s is: ',A_s%nza
    ALLOCATE(A_s%colind(A_s%nza), A_s%values(A_s%nza))
    A_s%values = 0.0
  
    DO nd = 1, TriMesh_nodenmb  
      DO i = 1, node_ngb_nodes(nd)%nmb
          A_s%colind(A_s%rowptr(nd) + i - 1) = node_ngb_nodes(nd)%addresses(i)
      END DO
    END DO
    CALL sparse_matrix_display_bandwidth(A_s)
  CASE(2)
    ! The dimension of the Matrix equals node2D
    A_s%dim = TriMesh_nodenmb
    ALLOCATE(A_s%rowptr(TriMesh_nodenmb + 1))
    A_s%rowptr(1) = 1
    DO nd = 1, TriMesh_nodenmb
        A_s%rowptr(nd+1) = A_s%rowptr(nd) + node_ngb_nodes_extended(nd)%nmb
    END DO
  
    A_s%nza = A_s%rowptr(TriMesh_nodenmb+1) - 1
    WRITE(*,*) 'The number of nonzeros in A_s is: ',A_s%nza
    ALLOCATE(A_s%colind(A_s%nza), A_s%values(A_s%nza))
    A_s%values = 0.0
  
    DO nd = 1, TriMesh_nodenmb  
      DO i = 1, node_ngb_nodes_extended(nd)%nmb
          A_s%colind(A_s%rowptr(nd) + i - 1) = node_ngb_nodes_extended(nd)%addresses(i)
      END DO
    END DO
    CALL sparse_matrix_display_bandwidth(A_s)  
  
  CASE DEFAULT
    WRITE(*,*) 'stencil_method =',stencil_method, 'is not specified!'
  END SELECT
END SUBROUTINE H2Ocean_Sparse_Matrix_Initialize


SUBROUTINE H2Ocean_Allocate_Nonhydro_Arrays
  USE Mod_H2Ocean
  USE Mod_TriMesh
  IMPLICIT NONE
  
  ALLOCATE(q(TriMesh_nodenmb))
  ALLOCATE(w_surf0(TriMesh_nodenmb),w_surf1(TriMesh_nodenmb))
  ALLOCATE(w_bot0(TriMesh_nodenmb),w_bot1(TriMesh_nodenmb))
  ALLOCATE(grad_d_dx(TriMesh_elementnmb),grad_d_dy(TriMesh_elementnmb))
  
  	  
  q = 0.0d0;
  w_surf0 = 0.0d0
  w_surf1 = 0.0d0
  w_bot0  = 0.0d0
  w_bot1  = 0.0d0
  grad_d_dx = 0.0d0
  grad_d_dy = 0.0d0
  
END SUBROUTINE H2Ocean_Allocate_Nonhydro_Arrays


SUBROUTINE H2Ocean_Nonhydro_Build_grad_d
  USE Mod_Precision 
  USE Mod_FEM
  USE Mod_H2Ocean
  USE Mod_TriMesh
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND)  :: el, elnodes(3)
  !The control volum of the water level will be 1/3 of all the neighbour elements
!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 

!$OMP   PARALLEL PRIVATE(el,elnodes,mythread) &
!$OMP&           SHARED(TriMesh_elementnmb,element,grad_d_dx,grad_d_dy, &
!$OMP&           bafuzeta_x,bafuzeta_y,Still_Depth,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_elementnmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
  DO el=1, TriMesh_elementnmb
      elnodes = element(el)%nodeinds(:)
      grad_d_dx(el) = DOT_PRODUCT(bafuzeta_x(:, el),Still_Depth(elnodes))
      grad_d_dy(el) = DOT_PRODUCT(bafuzeta_y(:, el),Still_Depth(elnodes))
  ENDDO
!$OMP END DO   
!$OMP END PARALLEL 

write(*,*)'maxval grad_d_dx grad_d_dy', maxval(grad_d_dx),maxval(grad_d_dy)
END SUBROUTINE H2Ocean_Nonhydro_Build_grad_d


SUBROUTINE H2Ocean_Nonhydro_Process
  USE Mod_H2Ocean
  IMPLICIT NONE
  
  
  
  CALL H2Ocean_NonHydro_Compute_w_bot
  
  
  SELECT CASE(stencil_method)
  CASE(1)  
    IF(enable_reduced_twolayer) THEN  
        CALL H2Ocean_NonHydro_Build_q_Spare_Matrix_ReducedTwoLayer
    ELSE
        CALL H2Ocean_NonHydro_Build_q_Spare_Matrix_1
    ENDIF
  CASE(2)
    CALL H2Ocean_NonHydro_Build_q_Spare_Matrix_2
  END SELECT
  
    
  CALL H2Ocean_NonHydro_Solve_q_Spare_Matrix
  
  IF(enable_reduced_twolayer) THEN
      CALL H2Ocean_NonHydro_Update_Velocity_ReducedTwoLayer
  ELSE
      CALL H2Ocean_NonHydro_Update_Velocity      
  ENDIF
  

END SUBROUTINE H2Ocean_Nonhydro_Process



SUBROUTINE H2Ocean_NonHydro_Compute_w_bot
    USE Mod_Precision
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_FEM
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nd,el, eledges(3), elnodes(3), j
    REAL(kind=REAL_DP)      :: temp, u_ave, v_ave, dep_dx, dep_dy
  !Varibles for OpenMP
!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 
 
!$OMP   PARALLEL PRIVATE(nd,temp,j,el,eledges,elnodes,u_ave,v_ave, &
!$OMP&           dep_dx,dep_dy,mythread)  &
!$OMP&           SHARED(w_bot1,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)        
    DO nd=1,TriMesh_nodenmb
        temp = 0.0d0;
        IF( zeta0(nd) + Still_Depth(nd) <= h_min) THEN
            w_bot1(nd) = 0.0d0
            CYCLE
        ENDIF 
        DO j=1, node_ngb_elements(nd)%nmb
            el = node_ngb_elements(nd)%addresses(j)
            eledges = element(el)%edgeinds(:) 
            elnodes = element(el)%nodeinds(:)  
            u_ave = SUM(velocity1(eledges)%u)/3.0d0
            v_ave = SUM(velocity1(eledges)%v)/3.0d0 
            temp = temp + (u_ave*grad_d_dx(el)+v_ave*grad_d_dy(el)) &
                         & *element(el)%area/3.0d0     
        ENDDO
        w_bot1(nd) = -temp/Zeta_Area(nd)
    ENDDO
 !$OMP END DO
 !$OMP END PARALLEL     
END SUBROUTINE H2Ocean_NonHydro_Compute_w_bot



SUBROUTINE H2Ocean_NonHydro_Build_q_Spare_Matrix_1
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_FEM
  USE Mod_Sparse_Matrix
  USE Mod_Dirichlet_Boundary
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND)  :: nd, i,el,eledges(3),elnodes(3)
  REAL(KIND=REAL_DP)      :: flux,u_ave, v_ave,discharge, h_flux
  INTEGER(KIND=INT_KIND)  :: j, ic,edg
  INTEGER(KIND=INT_KIND)  :: nd1, nd2 
  REAL(KIND=REAL_DP)      :: vol_sum
    
  INTEGER(KIND=INT_KIND)  :: offset
 
  !Varibles for OpenMP
!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 

  A_s%values = 0.0d0
  RHS_s = 0.0d0
 
!$OMP   PARALLEL PRIVATE(nd,i,el,elnodes,eledges,u_ave,v_ave,j,edg,offset,   &
!$OMP&                 discharge,h_flux,mythread)  &
!$OMP&           SHARED(A_s,RHS_s,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)   
  DO nd=1,TriMesh_nodenmb
       IF( zeta0(nd) +Still_Depth(nd) <=h_min) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) = 0.0d0
           CYCLE
       ENDIF
       IF( -Still_Depth(nd) >=bath_max) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) = 0.0d0
           CYCLE
       ENDIF       
       
       
       
        IF( node(nd)%ind == 3) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_3(nd)
           CYCLE
       ENDIF
       
       IF( node(nd)%ind == 2 .AND. .NOT. dirichlet_boundary) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2(nd)
           CYCLE
       ENDIF
       
       IF( node(nd)%ind == 2 .AND. dirichlet_boundary .AND. D_bo >=D_bo_Inf) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2(nd)
           CYCLE           
       ENDIF       
       IF( node(nd)%ind == 2 .AND. q_dirichlet_boundary .AND. q_D_bo <D_bo_Inf) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) = q_D_bo*Zeta_Area(nd)*dt
           CYCLE            
       ENDIF
 
       
      
       
       DO i = 1, node_ngb_elements(nd)%nmb
          el         = node_ngb_elements(nd)%addresses(i)
          elnodes    = element(el)%nodeinds(:)
          eledges    = element(el)%edgeinds(:)
          u_ave      = SUM(velocity1(eledges)%u)/3.0d0
          v_ave      = SUM(velocity1(eledges)%v)/3.0d0
          IF( MAXVAL(zeta0(elnodes))+MINVAL(Still_Depth(elnodes)) <=h_min) THEN
              CYCLE
          ENDIF 
          DO j=1,3  ! go through 3 edges
              edg = eledges(j)
              IF(ANY(edge(edg)%nodeinds(:)==nd)) THEN   
                 CALL H2Ocean_Compute_edge_discharge(nd,el,edg,u_ave,v_ave,discharge)
                 h_flux = Zeta0(nd)+Still_Depth(nd)
                 RHS_s(nd) = RHS_s(nd) + discharge*h_flux
                 CALL H2Ocean_NonHydro_Spare_Matrix_Fill_Edge(nd,el,edg,h_flux)
              ENDIF   
              !Vertical bottom velocity contribution
               !CALL H2Ocean_NonHydro_Spare_Matrix_Fill_W_bot(nd,el,edg)
          ENDDO  ! END DO j=1,3
       ENDDO     ! END DO i = 1, TriMesh.node_ngb_elements(nd)%nmb 
       
       ! vertical verlocity
       RHS_s(nd) = RHS_s(nd) - Zeta_Area(nd)*(w_surf0(nd)-2*w_bot1(nd)+w_bot0(nd)) 
       CALL sparse_matrix_locate(nd,nd,A_s,offset)
       A_s%values(offset) = A_s%values(offset) + 2.0*Zeta_Area(nd)*dt/(zeta0(nd)+Still_Depth(nd))        
  ENDDO  ! END DO nd=1,TriMesh_nodenmb      
 !$OMP END DO
 !$OMP END PARALLEL 
END SUBROUTINE H2Ocean_NonHydro_Build_q_Spare_Matrix_1

SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_W_bot(nd,el,edg)
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_Sparse_Matrix
  USE Mod_FEM
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND),INTENT(IN) :: nd,el,edg
  REAL(KIND=REAL_DP)     :: dhdx, dhdy
  
  INTEGER(KIND=INT_KIND) :: i, nd_el
  REAL(KIND=REAL_DP)     :: vol
  INTEGER(KIND=INT_KIND) :: offset
   
  INTEGER(KIND=INT_KIND) :: j,edge_el 
  
  
  
   
  
  dhdx = q_alpha*grad_zeta_dx(el) - (1.0d0-q_alpha)*grad_d_dx(el)
  dhdy = q_alpha*grad_zeta_dy(el) - (1.0d0-q_alpha)*grad_d_dy(el)
  vol = element(el)%area
  DO i=1,3
      nd_el = element(el)%nodeinds(i)
      CALL sparse_matrix_locate(nd,nd_el,A_s, offset)
          A_s%values(offset) = A_s%values(offset) &
          -  vol/3.0*1.0d0/3.0d0*dt  *(q_alpha* bafuzeta_x(i,el)  &
                                 +1.0/h_f(edg)*1.0d0/3.0d0*dhdx)*grad_zeta_dx(el) &
          -  vol/3.0*1.0d0/3.0d0*dt  *(q_alpha* bafuzeta_y(i,el)  &
                                 +1.0/h_f(edg)*1.0d0/3.0d0*dhdy)*grad_zeta_dy(el)                                  
 ENDDO  
   
END SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_W_bot


SUBROUTINE H2Ocean_NonHydro_Build_q_Spare_Matrix_2
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_FEM
  USE Mod_Sparse_Matrix
  USE Mod_Dirichlet_Boundary
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND)  :: nd, i,el,eledges(3),elnodes(3)
  REAL(KIND=REAL_DP)      :: flux,u_ave, v_ave,discharge, h_flux
  INTEGER(KIND=INT_KIND)  :: j, ic,edg
  INTEGER(KIND=INT_KIND)  :: nd1, nd2 
  REAL(KIND=REAL_DP)      :: vol_sum
    
  INTEGER(KIND=INT_KIND)  :: offset
  !Varibles for OpenMP
!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 

  A_s%values = 0.0d0
  RHS_s = 0.0d0
 
!$OMP   PARALLEL PRIVATE(nd,i,el,elnodes,eledges,u_ave,v_ave,j,edg,offset,   &
!$OMP&                 discharge,h_flux,mythread)  &
!$OMP&           SHARED(A_s,RHS_s,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)   
  DO nd=1,TriMesh_nodenmb
       IF( zeta0(nd) +Still_Depth(nd) <=h_min) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) = 0.0d0
           CYCLE
       ENDIF
        IF( node(nd)%ind == 3) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_3(nd)
           CYCLE
       ENDIF
       
       IF( node(nd)%ind == 2 .AND. .NOT. dirichlet_boundary) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2(nd)
           CYCLE
       ENDIF
       
       IF( node(nd)%ind == 2 .AND. dirichlet_boundary .AND. D_bo >=D_bo_Inf) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2(nd)
           CYCLE           
       ENDIF  
       
       IF( node(nd)%ind == 2 .AND. q_dirichlet_boundary .AND. q_D_bo <D_bo_Inf) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) = q_D_bo*Zeta_Area(nd)*dt
           CYCLE            
       ENDIF
                              
       DO i = 1, node_ngb_elements(nd)%nmb
          el         = node_ngb_elements(nd)%addresses(i)
          elnodes    = element(el)%nodeinds(:)
          eledges    = element(el)%edgeinds(:)
          u_ave      = SUM(velocity1(eledges)%u)/3.0d0
          v_ave      = SUM(velocity1(eledges)%v)/3.0d0
          IF( MAXVAL(zeta0(elnodes))+MINVAL(Still_Depth(elnodes)) <=h_min) THEN
              CYCLE
          ENDIF 
          DO j=1,3  ! go through 3 edges
              edg = eledges(j)
              IF(ANY(edge(edg)%nodeinds(:)==nd)) THEN   
                 CALL H2Ocean_Compute_edge_discharge(nd,el,edg,u_ave,v_ave,discharge)
                 h_flux = Zeta0(nd)+Still_Depth(nd)      
                 RHS_s(nd) = RHS_s(nd) + discharge*h_flux
                 CALL H2Ocean_NonHydro_Spare_Matrix_Fill_Edge_2(nd,el,edg,h_flux)
              ENDIF   
          ENDDO  ! END DO j=1,3
          
          !Vertical bottom velocity contribution
!          CALL H2Ocean_NonHydro_Spare_Matrix_Fill_W_bot(nd,el,edg)
              
       ENDDO     ! END DO i = 1, TriMesh.node_ngb_elements(nd)%nmb 
       
       ! vertical verlocity
       RHS_s(nd) = RHS_s(nd) - Zeta_Area(nd)*(w_surf0(nd)- w_bot1(nd)+w_bot0(nd)) 
       CALL sparse_matrix_locate(nd,nd,A_s,offset)
       A_s%values(offset) = A_s%values(offset) + 2.0*Zeta_Area(nd)*dt/(zeta0(nd)+Still_Depth(nd))        
  ENDDO  ! END DO nd=1,TriMesh_nodenmb      
 !$OMP END DO
 !$OMP END PARALLEL 
END SUBROUTINE H2Ocean_NonHydro_Build_q_Spare_Matrix_2

SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2(nd)
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_Sparse_Matrix
  USE Mod_FEM
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND),INTENT(IN) :: nd
  REAL(kind=REAL_DP)       :: nx_op,ny_op,vol_sum,vol 
  INTEGER                  :: i,j,k,edg_op, edg,el,nd_el
  INTEGER(KIND=INT_KIND)   :: offset
  REAL(kind=REAL_DP)       :: dhdx, dhdy
  edg_op = 0
  nx_op = 0.0d0
  ny_op = 0.0d0
  vol_sum = 0.0d0
  DO i=1, ob_edge_2%nmb
      edg = ob_edge_2%addresses(i)
	  IF(ANY(edge(edg)%nodeinds(:)  == nd)) THEN
	      edg_op = edg_op + 1
	      nx_op = nx_op + edge(edg)%ndx
	      ny_op = ny_op + edge(edg)%ndy
	  ENDIF
  ENDDO		 
  nx_op  = nx_op/edg_op
  ny_op  = ny_op/edg_op
	  
	  
  
  DO j=1, node_ngb_elements(nd)%nmb
	      el = node_ngb_elements(nd)%addresses(j)
	      vol = element(el)%area
	      vol_sum = vol_sum + vol
	      
	      dhdx = grad_zeta_dx(el) - grad_d_dx(el)
          dhdy = grad_zeta_dy(el) - grad_d_dy(el)
	      DO k=1,3
	          nd_el = element(el)%nodeinds(k)
	          CALL sparse_matrix_locate(nd, nd_el,A_s, offset)
	              A_s%values(offset) = A_s%values(offset) &
	              + vol*nx_op*bafuzeta_x(k, el) + vol*ny_op*bafuzeta_y(k, el) &
	              + vol*nx_op*1./3./(zeta0(nd)+Still_Depth(nd))*dhdx &
	              + vol*ny_op*1./3./(zeta0(nd)+Still_Depth(nd))*dhdy 
	      ENDDO
  ENDDO
  A_s%values(A_s%rowptr(nd) : A_s%rowptr(nd+1)-1 ) = &
  A_s%values(A_s%rowptr(nd) : A_s%rowptr(nd+1)-1 )/vol_sum
  RHS_s(nd) = 0.0d0
END SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2

SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_3(nd)
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_Sparse_Matrix
  USE Mod_FEM
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND),INTENT(IN) :: nd
  REAL(kind=REAL_DP)       :: nx_op,ny_op,vol_sum,vol 
  INTEGER                  :: i,j,k,edg_op, edg,el,nd_el
  INTEGER(KIND=INT_KIND)   :: offset
  REAL(kind=REAL_DP)       :: dhdx, dhdy
  edg_op = 0
  nx_op = 0.0d0
  ny_op = 0.0d0
  vol_sum = 0.0d0
  DO i=1, ob_edge_3%nmb
      edg = ob_edge_3%addresses(i)
	  IF(ANY(edge(edg)%nodeinds(:)  == nd)) THEN
	      edg_op = edg_op + 1
	      nx_op = nx_op + edge(edg)%ndx
	      ny_op = ny_op + edge(edg)%ndy
	  ENDIF
  ENDDO		 
  nx_op  = nx_op/edg_op
  ny_op  = ny_op/edg_op
	  
	  
  
  DO j=1, node_ngb_elements(nd)%nmb
	      el = node_ngb_elements(nd)%addresses(j)
	      vol = element(el)%area
	      vol_sum = vol_sum + vol
	      
	      dhdx = grad_zeta_dx(el) - grad_d_dx(el)
          dhdy = grad_zeta_dy(el) - grad_d_dy(el)
	      DO k=1,3
	          nd_el = element(el)%nodeinds(k)
	          CALL sparse_matrix_locate(nd, nd_el,A_s, offset)
	              A_s%values(offset) = A_s%values(offset) &
	              + vol*nx_op*bafuzeta_x(k, el) + vol*ny_op*bafuzeta_y(k, el) &
	              + vol*nx_op*1./3./(zeta0(nd)+Still_Depth(nd))*dhdx &
	              + vol*ny_op*1./3./(zeta0(nd)+Still_Depth(nd))*dhdy 
	      ENDDO
  ENDDO
  A_s%values(A_s%rowptr(nd) : A_s%rowptr(nd+1)-1 ) = &
  A_s%values(A_s%rowptr(nd) : A_s%rowptr(nd+1)-1 )/vol_sum
  RHS_s(nd) = 0.0d0
END SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_3 




SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_Edge(nd,el,edg,h_flux)
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_Sparse_Matrix
  USE Mod_FEM
  USE Mod_Dirichlet_Boundary
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND),INTENT(IN) :: nd,el,edg
  REAL(KIND=REAL_DP),INTENT(IN)     :: h_flux 
  
  REAL(KIND=REAL_DP)     :: dhdx, dhdy
  
  INTEGER(KIND=INT_KIND) :: i, nd_el
  REAL(KIND=REAL_DP)     :: dx, dy, L
  INTEGER(KIND=INT_KIND) :: offset
   
  INTEGER(KIND=INT_KIND) :: j,edge_el, eledges(3)
 
  
  
  eledges = element(el)%edgeinds(:)
  
  dhdx = q_alpha*grad_zeta_dx(el) - (1.0d0-q_alpha)*grad_d_dx(el)
  dhdy = q_alpha*grad_zeta_dy(el) - (1.0d0-q_alpha)*grad_d_dy(el)

  IF(edge(edg)%elementinds(1) == el) THEN
      dx = edge_vectors(edg,1)%dx
      dy = edge_vectors(edg,1)%dy
      L  = edge_vectors(edg,1)%length
      IF(edge(edg)%nodeinds(1) == nd) THEN
         dx = -edge_vectors(edg,1)%dx
         dy = -edge_vectors(edg,1)%dy
      ENDIF
  ENDIF
  
  IF(edge(edg)%elementinds(2) == el) THEN
      dx = edge_vectors(edg,2)%dx
      dy = edge_vectors(edg,2)%dy
      L  = edge_vectors(edg,2)%length
      IF(edge(edg)%nodeinds(1) == nd) THEN
         dx = -edge_vectors(edg,2)%dx
         dy = -edge_vectors(edg,2)%dy
      ENDIF      
  ENDIF   
  

  
  DO j=1,3
        edge_el = eledges(j)  
        DO i=1,3
          nd_el = element(el)%nodeinds(i)
          CALL sparse_matrix_locate(nd,nd_el,A_s, offset)

          A_s%values(offset) = A_s%values(offset) &
          + h_flux*L*dx*1.0d0/3.0d0*dt *(q_alpha* bafuzeta_y(i,el)  &
                                     +1.0/MAX(h_f(edge_el),h_f_min)*1.0d0/3.0d0*dhdy) &
          - h_flux*L*dy*1.0d0/3.0d0*dt *(q_alpha* bafuzeta_x(i,el)  &
                                     +1.0/MAX(h_f(edge_el),h_f_min)*1.0d0/3.0d0*dhdx)                                  
      ENDDO  
  ENDDO    
END SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_Edge




SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_Edge_2(nd,el,edg,h_flux)
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_Sparse_Matrix
  USE Mod_FEM
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND),INTENT(IN) :: nd,el,edg
  REAL(KIND=REAL_DP),INTENT(IN)     :: h_flux 
  
  REAL(KIND=REAL_DP)     :: dhdx, dhdy
  
  INTEGER(KIND=INT_KIND) :: i, nd_el
  REAL(KIND=REAL_DP)     :: dx, dy, L
  INTEGER(KIND=INT_KIND) :: offset
   
  INTEGER(KIND=INT_KIND) :: j,edge_el, eledges(3)
  INTEGER(KIND=INT_KIND) :: elin,el_edge
  REAL(KIND=REAL_DP)      :: alpha
  
  eledges = element(el)%edgeinds(:)
  IF(edge(edg)%elementinds(1) == el) THEN
      dx = edge_vectors(edg,1)%dx
      dy = edge_vectors(edg,1)%dy
      L  = edge_vectors(edg,1)%length
      IF(edge(edg)%nodeinds(1) == nd) THEN
         dx = -edge_vectors(edg,1)%dx
         dy = -edge_vectors(edg,1)%dy
      ENDIF
  ENDIF
  
  IF(edge(edg)%elementinds(2) == el) THEN
      dx = edge_vectors(edg,2)%dx
      dy = edge_vectors(edg,2)%dy
      L  = edge_vectors(edg,2)%length
      IF(edge(edg)%nodeinds(1) == nd) THEN
         dx = -edge_vectors(edg,2)%dx
         dy = -edge_vectors(edg,2)%dy
      ENDIF      
  ENDIF   
  
  DO j=1,3  ! Go through 3 vel
      edge_el = eledges(j)  
      DO elin=1,2 ! Go through two neighbour element
        el_edge = edge(edge_el)%elementinds(elin)
        IF(el_edge ==0) CYCLE
        alpha = element(el_edge)%area/(3.0*Vel_Area(edge_el))

        dhdx = q_alpha*grad_zeta_dx(el_edge) - (1.0d0-q_alpha)*grad_d_dx(el_edge)
        dhdy = q_alpha*grad_zeta_dy(el_edge) - (1.0d0-q_alpha)*grad_d_dy(el_edge)        
        DO i=1,3 !Go through 3 nodes of each el
          nd_el = element(el_edge)%nodeinds(i)
          CALL sparse_matrix_locate(nd,nd_el,A_s, offset)

          A_s%values(offset) = A_s%values(offset) &
          + alpha*h_flux*L*dx*1.0d0/3.0d0*dt *(q_alpha* bafuzeta_y(i,el_edge)  &
                                    +1.0/MAX(h_f(edge_el),h_f_min)*1.0d0/3.0d0*dhdy) &
          - alpha*h_flux*L*dy*1.0d0/3.0d0*dt *(q_alpha* bafuzeta_x(i,el_edge)  &
                                    +1.0/MAX(h_f(edge_el),h_f_min)*1.0d0/3.0d0*dhdx)                                  
        ENDDO
      ENDDO  
  ENDDO    
END SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_Edge_2

SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_Edge_ReducedTwoLayer(nd,el,edg,h_flux)
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_Sparse_Matrix
  USE Mod_FEM
  USE Mod_Dirichlet_Boundary
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND),INTENT(IN) :: nd,el,edg
  REAL(KIND=REAL_DP),INTENT(IN)     :: h_flux 
  
  REAL(KIND=REAL_DP)     :: dhdx, dhdy
  
  INTEGER(KIND=INT_KIND) :: i, nd_el
  REAL(KIND=REAL_DP)     :: dx, dy, L
  INTEGER(KIND=INT_KIND) :: offset
   
  INTEGER(KIND=INT_KIND) :: j,edge_el, eledges(3),elnodes(3)
 
  REAL(KIND=REAL_DP)     :: coef1,alpha_el
  
  eledges = element(el)%edgeinds(:)
  elnodes = element(el)%nodeinds(:)
  alpha_el = SUM(alpha_h(elnodes))/3.0d0
  coef1 =  (1.0+alpha_el)*(1.0+alpha_el)/2.0 +(alpha_el-alpha_el**2)/2.0
 
  dhdx = -grad_d_dx(el)
  dhdy = -grad_d_dy(el)

  IF(edge(edg)%elementinds(1) == el) THEN
      dx = edge_vectors(edg,1)%dx
      dy = edge_vectors(edg,1)%dy
      L  = edge_vectors(edg,1)%length
      IF(edge(edg)%nodeinds(1) == nd) THEN
         dx = -edge_vectors(edg,1)%dx
         dy = -edge_vectors(edg,1)%dy
      ENDIF
  ENDIF
  
  IF(edge(edg)%elementinds(2) == el) THEN
      dx = edge_vectors(edg,2)%dx
      dy = edge_vectors(edg,2)%dy
      L  = edge_vectors(edg,2)%length
      IF(edge(edg)%nodeinds(1) == nd) THEN
         dx = -edge_vectors(edg,2)%dx
         dy = -edge_vectors(edg,2)%dy
      ENDIF      
  ENDIF   
  
  DO j=1,3
        edge_el = eledges(j)  
        DO i=1,3
          nd_el = element(el)%nodeinds(i)
          CALL sparse_matrix_locate(nd,nd_el,A_s, offset)
          A_s%values(offset) = A_s%values(offset) &
          + coef1*h_flux*L*dx*1.0d0/3.0d0*dt *( bafuzeta_y(i,el))  &                         
          - coef1*h_flux*L*dy*1.0d0/3.0d0*dt *( bafuzeta_x(i,el))                                                            
      ENDDO  
  ENDDO    
END SUBROUTINE H2Ocean_NonHydro_Spare_Matrix_Fill_Edge_ReducedTwoLayer


SUBROUTINE H2Ocean_NonHydro_Solve_q_Spare_Matrix
   USE Mod_Sparse_Matrix
   USE Mod_H2Ocean
   USE Mod_Precision
   IMPLICIT NONE
   INTEGER            ::  maxiter
   REAL(kind=REAL_DP) :: tol1, tol2
   CHARACTER(LEN=200) :: method  
   
   maxiter  = 2000
   tol1     = 1.0e-7           ! The relative tolerance,
   tol2     = 1.0e-7           ! The absolute tolerance
   method = 'bcgstab'
   CALL SparseSolver(A_s,         &      ! The Matrix
 			  q,                  &      ! The solution
 			  RHS_s,              &      ! The right hand of the Ax=b
 			  maxiter,            &      ! Compute the vertical velocity
 			  tol1,               &      ! The relative tolerance,
 			  tol2,               &      ! The absolute tolerance
 			  method)                    ! The Iteration method		
END SUBROUTINE H2Ocean_NonHydro_Solve_q_Spare_Matrix


SUBROUTINE H2Ocean_NonHydro_Update_Velocity
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_Parameter
    USE Mod_FEM
    USE Mod_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nd, el, j, k, eledges(3), offset,edg,el1,el2
    INTEGER(KIND=INT_KIND)  :: nd1,nd2,elnodes(3),edg_el, nd_el
    REAL(kind=REAL_DP)      :: temp, u_ave,v_ave, zeta_dx, zeta_dy, discharge,h_ave
    REAL(kind=REAL_DP)      :: temp_u, temp_v, alpha
    REAL(kind=REAL_DP)      :: nx, ny, tmpu, tmpv
    REAL(kind=REAL_DP)      :: dhdx,dhdy
    REAL(kind=REAL_DP)      :: hu(3)
  !Varibles for OpenMP
!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 

 
 
!$OMP   PARALLEL PRIVATE(edg,temp_u,temp_v,el1,el2,nd1,nd2,eledges,elnodes, &
!$OMP&  temp,alpha,dhdx,dhdy,j,nd_el,nx,ny,tmpu,tmpv,nd,hu,mythread)  &
!$OMP&           SHARED(velocity1,w_surf1,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)      
    DO edg = 1, TriMesh_edgenmb    
      IF( h_wd(edg) <= h_min) THEN
          velocity1(edg)%u = 0.0d0
          velocity1(edg)%v = 0.0d0
          CYCLE
      ENDIF
      IF( h_f(edg) <= h_f_min) THEN
          velocity1(edg)%u = 0.0d0
          velocity1(edg)%v = 0.0d0
          CYCLE
      ENDIF   
      
      
      IF(ALL(node(edge(edg)%nodeinds)%ind == 2) .AND. u_dirichlet_boundary .AND. u_D_bo <D_bo_Inf) THEN
         velocity1(edg)%u = u_D_bo
         velocity1(edg)%v = 0.0d0
         CYCLE
      ENDIF
      
      
      temp_u = 0.0d0
      temp_v = 0.0d0
      el1 = edge(edg)%elementinds(1)
      el2 = edge(edg)%elementinds(2)
      nd1 = edge(edg)%nodeinds(1)
      nd2 = edge(edg)%nodeinds(2)
      eledges = element(el1)%edgeinds(:)
      elnodes = element(el1)%nodeinds(:)
 
      temp = 0.0d0
      hu =zeta0(elnodes)+Still_Depth(elnodes)
      IF(ALL(hu> h_min)) THEN
          alpha = element(el1)%area/(3.0*Vel_Area(edg))
          dhdx = q_alpha*grad_zeta_dx(el1) - (1.0d0-q_alpha)*grad_d_dx(el1)
          dhdy = q_alpha*grad_zeta_dy(el1) - (1.0d0-q_alpha)*grad_d_dy(el1)
          DO j=1,3
              nd_el = elnodes(j)
              temp_u = temp_u -  q_alpha*dt *q(nd_el)*bafuzeta_x(j, el1)*alpha  &
              - dt/(MAX(h_f(edg),h_f_min))*q(nd_el)/3.0*alpha *dhdx 
               
              temp_v = temp_v - q_alpha*dt*q(nd_el)*bafuzeta_y(j, el1)*alpha  &
              - dt/(MAX(h_f(edg),h_f_min))*q(nd_el)/3.0*alpha*dhdy  
          ENDDO                        
      ENDIF
      
      
      IF(el2 .NE. 0) THEN
           eledges = element(el2)%edgeinds(:)
           elnodes = element(el2)%nodeinds(:)
           hu =zeta0(elnodes)+Still_Depth(elnodes)
           IF(ALL(hu> h_min)) THEN
              alpha = element(el2)%area/(3.0*Vel_Area(edg))
              dhdx = q_alpha*grad_zeta_dx(el2) - (1.0d0-q_alpha)*grad_d_dx(el2)
              dhdy = q_alpha*grad_zeta_dy(el2) - (1.0d0-q_alpha)*grad_d_dy(el2)            
              DO j=1,3
                  nd_el = elnodes(j)
                  temp_u = temp_u -  q_alpha*dt*q(nd_el)*bafuzeta_x(j, el2)*alpha  &
                  - dt/(MAX(h_f(edg),h_f_min))*q(nd_el)/3.0*alpha*dhdx  
               
                  temp_v = temp_v - q_alpha*dt*q(nd_el)*bafuzeta_y(j, el2)*alpha  &
                  - dt/(MAX(h_f(edg),h_f_min))*q(nd_el)/3.0*alpha*dhdy 
              ENDDO               
          ENDIF
      ENDIF
         
      velocity1(edg)%u = velocity1(edg)%u + temp_u
      velocity1(edg)%v = velocity1(edg)%v + temp_v
      
     ! Solid Boundary Conditions
     !nx = edge(edg)%ndx
     !ny = edge(edg)%ndy     
     !IF (MINVAL(edge(edg)%elementinds(:)) == 0 .and. &
     !    COUNT(node(edge(edg)%nodeinds(:))%ind==1) == 2) THEN
     !   tmpu = ny*(ny*velocity1(edg)%u  - nx*velocity1(edg)%v)
     !   tmpv = nx*(nx*velocity1(edg)%v  - ny*velocity1(edg)%u)
     !   velocity1(edg)%u  = tmpu
     !   velocity1(edg)%v  = tmpv
     !ENDIF
     
  ENDDO    
 !$OMP END DO
  
  
  
!$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)  
  DO nd=1,TriMesh_nodenmb
  
      IF(zeta0(nd)+Still_Depth(nd) <= h_min) THEN
          w_surf1(nd) = 0.0d0
      ELSE
          w_surf1(nd) = w_surf0(nd)-(w_bot1(nd)-w_bot0(nd)) + 2.0* dt*q(nd)/(zeta1(nd)+Still_Depth(nd))
      ENDIF
  ENDDO
 !$OMP END DO
 !$OMP END PARALLEL 

  
  w_surf0 = w_surf1
  w_bot0  = w_bot1

END SUBROUTINE H2Ocean_NonHydro_Update_Velocity



SUBROUTINE H2Ocean_NonHydro_Update_Velocity_ReducedTwoLayer
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_Parameter
    USE Mod_FEM
    USE Mod_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nd, el, j, k, eledges(3), offset,edg,el1,el2
    INTEGER(KIND=INT_KIND)  :: nd1,nd2,elnodes(3),edg_el, nd_el
    REAL(kind=REAL_DP)      :: temp, u_ave,v_ave, zeta_dx, zeta_dy, discharge,h_ave
    REAL(kind=REAL_DP)      :: temp_u, temp_v, alpha
    REAL(kind=REAL_DP)      :: nx, ny, tmpu, tmpv
    REAL(kind=REAL_DP)      :: dhdx,dhdy
    REAL(kind=REAL_DP)      :: hu(3)
  !Varibles for OpenMP
!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 

 
 
!$OMP   PARALLEL PRIVATE(edg,temp_u,temp_v,el1,el2,nd1,nd2,eledges,elnodes, &
!$OMP&  temp,alpha,dhdx,dhdy,j,nd_el,nx,ny,tmpu,tmpv,nd,hu,mythread)  &
!$OMP&           SHARED(velocity1,w_surf1,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)      
    DO edg = 1, TriMesh_edgenmb    
      IF( h_wd(edg) <= h_min) THEN
          velocity1(edg)%u = 0.0d0
          velocity1(edg)%v = 0.0d0
          CYCLE
      ENDIF
      IF( h_f(edg) <= h_f_min) THEN
          velocity1(edg)%u = 0.0d0
          velocity1(edg)%v = 0.0d0
          CYCLE
      ENDIF   
      
      
      IF(ALL(node(edge(edg)%nodeinds)%ind == 2) .AND. u_dirichlet_boundary .AND. u_D_bo <D_bo_Inf) THEN
         velocity1(edg)%u = u_D_bo
         velocity1(edg)%v = 0.0d0
         CYCLE
      ENDIF
      
      
      temp_u = 0.0d0
      temp_v = 0.0d0
      el1 = edge(edg)%elementinds(1)
      el2 = edge(edg)%elementinds(2)
      nd1 = edge(edg)%nodeinds(1)
      nd2 = edge(edg)%nodeinds(2)
      eledges = element(el1)%edgeinds(:)
      elnodes = element(el1)%nodeinds(:)
 
      temp = 0.0d0
      hu =zeta0(elnodes)+Still_Depth(elnodes)
      IF(ALL(hu> h_min)) THEN
          alpha = element(el1)%area/(3.0*Vel_Area(edg))
          dhdx = q_alpha*grad_zeta_dx(el1) - (1.0d0-q_alpha)*grad_d_dx(el1)
          dhdy = q_alpha*grad_zeta_dy(el1) - (1.0d0-q_alpha)*grad_d_dy(el1)
          DO j=1,3
              nd_el = elnodes(j)
              temp_u = temp_u -  dt *q(nd_el)*bafuzeta_x(j, el1)*alpha  
               
              temp_v = temp_v -  dt*q(nd_el)*bafuzeta_y(j, el1)*alpha  
          ENDDO                        
      ENDIF
      
      
      IF(el2 .NE. 0) THEN
           eledges = element(el2)%edgeinds(:)
           elnodes = element(el2)%nodeinds(:)
           hu =zeta0(elnodes)+Still_Depth(elnodes)
           IF(ALL(hu> h_min)) THEN
              alpha = element(el2)%area/(3.0*Vel_Area(edg))
              dhdx = q_alpha*grad_zeta_dx(el2) - (1.0d0-q_alpha)*grad_d_dx(el2)
              dhdy = q_alpha*grad_zeta_dy(el2) - (1.0d0-q_alpha)*grad_d_dy(el2)            
              DO j=1,3
                  nd_el = elnodes(j)
                  temp_u = temp_u -  dt*q(nd_el)*bafuzeta_x(j, el2)*alpha 
               
                  temp_v = temp_v - dt*q(nd_el)*bafuzeta_y(j, el2)*alpha  
              ENDDO               
          ENDIF
      ENDIF
         
      velocity1(edg)%u = velocity1(edg)%u + (1.0+alpha_edge(edg))/2.0d0*temp_u
      velocity1(edg)%v = velocity1(edg)%v + (1.0+alpha_edge(edg))/2.0d0*temp_v
      
      velocity_delta(edg)%u = velocity_delta(edg)%u + 1.0/2.0d0*temp_u
      velocity_delta(edg)%v = velocity_delta(edg)%v + 1.0/2.0d0*temp_v
     
  ENDDO    
 !$OMP END DO
  
  
  
!$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)  
  DO nd=1,TriMesh_nodenmb
  
      IF(zeta0(nd)+Still_Depth(nd) <= h_min) THEN
          w_surf1(nd) = 0.0d0
      ELSE
          w_surf1(nd) = w_surf0(nd) + 2.0/(1.0-alpha_h(nd))*dt*q(nd)/(zeta1(nd)+Still_Depth(nd))
      ENDIF
  ENDDO
 !$OMP END DO
 !$OMP END PARALLEL 

  
  w_surf0 = w_surf1
  w_bot0  = w_bot1

END SUBROUTINE H2Ocean_NonHydro_Update_Velocity_ReducedTwoLayer


SUBROUTINE H2Ocean_NonHydro_Build_q_Spare_Matrix_ReducedTwoLayer
  USE Mod_Precision
  USE Mod_TriMesh
  USE Mod_H2Ocean
  USE Mod_FEM
  USE Mod_Sparse_Matrix
  USE Mod_Dirichlet_Boundary
  IMPLICIT NONE
  INTEGER(KIND=INT_KIND)  :: nd, i,el,eledges(3),elnodes(3)
  REAL(KIND=REAL_DP)      :: flux,u_ave, v_ave,discharge, h_flux
  INTEGER(KIND=INT_KIND)  :: j, ic,edg
  INTEGER(KIND=INT_KIND)  :: nd1, nd2 
  REAL(KIND=REAL_DP)      :: vol_sum
    
  INTEGER(KIND=INT_KIND)  :: offset
  REAL(KIND=REAL_DP)      :: alpha_el 
  !Varibles for OpenMP
!$  integer :: chunk, chunksize  
!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
!$  integer :: mythread 

  A_s%values = 0.0d0
  RHS_s = 0.0d0
 
!$OMP   PARALLEL PRIVATE(nd,i,el,elnodes,eledges,u_ave,v_ave,j,edg,offset,   &
!$OMP&                 discharge,h_flux,alpha_el,mythread)  &
!$OMP&           SHARED(A_s,RHS_s,CHUNK,nthreads)  


!$  NTHREADS = OMP_GET_NUM_THREADS()
!$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)   
  DO nd=1,TriMesh_nodenmb
       IF( zeta0(nd) +Still_Depth(nd) <=h_min) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) = 0.0d0
           CYCLE
       ENDIF
       IF( -Still_Depth(nd) >=bath_max) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) = 0.0d0
           CYCLE
       ENDIF       
       
       
       
        IF( node(nd)%ind == 3) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_3(nd)
           CYCLE
       ENDIF
       
       IF( node(nd)%ind == 2 .AND. .NOT. dirichlet_boundary) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2(nd)
           CYCLE
       ENDIF
       
       IF( node(nd)%ind == 2 .AND. dirichlet_boundary .AND. D_bo >=D_bo_Inf) THEN
           CALL H2Ocean_NonHydro_Spare_Matrix_Fill_OpenBoundary_2(nd)
           CYCLE           
       ENDIF       
       IF( node(nd)%ind == 2 .AND. q_dirichlet_boundary .AND. q_D_bo <D_bo_Inf) THEN
           CALL sparse_matrix_locate(nd,nd,A_s,offset)
           A_s%values(offset) = Zeta_Area(nd)*dt
           RHS_s(nd) =  1.0/(1.0+alpha_h(nd))*q_D_bo*Zeta_Area(nd)*dt
           CYCLE            
       ENDIF
 
       
      
       
       DO i = 1, node_ngb_elements(nd)%nmb
          el         = node_ngb_elements(nd)%addresses(i)
          elnodes    = element(el)%nodeinds(:)
          eledges    = element(el)%edgeinds(:)
          alpha_el   = SUM(alpha_h(elnodes))/3.0d0
          u_ave      = (1.0+alpha_el)*SUM(velocity1(eledges)%u)/3.0d0  &
                    +(alpha_el-alpha_el**2)*SUM(velocity_delta(eledges)%u)/3.0d0 
          v_ave      = (1.0+alpha_el)*SUM(velocity1(eledges)%v)/3.0d0  &
                    +(alpha_el-alpha_el**2)*SUM(velocity_delta(eledges)%v)/3.0d0 
          
          IF( MAXVAL(zeta0(elnodes))+MINVAL(Still_Depth(elnodes)) <=h_min) THEN
              CYCLE
          ENDIF 
          DO j=1,3  ! go through 3 edges
              edg = eledges(j)
              IF(ANY(edge(edg)%nodeinds(:)==nd)) THEN   
                 CALL H2Ocean_Compute_edge_discharge(nd,el,edg,u_ave,v_ave,discharge)
                 h_flux = Zeta0(nd)+Still_Depth(nd)
                 RHS_s(nd) = RHS_s(nd) + discharge*h_flux
                 CALL H2Ocean_NonHydro_Spare_Matrix_Fill_Edge_ReducedTwoLayer(nd,el,edg,h_flux)
              ENDIF   
              !Vertical bottom velocity contribution
               !CALL H2Ocean_NonHydro_Spare_Matrix_Fill_W_bot(nd,el,edg)
          ENDDO  ! END DO j=1,3
       ENDDO     ! END DO i = 1, TriMesh.node_ngb_elements(nd)%nmb 
       
       ! vertical verlocity
       RHS_s(nd) = RHS_s(nd) - Zeta_Area(nd)*(w_surf0(nd)) 
       CALL sparse_matrix_locate(nd,nd,A_s,offset)
       A_s%values(offset) = A_s%values(offset) + 2.0/(1.0-alpha_h(nd))*Zeta_Area(nd)*dt/(zeta0(nd)+Still_Depth(nd))        
  ENDDO  ! END DO nd=1,TriMesh_nodenmb      
 !$OMP END DO
 !$OMP END PARALLEL 
END SUBROUTINE H2Ocean_NonHydro_Build_q_Spare_Matrix_ReducedTwoLayer