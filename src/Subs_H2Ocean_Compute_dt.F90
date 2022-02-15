    SUBROUTINE H2Ocean_Compute_dt
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_Parameter
    USE Mod_TriMesh
    IMPLICIT NONE

    CALL H2Ocean_Compute_CFL_flux
    dt = MINVAL(CFL_flux)
    dt = MIN(dt,dt_save)
    
    dt = dt/ max(sqrt(maxval(velocity1(:)%u**2+velocity1(:)%v**2))/10.d0, 1.0d0)
    END SUBROUTINE H2Ocean_Compute_dt


    SUBROUTINE H2Ocean_Compute_CFL_flux
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_Parameter
    USE Mod_TriMesh
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nd, i,el,eledges(3),elnodes(3)
    REAL(KIND=REAL_DP)      :: flux,u_ave, v_ave,discharge, h_flux
    INTEGER(KIND=INT_KIND)  :: j, ic,edg
    INTEGER(KIND=INT_KIND)  :: nd1, nd2 
    REAL(KIND=REAL_DP)      :: depth_ave  
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  

    CFL_flux = 0.0d0
    !$OMP   PARALLEL PRIVATE(nd,i,el,elnodes,eledges,u_ave,v_ave,depth_ave, &
    !$OMP&         j,ic,nd1,nd2,edg,discharge,mythread)  &
    !$OMP&         SHARED(TriMesh_nodenmb,node_ngb_elements,element,edge,dt,&
    !$OMP&               h_min,edge_vectors, &
    !$OMP&         Zeta_Area,zeta0,zeta1,CFL_flux,CHUNK,nthreads)

    !$ NTHREADS = OMP_GET_NUM_THREADS()



    !$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO nd=1,TriMesh_nodenmb  
        DO i = 1, node_ngb_elements(nd)%nmb
            el         = node_ngb_elements(nd)%addresses(i)
            elnodes    = element(el)%nodeinds(:)
            eledges    = element(el)%edgeinds(:)
            depth_ave  = MAX(SUM(Still_Depth(elnodes))/3.0d0 + &
            SUM(zeta1(elnodes))/3.0d0,0.0d0)  

            u_ave      = SUM(velocity1(eledges)%u)/3.0d0  
            v_ave      = SUM(velocity1(eledges)%v)/3.0d0

            DO j=1,3
                IF(elnodes(j) == nd) THEN
                    ic = j;
                ENDIF
            ENDDO

            nd1 = elnodes(MOD(ic,3)+1 )
            nd2 = elnodes(MOD(ic+1,3)+1)

            DO j=1,3

                edg = eledges(j)
                ! edge1:
                IF(ANY(edge(edg)%nodeinds(:)==nd) .AND. &
                ANY(edge(edg)%nodeinds(:)==nd1) ) THEN          
                    IF(edge(edg)%elementinds(1) == el) THEN
                        CALL Cross_Product(u_ave,v_ave,edge_vectors(edg,1)%dx,  &
                        edge_vectors(edg,1)%dy,edge_vectors(edg,1)%length,discharge)
                        discharge = abs(discharge) + edge_vectors(edg,1)%length*SQRT(g*depth_ave)
                    ENDIF
                    IF(edge(edg)%elementinds(2) == el) THEN
                        CALL Cross_Product(u_ave,v_ave,edge_vectors(edg,2)%dx,  &
                        edge_vectors(edg,2)%dy,edge_vectors(edg,2)%length,discharge)
                        discharge = abs(discharge) + edge_vectors(edg,2)%length*SQRT(g*depth_ave)
                    ENDIF 
                    CFL_flux(nd) = CFL_flux(nd) + discharge 

                ENDIF   

                ! edge2
                IF(ANY(edge(edg)%nodeinds(:)==nd) .AND. &
                ANY(edge(edg)%nodeinds(:)==nd2) ) THEN

                    IF(edge(edg)%elementinds(1) == el) THEN
                        CALL Cross_Product(u_ave,v_ave,edge_vectors(edg,1)%dx,  &
                        edge_vectors(edg,1)%dy,edge_vectors(edg,1)%length,discharge)
                        discharge = abs(discharge) + edge_vectors(edg,1)%length*SQRT(g*depth_ave)
                    ENDIF

                    IF(edge(edg)%elementinds(2) == el) THEN
                        CALL Cross_Product(u_ave,v_ave,edge_vectors(edg,2)%dx,  &
                        edge_vectors(edg,2)%dy,edge_vectors(edg,2)%length,discharge)
                        discharge = abs(discharge) + edge_vectors(edg,2)%length*SQRT(g*depth_ave)
                    ENDIF  
                    CFL_flux(nd) = CFL_flux(nd) + discharge 
                ENDIF 
            ENDDO  ! END DO j=1,3
        ENDDO     ! END DO i = 1, TriMesh.node_ngb_elements(nd).nmb

        CFL_flux(nd) = Zeta_Area(nd) *max_CFL /CFL_flux(nd)

    ENDDO  ! END DO nd=1,TriMesh%nodenmb

    !$OMP END DO
    !$OMP END PARALLEL  
    END SUBROUTINE H2Ocean_Compute_CFL_flux