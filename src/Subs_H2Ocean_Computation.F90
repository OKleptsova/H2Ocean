    SUBROUTINE H2Ocean_Computation
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

    USE Mod_All_Variable
    USE Mod_Dirichlet_Boundary
    USE Mod_H2Ocean
    USE MOD_TriMesh
    IMPLICIT NONE
    INTEGER  :: max_loc(1)


    ! Step 1:
    CALL H2Ocean_Compute_WaterLevel_Gradient

    !Step 2:

    !Compute the water depth, h_wd, for wetting and drying of velocity points
    CALL H2Ocean_Compute_h_wd 


    ! Compute the water depth at velocity face 
    CALL H2Ocean_Compute_h_f

    IF(enable_advection) THEN   
        CALL H2Ocean_Compute_h_upwd 
    ENDIF  


    IF( (enable_advection .AND. advection_scheme ==2)  .OR. (enable_viscosity .AND. viscosity_scheme ==2) ) THEN
        CALL H2Ocean_Compute_Vel_Grad 
    ENDIF

    IF(enable_viscosity .AND. viscosity_scheme ==2) THEN
        CALL H2Ocean_Compute_Vel_2nd_Grad 
    ENDIF

    IF(Dirichlet_boundary) THEN
        CALL H2Ocean_Interp_D_bo(time_run)
    ENDIF

    IF(u_Dirichlet_boundary) THEN
        CALL H2Ocean_Interp_u_D_bo(time_run)
    ENDIF 


    !Step 3:

    CALL H2Ocean_Compute_New_Velocity_Explicit


    IF(enable_nonhydro) THEN
        IF(q_Dirichlet_boundary) THEN
            CALL H2Ocean_Interp_q_D_bo(time_run+dt)
        ENDIF
        CALL H2Ocean_Nonhydro_Process
    ENDIF

    CALL H2Ocean_Compute_Radiation


    CALL H2Ocean_Compute_New_WaterLevel_Explicit



    CALL H2Ocean_Compute_Max_WaterLevel


    !Step 4:
    CALL H2Ocean_Update

    END SUBROUTINE H2Ocean_Computation


    SUBROUTINE H2Ocean_Report
    USE Mod_All_Variable
    USE Mod_H2Ocean
    IMPLICIT NONE
    IF(MOD(i_run,print_int)==0) THEN
        WRITE(*,'(A3,I8.1,A4,f10.4,A12,e10.3,A9,e10.3,A5,e10.3)') 'i=',i_run,'t= ',time_run, "Max_Zeta = ", zetamax, &
        "Max_Vel=",   SQRT(MAXVAL(velocity1(:)%u**2+velocity1(:)%v**2)),'dt=',dt 
    END IF

    END SUBROUTINE H2Ocean_Report






    SUBROUTINE H2Ocean_Compute_h_upwd
    USE Mod_TriMESH
    USE Mod_H2Ocean 
    USE Mod_Precision
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: elnodes1(3),elnodes2(3),edges1(3),edges2(3)
    INTEGER(KIND=INT_KIND)  :: edg, el1, el2, nd1,nd2
    REAL(kind=REAL_DP)           :: hu1,hu2, q_star, hu3, hu4,hu(3)
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    !$OMP   PARALLEL PRIVATE(edg,el1,el2,elnodes1,elnodes2,q_star,hu3,hu4 , &
    !$OMP&                   mythread)  &
    !$OMP&           SHARED(h_upwd,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg=1, TriMesh_edgenmb
        el1 =  edge(edg)%elementinds(1)
        el2 =  edge(edg)%elementinds(2)
        elnodes1 = element(el1)%nodeinds(:)

        q_star = edge(edg)%dx*velocity1(edg)%v - &
        edge(edg)%dy*velocity1(edg)%u  

        IF(el2 .ne. 0 ) THEN
            elnodes2 =   element(el2)%nodeinds(:)
            edges2   =   element(el2)%edgeinds(:)
            hu3 = SUM(zeta1(elnodes1))/3.0d0  + SUM(Still_Depth(elnodes1))/3.0d0
            hu4 = SUM(zeta1(elnodes2))/3.0d0  + SUM(Still_Depth(elnodes2))/3.0d0
            IF(q_star>0.0) then
                h_upwd(edg) = hu4
            ELSEIF(q_star<0.0d0) THEN
                h_upwd(edg) = hu3
            ELSE
                h_upwd(edg) = (hu3 + hu4)/2.
            ENDIF
        ELSE
            hu3 = SUM(zeta0(elnodes1))/3.0d0 + SUM(Still_Depth(elnodes1))/3.0d0
            h_upwd(edg) = hu3
        ENDIF

        h_upwd(edg) = MAX(h_upwd(edg),0.0d0)
    ENDDO 
    !$OMP END DO   
    !$OMP END PARALLEL            
    END SUBROUTINE H2Ocean_Compute_h_upwd 


    !SUBROUTINE H2Ocean_Compute_h_upwd
    !    USE Mod_TriMESH
    !    USE Mod_H2Ocean 
    !    USE Mod_Precision
    !    IMPLICIT NONE
    !    INTEGER(KIND=INT_KIND)  :: elnodes1(3),elnodes2(3),edges1(3),edges2(3)
    !    INTEGER(KIND=INT_KIND)  :: edg, el1, el2, nd1,nd2
    !    REAL(kind=REAL_DP)           :: hu1,hu2, q_star, hu3, hu4,hu(3),h_ave
    !    !Varibles for OpenMP
    !!$  integer :: chunk, chunksize  
    !!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !!$  integer :: mythread  
    !
    !    
    !!$OMP   PARALLEL PRIVATE(edg,el1,el2,elnodes1,elnodes2,q_star,hu3,hu4 , &
    !!$OMP&                   mythread)  &
    !!$OMP&           SHARED(h_upwd,CHUNK,nthreads)  
    !
    !
    !!$  NTHREADS = OMP_GET_NUM_THREADS()
    !!$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    !
    !    DO edg=1, TriMesh_edgenmb
    !        el1 =  edge(edg)%elementinds(1)
    !        el2 =  edge(edg)%elementinds(2)
    !        elnodes1 = element(el1)%nodeinds(:)
    !
    !        q_star = edge(edg)%dx*velocity0(edg)%v - &
    !                 edge(edg)%dy*velocity0(edg)%u  
    !           
    !        IF(el2 .ne. 0 ) THEN
    !            elnodes2 =   element(el2)%nodeinds(:)
    !            edges2   =   element(el2)%edgeinds(:)
    !            hu3 = SUM(zeta1(elnodes1)) - SUM(zeta1(edge(edg)%nodeinds(:)))
    !            hu4 = SUM(zeta1(elnodes2)) - SUM(zeta1(edge(edg)%nodeinds(:)))  
    !            h_ave = element(el1)%area/3.0d0/Vel_Area(edg)*SUM(Still_Depth(elnodes1))/3.0d0 &
    !                 +  element(el2)%area/3.0d0/Vel_Area(edg)*SUM(Still_Depth(elnodes2))/3.0d0
    !            
    !            IF(q_star>0.0) then
    !                h_upwd(edg) = MAX(hu4 + h_ave,0.0d0)
    !            ELSEIF(q_star<0.0) THEN
    !                h_upwd(edg) = MAX(hu3 + h_ave,0.0d0)
    !            ELSE
    !                h_upwd(edg) = MAX((hu3 + hu4)/2. + h_ave,0.0d0)
    !            ENDIF
    !        ELSE
    !            hu3   = SUM(zeta1(elnodes1)) - SUM(zeta1(edge(edg)%nodeinds(:)))
    !            h_ave = SUM(Still_Depth(elnodes1))/3.0d0
    !            h_upwd(edg) = MAX(hu3 + h_ave,0.0d0)
    !        ENDIF
    !    ENDDO 
    !!$OMP END DO   
    !!$OMP END PARALLEL            
    !END SUBROUTINE H2Ocean_Compute_h_upwd 


    SUBROUTINE H2Ocean_Compute_New_WaterLevel_Explicit 
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_Parameter
    USE Mod_All_Variable
    USE Mod_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nd, i,el,eledges(3),elnodes(3)
    REAL(KIND=REAL_DP)      :: flux,u_ave, v_ave,discharge, h_flux
    INTEGER(KIND=INT_KIND)  :: j, ic,edg
    INTEGER(KIND=INT_KIND)  :: nd1, nd2 
    REAL(KIND=REAL_DP)      :: vol_sum



    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  



    !$OMP   PARALLEL PRIVATE(nd,flux,vol_sum,i,el,elnodes,eledges,u_ave,v_ave, &
    !$OMP&         j,ic,nd1,nd2,edg,discharge,h_flux,mythread)  &
    !$OMP&         SHARED(TriMesh_nodenmb,node_ngb_elements,element,edge,  &
    !$OMP&          dt,h_min,edge_vectors,&
    !$OMP&          flux_limiter_method, &
    !$OMP&          D_bo,Zeta_Area,Zeta_Wet_Area,zeta0,zeta1,CHUNK,nthreads)

    !$ NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO nd=1,TriMesh_nodenmb

        IF(dirichlet_boundary .AND. node(nd)%ind==2 .AND. D_bo<D_bo_Inf) THEN
            zeta1(nd) = D_bo
            CYCLE    
        ENDIF

        IF( -Still_Depth(nd)>= bath_max) THEN
            CYCLE
        ENDIF

        flux = 0.0d0
        vol_sum = 0.0d0
        DO i = 1, node_ngb_elements(nd)%nmb
            el         = node_ngb_elements(nd)%addresses(i)
            elnodes    = element(el)%nodeinds(:)
            eledges    = element(el)%edgeinds(:)
            u_ave      = SUM(velocity1(eledges)%u)/3.0d0
            v_ave      = SUM(velocity1(eledges)%v)/3.0d0
            IF( MAXVAL(zeta0(elnodes))+MINVAL(Still_Depth(elnodes)) <=h_min) THEN
                CYCLE
            ENDIF
            !vol_sum = vol_sum + element(el)%area/3.0d0

            DO j=1,3  ! go through 3 edges
                edg = eledges(j)
                IF(ANY(edge(edg)%nodeinds(:)==nd) ) THEN        
                    CALL H2Ocean_Compute_edge_discharge(nd,el,edg,u_ave,v_ave,discharge)
                    CALL H2Ocean_Flux_Limiter(nd,el,edg,discharge,flux_limiter_method,h_flux)         
                    flux = flux + discharge*dt*h_flux
                ENDIF   
            ENDDO  ! END DO j=1,3
        ENDDO     ! END DO i = 1, TriMesh.node_ngb_elements(nd).nmb

        IF(vol_sum>= Zeta_Area(nd)/1000.d0) THEN
            Zeta_Wet_Area(nd) = vol_sum
        ENDIF
        zeta1(nd) = (zeta0(nd)*Zeta_Area(nd) + flux)   &
        /(Zeta_Area(nd) + radiation(nd)*dt)

    ENDDO  ! END DO nd=1,TriMesh.nodenmb       
    !$OMP END DO
    !$OMP END PARALLEL  
    END SUBROUTINE H2Ocean_Compute_New_WaterLevel_Explicit






    SUBROUTINE H2Ocean_Compute_edge_discharge(nd,el,edg,u_ave,v_ave,discharge)
    USE Mod_Precision
    USE Mod_TriMesh
    USE Mod_H2Ocean
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND),INTENT(IN) :: nd,el,edg
    REAL(KIND=REAL_DP),INTENT(IN)     :: u_ave, v_ave
    REAL(KIND=REAL_DP),INTENT(OUT)    :: discharge
    REAL(KIND=REAL_DP)                :: dx, dy, L 

    IF(edge(edg)%elementinds(1) == el) THEN
        dx = edge_vectors(edg,1)%dx
        dy = edge_vectors(edg,1)%dy
        L  = edge_vectors(edg,1)%length

        IF(edge(edg)%nodeinds(1) == nd) THEN
            dx = -edge_vectors(edg,1)%dx
            dy = -edge_vectors(edg,1)%dy
        ENDIF     
        CALL Cross_Product(u_ave,v_ave,dx,dy,L,discharge)

    ENDIF
    IF(edge(edg)%elementinds(2) == el) THEN
        dx = edge_vectors(edg,2)%dx
        dy = edge_vectors(edg,2)%dy
        L  = edge_vectors(edg,2)%length
        IF(edge(edg)%nodeinds(1) == nd) THEN
            dx = -edge_vectors(edg,2)%dx
            dy = -edge_vectors(edg,2)%dy
        ENDIF     
        CALL Cross_Product(u_ave,v_ave,dx,dy,L,discharge)
    ENDIF   

    END SUBROUTINE H2Ocean_Compute_edge_discharge



    SUBROUTINE H2Ocean_Compute_Max_WaterLevel
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_Parameter
    USE Mod_All_Variable
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: nd

    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread    

    ! Compute zeta_max

    !$OMP   PARALLEL PRIVATE(nd,mythread)  &
    !$OMP&           SHARED(h_min,zetamax,zeta_min,zeta_max,CHUNK,nthreads)  

    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    DO nd=1,TriMesh_nodenmb 
        IF ( zeta1(nd)+ Still_Depth(nd)>= 1.0e-6) THEN    
            zeta_max(nd) = MAX(zeta_max(nd), zeta1(nd))
        ENDIF 
    ENDDO  
    !$OMP END DO 


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    DO nd=1,TriMesh_nodenmb
        IF ( zeta1(nd)+ Still_Depth(nd)>= h_min) THEN
            zeta_min(nd) = MIN(zeta_min(nd), zeta1(nd)) 
        ENDIF 
    ENDDO  
    !$OMP END DO 


    !$OMP END PARALLEL 

    !zetamax = MIN(MAXVAL(zeta1),MAXVAL(zeta_max))
    zetamax = -9999.d0
    DO nd=1,TriMesh_nodenmb
        IF ( zeta1(nd)+ Still_Depth(nd)>= h_min) THEN
            zetamax = MAX(zetamax,zeta1(nd))
        ENDIF 
    ENDDO



    END SUBROUTINE H2Ocean_Compute_Max_WaterLevel





    SUBROUTINE H2Ocean_Compute_h_wd
    USE Mod_TriMESH
    USE Mod_H2Ocean 
    USE Mod_Precision
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: edg, el1, el2, elnodes1(3),elnodes2(3)
    REAL(kind=REAL_DP)      :: zeta_el1, zeta_el2, depth_min,hu3,hu4
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    !$OMP   PARALLEL PRIVATE(edg,el1,el2,elnodes1,elnodes2,depth_min,mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,element,h_wd,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg = 1, TriMesh_edgenmb
        el1 = edge(edg)%elementinds(1)
        el2 = edge(edg)%elementinds(2)
        elnodes1 = element(el1)%nodeinds(:) 
        h_wd(edg) = MAXVAL(zeta1(elnodes1))      
        IF(el2 .ne. 0 ) THEN
            elnodes2   =  element(el2)%nodeinds(:)
            h_wd(edg)  = MAX(h_wd(edg),MAXVAL(zeta1(elnodes2)))
        ENDIF

        depth_min = MINVAL(Still_Depth(elnodes1))

        IF(el2 .NE. 0 ) THEN
            elnodes2 =  element(el2)%nodeinds(:)
            depth_min = MIN(MINVAL(Still_Depth(elnodes2)),depth_min)
        ENDIF   
        h_wd(edg) = MAX(h_wd(edg) +depth_min, 0.0d0)      
    ENDDO
    !$OMP END DO   
    !$OMP END PARALLEL            
    END SUBROUTINE H2Ocean_Compute_h_wd



    SUBROUTINE H2Ocean_Compute_h_f
    USE Mod_TriMESH
    USE Mod_H2Ocean 
    USE Mod_Precision
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: edg, el1, el2, elnodes1(3),elnodes2(3)
    REAL(kind=REAL_DP)      :: zeta_el1, zeta_el2, depth_min
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    !$OMP   PARALLEL PRIVATE(edg,el1,el2,elnodes1,elnodes2,zeta_el1,zeta_el2, &
    !$OMP&                   mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,element,h_f,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)


    DO edg = 1, TriMesh_edgenmb
        el1 = edge(edg)%elementinds(1)
        el2 = edge(edg)%elementinds(2)
        elnodes1 = element(el1)%nodeinds(:)
        IF(el2 .ne. 0 ) THEN
            elnodes2    =  element(el2)%nodeinds(:)
            zeta_el1    =  SUM(zeta1(elnodes1))/3.0d0 + SUM(Still_Depth(elnodes1))/3.0d0
            zeta_el2    =  SUM(zeta1(elnodes2))/3.0d0 + SUM(Still_Depth(elnodes2))/3.0d0
            h_f(edg)  = (element(el1)%area*zeta_el1+  element(el2)%area*zeta_el2)  &
            /(element(el1)%area+element(el2)%area)
        ELSE
            h_f(edg)   =  SUM(zeta1(elnodes1))/3.0d0 + SUM(Still_Depth(elnodes1))/3.0d0
        ENDIF     

        h_f(edg) = MAX(h_f(edg), 0.0d0)
    ENDDO 
    !$OMP END DO   
    !$OMP END PARALLEL          
    END SUBROUTINE H2Ocean_Compute_h_f


    ! SUBROUTINE H2Ocean_Compute_h_f
    !    USE Mod_TriMESH
    !    USE Mod_H2Ocean 
    !    USE Mod_Precision
    !    IMPLICIT NONE
    !    INTEGER(KIND=INT_KIND)  :: edg, el1, el2, elnodes1(3),elnodes2(3)
    !    REAL(kind=REAL_DP)      :: h_el1, h_el2, h_edg
    !    !Varibles for OpenMP
    !!$  integer :: chunk, chunksize  
    !!$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !!$  integer :: mythread  
    !
    !    
    !!$OMP   PARALLEL PRIVATE(edg,el1,el2,elnodes1,elnodes2,h_el1, h_el2, h_edg, &
    !!$OMP&                   mythread)  &
    !!$OMP&           SHARED(TriMesh_edgenmb,edge,element,h_f,CHUNK,nthreads)  
    !
    !
    !!$  NTHREADS = OMP_GET_NUM_THREADS()
    !!$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    !
    ! 
    !    DO edg = 1, TriMesh_edgenmb
    !        el1 = edge(edg).elementinds(1)
    !        el2 = edge(edg).elementinds(2)
    !        elnodes1 = element(el1).nodeinds(:)
    !        IF(el2 .ne. 0 ) THEN
    !            elnodes2    =  element(el2)%nodeinds(:)
    !            h_el1    =  SUM(zeta1(elnodes1))/3.0d0 + SUM(Still_Depth(elnodes1))/3.0d0
    !            h_el2    =  SUM(zeta1(elnodes2))/3.0d0 + SUM(Still_Depth(elnodes2))/3.0d0
    !            h_edg    =  SUM(zeta1(edge(edg)%nodeinds))/2.0d0  &
    !                      + SUM(Still_Depth(edge(edg)%nodeinds))/2.0d0
    !            h_f(edg)  = (element(el1)%area*(h_el1/3.+h_edg*2/3.)+  &
    !                         element(el2)%area*(h_el2/3.+h_edg*2/3.))  &
    !                         /(element(el1)%area+element(el2)%area)
    !        ELSE
    !            h_f(edg)   =   h_el1/3.+h_edg*2/3.
    !        ENDIF     
    !        
    !        h_f(edg) = MAX(h_f(edg), 0.0d0)
    !    ENDDO 
    !!$OMP END DO   
    !!$OMP END PARALLEL          
    !END SUBROUTINE H2Ocean_Compute_h_f


    SUBROUTINE H2Ocean_Compute_Vel_Grad
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND) :: edg, j,edg_ngb
    REAL(KIND=REAL_DP)     :: nx,ny, tmpu,tmpv

    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  
    grad_vel_dx = 0.0d0
    grad_vel_dy = 0.0d0

    !$OMP   PARALLEL PRIVATE(edg,j,edg_ngb, &
    !$OMP&                   mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,edge_ngb_edges,grad_vel_dx,&
    !$OMP&                  grad_vel_dy, &
    !$OMP&                  velocity1,vel_coef_dx,vel_coef_dy,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg=1,TriMesh_edgenmb
        grad_vel_dx(edg,1) = vel_coef_dx(edg,5)*velocity1(edg)%u
        grad_vel_dx(edg,2) = vel_coef_dx(edg,5)*velocity1(edg)%v
        grad_vel_dy(edg,1) = vel_coef_dy(edg,5)*velocity1(edg)%u
        grad_vel_dy(edg,2) = vel_coef_dy(edg,5)*velocity1(edg)%v      

        IF(MINVAL(edge(edg)%elementinds(:)) == 0) THEN
            ! With only one neighbour element, two neighbour edges
            DO j=1,2
                edg_ngb = edge_ngb_edges(edg,j) 
                grad_vel_dx(edg,1) = grad_vel_dx(edg,1)+vel_coef_dx(edg,j)*velocity1(edg_ngb)%u 
                grad_vel_dx(edg,2) = grad_vel_dx(edg,2)+vel_coef_dx(edg,j)*velocity1(edg_ngb)%v
                grad_vel_dy(edg,1) = grad_vel_dy(edg,1)+vel_coef_dy(edg,j)*velocity1(edg_ngb)%u
                grad_vel_dy(edg,2) = grad_vel_dy(edg,2)+vel_coef_dy(edg,j)*velocity1(edg_ngb)%v   
            ENDDO

        ELSE
            DO j=1,4
                edg_ngb = edge_ngb_edges(edg,j) 
                grad_vel_dx(edg,1) = grad_vel_dx(edg,1)+vel_coef_dx(edg,j)*velocity1(edg_ngb)%u 
                grad_vel_dx(edg,2) = grad_vel_dx(edg,2)+vel_coef_dx(edg,j)*velocity1(edg_ngb)%v
                grad_vel_dy(edg,1) = grad_vel_dy(edg,1)+vel_coef_dy(edg,j)*velocity1(edg_ngb)%u
                grad_vel_dy(edg,2) = grad_vel_dy(edg,2)+vel_coef_dy(edg,j)*velocity1(edg_ngb)%v   
            ENDDO      
        ENDIF   
    ENDDO
    !$OMP END DO   
    !$OMP END PARALLEL 
    END SUBROUTINE H2Ocean_Compute_Vel_Grad

   SUBROUTINE H2Ocean_Compute_Vel_2nd_Grad
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND) :: edg, j,edg_ngb
    REAL(KIND=REAL_DP)     :: nx,ny, tmpu,tmpv

    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  
    
    ! vel_2nd_grad(;,1)  du_dxx
    ! vel_2nd_grad(;,2)  du_dyy
    ! vel_2nd_grad(;,3)  dv_dxx 
    ! vel_2nd_grad(;,4)  dv_dyy
    vel_2nd_grad = 0.0d0

    !$OMP   PARALLEL PRIVATE(edg,j,edg_ngb, &
    !$OMP&                   mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,edge_ngb_edges,  &
    !$OMP&                   grad_vel_dx,grad_vel_dy,vel_2nd_grad,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg=1,TriMesh_edgenmb
        vel_2nd_grad(edg,1) = vel_coef_dx(edg,5)*grad_vel_dx(edg,1) 
        vel_2nd_grad(edg,2) = vel_coef_dy(edg,5)*grad_vel_dy(edg,1)       
        vel_2nd_grad(edg,3) = vel_coef_dx(edg,5)*grad_vel_dx(edg,2) 
        vel_2nd_grad(edg,4) = vel_coef_dy(edg,5)*grad_vel_dy(edg,2) 
        IF(MINVAL(edge(edg)%elementinds(:)) == 0) THEN
            ! With only one neighbour element, two neighbour edges
            DO j=1,2
                edg_ngb = edge_ngb_edges(edg,j) 
                vel_2nd_grad(edg,1) = vel_2nd_grad(edg,1) +vel_coef_dx(edg,j)*grad_vel_dx(edg_ngb,1) 
                vel_2nd_grad(edg,2) = vel_2nd_grad(edg,2) +vel_coef_dy(edg,j)*grad_vel_dy(edg_ngb,1)    
                vel_2nd_grad(edg,3) = vel_2nd_grad(edg,3) +vel_coef_dx(edg,j)*grad_vel_dx(edg_ngb,2)  
                vel_2nd_grad(edg,4) = vel_2nd_grad(edg,4) +vel_coef_dy(edg,j)*grad_vel_dy(edg_ngb,2)   
            ENDDO

        ELSE
            DO j=1,4
                edg_ngb = edge_ngb_edges(edg,j) 
                vel_2nd_grad(edg,1) = vel_2nd_grad(edg,1) +vel_coef_dx(edg,j)*grad_vel_dx(edg_ngb,1) 
                vel_2nd_grad(edg,2) = vel_2nd_grad(edg,2) +vel_coef_dy(edg,j)*grad_vel_dy(edg_ngb,1)    
                vel_2nd_grad(edg,3) = vel_2nd_grad(edg,3) +vel_coef_dx(edg,j)*grad_vel_dx(edg_ngb,2)  
                vel_2nd_grad(edg,4) = vel_2nd_grad(edg,4) +vel_coef_dy(edg,j)*grad_vel_dy(edg_ngb,2)  
            ENDDO      
        ENDIF   
    ENDDO
    !$OMP END DO   
    !$OMP END PARALLEL 
   END SUBROUTINE H2Ocean_Compute_Vel_2nd_Grad
   
   
   
   

    SUBROUTINE H2Ocean_Compute_WaterLevel_Gradient
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_FEM
    USE Mod_All_Variable
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)   :: el, elnodes(3), j, j_max_eta
    REAL(kind = REAL_DP )    :: hu(3), zeta_tmp(3), dep_level, zeta_level
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  

    grad_zeta_dx = 0.0d0
    grad_zeta_dy = 0.0d0

    !$OMP   PARALLEL PRIVATE(el,elnodes,hu,dep_level,j_max_eta,zeta_level,mythread)  &
    !$OMP&           SHARED(TriMesh_elementnmb,element,grad_zeta_dx,grad_zeta_dy, &
    !$OMP&           Still_Depth,zeta0,zeta1, &
    !$OMP&           h_min,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_elementnmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    DO el = 1,TriMesh_elementnmb
        elnodes = element(el)%nodeinds(:)
        hu = zeta1(elnodes)+Still_Depth(elnodes) 
        IF(ALL(hu<=h_min)) THEN
            grad_zeta_dx(el) = 0.0d0
            grad_zeta_dy(el) = 0.0d0
            CYCLE
        ELSEIF(ALL(hu>h_min)) THEN
            grad_zeta_dx(el) = SUM(bafuzeta_x(:, el)*zeta1(elnodes))
            grad_zeta_dy(el) = SUM(bafuzeta_y(:, el)*zeta1(elnodes))       
        ELSE
            zeta_tmp   =  0.0d0
            zeta_level =  MINVAL(zeta1(elnodes))        ! the lowest water  level
            dep_level  = -MINVAL(Still_Depth(elnodes))  ! the highest bottom level
            
            DO j=1,3
                IF(zeta1(elnodes(j))+Still_Depth(elnodes(j)) >h_min  ) THEN  ! if the node is wet and 
                    zeta_level = MAX(zeta_level,zeta1(elnodes(j)) )  ! Get the highest water level within the wet nodes
                ENDIF
            ENDDO
            
            DO j=1,3
               IF (ABS(zeta1(elnodes(j)) - zeta_level) <1.0e-5) THEN
                   j_max_eta = j
               ENDIF 
            ENDDO
            
            DO j=1,3
                zeta_tmp(j) =  MIN(MAX(- Still_Depth(elnodes(j_max_eta)),zeta1(elnodes(j))), zeta_level)
            ENDDO 
            grad_zeta_dx(el) = SUM(bafuzeta_x(:, el)*zeta_tmp)
            grad_zeta_dy(el) = SUM(bafuzeta_y(:, el)*zeta_tmp)         
        ENDIF
    ENDDO       
    !$OMP END DO   
    !$OMP END PARALLEL   
    
    !write(*,*) maxval(grad_zeta_dx),maxval(grad_zeta_dy)

    END SUBROUTINE H2Ocean_Compute_WaterLevel_Gradient
    
    
    SUBROUTINE H2Ocean_Compute_New_Velocity_Explicit
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_Parameter
    USE Mod_All_Variable
    USE Mod_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)       :: edg, elin,el, i,j, eledges(3),elnodes(3)
    REAL(kind=REAL_DP)           :: vol,vol_sum
    REAL(kind=REAL_DP)           :: Fx, Fy, adv_u,adv_v,coriolis_u,coriolis_v
    REAL(kind=REAL_DP)           :: vis_u,vis_v
    REAL(kind=REAL_DP)           :: r_vel
    !Solid boundary
    REAL(kind=REAL_DP)           :: nx, ny,tmpu, tmpv  

    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    !$OMP   PARALLEL PRIVATE(edg,Fx,Fy,adv_u,adv_v,vis_u,vis_v,elin,el,vol,vol_sum, &
    !$OMP&                   nx,ny,tmpu,tmpv,r_vel,coriolis_u,coriolis_v,mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,node,edge,element,velocity00,velocity0,velocity1, &
    !$OMP&                  h_wd,h_min,enable_advection, &
    !$OMP&                 enable_coriolis,enable_friction,dt,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg=1,TriMesh_edgenmb


        IF(ALL(node(edge(edg)%nodeinds)%ind == 2) .AND. u_dirichlet_boundary .AND. u_D_bo <D_bo_Inf) THEN
            velocity1(edg)%u = u_D_bo
            velocity1(edg)%v = 0.0d0
            CYCLE
        ENDIF

        Fx     = 0.0d0
        Fy     = 0.0d0
        adv_u  = 0.0d0
        adv_v  = 0.0d0
        vis_u  = 0.0d0
        vis_v  = 0.0d0

        r_vel  = 1.0d0
        coriolis_u = 0.0d0
        coriolis_v = 0.0d0
        ! 1. Water level gradient
        IF( h_wd(edg) <= h_min) THEN
            velocity1(edg)%u = 0.0d0
            velocity1(edg)%v = 0.0d0
            CYCLE
        ENDIF

        IF( h_f(edg) <= h_min) THEN
            velocity1(edg)%u = 0.0d0
            velocity1(edg)%v = 0.0d0
            CYCLE
        ENDIF


        !      vol_sum = 0.0d0
        DO elin = 1,2
            el  =  edge(edg)%elementinds(elin)
            IF(el == 0) CYCLE        
            IF(MAXVAL(zeta1(element(el)%nodeinds(:))) &
              +MINVAL(Still_Depth(element(el)%nodeinds(:))) < h_min) THEN  
                CYCLE
            ENDIF

            vol =  element(el)%area
            Fx  =  Fx - g*grad_zeta_dx(el)*vol/3.d0
            Fy  =  Fy - g*grad_zeta_dy(el)*vol/3.d0
        ENDDO


        Fx = Fx/Vel_Area(edg)
        Fy = Fy/Vel_Area(edg)    


        IF(enable_advection) THEN
            CALL H2Ocean_Compute_Advection(edg,adv_u,adv_v) 
        ENDIF

        IF(enable_coriolis) THEN
            CALL H2Ocean_Compute_Coriolis(edg,coriolis_u,coriolis_v) 
        ENDIF

        IF(enable_viscosity) THEN
            CALL H2Ocean_Compute_viscosity(edg,vis_u,vis_v) 
        ENDIF      

        IF(enable_friction) THEN
            SELECT CASE(friction_method)
            CASE(1)
                r_vel= 1.0d0 + dt*Cd*SQRT((velocity0(edg)%u)**2 + (velocity0(edg)%v)**2) &
                /MAX(h_f(edg), 1.0e-4)
            CASE(2)
                r_vel= 1.0d0 + dt*g*C_m**2*SQRT((velocity0(edg)%u)**2 + (velocity0(edg)%v)**2) &
                /MAX(h_f(edg), 1.0e-4)**(4./3.)
                
            CASE DEFAULT

            END SELECT
        ENDIF
 
        velocity1(edg)%u   = (velocity0(edg)%u +  Fx*dt +  adv_u*dt + coriolis_u*dt +vis_u*dt) &
        &     /r_vel 
        velocity1(edg)%v   = (velocity0(edg)%v +  Fy*dt +  adv_v*dt + coriolis_v*dt +vis_v*dt) &
        &     /r_vel


        !IF( SQRT((velocity1(edg)%u)**2 + (velocity1(edg)%v)**2)  > 10.0d0) THEN
        !    write(*,*) 'before'
        !    write(*,*) 'edg=',edg, SQRT((velocity0(edg)%u)**2 + (velocity0(edg)%v)**2),SQRT((velocity1(edg)%u)**2 + (velocity1(edg)%v)**2), &
        !    Fx,Fy,adv_u,adv_v,coriolis_u,coriolis_v,vis_u,vis_v,edge(edg)%nodeinds(1),edge(edg)%nodeinds(2)
        !    write(*,*) 'r_vel=',r_vel, h_f(edg),MAX(h_f(edg), 1.0e-4)**(4./3.)
        !ENDIF

        ! Solid Boundary Conditions
        IF (MINVAL(edge(edg)%elementinds(:)) == 0 .and. &
        COUNT(node(edge(edg)%nodeinds(:))%ind==1) == 2) THEN 
            nx = edge(edg)%ndx
            ny = edge(edg)%ndy 
            tmpu = ny*(ny*velocity1(edg)%u  - nx*velocity1(edg)%v)
            tmpv = nx*(nx*velocity1(edg)%v  - ny*velocity1(edg)%u)
            velocity1(edg)%u  = tmpu
            velocity1(edg)%v  = tmpv
        ENDIF      

    ENDDO
    !$OMP END DO   
    !$OMP END PARALLEL   
    END SUBROUTINE H2Ocean_Compute_New_Velocity_Explicit


    SUBROUTINE H2Ocean_Compute_Advection(edg,adv_u,adv_v)
    USE Mod_Precision
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_All_Variable
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND), INTENT(IN)   :: edg
    REAL(KIND=REAL_DP),INTENT(OUT)       :: adv_u, adv_v

    INTEGER(KIND=INT_KIND)               :: elin,el,elnodes(3),j, edg_ngb
    REAL(KIND=REAL_DP)                   :: vol,discharge,u_ave,v_ave,dx,dy,length
    REAL(KIND=REAL_DP)                   :: x_mid,y_mid, u_mid,v_mid
    REAL(KIND=REAL_DP)                   :: u_mid1,v_mid1,u_mid2,v_mid2,discharge1,discharge2
    adv_u = 0.0d0
    adv_v = 0.0d0

    SELECT CASE(advection_scheme)
    CASE(1)
        adv_u = 0.0d0
        adv_v = 0.0d0
        DO elin = 1,2
            el      = edge(edg)%elementinds(elin)
            IF (el == 0) CYCLE	
            elnodes =  element(el)%nodeinds(:)
            IF( h_wd(edg) <= h_min )  CYCLE
            IF( h_f(edg)  <= h_min )  CYCLE


            vol = element(el)%area 
            DO j= elin*2-1, elin*2              
                discharge = 0.0d0;  
                edg_ngb =  edge_ngb_edges(edg,j)

                IF( h_wd(edg_ngb) <= h_min )  CYCLE
                IF( h_f(edg_ngb)  <= h_min )  CYCLE             

                u_ave =  velocity0(edg_ngb)%u 
                v_ave =  velocity0(edg_ngb)%v    

                CALL Cross_Product(u_ave,v_ave,edge_vectors(edg,j+2)%dx,  &
                edge_vectors(edg,j+2)%dy,edge_vectors(edg,j+2)%length,&
                discharge)
                discharge = discharge*h_upwd(edg_ngb)
                IF(discharge >0.0d0) THEN
                    adv_u  =  adv_u + discharge*(velocity0(edg_ngb)%u  &
                    -velocity0(edg)%u) /vol  /MAX(h_f(edg),h_f_min)
                    adv_v  =  adv_v + discharge*(velocity0(edg_ngb)%v  &
                    -velocity0(edg)%v) /vol /MAX(h_f(edg),h_f_min)
                ENDIF	
            ENDDO

        END DO   ! elin  

    CASE(2)
        adv_u = 0.0d0
        adv_v = 0.0d0
        vol   = Vel_Area(edg)
        DO elin = 1,2

            el      = edge(edg)%elementinds(elin)
            IF (el == 0) CYCLE	
            elnodes = element(el)%nodeinds(:)
            IF( h_wd(edg) <= h_min )  CYCLE
            IF( h_f(edg)  <= h_min )  CYCLE
            DO j= elin*2-1, elin*2              
                discharge = 0.0d0;  
                edg_ngb =  edge_ngb_edges(edg,j)

                IF( h_wd(edg_ngb) <= h_min )  CYCLE
                IF( h_f(edg_ngb)  <= h_min )  CYCLE 

                SELECT CASE(j)
                CASE(1)
                    x_mid = 0.5*( node(edge(edg)%nodeinds(2))%x +      &
                    & element(edge(edg)%elementinds(1))%x)
                    y_mid = 0.5*( node(edge(edg)%nodeinds(2))%y +      &
                    & element(edge(edg)%elementinds(1))%y)
                CASE(2)
                    x_mid = 0.5*( node(edge(edg)%nodeinds(1))%x +      &
                    & element(edge(edg)%elementinds(1))%x)
                    y_mid = 0.5*( node(edge(edg)%nodeinds(1))%y +      &
                    & element(edge(edg)%elementinds(1))%y)
                CASE(3)
                    x_mid = 0.5*( node(edge(edg)%nodeinds(1))%x +      &
                    & element(edge(edg)%elementinds(2))%x)
                    y_mid = 0.5*( node(edge(edg)%nodeinds(1))%y +      &
                    & element(edge(edg)%elementinds(2))%y)
                CASE(4)
                    x_mid = 0.5*( node(edge(edg)%nodeinds(2))%x +      &
                    & element(edge(edg)%elementinds(2))%x)
                    y_mid = 0.5*( node(edge(edg)%nodeinds(2))%y +      &
                    & element(edge(edg)%elementinds(2))%y)
                END SELECT

                CALL H2Ocean_Vel_Interp(edg,x_mid,y_mid,u_mid1,v_mid1)
                CALL H2Ocean_Vel_Interp(edg_ngb,x_mid,y_mid,u_mid2,v_mid2) 

                u_ave = 0.5d0*(u_mid2+u_mid1)
                v_ave = 0.5d0*(v_mid2+v_mid1) 

                CALL Cross_Product(u_ave,v_ave,edge_vectors(edg,j+2)%dx,  &
                edge_vectors(edg,j+2)%dy,edge_vectors(edg,j+2)%length,&
                discharge)

                IF(discharge > 0.0d0) THEN            
                    discharge = discharge*h_f(edg_ngb)         
                    adv_u  =  adv_u + discharge*(velocity0(edg_ngb)%u  &
                    -velocity0(edg)%u) /vol  /MAX(h_f(edg),h_f_min)
                    adv_v  =  adv_v + discharge*(velocity0(edg_ngb)%v  &
                    -velocity0(edg)%v) /vol  /MAX(h_f(edg),h_f_min)
                ENDIF	
            ENDDO
        END DO   ! elin        

    CASE DEFAULT

    END SELECT
    END SUBROUTINE H2Ocean_Compute_Advection
    
    
    
    
    SUBROUTINE H2Ocean_Compute_Coriolis(edg,coriolis_u,coriolis_v)
    USE Mod_Precision
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_All_Variable
    USE Mod_Parameter
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND), INTENT(IN)   :: edg
    REAL(KIND=REAL_DP),INTENT(OUT)       :: coriolis_u,coriolis_v

    INTEGER(KIND=INT_KIND)               :: elin,el 
    REAL(KIND=REAL_DP)                   :: vol 
    REAL(KIND=REAL_DP)                   :: coriolis_u0,coriolis_v0,coriolis_u00,coriolis_v00
    REAL(KIND=REAL_DP)                   :: coriolis_param

    coriolis_u = 0.0d0
    coriolis_v = 0.0d0

    coriolis_u00 = 0.0d0
    coriolis_v00 = 0.0d0
    coriolis_u0  = 0.0d0
    coriolis_v0  = 0.0d0

    IF( Coordinate_type == 1) THEN
        coriolis_param =  2.0d0*omega*sintheta_coriolis;
    ELSEIF( Coordinate_type == 2) THEN
        coriolis_param = 2.0d0*omega*SIN(edge(edg)%y)
    ENDIF
    SELECT CASE(coriolis_method)

    CASE(1)
        DO elin = 1,2
            el = edge(edg)%elementinds(elin)
            IF (el == 0) CYCLE
            vol = element(el)%area/3.0d0
            coriolis_u0 =  coriolis_u0  + coriolis_param*velocity0(edg)%v*vol 
            coriolis_v0 =  coriolis_v0  - coriolis_param*velocity0(edg)%u*vol 
            coriolis_u00 = coriolis_u00 + coriolis_param*velocity00(edg)%v*vol 
            coriolis_v00 = coriolis_v00 - coriolis_param*velocity00(edg)%u*vol 
        END DO   ! elin   

        coriolis_u = 1.0/2.0*(3.0*coriolis_u0 - coriolis_u00)/Vel_Area(edg)
        coriolis_v = 1.0/2.0*(3.0*coriolis_v0 - coriolis_v00)/Vel_Area(edg)

    CASE DEFAULT

    END SELECT


    END SUBROUTINE H2Ocean_Compute_Coriolis


    SUBROUTINE H2Ocean_Compute_viscosity(edg,vis_u,vis_v)
    USE Mod_Precision
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_All_Variable
    USE Mod_Parameter
    USE Mod_FEM
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND), INTENT(IN)   :: edg
    REAL(KIND=REAL_DP),INTENT(OUT)       :: vis_u,vis_v

    INTEGER(KIND=INT_KIND)               :: elin,el,i 
    REAL(KIND=REAL_DP)                   :: vol,dux,dvx,duy,dvy,A_vel
    REAL(KIND=REAL_DP)                   :: vh
    vis_u = 0.0d0
    vis_v = 0.0d0
    SELECT CASE(viscosity_scheme)
    CASE(1)
    DO elin = 1,2
        el  = edge(edg)%elementinds(elin)
        IF (el == 0) CYCLE
        IF( h_wd(edg) <= h_min ) CYCLE	
        i = MINLOC(ABS(element(el)%edgeinds(:) - edg), 1)
        vol = element(el)%area  
        dux = SUM(bafuvel_x(:, el)*velocity1(element(el)%edgeinds(:))%u)
        duy = SUM(bafuvel_y(:, el)*velocity1(element(el)%edgeinds(:))%u)
        dvx = SUM(bafuvel_x(:, el)*velocity1(element(el)%edgeinds(:))%v)
        dvy = SUM(bafuvel_y(:, el)*velocity1(element(el)%edgeinds(:))%v)
        A_vel = 0.3*vol*SQRT(dux**2. + 0.5*(dux**2. + duy**2.) +  dvy**2.)
        vis_u = vis_u - A_vel*(bafuvel_x(i, el)*dux + bafuvel_y(i, el)*duy)*vol/3.
        vis_v = vis_v - A_vel*(bafuvel_x(i, el)*dvx + bafuvel_y(i, el)*dvy)*vol/3.

    ENDDO
    vis_u = vis_u /Vel_Area(edg)
    vis_v = vis_v /Vel_Area(edg)
    CASE(2)
        dux = grad_vel_dx(edg,1)
        duy = grad_vel_dy(edg,1)
        dvx = grad_vel_dx(edg,2)
        dvy = grad_vel_dy(edg,2)
        vh  = 0.2d0*Vel_Area(edg)*SQRT(dux**2 + dvy**2 + 0.5*(duy**2+dvx**2))
        vis_u = vis_u + vh*vel_2nd_grad(edg,1) + vh*vel_2nd_grad(edg,2)  
        vis_v = vis_v + vh*vel_2nd_grad(edg,3) + vh*vel_2nd_grad(edg,4)  
    CASE DEFAULT
        
    END SELECT
    
    END SUBROUTINE H2Ocean_Compute_viscosity



    SUBROUTINE H2Ocean_Flux_Limiter(nd,el,edg,discharge,method,h_flux)
    USE Mod_Precision
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_TriMesh_Interp
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND), INTENT(IN):: nd,el,edg,method
    REAL(KIND=REAL_DP)    , INTENT(IN) :: discharge
    REAL(KIND=REAL_DP)    , INTENT(OUT):: h_flux  

    INTEGER(KIND=INT_KIND)             :: nd_c,nd_d
    INTEGER(KIND=INT_KIND)             :: other_nd
    REAL(KIND=REAL_DP)                  :: dep_flux,zeta_flux

    IF(edge(edg)%nodeinds(1) == nd) THEN
        other_nd = edge(edg)%nodeinds(2)
    ELSE
        other_nd = edge(edg)%nodeinds(1)
    ENDIF


    IF(discharge >=0.0d0) THEN
        nd_c = other_nd
        nd_d = nd
    ELSE
        nd_c = nd
        nd_d = other_nd 
    ENDIF


    !  dep_flux = 0.5*(Still_Depth(nd)+Still_Depth(other_nd)) 
    SELECT CASE(method)

    CASE(1)
        dep_flux = MIN(Still_Depth(nd),Still_Depth(other_nd))
        h_flux   = MAX(dep_flux + zeta0(nd_c),0.0d0) 
    CASE(2)
        dep_flux  = Still_Depth(nd_c)
        zeta_flux = zeta0(nd_c)     
        h_flux    = MAX(dep_flux + zeta_flux,0.0d0)
    CASE(3)
        dep_flux  = 0.5*(Still_Depth(nd)+Still_Depth(other_nd))    
        h_flux    = MAX(dep_flux + zeta0(nd_c),0.0d0) 
    CASE DEFAULT
        WRITE(*,*) 'A flux limiter method need to be specified!'   
    END SELECT   
    RETURN
    END SUBROUTINE H2Ocean_Flux_Limiter


    SUBROUTINE H2Ocean_Vel_Interp(edg,x_mid,y_mid,u_mid,v_mid)
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_Parameter    
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND), INTENT(IN)  :: edg
    REAL(KIND=REAL_DP),     INTENT(IN)  :: x_mid,y_mid
    REAL(KIND=REAL_DP),     INTENT(OUT) :: u_mid,v_mid
    REAL(KIND=REAL_DP)  :: xc, yc, dx,dy, y_ave

    xc = edge(edg)%x
    yc = edge(edg)%y
    dx = x_mid - xc
    dy = y_mid - yc
    IF(Coordinate_type == 2) THEN
        y_ave = 0.5*(y_mid+yc)
        dx    = dx*cos(y_ave)
        dx    = dx * r_earth
        dy    = dy * r_earth
    ENDIF

    u_mid = velocity0(edg)%u + grad_vel_dx(edg,1)*dx + grad_vel_dy(edg,1)*dy
    v_mid = velocity0(edg)%v + grad_vel_dx(edg,2)*dx + grad_vel_dy(edg,2)*dy
    END SUBROUTINE H2Ocean_Vel_Interp
    !
    SUBROUTINE Cross_Product(u,v,dx,dy,length,sum)
    USE Mod_Precision
    IMPLICIT NONE 
    REAL(KIND=REAL_DP), INTENT(in)      :: u,v,dx,dy,length
    REAL(KIND=REAL_DP), INTENT(out)     :: sum

    sum = 0.0d0;
    sum = (dx*v - dy*u)*length
    RETURN
    END SUBROUTINE Cross_Product
    !
    !
    !
    SUBROUTINE H2Ocean_Update 
    USE Mod_H2Ocean
    IMPLICIT NONE

    ! Update all the variables
    zeta0      = zeta1
    velocity00 = velocity0
    velocity0  = velocity1
    END SUBROUTINE H2Ocean_Update




    SUBROUTINE H2Ocean_Compute_Radiation
    USE Mod_TriMesh
    USE Mod_H2Ocean  
    USE Mod_Parameter
    USE Mod_Dirichlet_Boundary
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  ::i,edg, nd1,nd2
    REAL(KIND=REAL_DP)           :: hu
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  
    radiation = 0.0d0
    !$OMP   PARALLEL PRIVATE(i,edg,hu,nd1,nd2,mythread)  &
    !$OMP&         SHARED(ob_edge_3,Still_Depth,Zeta1,edge,radiation,  &
    !$OMP&         CHUNK,nthreads)

    !$ NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(ob_edge_3%nmb/NTHREADS)+1  

    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO i=1,ob_edge_3%nmb
        edg = ob_edge_3%addresses(i)
        hu=max(0.5d0*SUM(Still_Depth(edge(edg)%nodeinds))  &
        + 0.5d0*SUM(Zeta1(edge(edg)%nodeinds)),0.0d0)
        nd1 = edge(edg)%nodeinds(1)
        nd2 = edge(edg)%nodeinds(2)
        radiation(nd1) = radiation(nd1) +   0.5d0*SQRT(g*hu)*edge(edg)%length
        radiation(nd2) = radiation(nd2) +   0.5d0*SQRT(g*hu)*edge(edg)%length
    ENDDO



    !$OMP END DO
    !$OMP END PARALLEL

    IF(.NOT. dirichlet_boundary .OR. D_bo>=D_bo_Inf) THEN
        DO i=1,ob_edge_2%nmb
            edg = ob_edge_2%addresses(i)
            hu=max(0.5d0*SUM(Still_Depth(edge(edg)%nodeinds))  &
            + 0.5d0*SUM(Zeta1(edge(edg)%nodeinds)),0.0d0)
            nd1 = edge(edg)%nodeinds(1)
            nd2 = edge(edg)%nodeinds(2)
            radiation(nd1) = radiation(nd1) +   0.5d0*SQRT(g*hu)*edge(edg)%length
            radiation(nd2) = radiation(nd2) +   0.5d0*SQRT(g*hu)*edge(edg)%length
        ENDDO  
    ENDIF
    END SUBROUTINE H2Ocean_Compute_Radiation
