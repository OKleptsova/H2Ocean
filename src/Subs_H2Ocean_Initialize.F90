    SUBROUTINE H2Ocean_Initialize 
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
    USE MOD_TriMesh
    USE Mod_H2Ocean
    USE Mod_All_Variable
    USE Mod_Dirichlet_Boundary
    IMPLICIT NONE

    ! Step1: Allocate arrays
    CALL H2Ocean_ALLOCATE_ARRAYS(TriMesh_nodenmb,TriMesh_edgenmb,TriMesh_elementnmb)

    !Step2:
    CALL H2Ocean_Create_Output_Folder


    !Step 3:
    CALL H2Ocean_Compute_Zeta_Area

    !Step 4:
    CALL H2Ocean_Compute_Vel_Area

    ! Step5: 
    CALL FEM_ShapeFunction_Initialize

    ! Step7:
    CALL H2Ocean_Build_edge_ngb_edges 


    ! Step 8:
    CALL H2Ocean_Compute_Vel_Coef


    ! Step 9:
    CALL H2Ocean_Build_Edge_Vectors







    ! Step 10:
    SELECT CASE(test_case)
    CASE(0) 
        CALL H2Ocean_Read_Initial_Files
    CASE(1)
        CALL H2Ocean_Solitary_Initialize
    CASE DEFAULT
        CALL H2Ocean_Read_Initial_Files
    END SELECT   



    IF(enable_nonhydro) THEN
        CALL H2Ocean_NonHydro_Initialize
        IF(enable_reduced_twolayer) THEN
            CALL H2Ocean_alpha_h_Initialize
        ENDIF
        
    ENDIF

    ! Step11:
    CALL H2Ocean_Output_Initialize


    IF (dirichlet_boundary) THEN
        CALL H2Ocean_Read_Dirichlet_Boundary
    ENDIF

    IF (u_dirichlet_boundary) THEN
        CALL H2Ocean_Read_u_Dirichlet_Boundary
    ENDIF 

    IF (q_dirichlet_boundary) THEN
        CALL H2Ocean_Read_q_Dirichlet_Boundary
    ENDIF     
    END SUBROUTINE H2Ocean_Initialize



    SUBROUTINE H2Ocean_Solitary_Initialize
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_Precision
    USE Mod_Parameter
    IMPLICIT NONE

    ! parameter to compute solitary wave profile
    REAL(kind=REAL_DP)    ::z0, dep0,eps,eps2,eps3,alph, L ,th,s,s2,s3,s4,s5,s6,th2,z,  &
    y2, y4, o0,o1,o2,o3,c,x,y
    INTEGER          :: i, nd                   
    ! This is only for solitary wave test case
    ! 1, Define the water depth
    Still_Depth(:) = 10.

    z0 = 2.0d0;
    dep0    = 10.d0;
    eps   = z0/dep0;
    eps2 = eps**2;
    eps3 = eps**3;
    alph =(1/dep0)* sqrt(3*eps/4)*(1-(5/8)*eps+(71/128)*eps2);

    L = 1/sqrt(3*z0/4/dep0**3)*acosh(sqrt(20.));

    DO i=1, TriMesh_nodenmb
        x=node(i)%x
        y=node(i)%y      
        s = 1./cosh(alph*(x-2.*L));
        th = tanh(alph*(x-2.*L));
        s2  = s**2;
        s3  = s**3;
        s4  = s**4;
        s5  = s**5;
        s6  = s**6;
        th2= th**2;  
        zeta0(i)    = dep0*(eps*s2-(3/4)*eps2*s2*th2+eps3*((5/8)*s2*th2-(101/80)*s4*th2));
    END DO
    zeta1 = zeta0

    DO i=1,TriMesh_edgenmb
        x = edge(i)%x

        s = 1./cosh(alph*(x-2.*L));
        th = tanh(alph*(x-2.*L));
        c = 1+0.5*eps-(3/20)*eps2+(3/56)*eps3;
        c=c*sqrt(g*dep0);
        s2  = s**2;
        s3  = s**3;
        s4  = s**4;
        s5  = s**5;
        s6  = s**6;
        th2= th**2;  
        z    = 1+eps*s2-(3/4)*eps2*s2*th2+eps3*((5/8)*s2*th2-(101/80)*s4*th2);
        y2 = (1./3.)*(z**3)/(z);
        y4 = (1./5.)*(z**5)/(z);   
        o0     = 1;
        o1     = 0.5-s2;
        o2     = -(3/20)-0.25*s2+s4+y2*((3/2)*s2-(9/4)*s4);
        o3     =(3/56)+(19/40)*s2+(1/5)*s4-(6/5)*s6+y2*(-(3/2)*s2-(15/4)*s4+(15/2)*s6)+y4*((-3/8)*s2+(45/16)*s4-(45/16)*s6);
        velocity0(i)%u = o0+eps*o1+eps2*o2+eps3*o3;
        velocity0(i)%u =-velocity0(i)%u*sqrt(g*dep0)+c;
        velocity0(i)%v = 0.0d0

    ENDDO

    velocity1(:)%u = velocity0(:)%u
    velocity1(:)%v = velocity0(:)%v


    IF(enable_nonhydro) THEN
        DO nd=1,TriMesh_nodenmb
            x = node(nd)%x
            s = 1./cosh(alph*(x-2.*L));
            th = tanh(alph*(x-2.*L));
            c = 1+0.5*eps-(3/20)*eps2+(3/56)*eps3;
            c=c*sqrt(g*dep0);
            s2  = s**2;
            s3  = s**3;
            s4  = s**4;
            s5  = s**5;
            s6  = s**6;
            th2= th**2;  
            z    = 1+eps*s2-(3/4)*eps2*s2*th2+eps3*((5/8)*s2*th2-(101/80)*s4*th2);
            y     = z
            y2    = y**2;
            y4    = y2**2;
            o0     = 0;
            o1     = -s2;
            o2     = (3/8)*s2+2*s4+y2*(0.5*s2-(3/2)*s4);
            o3     =(49/640)*s2-(17/20)*s4-(18/5)*s6+y2*(-(13/16)*s2-(25/16)*s4+(15/2)*s6)+y4*(-(3/40)*s2+(9/8)*s4-(27/16)*s6);
            w_surf0(nd) = sqrt(3*eps)*y*th*(o0+eps*o1+eps2*o2+eps3*o3);
            w_surf0(nd) = -w_surf0(nd)*sqrt(g*dep0);
        ENDDO
        w_surf1 = w_surf0 
    ENDIF

    END SUBROUTINE H2Ocean_Solitary_Initialize

    SUBROUTINE H2Ocean_Create_Output_Folder
    USE Mod_All_Variable
    USE Mod_IO
    IMPLICIT NONE

    CALL Create_Folder( output_folder )

    fid_time = get_file_unit();   
    OPEN(fid_time, file=trim(adjustl(output_folder))//'time.txt',status="replace")


    END SUBROUTINE H2Ocean_Create_Output_Folder

    SUBROUTINE FEM_ShapeFunction_Initialize
    USE Mod_FEM
    IMPLICIT NONE
    CALL FEM_Standard_Element_Definition 
    CALL FEM_Basis_Functions 
    END SUBROUTINE FEM_ShapeFunction_Initialize


    SUBROUTINE H2Ocean_Read_Initial_Files
    USE Mod_Control_File
    USE Mod_All_Variable
    USE Mod_H2Ocean
    USE Mod_IO
    USE MOD_TriMesh
    IMPLICIT NONE
    INTEGER         :: RES, STATUS
    LOGICAL         :: FEXIST
    INTEGER         :: file_id  
    INTEGER(KIND=8) :: nq, n
    REAL(KIND=REAL_DP)      :: zeta_in
    CHARACTER(len=400)      :: filename
    REAL(kind=REAL_DP), dimension(TriMesh_nodenmb) :: delta_b


    IF(enable_hot_start) THEN
        hot_start_bathymetry= hot_start_folder//'Hot_Start_Bathymetry.txt'
        ! Step1: Read Bathymetry
        CALL INQUIRE_File( TRIM(hot_start_bathymetry), FEXIST )
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The hot start bathymetryfile file doen not exist, the water depth will be set to zero!.'
            Still_Depth = 0.0d0
        ELSE
            ! READ bathymetry
            file_id=get_file_unit();
            OPEN (file_id, file=TRIM(hot_start_bathymetry), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_nodenmb) THEN
                WRITE(*,*) 'Wrong with the Hot_Start_Bathymetry!'
                STOP
            ENDIF
            DO n=1, TriMesh_nodenmb
                READ(file_id,*) Still_Depth(n)
            ENDDO
            CLOSE(file_id)   
        ENDIF


        ! Step2: Read Initial Water Levels
        hot_start_water_levels= hot_start_folder//'Hot_Start_Water_Levels.txt'
        CALL INQUIRE_File( TRIM(hot_start_water_levels), FEXIST )
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The init_zeta_file file doen not exist, the initial water levels will be set to zero!.'
            zeta0 = 0.0d0
        ELSE
            file_id=get_file_unit();
            ! READ init_zeta_file
            OPEN (file_id, file=TRIM(hot_start_water_levels), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_nodenmb) THEN
                WRITE(*,*) 'WRONG with the Hot_Start_Water_Levels!'
                STOP
            ENDIF
            DO n=1, TriMesh_nodenmb
                READ(file_id,*) zeta0(n)
            ENDDO    
            CLOSE(file_id)   
        ENDIF      

        zeta1 = zeta0

        hot_start_velocity00= hot_start_folder//'Hot_Start_Velocity00.txt'
        ! Step3: Read Initial Velocity
        CALL INQUIRE_File( TRIM(hot_start_velocity00), FEXIST )
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The hot_start_velocity00 file doen not exist, the initial velocity will be set to zero!.'
            velocity00%u = 0.0d0
            velocity00%v = 0.0d0
        ELSE
            file_id=get_file_unit();
            ! READ init_zeta_file
            OPEN (file_id, file=TRIM(hot_start_velocity00), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_edgenmb) THEN
                WRITE(*,*) 'WRONG with the Hot_Start_Velocity00!'
                STOP
            ENDIF
            DO n=1, TriMesh_edgenmb
                READ(file_id,*) velocity00(n)%u,velocity00(n)%v
            ENDDO    
            CLOSE(file_id)   
        ENDIF  

        hot_start_velocity0= hot_start_folder//'Hot_Start_Velocity0.txt'
        ! Step3: Read Initial Velocity
        CALL INQUIRE_File( TRIM(hot_start_velocity0), FEXIST )
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The hot_start_velocity0 file doen not exist, the initial velocity will be set to zero!.'
            velocity00%u = 0.0d0
            velocity00%v = 0.0d0
        ELSE
            file_id=get_file_unit();
            ! READ init_zeta_file
            OPEN (file_id, file=TRIM(hot_start_velocity0), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_edgenmb) THEN
                WRITE(*,*) 'WRONG with the Hot_Start_Velocity0!'
                STOP
            ENDIF
            DO n=1, TriMesh_edgenmb
                READ(file_id,*) velocity0(n)%u,velocity0(n)%v
            ENDDO    
            CLOSE(file_id)   
        ENDIF      




        hot_start_time= hot_start_folder//'Hot_Start_Time.txt'
        ! Step3: Read Initial Velocity
        CALL INQUIRE_File( TRIM(hot_start_time), FEXIST )
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The hot_start_time file doen not exist, the initial time will be set to zero!.'
            ini_time = 0.0d0
        ELSE
            file_id=get_file_unit();
            ! READ init_zeta_file
            OPEN (file_id, file=TRIM(hot_start_time), status='old')   
            READ(file_id,*) ini_time
        ENDIF    
    ELSE
        ! Step1: Read Bathymetry
      
        INQUIRE(file=TRIM(bathymetryfile), exist=FEXIST)
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The bathymetryfile file doen not exist, the water depth will be set to zero!.'
            Still_Depth = 0.0d0
            CLOSE(file_id)  
        ELSE
            ! READ bathymetry
            file_id=get_file_unit();
            OPEN (file_id, file=TRIM(bathymetryfile), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_nodenmb) THEN
                WRITE(*,*) 'wrong with the bathymetryfile!'
                STOP
            ENDIF
            DO n=1, TriMesh_nodenmb
                READ(file_id,*) Still_Depth(n)
            ENDDO
            CLOSE(file_id)   
        ENDIF


        ! Step2: Read Initial Water Levels
        INQUIRE(file=TRIM(init_zeta_file), exist=FEXIST)
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The init_zeta_file file doen not exist, the initial water levels will be set to zero!.'
            zeta0 = 0.0d0
        ELSE
            file_id=get_file_unit();
            ! READ init_zeta_file
            OPEN (file_id, file=TRIM(init_zeta_file), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_nodenmb) THEN
                WRITE(*,*) 'WRONG with the init_zeta_file!'
                STOP
            ENDIF
            DO n=1, TriMesh_nodenmb
                READ(file_id,*) zeta0(n)
            ENDDO    
            CLOSE(file_id)   
        ENDIF      


        IF(enable_bottom_change) THEN
        INQUIRE(file=TRIM(bottom_change_file), exist=FEXIST)
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The bottom_change_file file does not exist, the bottom levels will remain unchanged!.'
            delta_b = 0.0d0
        ELSE
            file_id=get_file_unit();
            ! READ bottom_change_file
            OPEN (file_id, file=TRIM(bottom_change_file), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_nodenmb) THEN
                WRITE(*,*) 'WRONG with the bottom_change_file!'
                STOP
            ENDIF
            DO n=1, TriMesh_nodenmb
                READ(file_id,*) delta_b(n)
                Still_Depth(n) = Still_Depth(n) - delta_b(n)
                zeta0(n)=zeta0(n)+delta_b(n)
            ENDDO    
            CLOSE(file_id)   
        ENDIF      
        ENDIF

        DO n=1, TriMesh_nodenmb
            IF( zeta0(n) + Still_Depth(n) < h_min) THEN
                zeta0(n) = -Still_Depth(n)
            ENDIF
        ENDDO
        zeta1 = zeta0




        ! Step3: Read Initial Velocity
        INQUIRE(file=TRIM(init_vel_file), exist=FEXIST)
        IF(.not.FEXIST) THEN
            WRITE(*,*) 'The init_vel_file file doen not exist, the initial velocity will be set to zero!.'
            velocity0%u = 0.0d0
            velocity0%v = 0.0d0
        ELSE
            file_id=get_file_unit();
            ! READ init_zeta_file
            OPEN (file_id, file=TRIM(init_vel_file), status='old')   
            READ(file_id,*) nq
            IF(nq .ne. TriMesh_edgenmb) THEN
                WRITE(*,*) 'WRONG with the init_vel_file!'
                STOP
            ENDIF
            DO n=1, TriMesh_edgenmb
                READ(file_id,*) velocity0(n)%u,velocity0(n)%v
            ENDDO    
            CLOSE(file_id)   
        ENDIF  

        velocity00 = velocity0
        velocity1  = velocity0        

    ENDIF   !IF(enable_hot_start) THEN       





    END SUBROUTINE H2Ocean_Read_Initial_Files



    SUBROUTINE H2Ocean_Compute_Zeta_Area
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)       :: nd, i, el
    REAL(kind=REAL_DP)           :: vol_sum 
    ! Compute the control area of the water level points,
    ! Here zeta locates at nodes
    ! The control volum of the water level will be 1/3 of all the neighbour elements


    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  



    !$OMP   PARALLEL PRIVATE(nd,i,vol_sum,el,mythread)  &
    !$OMP&           SHARED(Trimesh_nodenmb,node_ngb_elements,element,Zeta_Area,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_nodenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO nd=1, Trimesh_nodenmb
        vol_sum = 0.0d0
        DO i = 1, node_ngb_elements(nd)%nmb
            el         = node_ngb_elements(nd)%addresses(i)
            vol_sum    = vol_sum + element(el)%area/3.0d0
        ENDDO
        Zeta_Area(nd) = vol_sum
    ENDDO
    !$OMP END DO   
    !$OMP END PARALLEL   

    Zeta_Wet_Area = Zeta_Area
    END SUBROUTINE H2Ocean_Compute_Zeta_Area 


    SUBROUTINE H2Ocean_Compute_Vel_Area
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)       :: edg, i, el, elin
    REAL(kind=REAL_DP)                :: vol_sum
    !The control volum of the water level will be 1/3 of all the neighbour elements
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread 

    !$OMP   PARALLEL PRIVATE(edg,vol_sum,elin,el,mythread) &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,element,Vel_Area,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg=1, TriMesh_edgenmb
        vol_sum = 0.0d0
        DO elin = 1,2
            el   = edge(edg)%elementinds(elin)
            IF(el == 0) CYCLE
            vol_sum    = vol_sum + element(el)%area/3.0d0
        ENDDO
        Vel_Area(edg) = vol_sum  
    ENDDO
    !$OMP END DO   
    !$OMP END PARALLEL   
    END SUBROUTINE H2Ocean_Compute_Vel_Area 





    SUBROUTINE H2Ocean_Build_edge_ngb_edges
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_Precision
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND) :: edg, el1, el2,nd1, nd2, j, edg_j, eledges(3)
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    !$OMP   PARALLEL PRIVATE(edg,el1,el2,nd1,nd2,eledges,j,edg_j,mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,element,edge_ngb_edges,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg = 1, TriMesh_edgenmb
        el1 = edge(edg)%elementinds(1)
        el2 = edge(edg)%elementinds(2)
        nd1 = edge(edg)%nodeinds(1)
        nd2 = edge(edg)%nodeinds(2)
        eledges = element(el1)%edgeinds(:)

        DO j= 1,3
            edg_j = eledges(j)
            IF(ANY(edge(edg_j)%nodeinds(:)==nd2)  .and. edg_j .ne. edg) THEN
                edge_ngb_edges(edg,1) = edg_j
            ENDIF
            IF(ANY(edge(edg_j)%nodeinds(:)==nd1)  .and. edg_j .ne. edg) THEN
                edge_ngb_edges(edg,2) = edg_j
            ENDIF
        ENDDO
        IF(el2 .ne. 0) THEN
            eledges = element(el2)%edgeinds(:)
            DO j= 1,3
                edg_j = eledges(j)
                IF(ANY(edge(edg_j)%nodeinds(:)==nd1)  .and. edg_j .ne. edg) THEN
                    edge_ngb_edges(edg,3) = edg_j
                ENDIF
                IF(ANY(edge(edg_j)%nodeinds(:)==nd2)  .and. edg_j .ne. edg) THEN
                    edge_ngb_edges(edg,4) = edg_j
                ENDIF
            ENDDO
        ENDIF	
    ENDDO
    !$OMP END DO   
    !$OMP END PARALLEL     
    END SUBROUTINE H2Ocean_Build_edge_ngb_edges






    SUBROUTINE H2Ocean_Compute_Vel_Coef
    USE Mod_TriMesh
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_Parameter
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)   :: edg, j, edg_ngb
    REAL(KIND=REAL_DP)       :: A4(4,2),A2(2,2),A(2,2)
    REAL(KIND=REAL_DP)       :: x2(2),y2(2), x4(4),y4(4)
    REAL(KIND=REAL_DP)       :: xc,yc,x1,y1, y_ave
    REAL(KIND=REAL_DP)       :: delt

    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    !$OMP   PARALLEL PRIVATE(edg,xc,yc,x2,y2,x4,y4,j,x1,y1,edg_ngb, &
    !$OMP&                   A,y_ave,delt,mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,edge_ngb_edges, vel_coef_dx,vel_coef_dy,&
    !$OMP&                  Coordinate_type,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    DO edg = 1, TriMesh_edgenmb
        xc = edge(edg)%x
        yc = edge(edg)%y
        x2 = 0.0d0
        y2 = 0.0d0
        x4 = 0.0d0
        y4 = 0.0d0

        IF(MINVAL(edge(edg)%elementinds(:)) == 0) THEN
            ! With only one neighbour element, two neighbour edges
            DO j=1,2
                edg_ngb = edge_ngb_edges(edg,j)
                x1 = edge(edg_ngb)%x
                y1 = edge(edg_ngb)%y
                x2(j) = x1 -xc
                y2(j) = y1 -yc
                IF(Coordinate_type == 2) THEN
                    y_ave = 0.5d0*(y1+yc)
                    x4(j) = x4(j)*cos(y_ave)
                    x2(j) = r_earth *x2(j) 
                    y2(j) = r_earth *y2(j) 
                ENDIF
            ENDDO
            A = 0.0d0

            DO j=1,2
                A(1,1) = A(1,1) + x2(j)*x2(j)
                A(1,2) = A(1,2) + x2(j)*y2(j)
                A(2,1) = A(2,1) + x2(j)*y2(j)
                A(2,2) = A(2,2) + y2(j)*y2(j)
            ENDDO 

            delt = A(1,1)*A(2,2)-A(1,2)*A(2,1)

            vel_coef_dx(edg,:) = 0.0d0
            vel_coef_dy(edg,:) = 0.0d0
            DO j=1,2
                vel_coef_dx(edg,j) =  A(2,2)/delt*x2(j)  - A(1,2)/delt*y2(j)
                vel_coef_dy(edg,j) = -A(2,1)/delt*x2(j)  + A(1,1)/delt*y2(j)
            ENDDO
            DO j=1,2
                vel_coef_dx(edg,5) = vel_coef_dx(edg,5)-vel_coef_dx(edg,j) 
                vel_coef_dy(edg,5) = vel_coef_dy(edg,5)-vel_coef_dy(edg,j)
            ENDDO          

        ELSE
            DO j=1,4
                edg_ngb = edge_ngb_edges(edg,j)
                x1 = edge(edg_ngb)%x
                y1 = edge(edg_ngb)%y
                x4(j) = x1 - xc
                y4(j) = y1 - yc
                IF(Coordinate_type == 2) THEN
                    y_ave = 0.5d0*(y1+yc)
                    x4(j) = x4(j)*cos(y_ave)
                    x4(j) = r_earth *x4(j) 
                    y4(j) = r_earth *y4(j) 
                ENDIF
            ENDDO
            A = 0.0d0

            DO j=1,4
                A(1,1) = A(1,1) + x4(j)*x4(j)
                A(1,2) = A(1,2) + x4(j)*y4(j)
                A(2,1) = A(2,1) + x4(j)*y4(j)
                A(2,2) = A(2,2) + y4(j)*y4(j)

            ENDDO 

            delt = A(1,1)*A(2,2)-A(1,2)*A(2,1)


            vel_coef_dx(edg,:) = 0.0d0
            vel_coef_dy(edg,:) = 0.0d0
            DO j=1,4
                vel_coef_dx(edg,j) =  A(2,2)/delt*x4(j) - A(1,2)/delt*y4(j)
                vel_coef_dy(edg,j) = -A(2,1)/delt*x4(j) + A(1,1)/delt*y4(j)
            ENDDO


            DO j=1,4
                vel_coef_dx(edg,5) = vel_coef_dx(edg,5)-vel_coef_dx(edg,j) 
                vel_coef_dy(edg,5) = vel_coef_dy(edg,5)-vel_coef_dy(edg,j)
            ENDDO         
        ENDIF

    ENDDO  !DO edg = 1, TriMesh_edgenmb
    !$OMP END DO   
    !$OMP END PARALLEL 
    END SUBROUTINE H2Ocean_Compute_Vel_Coef 





    SUBROUTINE H2Ocean_Build_edg_upstream_el
    USE Mod_Precision
    USE Mod_H2Ocean
    USE Mod_TriMesh
    USE Mod_Control_File
    USE Mod_All_Variable
    USE Mod_IO
    IMPLICIT NONE

    REAL(KIND=REAL_DP)                       :: distance(50)
    INTEGER(KIND=INT_KIND)                   :: n
    INTEGER(KIND=INT_KIND)                   :: min_pos(1) 
    INTEGER(KIND=INT_KIND)                   :: edg,el, nd1,nd2
    REAL(KIND=REAL_DP)                       :: x1,y1,x2,y2, x0, y0
    INTEGER         :: RES, STATUS
    LOGICAL         :: FEXIST
    INTEGER         :: file_id  
    INTEGER(KIND=INT_KIND)   ::nq  
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    keyword = 'edg_upstream_el_file'
    CALL Read_control_file(ctrl_file,keyword,edg_upstream_el_file)



    CALL INQUIRE_File(TRIM(edg_upstream_el_file),FEXIST)

    IF(.not. FEXIST) THEN

        !$OMP   PARALLEL PRIVATE(edg,nd1,nd2,x1,y1,x2,y2,x0,y0,n,el, &
        !$OMP&                   min_pos,distance,mythread)  &
        !$OMP&           SHARED(TriMesh_edgenmb,edge,node,element,node_ngb_elements,&
        !$OMP&                  edge_upstream_el, &
        !$OMP&                  CHUNK,nthreads)  


        !$  NTHREADS = OMP_GET_NUM_THREADS()
        !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

        DO edg=1, TriMesh_edgenmb
            nd1 = edge(edg)%nodeinds(1) 
            nd2 = edge(edg)%nodeinds(2)
            x1  = node(nd1)%x
            y1  = node(nd1)%y
            x2  = node(nd2)%x
            y2  = node(nd2)%y

            x0 = 2.0d0*x1 - x2
            y0 = 2.0d0*y1 - y2


            distance = 9999.

            DO n =1, node_ngb_elements(nd1)%nmb
                el  = node_ngb_elements(nd1)%addresses(n)
                distance(n) = SQRT((x0-element(el)%x)**2+(y0-element(el)%y)**2) 
            ENDDO      

            min_pos = MINLOC(distance) 
            edge_upstream_el(edg,1) = node_ngb_elements(nd1)%addresses(min_pos(1))



            x0 = 2.0d0*x2 - x1
            y0 = 2.0d0*y2 - y1
            distance = 9999.
            DO n =1, node_ngb_elements(nd2)%nmb
                el  = node_ngb_elements(nd2)%addresses(n)
                distance(n) = SQRT((x0-element(el)%x)**2+(y0-element(el)%y)**2) 
            ENDDO    

            min_pos = MINLOC(distance) 
            edge_upstream_el(edg,2) = node_ngb_elements(nd2)%addresses(min_pos(1))
        ENDDO
        !$OMP END DO   
        !$OMP END PARALLEL    



        file_id=get_file_unit()
        WRITE(*,*) TRIM(edg_upstream_el_file)//' does not exist, now creat it!'
        OPEN (file_id, file=TRIM(edg_upstream_el_file), STATUS='Replace') 
        WRITE(file_id,*) TriMesh_edgenmb
        DO edg=1, TriMesh_edgenmb 
            WRITE(file_id,*) edge_upstream_el(edg,1),edge_upstream_el(edg,2)
        ENDDO    
        CLOSE(file_id)   
    ELSE
        file_id=get_file_unit();
        OPEN (file_id, file=TRIM(edg_upstream_el_file), status='old')   
        READ(file_id,*) nq
        IF(nq .ne. TriMesh_edgenmb) THEN
            write(*,*) 'wrong with the edg_upstream_el_file'
        ENDIF
        DO n=1, TriMesh_edgenmb
            READ(file_id,*) edge_upstream_el(n,1),edge_upstream_el(n,2)
        ENDDO
        CLOSE(file_id)   
    ENDIF     

    END SUBROUTINE H2Ocean_Build_edg_upstream_el





    SUBROUTINE H2Ocean_Build_Edge_Vectors
    USE Mod_Precision
    USE Mod_TriMesh
    USE Mod_H2Ocean
    USE Mod_Parameter
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)  :: edg, nd1, nd2, el1, el2, eledges(3), elnodes(3)
    INTEGER(KIND=INT_KIND)  :: i, j, nd_el, nd
    REAL(KIND=REAL_DP)      :: cosin,x1,x2,y1,y2
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread 


    ! 1-2, are the left and right sections, 
    ! which connect the centroid to the middle of the edge
    ! The 3-6 vectors, are the four neighbour vectors of the control volume of edge.
    ! used for the advection computation
    !  ALLOCATE(edge_vectors(edge2D,10))



    !$OMP   PARALLEL PRIVATE(edg,el1,el2,nd1,nd2, j,cosin,i,nd,x1,y1,x2,y2,&
    !$OMP&                   nd_el,mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,node,element,edge,edge_vectors, Coordinate_type,&
    !$OMP&                  CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

    DO edg= 1, TriMesh_edgenmb
        el1 = edge(edg)%elementinds(1)
        el2 = edge(edg)%elementinds(2)
        nd1 = edge(edg)%nodeinds(1)
        nd2 = edge(edg)%nodeinds(2)

        j = 1
        x1 = element(el1)%x
        y1 = element(el1)%y
        x2 = edge(edg)%x
        y2 = edge(edg)%y
        edge_vectors(edg,j)%dx = x2 - x1
        edge_vectors(edg,j)%dy = y2 - y1
        IF(Coordinate_type == 1) THEN
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

        ELSEIF(Coordinate_type == 2) THEN
            cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
        ENDIF


        j = 2
        x1 = 0.0d0
        y1 = 0.0d0
        x2 = 0.0d0
        y2 = 0.0d0
        edge_vectors(edg,j)%dx = 0.0d0
        edge_vectors(edg,j)%dy = 0.0d0
        edge_vectors(edg,j)%length = 0.0d0

        IF( el2 .NE. 0) THEN
            x2 = element(el2)%x
            y2 = element(el2)%y
            x1 = edge(edg)%x
            y1 = edge(edg)%y
            edge_vectors(edg,j)%dx = x2 - x1
            edge_vectors(edg,j)%dy = y2 - y1
            IF(Coordinate_type == 1) THEN
                edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
                + edge_vectors(edg,j)%dy**2) 
                edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
                edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

            ELSEIF(Coordinate_type == 2) THEN
            cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
            ENDIF
        ENDIF


        j = 3


        DO i= 1, 3
            nd= element(el1)%nodeinds(i)
            IF(COUNT(edge(edg)%nodeinds(:)==nd)==0) THEN
                nd_el = nd
                EXIT
            ENDIF
        ENDDO


        SELECT CASE(advection_scheme)
        CASE(1) 
            x1 = node(nd2)%x
            y1 = node(nd2)%y
            x2 = node(nd_el)%x
            y2 = node(nd_el)%y    
        CASE(2) 
            x1 = node(nd2)%x
            y1 = node(nd2)%y
            x2 = element(el1)%x
            y2 = element(el1)%y
        CASE DEFAULT

        END SELECT


        edge_vectors(edg,j)%dx = x2 - x1
        edge_vectors(edg,j)%dy = y2 - y1
        IF(Coordinate_type == 1) THEN
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

        ELSEIF(Coordinate_type == 2) THEN
            cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
        ENDIF





        j = 4

        SELECT CASE(advection_scheme)
        CASE(1) 
            x1 = node(nd_el)%x
            y1 = node(nd_el)%y
            x2 = node(nd1)%x
            y2 = node(nd1)%y    
        CASE(2) 
            x1 = element(el1)%x
            y1 = element(el1)%y
            x2 = node(nd1)%x
            y2 = node(nd1)%y
        CASE DEFAULT
        END SELECT
        edge_vectors(edg,j)%dx = x2 - x1
        edge_vectors(edg,j)%dy = y2 - y1
        IF(Coordinate_type == 1) THEN
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

        ELSEIF(Coordinate_type == 2) THEN
                       cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
        ENDIF



        j = 5
        x1 = 0.0d0
        y1 = 0.0d0
        x2 = 0.0d0
        y2 = 0.0d0
        edge_vectors(edg,j)%dx = 0.0d0
        edge_vectors(edg,j)%dy = 0.0d0
        edge_vectors(edg,j)%length = 0.0d0
        IF( el2 .NE. 0) THEN
            DO i= 1, 3
                nd= element(el2)%nodeinds(i)
                IF(COUNT(edge(edg)%nodeinds(:)==nd)==0) THEN
                    nd_el = nd
                    EXIT
                ENDIF
            ENDDO
        ENDIF    


        IF( el2 .NE. 0) THEN


            SELECT CASE(advection_scheme)
            CASE(1) 
                x2 = node(nd_el)%x
                y2 = node(nd_el)%y   
                x1 = node(nd1)%x
                y1 = node(nd1)%y   
            CASE(2) 
                x2 = element(el2)%x
                y2 = element(el2)%y
                x1 = node(nd1)%x
                y1 = node(nd1)%y
            CASE DEFAULT
            END SELECT



            edge_vectors(edg,j)%dx = x2 - x1
            edge_vectors(edg,j)%dy = y2 - y1
            IF(Coordinate_type == 1) THEN
                edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
                + edge_vectors(edg,j)%dy**2) 
                edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
                edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

            ELSEIF(Coordinate_type == 2) THEN
                      cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
            ENDIF
        ELSE
            x2 = edge(edg)%x
            y2 = edge(edg)%y
            x1 = node(nd1)%x
            y1 = node(nd1)%y
            edge_vectors(edg,j)%dx = x2 - x1
            edge_vectors(edg,j)%dy = y2 - y1
            IF(Coordinate_type == 1) THEN
                edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
                + edge_vectors(edg,j)%dy**2) 
                edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
                edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

            ELSEIF(Coordinate_type == 2) THEN
                        cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
            ENDIF            
        ENDIF   

        j = 6
        x1 = 0.0d0
        y1 = 0.0d0
        x2 = 0.0d0
        y2 = 0.0d0
        edge_vectors(edg,j)%dx = 0.0d0
        edge_vectors(edg,j)%dy = 0.0d0
        edge_vectors(edg,j)%length = 0.0d0

        IF( el2 .NE. 0) THEN


            SELECT CASE(advection_scheme)
            CASE(1) 
                x1 = node(nd_el)%x
                y1 = node(nd_el)%y
                x2 = node(nd2)%x
                y2 = node(nd2)%y  
            CASE(2) 
                x1 = element(el2)%x
                y1 = element(el2)%y
                x2 = node(nd2)%x
                y2 = node(nd2)%y
            CASE DEFAULT
            END SELECT



            edge_vectors(edg,j)%dx = x2 - x1
            edge_vectors(edg,j)%dy = y2 - y1
            IF(Coordinate_type == 1) THEN
                edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
                + edge_vectors(edg,j)%dy**2) 
                edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
                edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

            ELSEIF(Coordinate_type == 2) THEN
                            cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
            ENDIF

        ELSE
            x2 = node(nd2)%x
            y2 = node(nd2)%y
            x1 = edge(edg)%x
            y1 = edge(edg)%y
            edge_vectors(edg,j)%dx = x2 - x1
            edge_vectors(edg,j)%dy = y2 - y1
            IF(Coordinate_type == 1) THEN
                edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
                + edge_vectors(edg,j)%dy**2) 
                edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length 
                edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length

            ELSEIF(Coordinate_type == 2) THEN
                            cosin = cos(0.5*(y1+y2))
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx*cosin*r_earth
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy*r_earth
            edge_vectors(edg,j)%length = SQRT(edge_vectors(edg,j)%dx**2 &
            + edge_vectors(edg,j)%dy**2) 
            edge_vectors(edg,j)%dx = edge_vectors(edg,j)%dx/edge_vectors(edg,j)%length
            edge_vectors(edg,j)%dy = edge_vectors(edg,j)%dy/edge_vectors(edg,j)%length
            ENDIF
        ENDIF            
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL  

    END SUBROUTINE H2Ocean_Build_Edge_Vectors


