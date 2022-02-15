    SUBROUTINE H2Ocean_alpha_h_Initialize
    USE Mod_All_Variable
    USE Mod_H2Ocean
    USE MOD_TriMesh
    IMPLICIT NONE
    IF (k_wave_predefine < 1.0e-5 .AND. T_wave_predefine < 1.0e-5) THEN
        WRITE(*,*) 'You should define one of the variables k_wave_predefine or T_wave_predefine'
    ELSEIF (k_wave_predefine > 1.0e-5) THEN
        k_wave = k_wave_predefine
    ELSEIF (T_wave_predefine > 1.0e-5) THEN   
        CALL H2Ocean_Compute_k_wave(T_wave_predefine,TriMesh_nodenmb,Still_Depth,k_wave)
    ENDIF

    CALL H2Ocean_Compute_alpha_h

    CALL H2Ocean_Compute_alpha_edge

    END SUBROUTINE H2Ocean_alpha_h_Initialize

    SUBROUTINE H2Ocean_Compute_k_wave(T,nodenmb,dep,k)
    USE Mod_Precision
    USE Mod_Parameter
    IMPLICIT NONE
    REAL(KIND=REAL_DP),INTENT(IN)          :: T
    INTEGER(KIND=INT_KIND), INTENT(IN)     :: nodenmb
    REAL(KIND=REAL_DP),INTENT(IN)          :: dep(nodenmb)
    REAL(KIND=REAL_DP),INTENT(INOUT)       :: k(nodenmb)
    REAL(KIND=REAL_DP)                     :: L   ! wave length

    INTEGER(KIND=INT_KIND)      ::  nd

    DO nd=1,nodenmb
        CALL   ldis(T,dep(nd),L) 
        k(nd) = 2.0d0*pi/L
    ENDDO    

    END SUBROUTINE H2Ocean_Compute_k_wave

    SUBROUTINE ldis(T,d,L)
    USE Mod_Precision
    USE Mod_Parameter
    IMPLICIT NONE
    REAL(KIND=REAL_DP),INTENT(IN)             :: T,d
    REAL(KIND=REAL_DP),INTENT(INOUT)          :: L
    REAL(KIND=REAL_DP) :: sigma, y, d1,d2,d3,d4,d5,d6,Tem, kd,k,L1,tol,L0
    ! The first estimation is from Hunt's approximation, 
    ! Hunt, John N. 1979. "Direct Solution of Wave Dispersion Equation" Journal
    ! of the Waterways, Port, Coastal and Ocean Division, No. WW4, ASCE, pp 457-459
    sigma=2*pi/T;
    y=sigma*2*d/g;
    d1=0.6666666666;
    d2=0.3555555555;
    d3=0.1608465608;
    d4=0.0632098765;
    d5=0.0217540484;
    d6=0.0065407983;
    Tem=y/(1+d1*y+d2*y**2+d3*y**3+d4*y**4+d5*y**5*d6*y**6)+y**2;
    kd=SQRT(Tem);
    k=kd/d;
    L1=2*pi/k;
    tol=1.0e-12;
    L0=g*T**2/(2*pi);
    L=L0*TANH(2*pi*d/L1);
    DO WHILE (ABS(L/L1-1)>tol)
        L1=(L+L1)/2;
        L=L0*TANH(2*pi*d/L1);
    ENDDO
    END SUBROUTINE


    SUBROUTINE H2Ocean_Compute_alpha_h
    USE Mod_All_Variable
    USE MOD_TriMesh
    USE Mod_Precision
    USE Mod_H2Ocean
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)      ::  nd
    REAL(KIND=REAL_DP)         ::  kd
    alpha_h = fix_alpha
    
    IF( enable_fixed_alpha .NEQV. .TRUE. ) THEN
        DO nd=1,TriMesh_nodenmb
            kd = k_wave(nd) * Still_Depth(nd)
            CALL H2Ocean_alpha_h(kd,alpha_h(nd))
        ENDDO
    ENDIF
    END SUBROUTINE H2Ocean_Compute_alpha_h

    SUBROUTINE H2Ocean_alpha_h(kd,alpha)
    USE Mod_Precision
    IMPLICIT NONE
    REAL(KIND=REAL_DP),INTENT(IN)  :: kd
    REAL(KIND=REAL_DP),INTENT(OUT) :: alpha
    REAL(KIND=REAL_DP)   :: p1,p2,p3,p4,p5,p6,q1,q2,q3 
    REAL(KIND=REAL_DP)   :: a,b,c



    !p1 =  1.293e-009  ;
    !p2 =  -2.54e-006 ;
    !p3 =       1.001;
    !p4 =      -1.952 ;
    !p5 =      -0.177 ;
    !p6 =     -0.2825;
    !q1 =     0.09518 ;
    !q2 =     -0.6782 ;
    !q3 =        1.57 ;
    !
    !IF (kd<=4) THEN
    !    alpha = 0.5
    !ELSE
    !    alpha =    (p1*kd**5 + p2*kd**4 + p3*kd**3 + p4*kd**2 + p5*kd+ p6) /(kd**3 + q1*kd**2 + q2*kd + q3)
    !ENDIF


    a =      -2.042;
    b =      -1.007 ;
    c =      0.9996  ;
    IF (kd<=4) THEN
        alpha = 0.5
    ELSE
        alpha =    a*kd**b +c;
    ENDIF
    alpha = max(alpha,0.5d0)
    alpha = min(alpha,0.999)
    END SUBROUTINE H2Ocean_alpha_h


    SUBROUTINE H2Ocean_Compute_alpha_edge
    USE Mod_Parameter
    USE Mod_H2Ocean
    USE MOD_TriMesh
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND)     :: edg, el1, el2, elnodes1(3),elnodes2(3)
    REAL(kind=REAL_DP)      :: zeta_el1, zeta_el2
    !Varibles for OpenMP
    !$  integer :: chunk, chunksize  
    !$  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM     
    !$  integer :: mythread  


    !$OMP   PARALLEL PRIVATE(edg,el1,el2,elnodes1,elnodes2,zeta_el1,zeta_el2, &
    !$OMP&                   mythread)  &
    !$OMP&           SHARED(TriMesh_edgenmb,edge,element,alpha_edge,CHUNK,nthreads)  


    !$  NTHREADS = OMP_GET_NUM_THREADS()
    !$  CHUNK = int(TriMesh_edgenmb/NTHREADS)+1  
    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)


    DO edg = 1, TriMesh_edgenmb
        el1 = edge(edg)%elementinds(1)
        el2 = edge(edg)%elementinds(2)
        elnodes1 = element(el1)%nodeinds(:)
        IF(el2 .ne. 0 ) THEN
            elnodes2    =  element(el2)%nodeinds(:)
            zeta_el1    =  SUM(alpha_h(elnodes1))/3.0d0  
            zeta_el2    =  SUM(alpha_h(elnodes2))/3.0d0  
            alpha_edge(edg)  = (element(el1)%area*zeta_el1+  &
            element(el2)%area*zeta_el2)  &
            /(element(el1)%area+element(el2)%area)
        ELSE
            alpha_edge(edg)   =  SUM(alpha_h(elnodes1))/3.0d0  
        ENDIF     
    ENDDO 
    !$OMP END DO   
    !$OMP END PARALLEL  
    END SUBROUTINE H2Ocean_Compute_alpha_edge

