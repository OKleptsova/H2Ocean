    MODULE Mod_H2Ocean
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
    USE Mod_All_Variable
    USE Mod_Derived_DataType
    IMPLICIT NONE


    TYPE velocity_type
        REAL(KIND=REAL_DP)            :: u       ! u component
        REAL(KIND=REAL_DP)            :: v       ! v component
    END TYPE velocity_type

    !Zeta_Area --- the control area of water level zeta
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)      ::  Zeta_Area  

    !Zeta_Area --- the control area of water level zeta
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)      ::  Zeta_Wet_Area  



    !Vel_Area --- the control area of the velocities
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)      ::  Vel_Area

    ! Water levels, zeta1, at time n, zeta 2, at time n+1
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)        ::  zeta0, zeta1 

    ! Velocity, 0, at time n, 1, at time n+1
    TYPE(velocity_type),ALLOCATABLE, DIMENSION(:)        :: velocity0, velocity1   
    ! For Coriolis, at time n-1
    TYPE(velocity_type),ALLOCATABLE, DIMENSION(:)        :: velocity00  


    ! for reduced two layer method, the difference between velocities in two layers
    TYPE(velocity_type),ALLOCATABLE, DIMENSION(:)        :: velocity_delta 

    ! Water level gradient within each element
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)        :: grad_zeta_dx, grad_zeta_dy 

    ! the still water depth
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)        ::  Still_Depth    

    !! change in the bottom depth 
    !REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)        ::  delta_b    

    !h_wd   -- the water depth at velocity points, for wetting and drying 
    !h_f    -- the water depth at face, coefficient for water level gradient
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)         ::  h_wd
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)         ::  h_f
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)         ::  h_upwd  


    ! the four neighbour edges of each edge
    INTEGER(KIND=INT_KIND), ALLOCATABLE, DIMENSION(:,:) :: edge_ngb_edges
    ! the coefficients of least square method 
    REAL(kind=REAL_DP),     ALLOCATABLE, DIMENSION(:,:) :: vel_coef_dx,vel_coef_dy
    REAL(kind=REAL_DP),     ALLOCATABLE, DIMENSION(:,:) :: grad_vel_dx,grad_vel_dy
    REAL(kind=REAL_DP),     ALLOCATABLE, DIMENSION(:,:) :: vel_2nd_grad


    INTEGER(KIND=INT_KIND), ALLOCATABLE, DIMENSION(:,:) :: edge_upstream_el

    ! For advection
    !edge_vectors, store all the two vectors connected the edge
    TYPE(vector_type), ALLOCATABLE, DIMENSION(:,:)  ::  edge_vectors

    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)   :: CFL_flux
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)   :: zeta_max, zeta_min
    REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:)   :: radiation


    ! q-- the nonhydrostatic pressure 
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)   :: q 
    ! the vertical verlocity at the water surface
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)   :: w_surf0, w_surf1
    ! the vertical verlocity at the bottom
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)   :: w_bot0,  w_bot1

    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)   :: grad_d_dx,  grad_d_dy

    TYPE(address_type), ALLOCATABLE, DIMENSION(:)     :: node_ngb_nodes_extended  


    ! the coefficient of reduced two layers method
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)   :: alpha_h, alpha_edge
    REAL(KIND = REAL_DP), ALLOCATABLE, DIMENSION(:)   :: k_wave

    CONTAINS

    SUBROUTINE H2Ocean_ALLOCATE_ARRAYS(nodenum,edgenum,elementnum)
    IMPLICIT NONE
    INTEGER(KIND=INT_KIND),INTENT(IN) :: nodenum,edgenum,elementnum
    ALLOCATE(Zeta_Area(nodenum),Vel_Area(edgenum))
    ALLOCATE(zeta0(nodenum),zeta1(nodenum))
    ALLOCATE(velocity00(edgenum),velocity0(edgenum),velocity1(edgenum)) 
    ALLOCATE(Still_Depth(nodenum))

    ALLOCATE(Zeta_Wet_Area(nodenum))

    Zeta_Wet_Area = 0.0d0
    Zeta_Area   = 0.0d0
    Vel_Area    = 0.0d0
    zeta0       = 0.0d0
    zeta1       = 0.0d0 
    velocity00%u= 0.0d0
    velocity00%v= 0.0d0
    velocity0%u = 0.0d0
    velocity0%v = 0.0d0
    velocity1%u = 0.0d0
    velocity1%v = 0.0d0   
    Still_Depth = 0.0d0


    ALLOCATE(zeta_max(nodenum))
    zeta_max    = -999.d0

    ALLOCATE(zeta_min(nodenum))
    zeta_min    = 999.d0

    ALLOCATE(grad_zeta_dx(elementnum),grad_zeta_dy(elementnum))
    grad_zeta_dx = 0.0d0
    grad_zeta_dy = 0.0d0


    ALLOCATE(edge_ngb_edges(edgenum,4))
    edge_ngb_edges = 0

    ALLOCATE(vel_coef_dx(edgenum,5))
    ALLOCATE(vel_coef_dy(edgenum,5))
    vel_coef_dx = 0.0d0
    vel_coef_dy = 0.0d0

    ALLOCATE(grad_vel_dx(edgenum,2))
    grad_vel_dx = 0.0d0

    ALLOCATE(grad_vel_dy(edgenum,2))
    grad_vel_dy = 0.0d0    


    ALLOCATE(vel_2nd_grad(edgenum,4))
    vel_2nd_grad = 0.0d0    
    ! vel_2nd_grad(;,1)  du_dxx
    ! vel_2nd_grad(;,2)  du_dyy
    ! vel_2nd_grad(;,3)  dv_dxx 
    ! vel_2nd_grad(;,4)  dv_dyy


    ALLOCATE(edge_upstream_el(edgenum,2))


    ALLOCATE(h_wd(edgenum),h_f(edgenum))
    ALLOCATE(h_upwd(edgenum))


    ALLOCATE(edge_vectors(edgenum,6))

    IF(varied_time_step) THEN
        ALLOCATE(CFL_flux(nodenum))
    ENDIF

    ALLOCATE(radiation(nodenum))
    radiation = 0.0d0


    IF(enable_reduced_twolayer) THEN
        ALLOCATE(velocity_delta(edgenum))
        velocity_delta%u = 0.0d0
        velocity_delta%v = 0.0d0

        ALLOCATE(alpha_h(nodenum))
        alpha_h = 0.5d0

        ALLOCATE(alpha_edge(edgenum))
        alpha_edge = 0.5d0
        ALLOCATE(k_wave(nodenum))
        k_wave = 0.0d0
    ENDIF
    END SUBROUTINE H2Ocean_ALLOCATE_ARRAYS


    END MODULE Mod_H2Ocean
