MODULE Mod_All_Variable
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
  IMPLICIT NONE
  
  ! Recorde a processor-dependent approximation of the processor time in seconds
  REAL(kind=REAL_DP)            :: cputime_begin,  cputime_end
  ! Recorde a processor-dependent approximation of the wallclock time in seconds

  !REAL(KIND=REAL_DP)              ::  walltime_begin, walltime_end  
  REAL(KIND=4)              ::  walltime_begin, walltime_end  
  INTEGER(4)                :: computation_method
  
  LOGICAL                  :: enable_reduced_twolayer    = .FALSE.
  
  LOGICAL                  :: enable_nonhydro = .FALSE.
  REAL(kind=REAL_DP)        :: q_alpha = 0.5d0
  INTEGER(KIND=INT_KIND)   :: stencil_method  = 1
  
  INTEGER(KIND=INT_KIND)   :: test_case = 0
  ! Build in test case 
  ! 1 - Solitory wave propagation.
  
  
  CHARACTER(len=400)       :: bathymetryfile =''  
  CHARACTER(len=400)       :: bottom_change_file  =''
  CHARACTER(len=400)       :: init_zeta_file =''
  CHARACTER(len=400)       :: init_vel_file  =''
  LOGICAL                  :: enable_hot_start = .FALSE.
  CHARACTER(len=400)       :: hot_start_folder ='.\'
  
  CHARACTER(len=400)       :: hot_start_bathymetry,hot_start_water_levels
  CHARACTER(len=400)       :: hot_start_velocity00,hot_start_velocity0
  CHARACTER(len=400)       :: hot_start_time
  
  
  REAL(kind=REAL_DP)        :: h_min = 0.0001d0  
  REAL(kind=REAL_DP)        :: dt = 1.0d0,  dt_save = 1.0d0
  REAL(kind=REAL_DP)        :: time_run=0.0, final_time=10.0d0, ini_time=0.0d0    
  INTEGER(KIND=INT_KIND)   :: i_run = 0
  INTEGER(KIND=INT_KIND)   :: print_int = 1


! Output  
  INTEGER(KIND=INT_KIND)   :: output_method  = 1
                              !1  output every I steps
                              !2  output every I seconds
  INTEGER(KIND=INT_KIND)   :: output_int     = 1
  REAL(kind=REAL_DP)        :: output_real    = 1     
  
  REAL(kind=REAL_DP)       :: output_start_time = 0.0d0
  REAL(kind=REAL_DP)       :: output_end_time   = 0.0d0
  
  
  CHARACTER(len=400)       :: output_folder  = '.\';   
  LOGICAL                  :: output_to_bin  = .TRUE.
  LOGICAL                  :: output_to_txt  = .FALSE.
  LOGICAL                  :: enable_zeta_output  = .FALSE.     
  LOGICAL                  :: enable_vel_output   = .FALSE.
  LOGICAL                  :: enable_q_output     = .FALSE.
  LOGICAL                  :: enable_w_output     = .FALSE.
  INTEGER(KIND=INT_KIND)   :: bin_output_count    = 0 
  REAL(kind=REAL_DP)       :: output_real_count   = 0.0d0  
 
  INTEGER(KIND=INT_KIND)   :: maxmin_output_count    = 0 
  
  
  LOGICAL                  :: enable_advection = .FALSE.
  INTEGER(KIND=INT_KIND)   :: advection_scheme = 1
  
  
  LOGICAL                  :: enable_viscosity = .FALSE.
  INTEGER(KIND=INT_KIND)   :: viscosity_scheme = 2  
  
  INTEGER(KIND=INT_KIND)   :: flux_limiter_method = 1
    
  LOGICAL                  :: enable_friction = .FALSE.
  INTEGER(KIND=INT_KIND)   :: friction_method = 1
  REAL(kind=REAL_DP)       :: Cd  = 0.003
  REAL(kind=REAL_DP)       :: C_m = 0.012
  
  LOGICAL                  :: enable_coriolis = .FALSE.
  INTEGER(KIND=INT_KIND)   :: coriolis_method = 1 
  REAL(kind=REAL_DP)        :: sintheta_coriolis = 0.0d0
  
  
  
  LOGICAL                  ::  varied_time_step = .FALSE.
  REAL(kind=REAL_DP)       ::  max_CFL = 0.1
  REAL(kind=REAL_DP)       ::  h_f_min = 1.0e-3
  
  REAL(kind=REAL_DP)       ::  bath_max = 50.0
  
  
  
  ! Station output  
 
  LOGICAL                  :: enable_station_output = .FALSE.
  CHARACTER(len=400)       :: station_file  =''
  INTEGER                                              ::  n_station
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:)      ::  station_info
  INTEGER                                              ::  station_fid  
  INTEGER,            ALLOCATABLE, DIMENSION(:)        ::  station_el


  ! Station_vel output
  LOGICAL                  :: enable_station_vel_output = .FALSE.
  CHARACTER(len=400)       :: station_vel_file  =''
  INTEGER                                              ::  n_station_vel
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:)      ::  station_vel_info
  INTEGER                                              ::  station_vel_fid  
  INTEGER,            ALLOCATABLE, DIMENSION(:)        ::  station_vel_edg
  
  REAL(KIND=REAL_DP)                                   ::  station_output_real = 0.0
  REAL(KIND=REAL_DP)                                   ::  station_output_real_count= 0.0  
   
  
  
  ! Altim output  
  LOGICAL                  :: enable_altim_output = .FALSE.
  CHARACTER(len=400)       :: altim_file  =''
  INTEGER                                              ::  n_altim
  REAL(kind=REAL_DP), ALLOCATABLE, DIMENSION(:,:)      ::  altim_info
  INTEGER                                              ::  altim_fid  
  INTEGER,            ALLOCATABLE, DIMENSION(:)        ::  altim_el
  
  
  ! eta_max output
  LOGICAL                  :: enable_zeta_max_output  = .FALSE.
  CHARACTER(len=400)       :: zeta_max_file='' 
 
  LOGICAL                  :: enable_zeta_min_output  = .FALSE.
  CHARACTER(len=400)       :: zeta_min_file=''   
  
  LOGICAL                  :: enable_bottom_change    = .FALSE.
  
  
  INTEGER                  :: fid_time
  CHARACTER(len=400)       :: edg_upstream_el_file  ='edg_upstream_el.txt'
  
  REAL(kind=REAL_DP)       :: zetamax
  
  ! Two layer model parameters
  REAL(KIND = REAL_DP)        :: k_wave_predefine =0.0d0
  REAL(KIND = REAL_DP)        :: T_wave_predefine =0.0d0  
  REAL(KIND = REAL_DP)        :: fix_alpha = 0.15d0 
  LOGICAL                    :: enable_fixed_alpha = .TRUE. 
END  MODULE Mod_All_Variable   
