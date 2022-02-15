SUBROUTINE H2Ocean_Get_Ctrl_File
  USE Mod_Control_File
  IMPLICIT NONE
   
  ! Method 1: Specify the ctrl_file path here. 
!   ctrl_file='D:\Ocean\Programmes\H2Ocean_V2011_Final\workspace\standing_wave\control_file.txt';
  
  ! ctrl_file='D:\Ocean\Programmes\H2Ocean_V2011_Final\workspace\Japan_Tsunami_Pacific\control_file_01.txt';
  ! Method2: Pass the first command-line argument to the variable ctrl_file
  ! Method2 is used in most case. 
   CALL GET_COMMAND_ARGUMENT(1, ctrl_file)  
END SUBROUTINE H2Ocean_Get_Ctrl_File


SUBROUTINE H2Ocean_Read_Ctrl_File
  USE Mod_Control_File
  USE Mod_TriMesh
  USE Mod_All_Variable
  USE Mod_Dirichlet_Boundary
  IMPLICIT NONE
  
  WRITE(*,*)''
  WRITE(*,*)''
  WRITE(*,*)'+--------------------------------------------------------------------+'
  WRITE(*,*) 'H2Ocean are reading the following variables....'
  WRITE(*,*)'+--------------------------------------------------------------------+'
  
  WRITE(*,*) ''
  WRITE(*,*) '    +-----Reading Mesh Options....                 '
  ! Read the Coordinate_type, Cartesian Coordinate or in lat-lon Coordinate
  keyword = 'Coordinate_type'
  CALL Read_control_file(ctrl_file,keyword,Coordinate_type)   

  ! Read the node2Dfile, and element2Dfile
  keyword = 'node2Dfile'
  CALL Read_control_file(ctrl_file,keyword,node2Dfile)   

  keyword = 'element2Dfile'
  CALL Read_control_file(ctrl_file,keyword,element2Dfile) 
  
  
  
  WRITE(*,*) ''
  WRITE(*,*) '    +----- Reading Output Options....                 '
  keyword = 'output_folder'
  CALL Read_control_file(ctrl_file,keyword,output_folder)  

  keyword = 'output_method'
  CALL Read_control_file(ctrl_file,keyword,output_method)
  
  keyword = 'output_int'
  CALL Read_control_file(ctrl_file,keyword,output_int)
  
  keyword = 'output_real'
  CALL Read_control_file(ctrl_file,keyword,output_real)  

  keyword = 'output_start_time'
  CALL Read_control_file(ctrl_file,keyword,output_start_time)
    
  keyword = 'output_end_time'
  CALL Read_control_file(ctrl_file,keyword,output_end_time)  
 
  keyword = 'print_int'
  CALL Read_control_file(ctrl_file,keyword,print_int)

 
  keyword = 'enable_zeta_output'
  CALL Read_control_file(ctrl_file,keyword,enable_zeta_output)  

  keyword = 'enable_vel_output'
  CALL Read_control_file(ctrl_file,keyword,enable_vel_output) 

  keyword = 'enable_q_output'
  CALL Read_control_file(ctrl_file,keyword,enable_q_output) 
    
  keyword = 'enable_w_output'
  CALL Read_control_file(ctrl_file,keyword,enable_w_output)  
  
     
  keyword = 'output_to_txt'
  CALL Read_control_file(ctrl_file,keyword,output_to_txt)   

 
  keyword = 'output_to_bin'
  CALL Read_control_file(ctrl_file,keyword,output_to_bin) 
  
  keyword = 'enable_station_output'
  CALL Read_control_file(ctrl_file,keyword,enable_station_output) 

  keyword = 'station_file'
  CALL Read_control_file(ctrl_file,keyword,station_file)     
  
  keyword = 'enable_station_vel_output'
  CALL Read_control_file(ctrl_file,keyword,enable_station_vel_output) 

  keyword = 'station_vel_file'
  CALL Read_control_file(ctrl_file,keyword,station_vel_file)   
  
  keyword = 'station_output_real'
  CALL Read_control_file(ctrl_file,keyword,station_output_real)  
  
  
  keyword = 'enable_altim_output'
  CALL Read_control_file(ctrl_file,keyword,enable_altim_output) 

  keyword = 'altim_file'
  CALL Read_control_file(ctrl_file,keyword,altim_file)    
  
  
  keyword = 'enable_zeta_max_output'
  CALL Read_control_file(ctrl_file,keyword,enable_zeta_max_output)   
  
  keyword = 'zeta_max_file'
  CALL Read_control_file(ctrl_file,keyword,zeta_max_file)   
 
  
  WRITE(*,*) ''
  WRITE(*,*) '    +----- Reading Initial Options....'    
  
  
  keyword = 'enable_hot_start'
  CALL Read_control_file(ctrl_file,keyword,enable_hot_start)  
  
  keyword = 'bathymetryfile'
  CALL Read_control_file(ctrl_file,keyword,bathymetryfile) 
  
  keyword = 'init_zeta_file'
  CALL Read_control_file(ctrl_file,keyword,init_zeta_file) 
  
  keyword = 'init_vel_file'
  CALL Read_control_file(ctrl_file,keyword,init_vel_file) 

  keyword = 'hot_start_folder'
  CALL Read_control_file(ctrl_file,keyword,hot_start_folder)   
    
  
  WRITE(*,*) ''
  WRITE(*,*) '    +----- Reading Computation Options....                 ' 
  
  !Step1: Determine the Computation Method which will be used in the model.
 
  keyword = 'dt'
  CALL Read_control_file(ctrl_file,keyword,dt)   
  keyword = 'final_time'
  CALL Read_control_file(ctrl_file,keyword,final_time)   
  keyword = 'ini_time'
  CALL Read_control_file(ctrl_file,keyword,ini_time) 
    
  keyword = 'enable_reduced_twolayer'
  CALL Read_control_file(ctrl_file,keyword,enable_reduced_twolayer)   
  
  keyword = 'varied_time_step'
  CALL Read_control_file(ctrl_file,keyword,varied_time_step)  

  keyword = 'max_CFL'
  CALL Read_control_file(ctrl_file,keyword,max_CFL)  
    
 
  
  keyword = 'enable_nonhydro'
  CALL Read_control_file(ctrl_file,keyword,enable_nonhydro)  
  
  keyword = 'stencil_method'
  CALL Read_control_file(ctrl_file,keyword,stencil_method)   
  
  keyword = 'q_alpha'
  CALL Read_control_file(ctrl_file,keyword,q_alpha)  
  
  
  keyword = 'enable_advection'
  CALL Read_control_file(ctrl_file,keyword,enable_advection)  
  
  keyword = 'advection_scheme'
  CALL Read_control_file(ctrl_file,keyword,advection_scheme)    
  
  keyword = 'enable_viscosity'
  CALL Read_control_file(ctrl_file,keyword,enable_viscosity)   
  keyword = 'viscosity_scheme'
  CALL Read_control_file(ctrl_file,keyword,viscosity_scheme)   
  
  
  keyword = 'enable_coriolis'
  CALL Read_control_file(ctrl_file,keyword,enable_coriolis)  
  
  keyword = 'coriolis_method'
  CALL Read_control_file(ctrl_file,keyword,coriolis_method)
  
  keyword = 'sintheta_coriolis'
  CALL Read_control_file(ctrl_file,keyword,sintheta_coriolis)    
  
  
  keyword = 'flux_limiter_method'
  CALL Read_control_file(ctrl_file,keyword,flux_limiter_method) 
  
  
  keyword = 'enable_friction'
  CALL Read_control_file(ctrl_file,keyword,enable_friction)  
  
  keyword = 'friction_method'
  CALL Read_control_file(ctrl_file,keyword,friction_method)
  
  keyword = 'Cd'
  CALL Read_control_file(ctrl_file,keyword,Cd)    
  keyword = 'C_m'
  CALL Read_control_file(ctrl_file,keyword,C_m)
  
  keyword = 'h_min'
  CALL Read_control_file(ctrl_file,keyword,h_min)     
  
  keyword = 'h_f_min'
  CALL Read_control_file(ctrl_file,keyword,h_f_min)    
  
  
  keyword = 'enable_bottom_change'
  CALL Read_control_file(ctrl_file,keyword,enable_bottom_change)
  
  keyword = 'bottom_change_file'
  CALL Read_control_file(ctrl_file,keyword,bottom_change_file)
 
 
  keyword = 'test_case'
  CALL Read_control_file(ctrl_file,keyword,test_case)
  
  keyword = 'k_wave_predefine'
  CALL Read_control_file(ctrl_file,keyword,k_wave_predefine)

  keyword = 'T_wave_predefine'
  CALL Read_control_file(ctrl_file,keyword,T_wave_predefine)
  
  keyword = 'enable_fixed_alpha'
  CALL Read_control_file(ctrl_file,keyword,enable_fixed_alpha)
  
  keyword = 'fix_alpha'
  CALL Read_control_file(ctrl_file,keyword,fix_alpha)
  
      
  WRITE(*,*) ''
  WRITE(*,*) '    +----- Reading Boundary Condition Options....                 '
  keyword = 'dirichlet_boundary'
  CALL Read_control_file(ctrl_file,keyword,dirichlet_boundary)
  keyword = 'Dirichlet_file'
  CALL Read_control_file(ctrl_file,keyword,Dirichlet_file)
  
  keyword = 'u_dirichlet_boundary'
  CALL Read_control_file(ctrl_file,keyword,u_dirichlet_boundary)
  keyword = 'u_Dirichlet_file'
  CALL Read_control_file(ctrl_file,keyword,u_Dirichlet_file)  
  
  
  
  keyword = 'q_dirichlet_boundary'
  CALL Read_control_file(ctrl_file,keyword,q_dirichlet_boundary)
  keyword = 'q_Dirichlet_file'
  CALL Read_control_file(ctrl_file,keyword,q_Dirichlet_file)  
  
  
  

  WRITE(*,*)'+--------------------------------------------------------------------+'
  WRITE(*,*) 'H2Ocean finishes reading the ctrl file.'
  WRITE(*,*)'+--------------------------------------------------------------------+'
  WRITE(*,*)''
END SUBROUTINE H2Ocean_Read_Ctrl_File
