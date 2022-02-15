    PROGRAM H2Ocean
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
    USE Mod_All_Variable
    IMPLICIT NONE 


    CALL H2Ocean_Time_Stamp



    ! Display the version information of H2Ocean
    ! Such as the version number, the last update date
    CALL H2Ocean_Version_Info


    ! Read all the variables in the control file.
    CALL H2Ocean_Get_Ctrl_File
    CALL H2Ocean_Read_Ctrl_File


    !
    CALL H2Ocean_Read_Mesh_File


    ! The following steps will be taken over by different computation methods
    CALL H2Ocean_Initialize



    CALL H2Ocean_Time_Report


    WRITE(*,*) 'Staring Main Computation!!'
    dt_save = dt; 
    time_run  = ini_time
    i_run     = 0

    DO WHILE( final_time-time_run> 1.0e-7)
        i_run = i_run + 1


        IF(varied_time_step) THEN
            ! compute the new time step basing on CFL number    
            CALL H2Ocean_Compute_dt 
        ENDIF


        time_run = time_run + dt 

        CALL H2Ocean_Computation


        CALL H2Ocean_Output



        CALL H2Ocean_Report

    ENDDO      

    CALL H2Ocean_Finalize
    CALL H2Ocean_Time_Report       
    STOP
    END PROGRAM H2Ocean
