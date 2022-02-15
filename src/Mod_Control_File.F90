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

MODULE Mod_Control_File 
  USE Mod_Precision
  IMPLICIT NONE
  !keyword, used to read the control files.
  CHARACTER(len=400)         :: keyword 
  CHARACTER(len=400)         :: ctrl_file=''
 
  INTERFACE   Read_control_file
      MODULE PROCEDURE Read_control_file_Int
      MODULE PROCEDURE Read_control_file_Double
      MODULE PROCEDURE Read_control_file_Logical
      MODULE PROCEDURE Read_control_file_Character
  END INTERFACE
  
  CONTAINS
    SUBROUTINE Read_control_file_Int(infile,keyword,output)
        IMPLICIT NONE
        CHARACTER(len=400),INTENT(IN)      :: infile
        CHARACTER*(*)                      :: keyword
        INTEGER(KIND=INT_KIND), INTENT(OUT):: output
        CHARACTER(len=400) :: buffer, label
        CHARACTER(len=400) :: values
        INTEGER(KIND=INT_KIND) :: val 
        INTEGER :: pos
        INTEGER :: ios = 0, ios1 = 0
        INTEGER :: line = 0
        INTEGER :: file_ctrl_file = 100
        ios = 0
        ios1 = 0
        keyword =TRIM(ADJUSTL(keyword))
        line = 0;
        OPEN(file_ctrl_file, file=infile,status="old",action='read')
  
  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

        DO WHILE (ios == 0)
            READ(file_ctrl_file, '(A)', iostat=ios) buffer
    
            IF (ios == 0) THEN
                line = line + 1

                ! Find the first instance of whitespace.  Split label and data.
                pos = SCAN(buffer, '#')
                IF(pos .NE. 0) THEN
                    CYCLE;
                ENDIF
                pos    = SCAN(buffer, '=')
                label  = TRIM(ADJUSTL(buffer(1:pos-1)))
                values = TRIM(ADJUSTL(buffer(pos+1:)))
                IF (label == keyword) THEN
                    READ(values, *, iostat=ios1) val 
                    output = val 
                    GOTO 100
                ENDIF
            END IF
        END DO
    100    CLOSE(file_ctrl_file)    
    WRITE(*,'(a30,a3,i4)') TRIM(ADJUSTL(keyword)),' = ', output    
      
    END SUBROUTINE Read_control_file_Int

    SUBROUTINE Read_control_file_Double(infile,keyword,output)
        IMPLICIT NONE
        CHARACTER(len=400),INTENT(IN)      :: infile
        CHARACTER*(*)                      :: keyword
        REAL(kind=REAL_DP),INTENT(OUT)     :: output
        CHARACTER(len=400)                 :: buffer, label
        CHARACTER(len=400)                 :: values
        REAL(kind=REAL_DP)                 :: val 
        INTEGER :: pos
        INTEGER :: ios = 0, ios1 = 0
        INTEGER :: line = 0
        INTEGER :: file_ctrl_file = 100
        ios = 0
        ios1 = 0
        keyword =TRIM(ADJUSTL(keyword))
        line = 0;
        OPEN(file_ctrl_file, file=infile,status="old",action='read')
  
  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

        DO WHILE (ios == 0)
            READ(file_ctrl_file, '(A)', iostat=ios) buffer
    
            IF (ios == 0) THEN
                line = line + 1

                ! Find the first instance of whitespace.  Split label and data.
                pos = SCAN(buffer, '#')
                IF(pos .NE. 0) THEN
                    CYCLE;
                ENDIF
                pos    = SCAN(buffer, '=')
                label  = TRIM(ADJUSTL(buffer(1:pos-1)))
                values = TRIM(ADJUSTL(buffer(pos+1:)))
                IF (label == keyword) THEN
                    READ(values, *, iostat=ios1) val 
                    output = val 
                    GOTO 100
                ENDIF
            END IF
        END DO
    100    CLOSE(file_ctrl_file)      
    WRITE(*,'(a30,a3,e15.5)') TRIM(ADJUSTL(keyword)),' = ', output    
    
    END SUBROUTINE Read_control_file_Double
    
    
    SUBROUTINE Read_control_file_Logical(infile,keyword,output)
        IMPLICIT NONE
        CHARACTER(len=400),INTENT(IN)      :: infile
        CHARACTER*(*)                      :: keyword
        LOGICAL, INTENT(OUT)               :: output
        CHARACTER(len=400) :: buffer, label
        CHARACTER(len=400) :: values
        
        LOGICAL     :: val  
        INTEGER :: pos
        INTEGER :: ios = 0, ios1 = 0
        INTEGER :: line = 0
        INTEGER :: file_ctrl_file = 100
        ios = 0
        ios1 = 0
        keyword =TRIM(ADJUSTL(keyword))
        
        OPEN(file_ctrl_file, file=infile,status="old",action='read')
  
  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.
        line = 0;
        DO WHILE (ios == 0)
            READ(file_ctrl_file, '(A)', iostat=ios) buffer
            IF (ios == 0) THEN
                line = line + 1
                ! Find the first instance of whitespace.  Split label and data.
                pos = SCAN(buffer, '#')
                IF(pos .NE. 0) THEN
                    CYCLE;
                ENDIF
                pos    = SCAN(buffer, '=')

                label  = TRIM(ADJUSTL(buffer(1:pos-1)))
                values = TRIM(ADJUSTL(buffer(pos+1:)))
                IF (label == keyword) THEN
                    READ(values, *, iostat=ios1) val
                    output = val
                    GOTO 100
                ENDIF
            END IF
        END DO
     100   CLOSE(file_ctrl_file)   
     WRITE(*,'(a30,a3,l)') TRIM(ADJUSTL(keyword)),' = ', output   
       
    END SUBROUTINE Read_control_file_Logical    
    
    SUBROUTINE Read_control_file_Character(infile,keyword,output)
        IMPLICIT NONE
        CHARACTER(len=400),INTENT(IN)      :: infile
        CHARACTER*(*)                      :: keyword
        CHARACTER*(*) , INTENT(INOUT)      :: output
        CHARACTER(len=400) :: buffer, label
        CHARACTER(len=400) :: values
        CHARACTER(len=400) :: val 
        
        INTEGER :: pos
        INTEGER :: ios = 0, ios1 = 0
        INTEGER :: line = 0
        INTEGER :: file_ctrl_file = 100
        ios = 0
        ios1 = 0
        keyword =TRIM(ADJUSTL(keyword))
        
        OPEN(file_ctrl_file, file=infile,status="old",action='read')
  
  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.
        line = 0;
        DO WHILE (ios == 0)
            READ(file_ctrl_file, '(A)', iostat=ios) buffer
            IF(buffer == '') THEN
                CYCLE
            ENDIF
            IF (ios == 0) THEN
                line = line + 1

                ! Find the first instance of whitespace.  Split label and data.
                pos = SCAN(buffer, '#')
                IF(pos .NE. 0) THEN
                    CYCLE
                ENDIF
                pos    = SCAN(buffer, '=')
                
                label  = TRIM(ADJUSTL(buffer(1:pos-1)))
                values = TRIM(ADJUSTL(buffer(pos+1:)))
                IF (TRIM(ADJUSTL(label)) .eq. TRIM(ADJUSTL(keyword))) THEN
!                    READ(values, *, iostat=ios1) val
                    val = values
                    output = TRIM(ADJUSTL(val))
                    GOTO 100
                ENDIF
            END IF
        END DO
        output = TRIM(ADJUSTL(output))
    100 CLOSE(file_ctrl_file)
    WRITE(*,'(a30,a3,a)') TRIM(ADJUSTL(keyword)),' = ', TRIM(ADJUSTL(output))
    END SUBROUTINE Read_control_file_Character        
END MODULE Mod_Control_File 
