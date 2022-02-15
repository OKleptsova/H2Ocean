MODULE Mod_IO
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
  IMPLICIT NONE
  ! File Unit
  INTEGER, PARAMETER      :: STDERR=0, STDOUT=6, STDIN=5
  INTEGER, PARAMETER      :: LOCALUNIT=3
  INTEGER, DIMENSION(4), PARAMETER :: RESERVED_UNIT=(/STDIN, STDERR,&
       & STDOUT, LOCALUNIT/)      
  ! The level of verbosity of warning and information messages.
  INTEGER, SAVE :: VERBOSITY=1
  
  CONTAINS
  
  FUNCTION get_file_unit()
    !----------------------------------------------------------------------
    ! Return a free unit number to connect a file to. The checking process
    ! does not scale well but:
    !
    ! a) you can only have 100 unit numbers anyway.
    !
    ! b) opening files is slow so a bit of extra time here is not so
    ! important.
    !
    ! If no unit numbers are available then get_unit prints an error and
    ! returns -1.
    !----------------------------------------------------------------------
    INTEGER :: get_file_unit
    ! By saving this we cycle the starting point through the available unit
    ! numbers. Each unit therefore gets checked with the same frequency.
    INTEGER, SAVE :: i_file=10
    INTEGER :: iold
    LOGICAL :: used

    ! This tells us when we have looped around.
    iold=i_file
    
    unitloop: DO 
       ! Only unit numbers up to 99 are permitted.
       i_file=MOD(i_file+1,100)
       IF (i_file==iold) THEN
          WRITE(*,*) 'No available units'
          get_file_unit=-1
          RETURN
       END IF

       ! Check if this unit number has been used
       INQUIRE(unit=i_file, opened=used)
       IF (used) CYCLE unitloop

       ! Check against known list of preconnected unit numbers. This might
       ! not be necessary but, on at least one compiler (ifc), preconnected
       ! units do not show up through the inquire function.
       IF (ANY(i_file==RESERVED_UNIT)) CYCLE unitloop

       ! If we fall through to here then this unit number must be OK
       EXIT unitloop
    END DO unitloop

    get_file_unit=i_file

  END FUNCTION get_file_unit
  
  SUBROUTINE INQUIRE_File( FileName, FEXIST )
  IMPLICIT NONE
  CHARACTER(LEN=400)  :: FileName
  LOGICAL             :: FEXIST
  
  INQUIRE(FILE=TRIM(FileName),EXIST=FEXIST)
  END SUBROUTINE INQUIRE_File


  SUBROUTINE Create_Folder( FolderName )
  USE IFPORT
  IMPLICIT NONE
  CHARACTER(LEN=400)  :: FolderName
  INTEGER :: FEXIST, RES, STATUS
    ! 1. First check if the output folder exist, if not, creat it
    INQUIRE(DIRECTORY=TRIM( FolderName ),EXIST=FEXIST)
    IF(.not.FEXIST) THEN
        WRITE(*,*) 'The output folder does not exist, now creat it!'
        res=MAKEDIRQQ(TRIM( FolderName ))
    ENDIF 
  END SUBROUTINE Create_Folder
END MODULE Mod_IO


