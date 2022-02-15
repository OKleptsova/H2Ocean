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


MODULE Mod_VERSION
!  
! DESCRIPTION: This module contains noting else than the version string for H2Ocean.
!
! REVISION HISTORY:
  IMPLICIT NONE
  SAVE
  CHARACTER(len=7), PARAMETER  :: program_name='H2Ocean'
  INTEGER, PARAMETER           :: program_version=1
  INTEGER, PARAMETER           :: program_subversion=0
  INTEGER, PARAMETER           :: program_revision=0
  CHARACTER(len=10), PARAMETER :: program_date= '16-11-2011'
  
END MODULE Mod_VERSION
  
SUBROUTINE H2Ocean_Version_Info
  USE Mod_VERSION
  IMPLICIT NONE
  WRITE(*,*) '+--------------------------------------------------------------------+'
  WRITE(*,*) '+                               ',program_name
  WRITE(*,'(a13,i1,a1,i1)') ' + Version : ',program_version,'.',program_subversion
  WRITE(*,'(a13,i3.3)')     ' + Revision: ',program_revision
  WRITE(*,'(a13,a10)')      ' + Date    : ',program_date  
  WRITE(*,*) '+--------------------------------------------------------------------+'    
END SUBROUTINE H2Ocean_Version_Info