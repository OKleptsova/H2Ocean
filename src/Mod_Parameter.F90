!------------------------------------------------------------------------ 
! Copyright (C) 2009-2011 Technische Universiteit Delft,
!  Haiyang Cui, Guus Stelling and Julie Pietrzak
!  --|-----------------------------------------------------------|--
!    |  H.Cui@tudelft.nl	                                     |
!    |  Haiyang Cui                                              |
!    |  Environmental Fluid Mechanics section                    |
!    |  Faculty of Civil Engineering and Geosciences             |
!    |  Delft University of Technology                           |
!    |  Room 2.92,Building 23,Stevinweg 1,2628CN                 |
!    |  Delft, The Netherlands                                   |  
!    |  Tel:   +31 15 278 5433                                   |
!  --|-----------------------------------------------------------|--

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
!------------------------------------------------------------------------
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Haiyang Cui                               |
!   --|-----------------------------------------------------------|--
 
!---Authors--------------------
!
!   1.00: Haiyang Cui
!------------------------------

!---Updates--------------------
!
!------------------------------



!---Purpose--------------------
!   The TriMesh module. 
!   Main Subroutine is : TriMesh_Initialize
!   If you provide the node and element files, this subroutine will return
!   the derived data type Trimesh,  which inluces all the mesh information.
!-------------------------------
  
MODULE Mod_Parameter
  USE Mod_Precision
  IMPLICIT NONE
!--------------------------------------------------
!MODULE Mod_Parameter
!Some constant parametrs
!--------------------------------------------------  
  REAL(kind=REAL_DP), PARAMETER  :: pi      = 3.141592653589793
  REAL(kind=REAL_DP), PARAMETER  :: omega   = 2.0*pi/(24.0*60.0*60.0)
  REAL(kind=REAL_DP), PARAMETER  :: g       = 9.80616     ![m/s^2]
  REAL(kind=REAL_DP), PARAMETER  :: rad     = pi/180.0
  REAL(kind=REAL_DP), PARAMETER  :: r_earth = 6.3675e6  ![m]
END MODULE Mod_Parameter