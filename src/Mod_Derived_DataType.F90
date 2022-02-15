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
!------------------------------
  
MODULE Mod_Derived_DataType
 USE Mod_Precision
 IMPLICIT NONE
 
!------------------------------------------------------------------------
! MODULE Mod_Address_DataType
! Define the user data type: address_type
! it has two components: nmb, addresses
! nmb indicates the total number of entries.
! address(1:nmb) store the value of the entries.
!------------------------------------------------------------------------
 TYPE address_type
    INTEGER(KIND=INT_KIND)                           :: nmb
    INTEGER(KIND=INT_KIND), ALLOCATABLE,DIMENSION(:) :: addresses
 END TYPE address_type
 
 TYPE vector_type
    REAL(KIND=REAL_DP) :: dx      ! dx, point from point 1 to 2
    REAL(KIND=REAL_DP) :: dy      ! dy, point from point 1 to 2
    REAL(KIND=REAL_DP) :: length  ! length
 END TYPE vector_type
 
END MODULE Mod_Derived_DataType