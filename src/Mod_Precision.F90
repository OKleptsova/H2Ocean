!------------------------------------------------------------------------
! This module contains platform and problem specific parameters,
! and define floating point precision using kind intrinsic function.
!------------------------------------------------------------------------

!------------------------------------------------------------------------ 
! Copyright (C) 2009-2011 Technische Universiteit Delft,
!  Haiyang Cui, Guus Stelling and Julie Pietrzak
!
!  H.Cui@tudelft.nl	
!
!  Haiyang Cui
!  PhD student
!  Environmental Fluid Mechanics section
!  Faculty of Civil Engineering and Geosciences
!  Delft University of Technology
!  Delft, The Netherlands
!  Tel:   +31 15 278 5433
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
!------------------------------------------------------------------------
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
MODULE Mod_Precision
 IMPLICIT NONE

!--Single and Double Precision
   INTEGER, PARAMETER :: REAL_SP = SELECTED_REAL_KIND(6,30)  

!   INTEGER, PARAMETER :: REAL_DP = SELECTED_REAL_KIND(6,30)
   INTEGER, PARAMETER :: REAL_DP = SELECTED_REAL_KIND(6,70)
!   INTEGER, PARAMETER :: REAL_DP = SELECTED_REAL_KIND(15,307)
!--Kind of Integer for 8 digitals.
   INTEGER, PARAMETER :: INT_KIND = SELECTED_INT_KIND(8) 
END MODULE Mod_Precision
