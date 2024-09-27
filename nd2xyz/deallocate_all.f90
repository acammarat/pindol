!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nd2xyz. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Converts an ND trajectory into xyz format
! 
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of nd2xyz.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_all
  use var, only: natoms_UC, at_UC, mass_UC, &
       pos_eq_EC, at_EC, mass_pertype, at_pertype, pos_eq_UC, &
       eigen, vec, vec_red
  implicit none
  integer :: i
   
  deallocate ( mass_pertype, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for MASS_PERTYPE'
  deallocate ( at_pertype, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for AT_PERTYPE'
  deallocate ( natoms_UC, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for NATOMS_UC'
  deallocate ( pos_eq_UC, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for POS_EQ_UC'
  deallocate ( mass_UC, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for MASS_UC'
  deallocate ( at_UC, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for AT_UC'
  deallocate ( pos_eq_EC, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for POS_EQ_EC'
  deallocate ( at_EC, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for AT_EC'

  deallocate ( eigen, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for EIGEN'
  deallocate ( vec, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for VEC'
  deallocate ( vec_red, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for VEC_RED'
  
  return
end subroutine deallocate_all
