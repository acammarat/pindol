!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! phind. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Phi for Normal Dynamics
! 
! Calculates the first order anharmonic interaction strength Phi(l,l',l'')
! to be used for the PINDOL code.
! The implemented selection rule is Delta(q+q'+q''); 
! the selection rules described in RSC Advances, 9, 37491 (2019)
! https://doi.org/10.1039/C9RA08294H
! are not implemented.
! Phi(l,l',l'') is not written for any triplet involving at least 
! one of the skip modes as specified in the header of freq.nd
!
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of phind.
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
       pos_eq_EC, at_EC, delta, skipmod, &
       mass_pertype, at_pertype, pos_eq_UC, &
       eigen, vec, vec_red
  implicit none
  integer :: i
   

  write(*,'(a)',advance='no') ' Deallocating variables...'
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
  deallocate ( skipmod, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for skipmod'
  deallocate ( delta, stat = i )
  if ( i /= 0 ) stop ' ERROR: deallocation failed for DELTA'
  
  write(*,'(a)') ' done.'

  return
end subroutine deallocate_all
