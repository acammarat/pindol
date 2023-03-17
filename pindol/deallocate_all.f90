! pindol version 1.0, Copyright (C) 2023 P. Nicolini, A. Cammarata, M. Dašić
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
! This subroutine frees all memory
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine deallocate_all
  use pars, only: error_string
  use flags, only: distort_flag, nd_flag, readrestart_flag
  use input, only: mass_pertype, at_pertype, distort_mode, distort
  use refconf, only: natoms_UC, mass_UC, pos_eq_UC, at_UC
  use eigen, only: vec, vec_red, freq2, freq2_nd, eigen_real, eigen_imag
  use dvodemod, only: RWORK, IWORK
  use supercell, only: mass_EC, pos_eq_EC, pos_EC, vel_EC, at_EC
  use nd, only: Y, restart_Y, phi, phi_map
  implicit none
  integer :: i

  write(*,'(a)',advance='no') 'Deallocating variables...'

  ! read_input_atoms
  deallocate ( mass_pertype, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array MASS_PERTYPE.'
     stop
  end if
  deallocate ( at_pertype, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array AT_PERTYPE.'
     stop
  end if

  ! read_input_distort
  if ( distort_flag .and. distort_mode == 'mode' ) then
     deallocate ( distort, stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Deallocation failed for array DISTORT.'
        stop
     end if
  end if

  ! read_refconf
  deallocate ( mass_UC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array MASS_UC.'
     stop
  end if
  deallocate ( at_UC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array AT_UC.'
     stop
  end if
  deallocate ( natoms_UC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array NATOMS_UC.'
     stop
  end if
  deallocate ( pos_eq_UC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array POS_EQ_UC.'
     stop
  end if

  ! read_eigen
  deallocate ( vec, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array VEC.'
     stop
  end if
  deallocate ( vec_red, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array VEC_RED.'
     stop
  end if
  deallocate ( freq2, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array FREQ2.'
     stop
  end if
  deallocate ( freq2_nd, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array FREQ2_ND.'
     stop
  end if
  deallocate ( eigen_real, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array EIGEN_REAL.'
     stop
  end if
  deallocate ( eigen_imag, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array EIGEN_IMAG.'
     stop
  end if

  if ( nd_flag ) then
     ! read_phi
     deallocate ( phi, stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Deallocation failed for array PHI.'
        stop
     end if
     deallocate ( phi_map, stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Deallocation failed for array PHI_MAP.'
        stop
     end if
     ! run_nd
     deallocate ( Y, stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Deallocation failed for array Y.'
        stop
     end if
     ! initialize_dvode
     deallocate ( RWORK, stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Deallocation failed for array RWORK.'
        stop
     end if
     deallocate ( IWORK, stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Deallocation failed for array IWORK.'
        stop
     end if
  end if

  ! setup_supercell
  deallocate ( mass_EC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array MASS_EC.'
     stop
  end if
  deallocate ( at_EC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array AT_EC.'
     stop
  end if
  deallocate ( pos_eq_EC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array POS_EQ_EC.'
     stop
  end if
  deallocate ( pos_EC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array POS_EC.'
     stop
  end if
  deallocate ( vel_EC, stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Deallocation failed for array VEL_EC.'
     stop
  end if
  
  ! read_restart
  if ( readrestart_flag ) then
     deallocate ( restart_Y, stat = i )
     if ( i /= 0 ) then
        write(0,*) error_string, 'Deallocation failed for array RESTART_Y.'
        stop
     end if
  end if
  
  write(*,'(a)') 'DONE.'
  
  return
end subroutine deallocate_all
