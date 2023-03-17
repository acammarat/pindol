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
! This subroutine reads the matrix of the third-order phonon interaction strength (aka phi)
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Antonio Cammarata (Czech Technical University in Prague), cammaant@fel.cvut.cz
subroutine read_phi
  ! Note on units:
  !   [phi] = amu^1/2 Ang fs^-2
  use pars, only: error_string, eV_to_ieu
  use nd, only: nphi, phi, phi_map, jl1, jl2
  use int2str
  use io_units, only: phi_unit, tmp_unit
  implicit none
  integer :: i, j, jl(24)
  character(1) :: comment
  character(256) :: md5sum_phi, md5sum_file

  write(*,'(a)') 'Reading matrix of the phonon interaction strength...'

  ! open the phi input file
  open( unit=phi_unit, file='phi.nd', action='read', status='old' )
  read(phi_unit,*) !# phind v. 4.3 

  ! check md5sum qmatrix
  call system ( "md5sum qmatrix.nd | awk '{print $1}' > tmp.nd", i )
  if ( i /= 0 ) then
    write(0,*) error_string, 'md5sum "qmatrix.nd" exited with non-zero status.'
    stop
  end if
  read(phi_unit,*) comment, comment, comment, md5sum_phi
  open( unit=tmp_unit, file='tmp.nd', action='read', status='old' )
  read(tmp_unit,*) md5sum_file
  if ( trim(md5sum_phi) /= trim(md5sum_file) ) then
     write(0,*) error_string, '"qmatrix.nd" is not consistent with "phi.nd".'
     stop
  end if
  close( tmp_unit )

  ! check md5sum freq
  call system ( "md5sum freq.nd | awk '{print $1}' > tmp.nd", i )
  if ( i /= 0 ) then
    write(0,*) error_string, 'md5sum "freq.nd" exited with non-zero status.'
    stop
  end if
  read(phi_unit,*) comment, comment, comment, md5sum_phi
  open( unit=tmp_unit, file='tmp.nd', action='read', status='old' )
  read(tmp_unit,*) md5sum_file
  if ( trim(md5sum_phi) /= trim(md5sum_file) ) then
    write(0,*) error_string, '"freq.nd" is not consistent with "phi.nd".'
    stop
  end if
  close( tmp_unit )

  ! housekeeping
  call system ( "rm -f tmp.nd", i )
  if ( i /= 0 ) then
    write(0,*) error_string, 'Unable to remove tmp.nd.'
    stop
  end if

  ! read the number of non-null elements
  read(phi_unit,*) nphi
  if ( nphi < 1 ) then
     write(0,*) error_string, 'NPHI cannot be smaller than 1.'
     stop
  end if

  ! read the number of elements in each block
  read(phi_unit,*) jl(1:24)
  jl1(1) = 1
  jl2(1) = jl(1)
  do i = 2, 24
      jl1(i) = jl2(i-1) + 1
      jl2(i) = jl2(i-1) + jl(i)
  end do

  ! allocate the vectors storing non-null elements
  allocate ( phi(nphi,2), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array PHI.'
     stop
  end if
  allocate ( phi_map(nphi,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array PHI_MAP.'
     stop
  end if

  ! read the phi file
  do i = 1, nphi
     read(phi_unit,*) (phi_map(i,j),j=1,3), (phi(i,j),j=1,2)
  end do

  ! convert into internal energy units (phi is hard-coded in eV)
  phi(:,:) = phi(:,:) * eV_to_ieu
  
  ! close the file
  close( phi_unit )

  ! printout
  write(*,'(2x,2a)') i2a(nphi), ' phi elements read'

  write(*,'(a)') 'Reading matrix of the phonon interaction strength...DONE.'

  return
end subroutine read_phi
