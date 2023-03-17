!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nd2xyz. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Converts an ND trajectory into xyz format
! 
! If used for production, you should cite
! Phys. Rev. B XX, XXXXX (XXXX)
! https://doi.org/10.1103/xxx
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

subroutine write_last
  use functions, only: i2a, f2a
  use var, only: infile_q, infile_qd, atoms_UC, ncells_tot,  &
                 nptotq, side_eq_EC, natoms_UC, &
                 at_pertype
  implicit none
  integer :: i, j
  real(8) :: T
  real(8) :: pos_EC(atoms_UC,3,ncells_tot), Q(nptotq)
  character(256) :: command

  if ( (trim(infile_q(:)) == '0') .and. (trim(infile_qd(:)) == '0') ) then
    write(*,'(a)') ' ERROR: ND trajectory and velocity files must be provided'
    write(*,*)
    stop
  end if

  write(*,'(a)') ' Writing last configuration and velocities...'
  
  write(command,'(*(a))') 'tail -n 1 ', trim(infile_q), ' > last_q.tmp'
  call system(command, i)
  if ( i /= 0 ) then
    write(*,'(a)') ' ERROR: Unable to extract the last q set.'
    write(*,*)
    stop
  end if
  write(command,'(*(a))') 'tail -n 1 ', trim(infile_qd), ' > last_qd.tmp'
  call system(command, i)
  if ( i /= 0 ) then
    write(*,'(a)') ' ERROR: Unable to extract the last qd set.'
    write(*,*)
    stop
  end if

  open(20,file='POSCAR_last.vasp',action='WRITE') 

  open(10,file='last_q.tmp',action='READ') 
  read(10,*) T, Q(:)
  close(10)
  write(command,'(a)') 'rm last_q.tmp'
  call system(command, i)

  call normal_to_cartesian ( Q, pos_EC, 'pos' )
  write(20,'(a)') ' Last ND configuration'
  write(20,'(a)') '1.0'
  do i = 1, 3
    write(20,'(3(f20.16,1x))') side_eq_EC(i,:)
  end do
  write(20,'(*(a,1x))') at_pertype(:)
  write(20,'(*(i5,1x))') natoms_UC(:)*ncells_tot
  write(20,'(a)') 'Cartesian'
  do i = 1, atoms_UC
    do j = 1, ncells_tot
      write(20,'(3(f20.16,1x))') pos_EC(i,:,j)
    end do
  end do

  open(11,file='last_qd.tmp',action='READ') 
  read(11,*) T, Q(:)
  close(11)
  write(command,'(a)') 'rm last_qd.tmp'
  call system(command, i)
  call normal_to_cartesian ( Q, pos_EC, 'vel' )
  write(20,*)
  write(20,'(a)') 'Cartesian'
  do i = 1, atoms_UC
    do j = 1, ncells_tot
      write(20,'(3(f20.16,1x))') pos_EC(i,:,j)
    end do
  end do
  close(20)
  
  write(*,'(*(a))') ' Last configuration written into POSCAR_last.vasp'

!  deallocate ( Q, stat = i )
!  if ( i /= 0 ) stop ' ERROR: Dellocation failed for Q'
  
  return

end subroutine write_last
