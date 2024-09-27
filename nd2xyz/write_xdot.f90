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

subroutine write_xdot
  use functions, only: i2a, f2a
  use var, only: atoms_UC, ncells_tot, &
                 infile_qd, nptotq, t1, t2
  implicit none
  integer :: i, j, iconf
  real(8) :: T
  real(8) :: Qd(nptotq), pos_EC(atoms_UC,3,ncells_tot)

  if ( trim(infile_qd(:))/='0' ) then

    open(11,file=infile_qd,action='READ')
    write(*,'(a)') ' Converting ND velocities -> XYZ ...'
    open(21,file='ndvel.xyz',action='WRITE') 

    iconf = -1
    reading: do
      read(11,*,iostat=i) T
      if ( i /= 0 ) exit
      iconf = iconf + 1
      if ( T>=t1 ) then
        if ( T>t2 ) exit reading
        backspace(11)
        read(11,*,iostat=i) T, Qd(:)
        write(21,'(*(a))') 'Configuration ', i2a(iconf),' , time ', trim(f2a(T))
        call normal_to_cartesian ( Qd, pos_EC, 'vel' )
        do i = 1, atoms_UC
          do j = 1, ncells_tot
            write(21,'(*(a))') trim(f2a(pos_EC(i,1,j))), ' ', trim(f2a(pos_EC(i,2,j))), ' ', trim(f2a(pos_EC(i,3,j)))
          end do
        end do
      end if
    end do reading
 
    write(*,'(*(a))') ' Velocities written into ndvel.xyz'
    close(11)
    close(21)

  end if

  return

end subroutine write_xdot
