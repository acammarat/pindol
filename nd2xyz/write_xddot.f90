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

subroutine write_xddot
  use functions, only: i2a, f2a
  use var, only: atoms_UC, ncells_tot, &
                 infile_qdd, nptotq, t1, t2
  implicit none
  integer :: i, j, iconf
  real(8) :: T
  real(8) :: Qdd(nptotq), pos_EC(atoms_UC,3,ncells_tot)

  if ( trim(infile_qdd(:))/='0' ) then

    open(12,file=infile_qdd,action='READ')
    write(*,'(a)') ' Converting ND accelerations -> XYZ ...'
    open(21,file='ndacc.xyz',action='WRITE') 

    iconf = -1
    reading: do
      read(12,*,iostat=i) T
      if ( i /= 0 ) exit
      iconf = iconf + 1
      if ( T>=t1 ) then
        if ( T>t2 ) exit reading
        backspace(12)
        read(12,*,iostat=i) T, Qdd(:)
        write(21,'(*(a))') 'Configuration ', i2a(iconf),' , time ', trim(f2a(T))
        call normal_to_cartesian ( Qdd, pos_EC, 'acc' )
        do i = 1, atoms_UC
          do j = 1, ncells_tot
            write(21,'(3(f20.16,1x))') pos_EC(i,:,j)
          end do
        end do
      end if
    end do reading
 
    write(*,'(*(a))') ' Accelerations written in ndacc.xyz'
    close(12)
    close(21)

  end if

  return

end subroutine write_xddot
