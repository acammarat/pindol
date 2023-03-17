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

subroutine write_xyz
  use functions, only: i2a, f2a
  use var, only: infile_q, atoms_UC, ncells_tot, at_EC, atoms_EC, &
                 nptotq, t1, t2
  implicit none
  integer :: i, j, iconf, iave
  real(8) :: T
  real(8) :: pos_ave_EC(atoms_UC,3,ncells_tot)
  real(8) :: Q(nptotq), pos_EC(atoms_UC,3,ncells_tot)
 
  if ( trim(infile_q(:))/='0' ) then
  
    write(*,'(a)') ' Converting trajectory ND -> XYZ ...'
  
    open(20,file='ndtraj.xyz',action='WRITE') 
    open(10,file=infile_q,action='READ') 
  
    iconf = -1
    iave = 0
    pos_ave_EC(:,:,:) = 0.d0
    reading: do
      read(10,*,iostat=i) T
      if ( i /= 0 ) exit
      iconf = iconf + 1
      if ( T>=t1 ) then
        if ( T>t2 ) exit reading
        backspace(10)
        read(10,*,iostat=i) T, Q(:)
        iave = iave + 1
        write(20,'(a)') i2a(atoms_EC)
        write(20,'(*(a))') 'Configuration ', i2a(iconf),' , time ', trim(f2a(T))
        call normal_to_cartesian ( Q, pos_EC, 'pos' )
        do i = 1, atoms_UC
          do j = 1, ncells_tot
            write(20,'(*(a))') at_EC(i,j), ' ', trim(f2a(pos_EC(i,1,j))), ' ', trim(f2a(pos_EC(i,2,j))), ' ', trim(f2a(pos_EC(i,3,j)))
            pos_ave_EC(i,:,j) = pos_ave_EC(i,:,j) + pos_EC(i,:,j)
          end do
        end do
      end if
    end do reading
  
    if ( iave == 0 ) then
      write(*,'(a)') ' ERROR: no configurations found within the specified time interval.'
      stop
    end if
  
    write(*,'(*(a))') ' Trajectory written into ndtrj.xyz'
    close(20)
    close(10)
  
    pos_ave_EC(:,:,:) = pos_ave_EC(:,:,:) / dble(iave)
  
    open(20,file='ndtraj_ave.xyz',action='WRITE')
    write(20,'(a)') i2a(atoms_EC)
    write(20,'(*(a))') 'Ave. conf. time window: ', trim(f2a(t2-t1)), ' num. conf: ', i2a(iave)
    do i = 1, atoms_UC
      do j = 1, ncells_tot
        write(20,'(*(a))') at_EC(i,j), ' ', trim(f2a(pos_ave_EC(i,1,j))), ' ', trim(f2a(pos_ave_EC(i,2,j))), ' ', trim(f2a(pos_ave_EC(i,3,j)))
      end do
    end do
    close(20)

  end if

  return

end subroutine write_xyz
