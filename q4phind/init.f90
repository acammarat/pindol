!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! q4phind. Copyright (C) Antonio Cammarata                            
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Program to calculate the atomic character
! Order a q-point set as GM, H, S and remove possible complex-conjugated couples
! to generate qmatrix.nd and freq.nd compatible with phind v >= 4.3
!
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of q4phind.
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

subroutine init
use functions, only: i2a
use var, only: npinp, vec_inp

  integer :: i
  real(8) :: dum
  character(256) :: infile
  logical :: file_exists
  
  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '') ) then
    write(*,'(a)') ' Syntax: q4phind <q-point set>'
    write(*,*)
    stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  write(*,'(*(a))') ' Reading file: ', trim(infile)

  open(unit=10,file=infile,action='read')

  npinp = 0
  do
    read(10,*,iostat=i) dum
    if ( i == 0 ) then
      npinp = npinp + 1
    else if ( i < 0 ) then
      exit
    end if
  end do
  rewind(10)
  write(*,'(*(a))') ' Number of q-points found in input file: ', i2a(npinp)

  allocate (vec_inp(npinp,3), stat=i)
  if (i/=0) stop 'Allocation failed for vec_inp'
  do i = 1, npinp
    read(10,*) vec_inp(i,:)
  end do
  close(10)

end subroutine init 
