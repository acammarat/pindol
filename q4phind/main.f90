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
!
! v 0.3
! - corrected a bug when treating q-points at the border of the Brillouin zone
!
! v 0.2
! - added check on complex conjugated couples and identical vectors
!   for the border zone cases
! - added the generation of optimal mesh for maximizing the ph-ph interactions
!
! v 0.1
! 
! Antonio Cammarata
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Format of the input file
!
! double double double  reduced components of the first q-point
! ...
! double double double  reduced components of the last q-point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module var

  ! global parameters
  character(3), parameter :: version='0.3'
  character(7), parameter :: progname = 'Q4PHIND'
  ! init
  integer,save  :: npG, npH, npS, nptot, npinp, npopt
  real(8), save, allocatable :: vec_red(:,:), vec_inp(:,:)
  ! output
  real(8), save, allocatable :: vec_opt(:,:), vec_red_opt(:,:)

end module var

module functions
contains
  
  function i2a(i) result(out)
    character(:), allocatable :: out
    integer, intent(in) :: i
    character(range(i)+2) :: x

    write(x,'(i0)') i

    out = trim(x)
    
  end function i2a

end module functions

program q4phind

  call show_logo

  call init

  call check_qpt

  call write_qordered

  call write_qopt

  call credits

  stop
end program q4phind
