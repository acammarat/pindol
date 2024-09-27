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

subroutine credits
  use var, only: progname

  write(*,*)
  write(*,'(*(a))') ' Suggested reference for the acknowledgment of ',progname,' usage.'
  write(*,*)
  write(*,'(*(a))') ' The users of ',progname, ' have little formal obligations'
  write(*,'(a)') ' specified in the GNU General Public License, http://www.gnu.org/copyleft/gpl.txt .'
  write(*,'(a)') ' However, it is common practice in the scientific literature,'
  write(*,'(a)') ' to acknowledge the efforts of people that have made the research possible.'
  write(*,'(a)') ' In this spirit, please cite '
  write(*,*)
  write(*,'(a)') ' A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)'
  write(*,'(a)') ' https://doi.org/10.1063/5.0224108'
  write(*,*)
  write(*,*)

end subroutine credits
