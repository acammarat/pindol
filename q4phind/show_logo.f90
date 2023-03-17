!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! q4phind. Copyright (C) Antonio Cammarata                            
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Order a q-point set as GM, H, S and remove possible complex-conjugated couples
! to generate qmatrix.nd and freq.nd compatible with phind v >= 4.3
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

subroutine show_logo
  use var, only: version

  write(*,'(a)')     "       _  _       _                   " 
  write(*,'(a)')     "  __ _| || |    _| |_              _  "
  write(*,'(a)')     " / _` | || |_  /     \   _ __   __| | "
  write(*,'(a)')     "| (_| |__   _|( (| |) ) | '_ \ / _` | "
  write(*,'(a)')     " \__, |  |_|   \_   _/  | | | | (_| | "
  write(*,'(a)')     "    |_|          |_|    |_| |_|\__,_| "
  write(*,'(*(a))')  "                              ",version
  write(*,*)
         
  
  return
end subroutine show_logo
