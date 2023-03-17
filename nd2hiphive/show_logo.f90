!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nd2hiphive. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Processes an ND trajectory and Creates input files for hiPhive 
! (https://hiphive.materialsmodeling.org/)
! to calculate the effective interatomic force constants
! 
! If used for production, you should cite
! Phys. Rev. B XX, XXXXX (XXXX)
! https://doi.org/10.1103/xxx
!
!    This file is part of nd2hiphive.
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

  write(*,'(a)')     "            _ ____  _     _       _     _           "
  write(*,'(a)')     "  _ __   __| |___ \| |__ (_)_ __ | |__ (_)_   _____ "
  write(*,'(a)')     " | '_ \ / _` | __) | '_ \| | '_ \| '_ \| \ \ / / _ \"
  write(*,'(a)')     " | | | | (_| |/ __/| | | | | |_) | | | | |\ V /  __/"
  write(*,'(a)')     " |_| |_|\__,_|_____|_| |_|_| .__/|_| |_|_| \_/ \___|"
  write(*,'(a)')     "                           |_|                      "
  write(*,'(*(a))')  "                                            ",version
  write(*,*)
  
  return
end subroutine show_logo
