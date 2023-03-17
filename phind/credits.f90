!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! phind. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Phi for Normal Dynamics
! 
! Calculates the first order anharmonic interaction strength Phi(l,l',l'')
! to be used for the PINDOL code.
! The implemented selection rule is Delta(q+q'+q''); 
! the selection rules described in RSC Advances, 9, 37491 (2019)
! https://doi.org/10.1039/C9RA08294H
! are not implemented.
! Phi(l,l',l'') is not written for any triplet involving at least 
! one of the skip modes as specified in the header of freq.nd
!
! If used for production, you should cite
! Phys. Rev. B XX, XXXXX (XXXX)
! https://doi.org/10.1103/xxx
!
!    This file is part of phind.
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
  write(*,'(a)') ' Normal Dynamics: solving Newton’s equations in the reciprocal space'
  write(*,'(a)') ' A. Cammarata et al., Phys. Rev. Lett XX, XXXXX (XXXX)'
  write(*,'(a)') ' https://doi.org/10.1103/xxx'
  write(*,*)
  write(*,'(a)') ' Sampling dynamical trajectories in the reciprocal space'
  write(*,'(a)') ' A. Cammarata et al., Phys. Rev. B XX, XXXXX (XXXX)'
  write(*,'(a)') ' https://doi.org/10.1103/xxx'
  write(*,*)
  write(*,*)

end subroutine credits
