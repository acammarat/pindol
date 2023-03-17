! pindol version 1.0, Copyright (C) 2023 P. Nicolini, A. Cammarata, M. Dašić
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
! This subroutine prints the suggested acknowledgements
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Antonio Cammarata (Czech Technical University in Prague), cammaant@fel.cvut.cz
subroutine credits
  use pars, only: version

  write(*,*)
  write(*,'(a)') ' Suggested reference for the acknowledgment of PINDOL '//version//' usage.'
  write(*,*)
  write(*,'(a)') ' The users of PINDOL '//version//' have little formal obligations'
  write(*,'(a)') ' specified in the GNU General Public License, https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt.'
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

end subroutine credits
