! pindol version 1.0.1, Copyright (C) 2023 P. Nicolini, A. Cammarata, M. Dašić
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
! This subroutine calculates the kinetic energy of the system according to
! 1/2 sum_l Qd_l (Qd_l)*
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine calc_kinetic_energy ( Qd, Ekin )
  ! Note on units:
  !   [Qd] = amu^1/2 Ang fs^-1
  !   [Ekin] = amu Ang^2 fs^-2
  use eigen, only: nptotq, npdiffq, npGq, npHq, npuniqueq
  implicit none
  real(8), intent(in) :: Qd(nptotq)
  real(8), intent(out) :: Ekin
  integer :: l

  ! init
  Ekin = 0.d0

  ! Gamma-point, if present
  do l = 1, npGq
     Ekin = Ekin + 0.5d0 * Qd(l) * Qd(l)
  end do

  ! set H, if present
  do l = npGq+1, npGq+npHq
     Ekin = Ekin + 0.5d0 * ( Qd(l) * Qd(l) + Qd(l+npdiffq) * Qd(l+npdiffq) )
  end do

  ! set S, if present (the factor 1/2 is omitted because the sum should run over both S and S*)
  do l = npGq+npHq+1, npuniqueq
     Ekin = Ekin + Qd(l) * Qd(l) + Qd(l+npdiffq) * Qd(l+npdiffq)
  end do

  return
end subroutine calc_kinetic_energy
