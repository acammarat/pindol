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
! This subroutine calculates the harmonic potential energy of the system according to:
! 1/2 sum_l omega_l^2 Q_l (Q_l)*
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine calc_harmonic_potential_energy ( Q, Epot_harm )
  ! Note on units:
  !   [Q] = amu^1/2 Ang
  !   [freq2] = fs^-2
  !   [Epot_harm] = amu Ang^2 fs^-2
  use eigen, only: freq2, nptotq, npdiffq, npGq, npHq, npuniqueq
  implicit none
  real(8), intent(in) :: Q(nptotq)
  real(8), intent(out) :: Epot_harm 
  integer :: l

  ! init
  Epot_harm = 0.d0

  ! Gamma-point, if present
  do l = 1, npGq
     Epot_harm = Epot_harm + 0.5d0 * freq2(l) * Q(l) * Q(l)
  end do

  ! set H, if present
  do l = npGq+1, npGq+npHq
     Epot_harm = Epot_harm + 0.5d0 * freq2(l) * ( Q(l) * Q(l) + Q(l+npdiffq) * Q(l+npdiffq) )
  end do

  ! set S, if present  (the factor 1/2 is omitted because the sum should run over both S and S*)
  do l = npGq+npHq+1, npuniqueq
     Epot_harm = Epot_harm + freq2(l) * ( Q(l) * Q(l) + Q(l+npdiffq) * Q(l+npdiffq) )
  end do

  return
end subroutine calc_harmonic_potential_energy

! This subroutine calculates the anharmonic potential energy of the system according to:
! 1/6 sum_l sum_l1 sum_l2 phi_{l,l1,l2} Q_l Q_l1 Q_l2
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine calc_anharmonic_potential_energy ( Q, Epot_anh )
    ! Note on units:
  !   [Q] = amu^1/2 Ang
  !   [Qdd] = amu^1/2 Ang fs^-2
  !   [Epot_anh] = amu Ang^2 fs^-2
  use eigen, only: nptotq, npGq, npHq, npdiffq, npuniqueq
  implicit none
  real(8), intent(in) :: Q(nptotq)
  real(8), intent(out) :: Epot_anh
  integer :: l
  real(8) :: Qdd(nptotq)
  
  ! init
  Epot_anh = 0.d0

  ! compute the third order part: -1/2 ( sum_l1 sum_l2 phi_{l,l1,l2} Q_l1 Q_l2 )*
  call calc_double_sum_phi ( Q, Qdd )

  ! remove the minux sign and the complex-conjugate
  do l = 1, nptotq
     Qdd(l) = -Qdd(l)
  end do
  
  ! Gamma-point, if present
  do l = 1, npGq
     Epot_anh = Epot_anh + Q(l) * Qdd(l)
  end do

  ! set H, if present
  do l = npGq+1, npGq+npHq
     Epot_anh = Epot_anh + Q(l) * Qdd(l) - Q(l+npdiffq) * Qdd(l+npdiffq)
  end do

  ! set S, if present (the factor 2 comes from the fact that the sum should run over both S and S*)
  do l = npGq+npHq+1, npuniqueq
     Epot_anh = Epot_anh + 2.d0 * ( Q(l) * Qdd(l) - Q(l+npdiffq) * Qdd(l+npdiffq) )
  end do

  ! divide by the missing 1/3 factor
  Epot_anh = Epot_anh / 3.d0
  
  return
end subroutine calc_anharmonic_potential_energy
