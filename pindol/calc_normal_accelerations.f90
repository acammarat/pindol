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
! This subroutine computes the second-order time-derivatives of the normal coordinates
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine calc_normal_accelerations ( Y, Qdd )
  ! Note on units:
  !   [Y(1:nptotq)] = amu^1/2 Ang
  !   [Y(nptotq+1:2*nptotq)] = amu^1/2 Ang fs^-1
  !   [Y(2*nptotq+1)] = .
  !   [Y(2*nptotq+2)] = fs^-1
  !   [freq2] = fs^-2
  !   [Qdd] = amu^1/2 Ang fs^-2
  use flags, only: nvt_flag
  use eigen, only: freq2, ndof, nptotq, npdiffq, npGq, npHq, npuniqueq
  implicit none
  real(8), intent(in) :: Y(ndof)
  real(8), intent(out) :: Qdd(nptotq)
  integer :: l
  
  ! compute the third order part
  ! i.e. -1/2 ( sum_l1 sum_l2 phi_{l,l1,l2} Q_l1 Q_l2 )*
  call calc_double_sum_phi ( Y(1:nptotq), Qdd )

  ! compute the harmonic part, 
  ! i.e. -omega_l^2 Q_l
  
  ! Gamma-point, if present
  do l = 1, npGq
     Qdd(l) = Qdd(l) - freq2(l) * Y(l)
  end do

  ! set H, if present
  do l = npGq+1, npGq+npHq
     Qdd(l)         = Qdd(l)         - freq2(l) * Y(l)         ! real part
     Qdd(l+npdiffq) = Qdd(l+npdiffq) - freq2(l) * Y(l+npdiffq) ! imaginary part
  end do
  
  ! set S, if present
  do l = npGq+npHq+1, npuniqueq
     Qdd(l)         = Qdd(l)         - freq2(l) * Y(l)         ! real part
     Qdd(l+npdiffq) = Qdd(l+npdiffq) - freq2(l) * Y(l+npdiffq) ! imaginary part
  end do

  ! NVT dynamics
  if ( nvt_flag ) then

     ! system-thermostat coupling
     do l = 1, nptotq
        Qdd(l) = Qdd(l) - Y(2*nptotq+2) * Y(nptotq+l)
     end do
     
  end if

  return
end subroutine calc_normal_accelerations
