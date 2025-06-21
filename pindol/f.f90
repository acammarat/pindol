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
! This subroutine computes the time-derivative of the dynamical variables according to:
! d/dt Q_l = dot Q_l
! d/dt dot Q_l = -freq2 (Q_l)* - 1/2 ( sum_l1 sum_l2 phi_{l,l1,l2} Q_l1 Q_l2 )*
! and it is called by the dvode solver
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Antonio Cammarata (Czech Technical University in Prague), cammaant@fel.cvut.cz
!   - Miljan Dašić (Czech Technical University in Prague, Czech Republic & Institute of Physics Belgrade, Serbia), miljan.dasic@scl.rs
subroutine F (NEQ, time, Y, YDOT, RPAR, IPAR)
  ! Note on units:
  !   [Y(1:nptotq)] = amu^1/2 Ang
  !   [Y(nptotq+1:2*nptotq)] = amu^1/2 Ang fs^-1
  !   [Y(2*nptotq+1)] = .
  !   [Y(2*nptotq+2)] = fs^-1
  !   [YDOT(1:nptotq)] = amu^1/2 Ang fs^-1
  !   [YDOT(nptotq+1:2*nptotq)] = amu^1/2 Ang fs^-2
  !   [YDOT(2*nptotq+1)] = fs^-1
  !   [YDOT(2*nptotq+2)] = fs^-2
  !   [freq2_nd] = fs^-2
  !   [Ekin2] = amu Ang^2 fs^-2
  !   [kB_ieu] = amu Ang^2 fs^2 K^-1
  !   [*temperature*] = K  
  !   [thermomass] = amu Ang^2 K^-1
  use pars, only: kB_ieu
  use flags, only: nvt_flag
  use eigen, only: ndof_sub, freq2_nd, npdiffq, npGq, npuniqueq, nptotq, ndof_temperature
  use nd, only: thermomass, temperature_begin, temperature_end, time_total, time_init, temperature, temperature_ramp
  implicit none
  integer, intent(in) :: NEQ, IPAR
  real(8), intent(in) :: time, RPAR
  real(8), intent(in) :: Y(NEQ)
  real(8), intent(out) :: YDOT(NEQ)
  integer :: l
  real(8) :: Ekin2, target_temperature
  
  ! populating the first half (1:nptotq) of the YDOT vector
  do l = 1, nptotq
     YDOT(l) = Y(nptotq+l)
  end do

  ! compute the third order part
  call calc_double_sum_phi ( Y(1:nptotq), YDOT(nptotq+1:2*nptotq) )
  
  ! compute the harmonic part
  do l = 1, npGq ! Gamma-point
     YDOT(nptotq+l) = YDOT(nptotq+l) - freq2_nd(l) * Y(l)
  end do
  do l = npGq+1, npuniqueq ! set H + set S
     YDOT(nptotq+l) = YDOT(nptotq+l) - freq2_nd(l) * Y(l)
     YDOT(nptotq+npdiffq+l) = YDOT(nptotq+npdiffq+l) - freq2_nd(l) * Y(npdiffq+l)
  end do

  ! NVT dynamics
  if ( nvt_flag ) then

     ! system-thermostat coupling
     do l = 1, nptotq
        YDOT(nptotq+l) = YDOT(nptotq+l) - Y(2*nptotq+2) * Y(nptotq+l)
     end do
     
     ! computing twice the kinetic energy
     call calc_kinetic_energy ( Y(nptotq+1:2*nptotq), Ekin2 )
     Ekin2 = 2.d0 * Ekin2
     
     ! computing the current target temperature
     if ( temperature_ramp ) then
        target_temperature = temperature_begin + ( temperature_end - temperature_begin ) / time_total * ( time - time_init )
     else
        target_temperature = temperature
     end if

     ! updating thermostat variables
     YDOT(2*nptotq+1) = Y(2*nptotq+2)
     YDOT(2*nptotq+2) = ( Ekin2 - dble(ndof_temperature-ndof_sub) * kB_ieu * target_temperature ) / ( thermomass * target_temperature )
     
  end if
     
  return
end subroutine F

! Jacobian matrix not implemented
subroutine jac
end subroutine jac
