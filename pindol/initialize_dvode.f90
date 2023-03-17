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
! This subroutine initialize all variables needed by the dvode solver
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
!   - Miljan Dašić (Czech Technical University in Prague, Czech Republic & Institute of Physics Belgrade, Serbia), miljan.dasic@scl.rs
subroutine initialize_dvode
  use pars, only: error_string
  use eigen, only: ndof
  use dvodemod
  implicit none
  integer :: i

  ! real work array (LRW+RWORK)
  if ( dvode_method == 10 ) then
     LRW = 20 + 16 * ndof
  else if ( dvode_method == 21 .or. dvode_method == 22 ) then
     LRW = 22 + 9 * ndof + 2 * ndof*ndof
  end if
  allocate ( RWORK(LRW), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array RWORK.'
     stop
  end if
  
  ! integer work array (LIW+IWORK)
  if ( dvode_method == 10 ) then
     LIW = 30
  else if ( dvode_method == 21 .or. dvode_method == 22 ) then
     LIW = 30 + ndof
  end if
  allocate ( IWORK(LIW), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array IWORK.'
     stop
  end if

  ! array of tolerance not implemented
  ITOL = 1

  ! tolerance parameters
  if ( dvode_tol_flag == 'rel' ) then
     RTOL = dvode_tol
     ATOL = 0.d0
  else if ( dvode_tol_flag == 'abs' ) then
     RTOL = 0.d0
     ATOL = dvode_tol
  else if ( dvode_tol_flag == 'both' ) then
     RTOL = dvode_tol
     ATOL = dvode_tol
  end if

  ! ITASK = 1 for normal computation of output values of Y at t = TOUT.
  ITASK = 1 

  ! IOPT = 0 to indicate no optional input used.
  IOPT = 0 

  ! ISTATE = 1 initial value
  ISTATE = 1
  
  return
end subroutine initialize_dvode
