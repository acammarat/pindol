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
! This subroutine performs the projection from direct space to normal coordinates
! it takes as an input the current atomic positions or velocities or accelerations
! (together with the equilibrium positions in order to calculate displacements and phase factors)
! and it projects onto the subset of unique points in the reciprocal space
! Contributors:
!   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com
subroutine cartesian_to_normal ( direct, normal, tag )
  ! Note on units:
  !   [pos_eq_EC] = Ang
  !   [eigen_real] = .
  !   [eigen_imag] = .
  !   [mass_UC] = amu
  !   [vec] = Ang^-1
  ! if tag == 'pos':
  !   [direct] = Ang
  !   [normal] = amu^1/2 Ang
  ! if tag == 'vel':
  !   [direct] = Ang fs^-1
  !   [normal] = amu^1/2 Ang fs^-1
  ! if tag == 'acc':
  !   [direct] = Ang fs^-2
  !   [normal] = amu^1/2 Ang fs^-2
  use eigen, only: eigen_real, eigen_imag, vec, nq, npGq, nptotq, npuniqueq, npdiffq
  use refconf, only: atoms_UC, mass_UC
  use supercell, only: ncells_tot, pos_eq_EC
  implicit none
  real(8), intent(in) :: direct(atoms_UC,3,ncells_tot)
  character(3), intent(in) :: tag
  real(8), intent(out) :: normal(nptotq)
  integer :: i, l, k, n, ib, ie
  real(8) :: dotp, rexp, iexp, rdot, idot, sum_real, sum_imag

  ! initialize everything to zero
  normal(:) = 0.d0

  ! project cartesian components onto normal modes
  ! e.g. Q(l) = 1/sqrt(N0) * sum_(i runs over atoms in the unit cell, n runs over unit cells in the supercell)
  !                              [ sqrt(m_i) * exp( -i * k(l) . r^eq_i,n ) * displ_i,n . eigen*_i(l) ]

  ! at the Gamma-point (if present)
  do l = 1, npGq
     do i = 1, atoms_UC
        ib = (i-1) * 3 + 1
        ie = ib + 2
        sum_real = 0.d0
        do n = 1, ncells_tot
           if ( tag == 'pos' ) then
              sum_real = sum_real + dot_product( direct(i,:,n)-pos_eq_EC(i,:,n), eigen_real(l,ib:ie) )
           else
              sum_real = sum_real + dot_product( direct(i,:,n), eigen_real(l,ib:ie) )
           end if
        end do
        normal(l) = normal(l) + sqrt(mass_UC(i)) * sum_real
     end do
  end do

  ! the same for the others unique points, but here we need to consider real and imaginary parts
  do l = npGq+1, npuniqueq
     k = ceiling( dble(l) / dble(nq) ) 
     do i = 1, atoms_UC
        ib = (i-1) * 3 + 1
        ie = ib + 2
        sum_real = 0.d0
        sum_imag = 0.d0
        do n = 1, ncells_tot
           dotp = dot_product( vec(k,:), pos_eq_EC(i,:,n) )
           rexp = cos( dotp )
           iexp = -sin( dotp )
           if ( tag == 'pos' ) then
              rdot = dot_product( direct(i,:,n)-pos_eq_EC(i,:,n),  eigen_real(l,ib:ie) )
              idot = dot_product( direct(i,:,n)-pos_eq_EC(i,:,n), -eigen_imag(l,ib:ie) )
           else
              rdot = dot_product( direct(i,:,n),  eigen_real(l,ib:ie) )
              idot = dot_product( direct(i,:,n), -eigen_imag(l,ib:ie) )
           end if
           sum_real = sum_real + rexp * rdot - iexp * idot
           sum_imag = sum_imag + rexp * idot + iexp * rdot
        end do
        normal(l) = normal(l) + sqrt(mass_UC(i)) * sum_real
        normal(l+npdiffq) = normal(l+npdiffq) + sqrt(mass_UC(i)) * sum_imag
     end do
  end do

  ! divide by the factor sqrt(N0)
  do l = 1, nptotq
     normal(l) = normal(l) / sqrt(dble(ncells_tot))
  end do

  return
end subroutine cartesian_to_normal
