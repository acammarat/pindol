!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nd2xyz. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Converts an ND trajectory into xyz format
! 
! If used for production, you should cite
! A. Cammarata, M. Dasic, P. Nicolini, J. Chem. Phys. 161, 084111 (2024)
! https://doi.org/10.1063/5.0224108
!
!    This file is part of nd2xyz.
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

! this subroutine performs the backprojection from normal coordinates to direct space 
! it takes as an input the current components of the subset of unique points in the reciprocal space
! (together with the equilibrium positions in order to calculate phase factors and absolute positions from displacements)
! and it backprojects onto atomic positions and velocities
! note on units:
!   [pos_eq_EC] = Ang
!   [eigen] = .
!   [mass_UC] = amu
!   [vec] = Ang^-1
! if tag == pos
!   [normal] = amu^1/2 Ang
!   [direct] = Ang
! if tag == vel
!   [normal] = amu^1/2 Ang fs^-1
!   [direct] = Ang fs^-1
! if tag == acc
!   [normal] = amu^1/2 Ang fs^-2
!   [direct] = Ang fs^-2
subroutine normal_to_cartesian ( normal, direct, tag )
  use var, only: atoms_UC, ncells_tot, mass_UC, pos_eq_EC, eigen, vec, nq, npGq, npHq, nptotq, npdiffq, npuniqueq
  implicit none
  real(8), intent(in) :: normal(nptotq)
  character(3), intent(in) :: tag
  real(8), intent(out) :: direct(atoms_UC,3,ncells_tot)
  integer :: i, l, k, n, ib, ie
  real(8) :: dotp, rexp, iexp, rQ, iQ, reig(3), ieig(3)

  
  ! initialize everything to zero
  direct(:,:,:) = 0.d0

  ! project normal modes onto cartesian components
  ! e.g. displ_i,n = 1 / sqrt(N0*m_i) * sum_(l runs over normal modes) Q(l) * exp( i * k . r_in,0 ) * eigen*_i(l)

  ! at the Gamma-point (if present)
  do i = 1, atoms_UC
     ib = (i-1) * 3 + 1
     ie = ib + 2
     do n = 1, ncells_tot
        do l = 1, npGq  
           direct(i,:,n) = direct(i,:,n) + normal(l) * dble(eigen(l,ib:ie))
        end do
     end do
  end do

  ! at q-points belonging to the set H
  do i = 1, atoms_UC
     ib = (i-1) * 3 + 1
     ie = ib + 2
     do l = npGq+1, npGq+npHq
        reig(:) = dble( eigen(l,ib:ie) )
        ieig(:) = imag( eigen(l,ib:ie) )
        rQ = normal(l)
        iQ = normal(l+npdiffq)
        do n = 1, ncells_tot
           k = ceiling( dble(l) / dble(nq) ) 
           dotp = dot_product( vec(k,:), pos_eq_EC(i,:,n) )
           rexp = cos( dotp )
           iexp = sin( dotp )
           direct(i,:,n) = direct(i,:,n) + rexp*rQ*reig(:) - rexp*iQ*ieig(:) - iexp*rQ*ieig(:) - iexp*iQ*reig(:)
        end do
     end do
  end do
  
  ! at q-points belonging to the set S (note that the factor 2 accounts for the set S*)
  do i = 1, atoms_UC
     ib = (i-1) * 3 + 1
     ie = ib + 2
     do l = npGq+npHq+1, npuniqueq
        reig(:) = dble( eigen(l,ib:ie) )
        ieig(:) = imag( eigen(l,ib:ie) )
        rQ = normal(l)
        iQ = normal(l+npdiffq)
        do n = 1, ncells_tot
           k = ceiling( dble(l) / dble(nq) ) 
           dotp = dot_product( vec(k,:), pos_eq_EC(i,:,n) )
           rexp = cos( dotp )
           iexp = sin( dotp )
           direct(i,:,n) = direct(i,:,n) + 2.d0 * ( rexp*rQ*reig(:) - rexp*iQ*ieig(:) - iexp*rQ*ieig(:) - iexp*iQ*reig(:) )
        end do
     end do
  end do

  ! TBAC - divide by the factor sqrt(N0*m_i)
  do i = 1, atoms_UC
     direct(i,:,:) = direct(i,:,:) / sqrt(mass_UC(i)*dble(ncells_tot))
  end do

  ! in case of positions, add the equilibrium configuration
  if ( tag == 'pos' ) then
     direct(:,:,:) = direct(:,:,:) + pos_eq_EC(:,:,:)
  end if
  
  return
end subroutine normal_to_cartesian
