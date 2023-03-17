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

subroutine calc_phiun
  use var, only: atoms_UC, pos_eq_EC, pos_eq_UC, &
       nq, vec, eigen, mass_UC, atoms_EC, &
       ncells_tot, delta, ll, phi

  implicit none
  integer :: i, j, k, ix, iy, iz
  integer :: n2, n3, ix1, ix2, ix3 
  integer :: k1, k2, k3, l, l1, l2 
  integer(8) :: il
  real(8) :: t1, t2
  real(8) :: fc3(atoms_EC,atoms_EC,atoms_EC,3,3,3)
  complex(8) :: s1, s2, s3, s4
  logical :: file_exists


! The formula implemented here is a modified version of (10.74)
! from Wallace, Thermodynamics of Crystals,
! Chap. 3 "Lattice Dynamics", Section "Harmonic Phonons".
! The factor sqrt(hbar^3/(8 omega_l omega_l' omega_ll') has been removed
! according to the formulation in the supporting information of ND paper.
! The factor sqrt(1/N0^3) has been included 
! according to the formulation in the supporting information of ND paper.

  inquire(file='fc3.dat',exist=file_exists)
  if ( file_exists .eqv. .false. ) then
    write(*,'(a)') ' ERROR: force constant file fc3.dat not found.'
    write(*,*)
    stop
  end if

  write(*,'(a)',advance='no') ' Reading the force constant file fc3.dat...'
  open(unit=15,file='fc3.dat',action='READ')

  do i = 1, atoms_EC
    do j = 1, atoms_EC
      do k = 1, atoms_EC
        do ix = 1, 3
          do iy = 1, 3
            do iz = 1, 3
              read(15,*) fc3(i, j, k, ix, iy, iz) ! eV/Ang^3
            end do
          end do
        end do
      end do
    end do
  end do

  close(15)
  write(*,'(a)') ' done.'

  allocate ( phi(ll), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for PHI'

  write(*,'(a)',advance='no') ' Calculating phonon-phonon scattering tensor...'

  call cpu_time(t1)
  !$omp parallel do schedule ( dynamic ) private ( l, l1, l2, k1, k2, k3, i, j, k, n2, n3, s1, s2, s3, s4, ix1, ix2, ix3 )
  do il = 1, ll
    l  = delta(il,1)
    l1 = delta(il,2)
    l2 = delta(il,3)

    k1 = ceiling(dble(l)/dble(nq))
    k2 = ceiling(dble(l1)/dble(nq))
    k3 = ceiling(dble(l2)/dble(nq))

    s1 = dcmplx(0.d0,0.d0)
    do i = 1, atoms_UC
      s2 = dcmplx(0.d0,0.d0)
      do j = 1, atoms_UC
        do k = 1, atoms_UC
          s3 = dcmplx(0.d0,0.d0)
          do ix1 = 1, 3
            do ix2 = 1, 3
              do ix3 = 1, 3
                s4 = dcmplx(0.d0,0.d0)
    
                do n2 = 1, ncells_tot
                  do n3 = 1, ncells_tot
                    s4 = s4 + fc3((i-1)*ncells_tot+1,(j-1)*ncells_tot+n2,(k-1)*ncells_tot+n3,ix1,ix2,ix3) * &
                          zexp( dcmplx( 0.d0, dot_product( vec(k2,1:3), pos_eq_EC(j,1:3,n2) ) ) ) * &
                          zexp( dcmplx( 0.d0, dot_product( vec(k3,1:3), pos_eq_EC(k,1:3,n3) ) ) ) 
                  end do
                end do
                s3 = s3 + eigen(l,(i-1)*3+ix1) * eigen(l1,(j-1)*3+ix2) * &
                          eigen(l2,(k-1)*3+ix3) * s4
              end do
            end do
          end do
          s2 = s2 + mass_UC(j) * mass_UC(k) * s3
        end do
      end do
      s1 = s1 + zexp( dcmplx( 0.d0, dot_product( vec(k1,1:3), pos_eq_UC(i,1:3) ) ) ) * mass_UC(i) * s2 
      s1 = s1 / sqrt(dble(ncells_tot**3)) ! [s1]=[uma^(-3/2).eV.Ang^-3]
    end do
    phi(il)%kk(1) = l
    phi(il)%kk(2) = l1
    phi(il)%kk(3) = l2
    phi(il)%el(1) = dble(s1)
    phi(il)%el(2) = dimag(s1)

  end do
  !$omp end parallel do
  close(10)
  call cpu_time(t2)
  write(*,'(a)') ' done.'
  write(*,'(a,f18.3)') '  Cumulative (all threads) execution time of this step [s]: ', t2-t1

  return
end subroutine calc_phiun


